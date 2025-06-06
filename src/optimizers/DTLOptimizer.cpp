#include <IO/Logger.hpp>
#include <algorithm>
#include <cmath>
#include <corax/optimize/opt_generic.h>
#include <iomanip>
#include <iostream>
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <limits>
#include <optimizers/DTLOptimizer.hpp>
#include <parallelization/ParallelContext.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
#ifdef WITH_GSL
#include <gsl/gsl_multimin.h>
#endif

static bool isValidLikelihood(double ll) {
  return std::isnormal(ll) && ll < -0.0000001;
}

static bool lineSearchParameters(FunctionToOptimize &function,
                                 Parameters &currentRates,
                                 const Parameters &gradient,
                                 unsigned int &llComputationsLine,
                                 OptimizationSettings &settings) {
  // alpha controls the current size step
  double alpha = settings.startingAlpha;
  // we stop the search when alpha < minAlpha
  const double minAlpha = settings.minAlpha;
  bool noImprovement = true;
  if (settings.verbose) {
    Logger::info << "lineSearch from ll=" << currentRates.getScore()
                 << std::endl;
    Logger::info << "gradient=" << gradient << std::endl;
  }
  while (alpha > minAlpha) {
    Parameters normalizedGradient(gradient);
    normalizedGradient.normalize(alpha);
    Parameters proposal = currentRates + (normalizedGradient * alpha);
    function.evaluate(proposal);
    llComputationsLine++;
    auto improvement = proposal.getScore() - currentRates.getScore();
    if (improvement > 0) {
      if (settings.verbose) {
        Logger::info << "improv: alpha=" << alpha << ", params=" << proposal
                     << std::endl;
      }
      currentRates = proposal;
      alpha *= 1.5;
      if (improvement > settings.lineSearchMinImprovement) {
        noImprovement = false;
      }
    } else {
      if (settings.verbose) {
        Logger::info << "no improv: alpha=" << alpha << ", params=" << proposal
                     << std::endl;
      }
      alpha *= 0.5;
      if (!noImprovement && proposal.dimensions() > 1) {
        // it's time to recompute the gradient, unless we only have one
        // dimension
        return true;
      }
    }
  }
  return !noImprovement;
}

struct TargetParam {
  FunctionToOptimize *function;
  unsigned int n;
  bool verbose;
};

double myTargetFunction(void *function, double *value) {
  auto targetParam = (TargetParam *)(function);
  auto f = (FunctionToOptimize *)(targetParam->function);
  unsigned int n = targetParam->n;
  bool verbose = targetParam->verbose;
  Parameters param(n);
  for (unsigned int i = 0; i < n; ++i) {
    param[i] = value[i];
  }
  auto v = f->evaluate(param);
  if (verbose) {
    Logger::info << "params=" << param << std::endl;
  }
  return -v;
}

Parameters optimizeParametersLBFGSB(FunctionToOptimize &function,
                                    const Parameters &startingParameters,
                                    OptimizationSettings settings) {
  unsigned int n = startingParameters.dimensions();
  float lb = 1.0e-10;
  float ub = 2.0;
  std::vector<double> xmin(n, lb);
  std::vector<double> xmax(n, ub);
  std::vector<double> x(n, 0.5);
  for (unsigned int i = 0; i < n; ++i) {
    x[i] = startingParameters[i];
  }
  std::vector<int> bound(n, CORAX_OPT_LBFGSB_BOUND_BOTH);
  TargetParam targetFunction;
  targetFunction.function = &function;
  targetFunction.n = startingParameters.dimensions();
  targetFunction.verbose = settings.verbose;
  void *params = &targetFunction;
  // float factr = parser.getValue("lbfgsb.factr");
  // float pgtol = parser.getValue("lbfgsb.pgtol");
  float factr = static_cast<float>(settings.factr);
  float pgtol = 0.001;
  if (settings.verbose) {
    Logger::timed << "Starting LBFGSB search" << std::endl;
  }
  corax_opt_minimize_lbfgsb(&x[0], &xmin[0], &xmax[0], &bound[0], n,
                            factr, // convergence tolerance
                            pgtol, // gradient epsilon
                            params, myTargetFunction);
  Parameters res(n);
  for (unsigned int i = 0; i < n; ++i) {
    res[i] = x[i];
  }
  function.evaluate(res);
  if (settings.verbose) {
    Logger::timed << "opt_params=" << res << std::endl;
  }
  return res;
}

class FunctionOneDim : public FunctionToOptimize {
public:
  FunctionOneDim(const Parameters &parameters, unsigned int index,
                 FunctionToOptimize &fun)
      : _parameters(parameters), _index(index), _fun(fun) {
    assert(_parameters.dimensions() != 1);
  }

  virtual ~FunctionOneDim() {};
  virtual double evaluate(Parameters &parameters) {
    assert(parameters.dimensions() == 1);
    _parameters[_index] = parameters[0];
    auto res = _fun.evaluate(_parameters);
    parameters[0] = _parameters[_index];
    parameters.setScore(res);
    return res;
  }

private:
  Parameters _parameters;
  unsigned int _index;
  FunctionToOptimize &_fun;
};

static Parameters
optimizeParametersGradient(FunctionToOptimize &function,
                           const Parameters &startingParameters,
                           OptimizationSettings settings) {
  if (startingParameters.dimensions() == 0) {
    return Parameters();
  }
  double epsilon = settings.epsilon;
  Parameters currentRates = startingParameters;
  function.evaluate(currentRates);
  unsigned int llComputationsGrad = 0;
  unsigned int llComputationsLine = 0;
  unsigned int dimensions = startingParameters.dimensions();
  Parameters gradient(dimensions);
  if (settings.verbose) {
    Logger::timed << "Starting gradient descent search" << std::endl;
    Logger::info << "gradient epsilon=" << epsilon << std::endl;
  }
  bool stop = false;
  while (!stop) {
    std::vector<Parameters> closeRates(dimensions, currentRates);
    for (unsigned int i = 0; i < dimensions; ++i) {
      Parameters closeRates = currentRates;
      closeRates[i] += epsilon;
      function.evaluate(closeRates);
      llComputationsGrad++;
      gradient[i] =
          (currentRates.getScore() - closeRates.getScore()) / (-epsilon);
    }
    double oldScore = currentRates.getScore();
    stop |= !lineSearchParameters(function, currentRates, gradient,
                                  llComputationsLine, settings);
    stop |= (currentRates.getScore() - oldScore) <
            settings.optimizationMinImprovement;
    if (!stop) {
      settings.onBetterParametersFoundCallback();
    }
  }
  auto res = currentRates;
  function.evaluate(res);
  if (settings.verbose) {
    Logger::timed << "opt_params=" << res << std::endl;
  }
  return res;
}

static Parameters
optimizeParametersIndividually(FunctionToOptimize &function,
                               const Parameters &startingParameters,
                               OptimizationSettings settings) {
  if (startingParameters.dimensions() <= 1) {
    return startingParameters;
  }
  settings.verbose = false;
  const unsigned int N = startingParameters.dimensions();
  Parameters currentParameters(startingParameters);
  settings.optimizationMinImprovement = 1000000.0; // we only want one iteration
  if (settings.verbose) {
    Logger::timed << "Starting individual parameter optimization on " << N
                  << " parameters" << std::endl;
  }
  for (unsigned int i = 0; i < N; ++i) {
    FunctionOneDim fun(currentParameters, i, function);
    Parameters individualParam(1);
    assert(individualParam.dimensions() == 1);
    individualParam[0] = currentParameters[i];
    individualParam.setScore(currentParameters.getScore());
    individualParam =
        DTLOptimizer::optimizeParameters(fun, individualParam, settings);
    if (settings.verbose) {
      Logger::info << "Individual opt i=" << i
                   << ": pbefore=" << currentParameters[i]
                   << ", pafter=" << individualParam[0] << ", llDiff="
                   << individualParam.getScore() - currentParameters.getScore()
                   << std::endl;
    }
    currentParameters[i] = individualParam[0];
    currentParameters.setScore(individualParam.getScore());
    settings.onBetterParametersFoundCallback();
  }
  return currentParameters;
}

class PerCoreFunction : public FunctionToOptimize {
public:
  PerCoreFunction(PerCoreEvaluations &evaluations)
      : _evaluations(evaluations) {}
  virtual double evaluate(Parameters &parameters) {
    parameters.ensurePositivity();
    double ll = 0.0;
    for (auto evaluation : _evaluations) {
      evaluation->setRates(parameters);
      ll += evaluation->evaluate();
    }
    ParallelContext::sumDouble(ll);
    if (!isValidLikelihood(ll)) {
      ll = -std::numeric_limits<double>::infinity();
    }
    parameters.setScore(ll);
    return ll;
  }

private:
  PerCoreEvaluations &_evaluations;
};

Parameters
DTLOptimizer::optimizeParameters(PerCoreEvaluations &evaluations,
                                 const Parameters &startingParameters,
                                 OptimizationSettings settings) {
  PerCoreFunction function(evaluations);
  return optimizeParameters(function, startingParameters, settings);
}

ModelParameters DTLOptimizer::optimizeModelParameters(
    PerCoreEvaluations &evaluations, bool optimizeFromStartingParameters,
    const ModelParameters &startingParameters, OptimizationSettings settings) {

  ModelParameters res = startingParameters;
  if (!startingParameters.info.perFamilyRates) {
    const Parameters *startingRates =
        optimizeFromStartingParameters ? &startingParameters.rates : nullptr;
    res.rates = DTLOptimizer::optimizeParametersGlobalDTL(
        evaluations, startingRates, settings);
  } else {
    ParallelContext::pushSequentialContext(); // work locally
    for (unsigned int i = 0; i < evaluations.size(); ++i) {
      Parameters localRates = startingParameters.getRates(i);
      const Parameters *startingRates =
          optimizeFromStartingParameters ? &localRates : nullptr;
      PerCoreEvaluations localEvaluation;
      localEvaluation.push_back(evaluations[i]);
      localRates = DTLOptimizer::optimizeParametersGlobalDTL(
          localEvaluation, startingRates, settings);
      res.setRates(i, localRates);
    }
    ParallelContext::popContext();
  }
  return res;
}

Parameters
DTLOptimizer::optimizeParametersGlobalDTL(PerCoreEvaluations &evaluations,
                                          const Parameters *startingParameters,
                                          OptimizationSettings settings) {
  unsigned int freeParameters = 0;
  if (evaluations.size()) {
    freeParameters = evaluations[0]->getRecModelInfo().modelFreeParameters();
  }
  ParallelContext::maxUInt(freeParameters);
  if (freeParameters == 0) {
    return Parameters();
  }
  std::vector<Parameters> startingRates;
  if (startingParameters) {
    startingRates.push_back(*startingParameters);
  }
  if (freeParameters == 1) {
    Parameters p(1);
    p[0] = 0.1;
    startingRates.push_back(p);
    p[0] = 0.3;
    startingRates.push_back(p);
    p[0] = 1.0;
    startingRates.push_back(p);
    p[0] = 10.0;
    startingRates.push_back(p);
  } else if (freeParameters == 2) {
    startingRates.push_back(Parameters(0.1, 0.2));
    startingRates.push_back(Parameters(0.2, 0.2));
    startingRates.push_back(Parameters(0.5, 0.5));
    startingRates.push_back(Parameters(0.5, 1.0));
    startingRates.push_back(Parameters(0.01, 0.01));
  } else if (freeParameters == 3) {
    startingRates.push_back(Parameters(0.1, 0.2, 0.1));
    startingRates.push_back(Parameters(0.01, 0.01, 0.01));
  } else {
    startingRates.push_back(Parameters(0.5, 0.5, 0.2, 0.01));
    startingRates.push_back(Parameters(0.1, 0.2, 0.1, 0.1));
    startingRates.push_back(Parameters(0.2, 0.2, 0.0, 0.1));
    startingRates.push_back(Parameters(0.01, 0.01, 0.01, 0.01));
  }
  ParallelContext::barrier();
  Parameters best;
  best.setScore(-10000000000);
  for (auto rates : startingRates) {
    Parameters newRates = optimizeParameters(evaluations, rates, settings);
    bool stop = (fabs(newRates.getScore() - best.getScore()) <
                 settings.optimizationMinImprovement);
    stop = false;
    if (newRates.getScore() > best.getScore()) {
      best = newRates;
    }
    if (stop) {
      break;
    }
  }
  return best;
}

Parameters
DTLOptimizer::optimizeParametersPerSpecies(PerCoreEvaluations &evaluations,
                                           unsigned int speciesNodesNumber) {
  Parameters globalRates = optimizeParametersGlobalDTL(evaluations);
  Parameters startingSpeciesRates(speciesNodesNumber, globalRates);
  Parameters rates =
      DTLOptimizer::optimizeParameters(evaluations, startingSpeciesRates);
  return rates;
}

static Parameters findBestPointNelderMear(Parameters r1, Parameters r2,
                                          unsigned int iterations,
                                          FunctionToOptimize &function) {
  Parameters best = r1;
  best.setScore(-100000000000);
  for (unsigned int i = 0; i < iterations; ++i) {
    Parameters current =
        r1 + ((r2 - r1) * (double(i) / double(iterations - 1)));
    function.evaluate(current);
    if (current < best) {
      best = current;
    }
  }
  return best;
}

static Parameters
optimizeParametersNelderMear(FunctionToOptimize &function,
                             const Parameters &startingParameters) {
  std::vector<Parameters> rates;
  rates.push_back(startingParameters);
  auto N = startingParameters.dimensions();
  for (unsigned int r = 0; r < N; ++r) {
    auto p = startingParameters;
    p[r] -= 0.09;
    rates.push_back(p);
  }
  // n + 1 points for n dimensions
  assert(rates.size() == N + 1);
  for (auto &r : rates) {
    function.evaluate(r);
  }
  Parameters worstRate = startingParameters;
  unsigned int currentIt = 0;
  while (worstRate.distance(rates.back()) > 0.005) {
    std::sort(rates.begin(), rates.end());
    worstRate = rates.back();
    // centroid
    Parameters x0(worstRate.dimensions());
    for (unsigned int i = 0; i < rates.size() - 1; ++i) {
      x0 = x0 + rates[i];
    }
    x0 = x0 / double(rates.size() - 1);
    // reflexion, exansion and contraction at the same time
    Parameters x1 = x0 - (x0 - rates.back()) * 0.5;
    Parameters x2 = x0 + (x0 - rates.back()) * 1.5;
    unsigned int iterations = 8;
    Parameters xr = findBestPointNelderMear(x1, x2, iterations, function);
    if (xr < rates[rates.size() - 1]) {
      rates.back() = xr;
    }
    currentIt++;
  }
  std::sort(rates.begin(), rates.end());
  function.evaluate(rates[0]);
  Logger::timed << "Simplex converged after " << currentIt << " iterations"
                << std::endl;
  return rates[0];
}

#ifdef WITH_GSL

static double gslEval(const gsl_vector *v, void *params) {
  auto targetParam = (TargetParam *)(params);
  auto f = (FunctionToOptimize *)(targetParam->function);
  unsigned int n = v->size;
  Parameters values(n);
  for (unsigned int i = 0; i < n; ++i) {
    values[i] = gsl_vector_get(v, i);
  }
  auto res = f->evaluate(values);
  Logger::timed << values << std::endl;
  return -res;
}

static Parameters
optimizeParametersGSLSimplex(FunctionToOptimize &function,
                             const Parameters &startingParameters,
                             OptimizationSettings settings) {
  unsigned int n = startingParameters.dimensions();
  TargetParam targetFunction;
  targetFunction.function = &function;
  targetFunction.n = startingParameters.dimensions();

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  /* Starting point */
  auto x = gsl_vector_alloc(n);
  for (unsigned int i = 0; i < n; ++i) {
    gsl_vector_set(x, i, startingParameters[i]);
  }

  /* Set initial step sizes to 1 */
  auto ss = gsl_vector_alloc(n);
  gsl_vector_set_all(ss, 1.0);

  /* Initialize method and iterate */
  minex_func.n = n;
  minex_func.f = &gslEval;
  minex_func.params = &targetFunction;

  auto s = gsl_multimin_fminimizer_alloc(T, n);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    if (status) {
      break;
    }
    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, 1e-2);
  } while (status == GSL_CONTINUE && iter < 1000);

  Parameters res(n);
  for (unsigned int i = 0; i < n; ++i) {
    res[i] = gsl_vector_get(s->x, i);
  }
  res.setScore(s->fval);
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
  return res;
}

#endif

Parameters DTLOptimizer::optimizeParameters(FunctionToOptimize &function,
                                            Parameters startingParameters,
                                            OptimizationSettings settings) {
  auto res = startingParameters;
  switch (settings.strategy) {
  case RecOpt::Gradient:
    res = optimizeParametersGradient(function, startingParameters, settings);
    break;
  case RecOpt::LBFGSB:
    res = optimizeParametersLBFGSB(function, startingParameters, settings);
    break;
  case RecOpt::GSL_SIMPLEX:
#ifdef WITH_GSL
    res = optimizeParametersGSLSimplex(function, startingParameters, settings);
#else
    std::cerr
        << "Error, GSL routine not available, please install GSL and recompile"
        << std::endl;
#endif
    break;
  case RecOpt::Simplex:
    res = optimizeParametersNelderMear(function, startingParameters);
    break;
  default:
    assert(false);
    break;
  }
  if (settings.individualParamOpt && startingParameters.dimensions() > 1) {
    double llDiff = 0.0;
    unsigned int it = 0;
    do {
      auto ll = res.getScore();
      res = optimizeParametersIndividually(function, res, settings);
      llDiff = res.getScore() - ll;
      ll = res.getScore();
      ++it;
      if (settings.verbose) {
        Logger::timed << "llDiff after one round of individual opt: " << llDiff
                      << std::endl;
      }
    } while (llDiff > settings.individualParamOptMinImprovement &&
             it < settings.individualParamOptMaxIt);
  }
  return res;
}
