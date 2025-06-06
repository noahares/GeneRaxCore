#pragma once

class JointTree;
class PerCoreGeneTrees;
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <maths/ModelParameters.hpp>
#include <maths/Parameters.hpp>
#include <memory>
#include <string>
#include <util/enums.hpp>

class RootedTree;

class DTLOptimizerListener {
public:
  /**
   *  Callback called when we find a better set of parameters
   *  (e.g. to save checkpoint)
   */
  virtual void onBetterParametersFoundCallback() = 0;
};

enum class LBFGSBPrecision : int64_t {
  HIGH = 1,
  MEDIUM = (int64_t)1e7,
  LOW = (int64_t)1e12,
};

struct OptimizationSettings {
  OptimizationSettings()
      : strategy(RecOpt::Gradient), lineSearchMinImprovement(0.1),
        optimizationMinImprovement(3.0), minAlpha(0.0000001),
        startingAlpha(0.1), epsilon(0.0000001), verbose(false),
        individualParamOpt(false), individualParamOptMinImprovement(10.0),
        individualParamOptMaxIt(3), factr(LBFGSBPrecision::HIGH) {}

  RecOpt strategy;
  double lineSearchMinImprovement;
  double optimizationMinImprovement;
  double minAlpha;
  double startingAlpha;
  double epsilon;
  bool verbose;
  bool individualParamOpt;
  double individualParamOptMinImprovement;
  unsigned int individualParamOptMaxIt;
  std::vector<DTLOptimizerListener *> listeners;
  LBFGSBPrecision factr;

  void onBetterParametersFoundCallback() {
    for (auto listener : listeners) {
      listener->onBetterParametersFoundCallback();
    }
  }
};

class FunctionToOptimize {
public:
  virtual ~FunctionToOptimize() {};
  virtual double evaluate(Parameters &parameters) = 0;
};

class DTLOptimizer {
public:
  DTLOptimizer() = delete;

  /**
   *  Generic parallel method for parameter optimization. Finds the
   *  parameters (starting from startingParameters) that optimizes the
   *  function equal to the sum (over parallel ranks) of the evaluations
   *  functions.
   *  @param evaluations the subset of functions allocated to the
   *                     current core
   *  @param startingParameters starting parameters
   *  @return The parameters that maximize the function
   */
  static Parameters
  optimizeParameters(FunctionToOptimize &function,
                     Parameters startingParameters,
                     OptimizationSettings settings = OptimizationSettings());

  static Parameters
  optimizeParameters(PerCoreEvaluations &evaluations,
                     const Parameters &startingParameters,
                     OptimizationSettings settings = OptimizationSettings());

  /**
   *  Finds the global parameters that maximize evaluations. Global
   *  parameters means that they are not per-species nor per-families
   *  @param evaluations the subset of functions allocated to the
   *                     current core
   *  @param startingParameters if not set, several preselected starting
   *                            parameters will be tried
   *  @return The parameters that maximize the function
   */
  static Parameters optimizeParametersGlobalDTL(
      PerCoreEvaluations &evaluations,
      const Parameters *startingParameters = nullptr,
      OptimizationSettings settings = OptimizationSettings());

  /**
   * Same as optimizeParameters, but with a ModelParameters as input.
   */
  static ModelParameters optimizeModelParameters(
      PerCoreEvaluations &evaluations, bool optimizeFromStartingParameters,
      const ModelParameters &startingParameters,
      OptimizationSettings settings = OptimizationSettings());

  /**
   * Finds the per-species parameters that maximize  evaluations
   *  @param evaluations the subset of functions allocated to the
   *                     current core
   *  @param speciesNodesNumber number of species nodes
   *  @return The parameters that maximize the function
   */
  static Parameters
  optimizeParametersPerSpecies(PerCoreEvaluations &evaluations,
                               unsigned int speciesNodesNumber);
};
