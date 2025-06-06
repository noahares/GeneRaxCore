#include "SpeciesTreeOptimizer.hpp"

#include "DTLOptimizer.hpp"
#include <IO/FileSystem.hpp>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <routines/Routines.hpp>
#include <search/SpeciesSPRSearch.hpp>
#include <search/SpeciesTransferSearch.hpp>
#include <support/ICCalculator.hpp>
#include <util/Paths.hpp>

static std::unique_ptr<SpeciesTree>
makeSpeciesTree(const std::string &speciesTreeFile,
                const Families &initialFamilies) {
  if (speciesTreeFile == "random") {
    return std::make_unique<SpeciesTree>(initialFamilies);
  } else {
    return std::make_unique<SpeciesTree>(speciesTreeFile);
  }
}

SpeciesTreeOptimizer::SpeciesTreeOptimizer(
    const std::string speciesTreeFile, const Families &initialFamilies,
    const RecModelInfo &recModelInfo, const Parameters &startingRates,
    bool userDTLRates, const std::string &outputDir,
    const SpeciesTreeSearchParams &searchParams)
    : _speciesTree(makeSpeciesTree(speciesTreeFile, initialFamilies)),
      _geneTrees(std::make_unique<PerCoreGeneTrees>(initialFamilies, true)),
      _initialFamilies(initialFamilies), _outputDir(outputDir),
      _firstOptimizeRatesCall(true), _userDTLRates(userDTLRates),
      _modelRates(startingRates, 1, recModelInfo), _searchParams(searchParams),
      _okForClades(0), _koForClades(0),
      _searchState(
          *_speciesTree,
          Paths::getSpeciesTreeFile(_outputDir, "inferred_species_tree.newick"),
          _geneTrees->getTrees().size()) {

  _modelRates.info.perFamilyRates = false; // we set it back a few
                                           // lines later
  updateEvaluations();
  _modelRates = ModelParameters(startingRates, _geneTrees->getTrees().size(),
                                recModelInfo);
  _speciesTree->addListener(this);
  saveCurrentSpeciesTreeId();
  _computeAllGeneClades();
  _searchState.farFromPlausible &=
      _unsupportedCladesNumber() >
      std::max<unsigned int>(_speciesTree->getTree().getNodeNumber() / 4, 1);
}

static bool testAndSwap(size_t &hash1, size_t &hash2) {
  std::swap(hash1, hash2);
  return hash1 != hash2;
}

void SpeciesTreeOptimizer::optimize(SpeciesSearchStrategy strategy) {
  computeRecLikelihood();
  size_t hash1 = 0;
  size_t hash2 = 0;
  unsigned int index = 0;
  switch (strategy) {
  case SpeciesSearchStrategy::SPR:
    for (unsigned int radius = 1; radius <= _searchParams.sprRadius; ++radius) {
      _searchState.bestLL = _evaluator.optimizeModelRates();
      sprSearch(radius);
    }
    break;
  case SpeciesSearchStrategy::TRANSFERS:
    transferSearch();
    rootSearch(_searchParams.rootBigRadius, false);
    transferSearch();
    rootSearch(_searchParams.rootBigRadius, true);
    break;
  case SpeciesSearchStrategy::HYBRID:
    /**
     *  Alternate transfer search and normal
     *  SPR search, until one does not find
     *  a better tree. Run each at least once.
     */
    if (!_searchState.farFromPlausible) {
      _searchState.bestLL = _evaluator.optimizeModelRates();
      computeRecLikelihood();
      rootSearch(_searchParams.rootSmallRadius, false);
    }
    do {
      if (index++ % 2 == 0) {
        transferSearch();
      } else {
        sprSearch(_searchParams.sprRadius);
      }
      if (!_searchState.farFromPlausible) {
        rootSearch(_searchParams.rootSmallRadius, false);
      }
      hash1 = _speciesTree->getHash();
    } while (testAndSwap(hash1, hash2));
    _searchState.bestLL = _evaluator.optimizeModelRates(true);
    rootSearch(_searchParams.rootBigRadius, true);
    break;
  case SpeciesSearchStrategy::REROOT:
    rootSearch(_searchParams.rootBigRadius, true);
    break;
  case SpeciesSearchStrategy::EVAL:
    _searchState.bestLL = _evaluator.optimizeModelRates(true);
    Logger::info << "Reconciliation likelihood: " << computeRecLikelihood()
                 << std::endl;
    break;
  case SpeciesSearchStrategy::SKIP:
    assert(false);
  }
}

SpeciesTreeOptimizer::~SpeciesTreeOptimizer() {
  _speciesTree->removeListener(this);
}

double SpeciesTreeOptimizer::rootSearch(unsigned int maxDepth,
                                        bool outputConsel) {
  TreePerFamLLVec treePerFamLLVec;
  RootLikelihoods rootLikelihoods(_evaluations.size());
  Logger::info << std::endl;
  SpeciesRootSearch::rootSearch(*_speciesTree, _evaluator, _searchState,
                                maxDepth, &rootLikelihoods,
                                (outputConsel ? &treePerFamLLVec : nullptr));
  saveCurrentSpeciesTreeId();
  {
    auto newick = _speciesTree->getTree().getNewickString();
    PLLRootedTree tree(newick, false);
    rootLikelihoods.fillTree(tree);
    auto out = Paths::getSpeciesTreeFile(_outputDir, "species_tree_llr.newick");
    tree.save(out);
  }
  {
    auto newick = _speciesTree->getTree().getNewickString();
    PLLRootedTree tree(newick, false);
    rootLikelihoods.fillTreeBootstraps(tree);
    auto out = Paths::getSpeciesTreeFile(_outputDir,
                                         "species_tree_root_support.newick");
    tree.save(out);
  }
  if (outputConsel) {
    std::string treesOutput = Paths::getConselTreeList(_outputDir, "roots");
    std::string llOutput = Paths::getConselLikelihoods(_outputDir, "roots");
    savePerFamilyLikelihoods(treePerFamLLVec, treesOutput, llOutput);
  }
  auto ll = computeRecLikelihood();
  return ll;
}

std::vector<double> SpeciesTreeOptimizer::_getSupport() {
  std::string temp = FileSystem::joinPaths(_outputDir, "tmp");
  std::vector<double> idToSupport;
  bool paralogyAware = true;
  unsigned int eqpicRadius = 3;
  ICCalculator::computeScores(_speciesTree->getTree(), _initialFamilies,
                              paralogyAware, eqpicRadius, temp, idToSupport);
  for (auto node : _speciesTree->getTree().getLeaves()) {
    idToSupport[node->node_index] = 1.0;
  }
  return idToSupport;
}

double SpeciesTreeOptimizer::transferSearch() {
  Logger::info << std::endl;
  SpeciesTransferSearch::transferSearch(*_speciesTree, _evaluator,
                                        _searchState);
  return _searchState.bestLL;
}

double SpeciesTreeOptimizer::sprSearch(unsigned int radius) {
  Logger::info << std::endl;
  SpeciesSPRSearch::SPRSearch(*_speciesTree, _evaluator, _searchState, radius);
  return _searchState.bestLL;
}

std::string
SpeciesTreeOptimizer::saveCurrentSpeciesTreeId(std::string name,
                                               bool masterRankOnly) {
  std::string res = Paths::getSpeciesTreeFile(_outputDir, name);
  saveCurrentSpeciesTreePath(res, masterRankOnly);
  return res;
}

void SpeciesTreeOptimizer::saveCurrentSpeciesTreePath(const std::string &str,
                                                      bool masterRankOnly) {
  _speciesTree->saveToFile(str, masterRankOnly);
  ParallelContext::barrier();
}

double SpeciesTreeOptimizer::computeRecLikelihood() {
  return _evaluator.computeLikelihood();
}

void SpeciesTreeOptimizer::updateEvaluations() {
  assert(_geneTrees);
  auto &trees = _geneTrees->getTrees();
  _evaluations.resize(trees.size());
  for (unsigned int i = 0; i < trees.size(); ++i) {
    auto &tree = trees[i];
    std::string enforcedRootedGeneTree;
    if (_modelRates.info.forceGeneTreeRoot) {
      enforcedRootedGeneTree = tree.startingGeneTreeFile;
    }
    _evaluations[i] = std::make_shared<ReconciliationEvaluation>(
        _speciesTree->getTree(), *tree.geneTree, tree.mapping, _modelRates.info,
        enforcedRootedGeneTree);
    _evaluations[i]->setRates(_modelRates.getRates(i));
    _evaluations[i]->setPartialLikelihoodMode(
        PartialLikelihoodMode::PartialSpecies);
  }
  _previousGeneRoots.resize(_evaluations.size());
  std::fill(_previousGeneRoots.begin(), _previousGeneRoots.end(), nullptr);
  _evaluator.init(_evaluations, *_geneTrees, _modelRates,
                  _modelRates.info.rootedGeneTree,
                  _modelRates.info.pruneSpeciesTree, _userDTLRates);
}

void SpeciesTreeOptimizer::onSpeciesTreeChange(
    const std::unordered_set<corax_rnode_t *> *nodesToInvalidate) {
  for (auto &evaluation : _evaluations) {
    evaluation->onSpeciesTreeChange(nodesToInvalidate);
  }
}

std::string getCladesSetPath(const std::string &outputDir, int rank) {
  std::string basePath = "clades_" + std::to_string(rank) + ".txt";
  return FileSystem::joinPaths(outputDir, basePath);
}

void SpeciesTreeOptimizer::_computeAllGeneClades() {
  ParallelContext::barrier();
  // Compute local clades
  auto speciesLabelToInt = _speciesTree->getTree().getLabelToIntMap();
  CladeSet allClades;
  for (auto &tree : _geneTrees->getTrees()) {
    auto cladesSet =
        Clade::buildCladeSet(*tree.geneTree, tree.mapping, speciesLabelToInt);
    allClades.insert(cladesSet.begin(), cladesSet.end());
  }
  // Write local clades
  std::ofstream os(getCladesSetPath(_outputDir, ParallelContext::getRank()));
  for (auto clade : allClades) {
    os << clade << " ";
  }
  os.close();
  ParallelContext::barrier();
  // Load all clades
  _geneClades.clear();
  for (unsigned int rank = 0; rank < ParallelContext::getSize(); ++rank) {
    std::ifstream is(getCladesSetPath(_outputDir, rank));
    unsigned int clade = 0;
    while (is >> clade) {
      _geneClades.insert(clade);
    }
  }
  assert(ParallelContext::isIntEqual(_geneClades.size()));
  ParallelContext::barrier();
  std::remove(getCladesSetPath(_outputDir, ParallelContext::getRank()).c_str());
}

unsigned int SpeciesTreeOptimizer::_unsupportedCladesNumber() {
  auto speciesClades = Clade::buildCladeSet(_speciesTree->getTree());
  unsigned int intersectionSize = intersection_size(speciesClades, _geneClades);
  return speciesClades.size() - intersectionSize;
}

static std::string getSubtreeID(corax_rnode_t *subtree) {
  if (!subtree->left) {
    return std::string(subtree->label);
  }
  std::string res("(");
  std::string id1 = getSubtreeID(subtree->left);
  std::string id2 = getSubtreeID(subtree->right);
  if (id1 > id2) {
    std::swap(id1, id2);
  }
  return std::string("(") + id1 + "," + id2 + ")";
}

void SpeciesTreeOptimizer::savePerFamilyLikelihoods(
    const TreePerFamLLVec &treePerFamLLVec, const std::string &treesOutput,
    const std::string &llOutput) {
  ParallelContext::barrier();
  if (ParallelContext::getRank() == 0 && treePerFamLLVec.size() != 0) {
    std::ofstream osLL(llOutput);
    std::ofstream osTrees(treesOutput);
    auto treesNumber = treePerFamLLVec.size();
    auto familiesNumber = treePerFamLLVec[0].second.size();
    osLL << treesNumber << " " << familiesNumber << std::endl;
    unsigned int index = 1;
    for (const auto &treePerFamLL : treePerFamLLVec) {
      const auto &tree = treePerFamLL.first;
      const auto perFamLL = treePerFamLL.second;
      std::string treeName = "tree";
      treeName += std::to_string(index++);
      osTrees << tree << std::endl;
      osLL << treeName;
      for (auto ll : perFamLL) {
        osLL << " " << ll;
      }
      osLL << std::endl;
    }
  }
  ParallelContext::barrier();
}

double SpeciesTreeLikelihoodEvaluator::computeLikelihood(PerFamLL *perFamLL) {
  if (_rootedGeneTrees) {
    for (auto evaluation : *_evaluations) {
      evaluation->setRoot(nullptr);
    }
  }
  if (perFamLL) {
    perFamLL->clear();
  }
  double sumLL = 0.0;
  for (auto &evaluation : *_evaluations) {
    auto ll = evaluation->evaluate();
    if (perFamLL) {
      perFamLL->push_back(ll);
    }
    sumLL += ll;
  }
  ParallelContext::sumDouble(sumLL);
  return sumLL;
}

double SpeciesTreeLikelihoodEvaluator::computeLikelihoodFast() {
  double sumLL = 0.0;
  for (auto &evaluation : *_evaluations) {
    auto ll = evaluation->evaluate();
    sumLL += ll;
  }
  ParallelContext::sumDouble(sumLL);
  return sumLL;
}

bool SpeciesTreeLikelihoodEvaluator::providesFastLikelihoodImpl() const {
  return _rootedGeneTrees;
}

double SpeciesTreeLikelihoodEvaluator::optimizeModelRates(bool thorough) {
  if (_userDTLRates) {
    return computeLikelihood();
  }
  auto rates = *_modelRates;
  OptimizationSettings settings;
  double ll = computeLikelihood();
  if (!thorough) {
    settings.lineSearchMinImprovement = 10.0;
    settings.minAlpha = 0.01;
    settings.optimizationMinImprovement = std::max(3.0, ll / 1000.0);
  }
  bool _firstOptimizeRatesCall = false;
  *_modelRates = DTLOptimizer::optimizeModelParameters(
      *_evaluations, !_firstOptimizeRatesCall, *_modelRates, settings);
  _firstOptimizeRatesCall = false;
  unsigned int i = 0;
  for (auto &evaluation : *_evaluations) {
    evaluation->setRates(_modelRates->getRates(i++));
  }
  if (!_modelRates->info.perFamilyRates) {
    Logger::timed << "[Species search] Best rates: " << _modelRates->rates
                  << std::endl;
  }
  return computeLikelihood();
}

void SpeciesTreeLikelihoodEvaluator::getTransferInformation(
    SpeciesTree &speciesTree, TransferFrequencies &frequencies,
    PerSpeciesEvents &perSpeciesEvents,
    PerCorePotentialTransfers &potentialTransfers) {
  ParallelContext::barrier();
  unsigned int reconciliationSamples = 0; // use ML reconciliation
  Routines::getTransfersFrequencies(speciesTree.getTree(), *_geneTrees,
                                    *_modelRates, reconciliationSamples,
                                    frequencies, potentialTransfers);
  const bool forceTransfers = true;
  Routines::getPerSpeciesEvents(speciesTree.getTree(), *_geneTrees,
                                *_modelRates, reconciliationSamples,
                                perSpeciesEvents, forceTransfers);
}

void SpeciesTreeLikelihoodEvaluator::pushRollback() {
  if (_rootedGeneTrees) {
    _previousGeneRoots.push(std::vector<corax_unode_t *>());
    for (auto evaluation : *_evaluations) {
      _previousGeneRoots.top().push_back(evaluation->getRoot());
    }
  }
}

void SpeciesTreeLikelihoodEvaluator::popAndApplyRollback() {
  if (_rootedGeneTrees) {
    for (unsigned int i = 0; i < _evaluations->size(); ++i) {
      (*_evaluations)[i]->setRoot(_previousGeneRoots.top()[i]);
    }
    _previousGeneRoots.pop();
  }
}
