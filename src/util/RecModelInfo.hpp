#pragma once

#include <IO/ArgumentsHelper.hpp>
#include <maths/Parameters.hpp>

#include "enums.hpp"

struct RecModelInfo {
  // reconciliation model (UndatedDTL, UndatedDL, etc)
  RecModel model;
  // optimizer for rate optimization (Gradient, LBFGSB, etc)
  RecOpt recOpt;
  // if set to true, each family can have different set of rates
  bool perFamilyRates;
  // number of gamma categories for rate heterogeneity among families
  size_t gammaCategories;
  // at which ancestral species do we consider that originations
  // are possible and with which probability
  OriginationStrategy originationStrategy;
  // if set to true, for each family we prune from the species
  // tree the taxa not covered by this family
  bool pruneSpeciesTree;
  // if set to true, the reconciliation likelihood is calculated
  // as sum over the gene tree root ML position and its neighbouring
  // positions (specific to GeneRax)
  bool rootedGeneTree;
  // if set to true, only the gene root position from the starting
  // gene tree is considered (specific to GeneRax)
  bool forceGeneTreeRoot;
  // if set to true, different gene root positions are weighted
  // according to their MAD scores (specific to GeneRax)
  bool madRooting;
  // if the reconciliation model accounts for polytomies, branches
  // with lengths <= branchLengthThreshold will be contracted
  double branchLengthThreshold;
  // horizontal gene transfer constraint
  TransferConstaint transferConstraint;
  // disable duplications
  bool noDup;
  // disable DL inference
  bool noDL;
  // disable TL inference to speed up likelihood evaluation
  bool noTL;
  // path to a file which sets for each species the probability of
  // a gene copy to be not lost if not observed in the data
  std::string fractionMissingFile;
  // use less RAM, but likelihood evaluation might be slower
  // (specific to AleRax)
  bool memorySavings;

  /**
   *  Default constructor
   */
  RecModelInfo()
      : model(RecModel::UndatedDTL), recOpt(RecOpt::Gradient),
        perFamilyRates(true), gammaCategories(1),
        originationStrategy(OriginationStrategy::ROOT), pruneSpeciesTree(true),
        rootedGeneTree(true), forceGeneTreeRoot(false), madRooting(false),
        branchLengthThreshold(-1.0),
        transferConstraint(TransferConstaint::PARENTS), noDup(false),
        noDL(false), noTL(false), memorySavings(false) {}

  /**
   *  Constructor
   */
  RecModelInfo(RecModel model, RecOpt recOpt, bool perFamilyRates,
               unsigned int gammaCategories,
               OriginationStrategy originationStrategy, bool pruneSpeciesTree,
               bool rootedGeneTree, bool forceGeneTreeRoot, bool madRooting,
               double branchLengthThreshold,
               TransferConstaint transferConstraint, bool noDup, bool noDL,
               bool noTL, const std::string &fractionMissingFile,
               bool memorySavings)
      : model(model), recOpt(recOpt), perFamilyRates(perFamilyRates),
        gammaCategories(gammaCategories),
        originationStrategy(originationStrategy),
        pruneSpeciesTree(pruneSpeciesTree), rootedGeneTree(rootedGeneTree),
        forceGeneTreeRoot(forceGeneTreeRoot), madRooting(madRooting),
        branchLengthThreshold(branchLengthThreshold),
        transferConstraint(transferConstraint), noDup(noDup), noDL(noDL),
        noTL(noTL), fractionMissingFile(fractionMissingFile),
        memorySavings(memorySavings) {}

  void readFromArgv(char **argv, int &i) {
    model = RecModel(atoi(argv[i++]));
    recOpt = RecOpt(atoi(argv[i++]));
    perFamilyRates = bool(atoi(argv[i++]));
    gammaCategories = atoi(argv[i++]);
    originationStrategy = Enums::strToOrigination(argv[i++]);
    pruneSpeciesTree = bool(atoi(argv[i++]));
    rootedGeneTree = bool(atoi(argv[i++]));
    forceGeneTreeRoot = bool(atoi(argv[i++]));
    madRooting = bool(atoi(argv[i++]));
    transferConstraint = ArgumentsHelper::strToTransferConstraint(argv[i++]);
    noDup = bool(atoi(argv[i++]));
    noDL = bool(atoi(argv[i++]));
    noTL = bool(atoi(argv[i++]));
    branchLengthThreshold = double(atof(argv[i++]));
    fractionMissingFile = std::string(argv[i++]);
    if (fractionMissingFile == "NONE") {
      fractionMissingFile = std::string();
    }
    memorySavings = bool(atoi(argv[i++]));
  }

  std::vector<std::string> getArgv() const {
    std::vector<std::string> argv;
    argv.push_back(std::to_string(static_cast<int>(model)));
    argv.push_back(std::to_string(static_cast<int>(recOpt)));
    argv.push_back(std::to_string(static_cast<int>(perFamilyRates)));
    argv.push_back(std::to_string(static_cast<int>(gammaCategories)));
    argv.push_back(Enums::originationToStr(originationStrategy));
    argv.push_back(std::to_string(static_cast<int>(pruneSpeciesTree)));
    argv.push_back(std::to_string(static_cast<int>(rootedGeneTree)));
    argv.push_back(std::to_string(static_cast<int>(forceGeneTreeRoot)));
    argv.push_back(std::to_string(static_cast<int>(madRooting)));
    argv.push_back(
        ArgumentsHelper::transferConstraintToStr(transferConstraint));
    argv.push_back(std::to_string(static_cast<int>(noDup)));
    argv.push_back(std::to_string(static_cast<int>(noDL)));
    argv.push_back(std::to_string(static_cast<int>(noTL)));
    argv.push_back(std::to_string(branchLengthThreshold));
    if (fractionMissingFile.size()) {
      argv.push_back(fractionMissingFile);
    } else {
      argv.push_back(std::string("NONE"));
    }
    argv.push_back(std::to_string(static_cast<int>(memorySavings)));
    return argv;
  }

  static int getArgc() { return 15; }

  std::vector<char> getParamTypes() const {
    std::vector<char> res;
    for (auto str : Enums::parameterNames(model)) {
      assert(str.size() == 1);
      res.push_back(str[0]);
    }
    if (originationStrategy == OriginationStrategy::OPTIMIZE) {
      res.push_back('O');
    }
    return res;
  }

  unsigned int modelFreeParameters() const {
    return Enums::freeParameters(model) +
           (originationStrategy == OriginationStrategy::OPTIMIZE ? 1 : 0);
  }

  /**
   *  Return global parameters with the appropriate
   *  number of values (all set to 0.1)
   */
  Parameters getDefaultGlobalParameters() const {
    Parameters res(modelFreeParameters());
    for (unsigned int i = 0; i < res.dimensions(); ++i) {
      res[i] = 0.1;
    }
    if (noDup) {
      res[0] = 0.0;
    }
    return res;
  }

  /**
   *  Take user-defined parameters and return parameters
   *  with the appropriate dimensions. If the input parameters
   *  have too many values, the last ones are discarded, and if
   *  they do not have enough values, they are completed with 0.1
   */
  Parameters getParametersFromUser(const Parameters &user) const {
    Parameters res(modelFreeParameters());
    for (unsigned int i = 0; i < res.dimensions(); ++i) {
      if (user.dimensions() > i) {
        res[i] = user[i];
      } else {
        res[i] = 0.1;
      }
    }
    return res;
  }

  bool isDated() const {
    return transferConstraint == TransferConstaint::RELDATED;
  }
};
