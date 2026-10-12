#pragma once

#include <IO/Families.hpp>
#include <IO/GeneSpeciesMapping.hpp>

class PLLUnrootedTree;

/**
 * Holds the gene trees (with there mappings to the species tree)
 * allocated to the current parallel core.
 *
 * The sets of PerCoreGeneTrees over all parallel cores are a
 * partition of all gene trees.
 */
class PerCoreGeneTrees {
public:
  struct GeneTree {
    std::string name;
    std::string startingGeneTreeFile;
    unsigned int familyIndex;
    GeneSpeciesMapping mapping;
    PLLUnrootedTree *geneTree;
    bool ownTree; // If true, I am responsible for destroying the tree
    ~GeneTree();
  };

  /**
   *  Parse and allocate to the current core the gene trees
   *  from the description of the gene families.
   *  @param families families description
   */
  PerCoreGeneTrees(const Families &families, bool acceptMultipleTrees = false,
                   bool ccpMode = false);
  /**
   * Create an instance with a unique gene tree, without
   * accouting for parallelization.
   *
   * This is useful for applying some methods, whose interface
   * require a PerCoreGeneTrees object, to a unique gene tree
   * without parallelization.
   *  @param mapping The gene to species mapping
   *  @param geneTree The gene tree
   */
  PerCoreGeneTrees(const GeneSpeciesMapping &mapping,
                   PLLUnrootedTree &geneTree);

  // forbid copy and move
  PerCoreGeneTrees(const PerCoreGeneTrees &) = delete;
  PerCoreGeneTrees &operator=(const PerCoreGeneTrees &) = delete;
  PerCoreGeneTrees(PerCoreGeneTrees &&) = delete;
  PerCoreGeneTrees &operator=(PerCoreGeneTrees &&) = delete;
  ~PerCoreGeneTrees() = default;

  /**
   *  @return Trees allocated to the current core
   */
  std::vector<GeneTree> &getTrees() { return _geneTrees; }
  const std::vector<GeneTree> &getTrees() const { return _geneTrees; }

  /**
   *  @param speciesTreeFile path to the species tree file
   *  @return true if the mappings are valid.
   */
  bool checkMappings(const std::string &speciesTreeFile);

  /**
   *  @param families All families
   *  @param[out] perCoreFamilies Families allocated to the current core
   */
  static void getPerCoreFamilies(const Families &families,
                                 Families &perCoreFamilies);

  /**
   *  @param families All families
   *  @param ccpMode To count clades in CCP files or leaves in gene trees
   *  @return The number of observed clades or gene leaves per family
   */
  static std::vector<unsigned int>
  parallelGetTreeSizes(const Families &families, bool ccpMode = false);

private:
  std::vector<GeneTree> _geneTrees;
  std::vector<unsigned int> _treeSizes;
};
