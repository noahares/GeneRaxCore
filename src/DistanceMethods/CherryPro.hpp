#pragma once

#include <IO/Families.hpp>
#include <memory>
#include <trees/PLLRootedTree.hpp>
#include <vector>

class CherryPro {
public:
  CherryPro() = delete;

  /**
   *
   */
  static std::unique_ptr<PLLRootedTree>
  geneTreeCherryPro(const Families &families);
};
