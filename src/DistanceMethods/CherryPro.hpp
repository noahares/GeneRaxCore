#pragma once

#include <memory>

#include <IO/Families.hpp>
#include <trees/PLLRootedTree.hpp>

class CherryPro {
public:
  CherryPro() = delete;

  /**
   *
   */
  static std::unique_ptr<PLLRootedTree>
  geneTreeCherryPro(const Families &families);
};
