#pragma once

#include <IO/Families.hpp>
#include <memory>
#include <trees/PLLRootedTree.hpp>
#include <vector>

class Cherry {
public:
  Cherry() = delete;

  /**
   *
   */
  static std::unique_ptr<PLLRootedTree>
  geneTreeCherry(const Families &families);
};
