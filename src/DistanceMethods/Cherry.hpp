#pragma once

#include <memory>

#include <IO/Families.hpp>
#include <trees/PLLRootedTree.hpp>

class Cherry {
public:
  Cherry() = delete;

  /**
   *
   */
  static std::unique_ptr<PLLRootedTree>
  geneTreeCherry(const Families &families);
};
