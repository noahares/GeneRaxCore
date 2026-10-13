#pragma once

#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include <parallelization/ParallelContext.hpp>

class ParallelOfstream {
public:
  ParallelOfstream(const std::string &fileName, bool masterRankOnly = true)
      : _os(nullptr) {
    if (!ParallelContext::getRank() || !masterRankOnly) {
      _os = std::make_unique<std::ofstream>(fileName);
    } else {
      _os = std::make_unique<std::ostream>(nullptr);
    }
  }

  ~ParallelOfstream() { close(); }

  // forbid copy, allow move
  ParallelOfstream(const ParallelOfstream &) = delete;
  ParallelOfstream &operator=(const ParallelOfstream &) = delete;
  ParallelOfstream(ParallelOfstream &&) = default;
  ParallelOfstream &operator=(ParallelOfstream &&) = default;

  void close() {
    // force the _os destruction if it holds an ofstream
    _os = std::make_unique<std::ostream>(nullptr);
  }

  template <typename T> std::ostream &operator<<(const T &input) {
    *_os << input;
    return *_os;
  }

  std::ostream &operator<<(std::ostream &(*manip)(std::ostream &)) {
    manip(*_os);
    return *_os;
  }

private:
  std::unique_ptr<std::ostream> _os;
};
