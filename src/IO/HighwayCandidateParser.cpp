#include "HighwayCandidateParser.hpp"
#include <IO/IO.hpp>
#include <sstream>
#include <string>

bool isParent(corax_rnode_t *parent, corax_rnode_t *child) {
  while (child) {
    if (parent == child) {
      return true;
    }
    child = child->parent;
  }
  return false;
}

bool readTaxa(
    std::stringstream &iss, std::vector<Highway> &highways,
    const std::unordered_map<std::string, corax_rnode_t *> &labelToNode) {
  std::string source, target, rate_s;
  std::vector<corax_rnode_t *> source_nodes;
  std::vector<corax_rnode_t *> target_nodes;

  if (!std::getline(iss, source, ',') || !std::getline(iss, target, ',')) {
    return false;
  }
  auto parse_node = [&](std::string &str, std::vector<corax_rnode_t *> &nodes) {
    IO::removeSpaces(str);
    if (str == "*") {
      for (auto it : labelToNode) {
        nodes.push_back(it.second);
      }
    } else {
      auto node = labelToNode.find(str);
      if (node == labelToNode.end()) {
        return false;
      }
      nodes.push_back(node->second);
    }
    return true;
  };
  if (!parse_node(source, source_nodes) || !parse_node(target, target_nodes))
    return false;
  for (auto from : source_nodes) {
    for (auto to : target_nodes) {
      if (!isParent(to, from)) {
        highways.push_back(Highway(from, to));
      }
    }
  }
  if (std::getline(iss, rate_s, '\n')) {
    try {
      double rate = std::stod(rate_s);
      for (auto &hw : highways) {
        hw.proba = rate;
      }
    } catch (std::exception &err) {
      return false;
    }
  }
  return true;
}

std::vector<Highway>
HighwayCandidateParser::parse(const std::string &candidateFile,
                              PLLRootedTree &speciesTree) {
  std::vector<Highway> candidates;
  auto labelToNode = speciesTree.getLabelToNode(false);
  std::ifstream is(candidateFile);
  std::string line;
  size_t lineNumber = 0;
  while (std::getline(is, line)) {
    lineNumber++;
    if (IO::isBlanck(line)) {
      continue;
    }
    IO::removeSpaces(line);
    if (line[0] == '#') {
      continue;
    }
    std::stringstream iss(line);
    if (!readTaxa(iss, candidates, labelToNode)) {
      Logger::info << "Failed to parse line " << lineNumber
                   << " of highways file " << candidateFile
                   << "! Continuing with next line" << std::endl;
    }
  }
  return candidates;
}
