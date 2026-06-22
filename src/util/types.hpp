#pragma once

#include <array>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

typedef struct corax_unode_s corax_unode_t;
typedef struct corax_utree_s corax_utree_t;
typedef struct corax_rnode_s corax_rnode_t;
typedef struct corax_rtree_s corax_rtree_t;

using VectorDouble = std::vector<double>;
using MatrixDouble = std::vector<VectorDouble>;
using DistanceMatrix = MatrixDouble;
using VectorUint = std::vector<unsigned int>;
using MatrixUint = std::vector<VectorUint>;
using VectorString = std::vector<std::string>;
using StringToUint = std::unordered_map<std::string, unsigned int>;
using StringToInt = std::unordered_map<std::string, int>;

using SPID = unsigned int; // species ID
using BID = unsigned int;  // branch ID
using CID = unsigned int;  // clade ID
using GID = unsigned int;  // gene ID
using TaxaSet = std::unordered_set<SPID>;
using MetaQuartet = std::array<TaxaSet, 4>;
using UInt3 = std::array<unsigned int, 3>;
using UInt4 = std::array<unsigned int, 4>;
using PairUInt = std::pair<unsigned int, unsigned int>;
using PerFamLL = std::vector<double>;
using DatedBackup = std::vector<unsigned int>;
using RatesVector = std::vector<std::vector<double>>;

struct TransferFrequencies {
  MatrixUint count;
  VectorString idToLabel;
};
