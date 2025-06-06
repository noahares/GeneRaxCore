#pragma once

/*
 * This file was copied from RAxML-NG
 */

#include <array>
#include <map>
#include <random>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <corax/core/partition.h>

enum class StartingTree { random, parsimony, user };

enum class Command {
  none = 0,
  help,
  version,
  evaluate,
  search,
  bootstrap,
  all,
  support,
  bsconverge,
  terrace,
  check,
  parse,
  start
};

enum class FileFormat {
  autodetect = 0,
  fasta,
  phylip,
  iphylip,
  vcf,
  catg,
  binary
};

enum class DataType {
  autodetect = 0,
  dna,
  protein,
  binary,
  multistate,
  diploid10
};

enum class ParamValue {
  undefined = 0,
  equal = 1,
  user = 2,
  model = 3,
  empirical = 4,
  ML = 5
};

enum class AscBiasCorrection {
  none = 0,
  lewis = CORAX_ATTRIB_AB_LEWIS,
  felsenstein = CORAX_ATTRIB_AB_FELSENSTEIN,
  stamatakis = CORAX_ATTRIB_AB_STAMATAKIS,
};

enum class BootstopCriterion { none = 0, autoMRE, autoMR, autoFC };

enum class LoadBalancing { naive = 0, kassian, benoit };

const std::string ParamValueNames[] = {"undefined", "equal",     "user",
                                       "model",     "empirical", "ML"};

typedef std::vector<double> doubleVector;
typedef std::vector<int> intVector;
typedef std::vector<unsigned int> uintVector;
typedef std::vector<std::string> NameList;
typedef std::pair<size_t, std::string> IdNamePair;
typedef std::vector<IdNamePair> IdNameVector;
typedef std::unordered_map<size_t, std::string> IdNameMap;
typedef std::unordered_map<std::string, size_t> NameIdMap;
typedef std::unordered_map<std::string, std::string> NameMap;
typedef std::set<size_t> IDSet;
typedef std::vector<size_t> IDVector;
typedef std::map<StartingTree, size_t> StartingTreeMap;

typedef unsigned int WeightType;
typedef std::vector<WeightType> WeightVector;
typedef std::vector<WeightVector> WeightVectorList;
typedef std::unordered_map<size_t, WeightVector> WeightVectorMap;

typedef std::default_random_engine RandomGenerator;

/*
 * workaround needed for using enum as std::map key
 * code from: http://stackoverflow.com/a/24847480
 * */
struct EnumClassHash {
  template <typename T> std::size_t operator()(T t) const {
    return static_cast<std::size_t>(t);
  }
};

/* generic std::exception class */
class RaxmlException : public std::exception {
public:
  RaxmlException(const std::string &message) : _message(message) {}

  virtual const char *what() const noexcept { return message().c_str(); }

  virtual const std::string &message() const {
    update_message();
    return _message;
  }

  virtual ~RaxmlException() {}

protected:
  mutable std::string _message;

  virtual void update_message() const {}

  template <typename... Args>
  std::string format_message(const std::string &fmt, Args... args) const {
    size_t size = std::snprintf(nullptr, 0, fmt.c_str(), args...) + 1;
    std::string msg;
    msg.resize(size);
    std::snprintf(&msg[0], size, fmt.c_str(), args...);
    return msg;
  }
};
