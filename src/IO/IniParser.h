#ifndef INI_PARSER_H
#define INI_PARSER_H

#include <map>
#include <string>

class IniParser {
public:
  static IniParser &getInstance(); // Singleton instance getter
  bool load(const std::string &filename);
  float getValue(const std::string &key) const;

private:
  IniParser() : is_loaded(false) {}                         // Private constructor for singleton
  IniParser(const IniParser &) = delete; // Delete copy constructor
  IniParser &
  operator=(const IniParser &) = delete; // Delete assignment operator
  std::map<std::string, float> data;
  bool is_loaded;
};

#endif // INI_PARSER_H
