#include "IniParser.h"
#include <fstream>
#include <iostream>
#include <sstream>

IniParser& IniParser::getInstance() {
    static IniParser instance; // Static instance ensures single creation
    instance.load("optimizer_params.ini");
    return instance;
}

bool IniParser::load(const std::string& filename) {
  if (is_loaded) {
    return true;
  }
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Failed to open INI file: " << filename << std::endl;
    return false;
  }

  std::string line;
  std::string currentSection;
  while (std::getline(file, line)) {
    if (line.empty() || line[0] == ';' || line[0] == '#') {
      // Skip empty lines and comments
      continue;
    }
    if (line[0] == '[' && line[line.length() - 1] == ']') {
      // Section line
      currentSection = line.substr(1, line.length() - 2);
    } else {
      // Key-value pair
      std::istringstream iss(line);
      std::string key, value;
      if (std::getline(iss, key, '=')) {
        if (std::getline(iss, value)) {
          // Remove leading/trailing whitespace
          key.erase(0, key.find_first_not_of(" \t"));
          key.erase(key.find_last_not_of(" \t") + 1);
          value.erase(0, value.find_first_not_of(" \t"));
          value.erase(value.find_last_not_of(" \t") + 1);
          // Store key-value pair
          float float_value = std::stof(value);
          if (!currentSection.empty()) {
            data[currentSection + "." + key] = float_value;
          } else {
            data[key] = float_value;
          }
        }
      }
    }
  }
  file.close();
  is_loaded = true;
  return true;
}

// TODO: make this generic over the value type
float IniParser::getValue(const std::string& key, const double default_value) const {
  auto it = data.find(key);
  if (it != data.end()) {
    return it->second;
  }
  return default_value;
}
