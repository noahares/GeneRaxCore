#pragma once

#include <fstream>
#include <iterator>
#include <string>
#include <parallelization/ParallelContext.hpp>
#include <sys/stat.h>
#include <sys/types.h>



class FileSystem {
public:
  FileSystem() = delete;

  static std::string joinPaths(const std::string &p1, const std::string &p2)
  {
    std::string sep;
#if defined(_WIN32)
    sep = "\\";
#else
    sep = "/";
#endif
    return p1 + sep + p2;
  }

  static void mkdir(const std::string &dirPath, bool masterRankOnly)
  {
    if (masterRankOnly && ParallelContext::getRank() != 0) {
      return;
    }
#if defined(_WIN32)
    _mkdir(dirPath.c_str()); // can be used on Windows
#else
    mode_t nMode = 0733; // UNIX style permissions
    ::mkdir(dirPath.c_str(), nMode); // can be used on non-Windows
#endif
  }

  static bool dirExists(const std::string &dirPath)
  {
    struct stat info;
    if(stat(dirPath.c_str(), &info) != 0) {
      return false;
    } else if (info.st_mode & S_IFDIR) {
      return true;
    } else {
      return false;
    }
  }

  static bool exists(const std::string &filePath)
  {
    std::ifstream f(filePath);
    return f.good();
  }

  static void getFileContent(const std::string &filePath, std::string &content)
  {
    std::ifstream ifs(filePath);
    content.assign((std::istreambuf_iterator<char>(ifs)),
        (std::istreambuf_iterator<char>()) );
  }

  static void replaceWithContentIfFile(std::string &str)
  {
    std::ifstream ifs(str);
    if (ifs.good()) {
      str.assign((std::istreambuf_iterator<char>(ifs)),
          (std::istreambuf_iterator<char>()) );
    }
  }

  static void copy(const std::string &f1, const std::string &f2, bool masterRankOnly)
  {
    if (masterRankOnly && ParallelContext::getRank() != 0) {
      return;
    }
    std::ifstream src(f1, std::ios::binary);
    std::ofstream dst(f2, std::ios::binary);
    dst << src.rdbuf();
  }

};


