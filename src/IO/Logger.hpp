#pragma once

#include <chrono>
#include <fstream>
#include <iostream>
#include "parallelization/ParallelContext.hpp"

using TimePoint = std::chrono::high_resolution_clock::time_point;



/**
 *  All logs are printed to std::cout.
 *  In parallel runs, only the master can log.
 *  init() must be called before using the Logger.
 *  If initFileOutput is called, logs are also printed
 *  to the given file
 *
 *  Logger::info is used for normal logs
 *  Logger::timed also prints the elapsed time since the start of the program
 *  Logger::error is used for errors and prefixes messages with [Error]
 */
class Logger: public std::ofstream
{
private:
  enum LoggerType {
    lt_info, lt_error, lt_timed, lt_perrank
  };
  LoggerType _type;
  std::ostream *_os; // I do not own this one
  bool _silent;

  Logger();
  void setType(LoggerType type) {_type = type;}
  void setStream(std::ostream &os) {_os = &os;}

public:
  static void init();
  static void close();

  static void initFileOutput(const std::string &output);

  static void mute() {info._silent = true;}
  static void unmute() {info._silent = false;}

  // if true, the MPI rank can't write logs
  bool isSilent() {
    if (!_silent) {
      return false;
    }
    if (_type == lt_perrank) {
      return false;
    }
    return (_type == lt_timed || _type == lt_info) && ParallelContext::getRank();
  }

  static void enableLogFile(bool enable) {
    logFile = enable ? saveLogFile : 0;
  }

  static long getElapsedSec() {
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    return std::chrono::duration_cast<std::chrono::seconds>(elapsed).count();
  }

  template <typename T>
  Logger& operator << (T&& message)
  {
    if (isSilent()) {
      return *this;
    }
    if (_type == lt_timed) {
      auto seconds = getElapsedSec();
      auto hours  = seconds / 3600;
      auto minutes = (seconds % 3600) / 60;
      seconds = seconds % 60;
      char time[30];
      sprintf(time, "%02ld:%02ld:%02ld", hours, minutes, seconds);
      *_os << "[" << time << "] " << message;
      if (logFile) {
        *logFile << "[" << time << "] " << message;
      }
      return Logger::info;
    } else if (_type == lt_error) {
      *_os << "[Error] " << message;
      if (logFile) {
        *logFile << "[Error] " << message;
      }
      // Let a Logger::info continue the error message
      // even for a non-master rank
      unmute();
      return Logger::info;
    } else if (_type == lt_perrank) {
      initRankFileOutput();
      *rankLogFile << message;
      return *this;
    } else {
      *_os << message;
      _os->flush();
      if (logFile) {
        *logFile << message;
      }
      return *this;
    }
  }

  Logger& operator << (std::ostream& (*manip)(std::ostream&))
  {
    if (isSilent()) {
      return *this;
    }
    if (_type != lt_perrank) {
      manip(*_os);
      if (logFile) {
        manip(*logFile);
      }
      // After printing std::endl ensure that for
      // all non-master ranks Logger::info is silent
      mute();
    } else {
      manip(*rankLogFile);
    }
    return *this;
  }

  static void initRankFileOutput();

  static Logger info;
  static Logger error;
  static Logger timed;
  static Logger perrank;
  static TimePoint start;
  static std::string outputdir;
  static std::ofstream *logFile;
  static std::ofstream *rankLogFile;
  static std::ofstream *saveLogFile;
  static bool inited;

};


