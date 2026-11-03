#pragma once

#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

using TimePoint = std::chrono::high_resolution_clock::time_point;

/**
 *  All logs are printed to std::cout.
 *  All parallel ranks can write error logs.
 *  Only the master rank can write info and timed logs.
 *  init() must be called before using Logger.
 *  If initFileOutput(file) is called, logs are also printed to the file.
 *  close() must be called before reiniting or destructing Logger.
 *
 *  Logger::info is used for normal logs
 *  Logger::timed also prints the time elapsed since the start of the program
 *  Logger::error is used for errors and prefixes messages with [Error]
 */
class Logger {
private:
  enum LoggerType { lt_info, lt_error, lt_timed, lt_perrank };
  LoggerType _type;
  std::ostream *_os; // not owned by this class
  bool _muted; // if true, non-master ranks can't write info and timed logs

  Logger() : _type(lt_perrank), _os(nullptr), _muted(true) {}
  ~Logger() = default; // forbid explicit destruction (called only on exit)

  void setType(LoggerType type) { _type = type; }
  void setStream(std::ostream *os) { _os = os; }
  bool isSilent() const; // if true, the MPI rank can't write logs

  // disable/enable non-master ranks writing info logs
  // NOTE: these two functions as well as the _muted member are only for
  // internal use by the stream operators and should not be touched elsewhere
  static void mute() { info._muted = true; }
  static void unmute() { info._muted = false; }

public:
  // forbid copy and move
  Logger(const Logger &) = delete;
  Logger &operator=(const Logger &) = delete;
  Logger(Logger &&) = delete;
  Logger &operator=(Logger &&) = delete;

  static void init();
  static void close();

  static void initRankFileOutput();
  static void initFileOutput(const std::string &output);

  static void enableLogFile(bool enable) {
    logFileOs = (enable) ? saveLogFileOs : nullptr;
  }

  static long getElapsedSec() {
    if (!inited) {
      return 0;
    }
    TimePoint finish = std::chrono::high_resolution_clock::now();
    auto elapsed =
        std::chrono::duration_cast<std::chrono::seconds>(finish - start);
    return static_cast<long>(elapsed.count());
  }

  template <typename T> Logger &operator<<(const T &message) {
    if (!inited || isSilent()) {
      return *this;
    }
    if (_type == lt_timed) { // Logger::timed call
      auto seconds = getElapsedSec();
      auto hours = seconds / 3600;
      auto minutes = (seconds % 3600) / 60;
      seconds = seconds % 60;
      char time[30];
      sprintf(time, "%02ld:%02ld:%02ld", hours, minutes, seconds);
      *_os << "[" << time << "] " << message;
      if (logFileOs) {
        *logFileOs << "[" << time << "] " << message;
      }
      // Let Logger::info continue the timed message
      // only for the master rank
      mute();
      return Logger::info;
    } else if (_type == lt_error) { // Logger::error call
      *_os << "[Error] " << message;
      if (logFileOs) {
        *logFileOs << "[Error] " << message;
      }
      // Let Logger::info continue the error message
      // even for a non-master rank
      unmute();
      return Logger::info;
    } else if (_type == lt_info) { // Logger::info call
      *_os << message;
      if (logFileOs) {
        *logFileOs << message;
      }
    } else { // Logger::perrank call
      if (rankLogFileOs) {
        *rankLogFileOs << message;
      }
    }
    // Continue the message using either Logger::info
    // or Logger::perrank
    return *this;
  }

  Logger &operator<<(std::ostream &(*manip)(std::ostream &)) {
    if (!inited || isSilent()) {
      return *this;
    }
    if (_type != lt_perrank) {
      manip(*_os);
      if (logFileOs) {
        manip(*logFileOs);
      }
      // After printing std::endl, ensure that Logger::info
      // is silent for all non-master ranks (the default state)
      // NOTE: any Logger statement is expected to end as std::endl;
      if (manip == static_cast<decltype(manip)>(std::endl)) {
        mute();
      }
    } else {
      if (rankLogFileOs) {
        manip(*rankLogFileOs);
      }
    }
    return *this;
  }

  static Logger info;
  static Logger error;
  static Logger timed;
  static Logger perrank;
  static TimePoint start;
  static std::string outputPref;
  static std::ofstream *logFileOs;
  static std::ofstream *saveLogFileOs;
  static std::ofstream *rankLogFileOs;
  static bool inited;
  static bool fileInited;
};
