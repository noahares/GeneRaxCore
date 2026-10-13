#include "Logger.hpp"

#include <parallelization/ParallelContext.hpp>

Logger Logger::info;
Logger Logger::error;
Logger Logger::timed;
Logger Logger::perrank;
TimePoint Logger::start;
std::string Logger::outputPref;
std::ofstream *Logger::logFileOs = nullptr;
std::ofstream *Logger::saveLogFileOs = nullptr;
std::ofstream *Logger::rankLogFileOs = nullptr;
bool Logger::inited = false;
bool Logger::fileInited = false;

bool Logger::isSilent() const {
  return _muted && (_type == lt_timed || _type == lt_info) &&
         ParallelContext::getRank();
}

void Logger::init() {
  if (!inited) {
    inited = true;
    info.setType(lt_info);
    info.setStream(&std::cout);
    error.setType(lt_error);
    error.setStream(&std::cout);
    timed.setType(lt_timed);
    timed.setStream(&std::cout);
    perrank.setType(lt_perrank);
    perrank.setStream(nullptr);
    start = std::chrono::high_resolution_clock::now();
  }
  ParallelContext::barrier();
}

void Logger::initRankFileOutput() {
  if (inited && fileInited && !rankLogFileOs) {
    auto rankLogFile = Logger::outputPref + "_rank_" +
                       std::to_string(ParallelContext::getRank()) + ".log";
    rankLogFileOs = new std::ofstream(rankLogFile);
  }
}

void Logger::initFileOutput(const std::string &output) {
  if (inited && !fileInited) {
    outputPref = output;
    fileInited = true;
    // only the master rank can write to the log file
    if (ParallelContext::getRank() == 0) {
      auto logFile = Logger::outputPref + ".log";
      std::cout << "Logs will also be printed into " << logFile << std::endl;
      // init the log file, it will be appended on a restart
      logFileOs = new std::ofstream(logFile, std::ios::app);
      saveLogFileOs = logFileOs;
    }
  }
}

void Logger::close() {
  outputPref.clear();
  delete saveLogFileOs;
  delete rankLogFileOs;
  logFileOs = saveLogFileOs = rankLogFileOs = nullptr;
  inited = fileInited = false;
  ParallelContext::barrier();
}
