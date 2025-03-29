#pragma once

#include <string>



/**
 *  Schedules jobs with the external dependency MPIScheduler,
 *  allowing to schedule independent jobs in parallel with a
 *  specified number of cores per job.
 *  @param outputDir The GeneRax run output directory
 *  @param commandFile The path to the command file
 *                     (see the MPIScheduler documentation)
 *  @param splitImplem Split or fork implementation
 *                     (see the MPIScheduler documentation)
 *  @param execPath The path to the executable to schedule
 *                  (in practice, it will be the same as the
 *                  main executable, but with specific arguments)
 */
class Scheduler {
public:
  Scheduler() = delete;

  static void schedule(const std::string &outputDir,
      const std::string &commandFile,
      bool splitImplem,
      const std::string &execPath);

};


