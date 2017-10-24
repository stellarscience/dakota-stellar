//==============================================================================
// This software developed by Stellar Science Ltd Co.
// Copyright (C) 2016 Stellar Science.
// See README-Stellar for license
//------------------------------------------------------------------------------
#ifndef POLLING_SCHEDULE_HPP
#define POLLING_SCHEDULE_HPP

#include <vector>

namespace Dakota {

  /// Encapsulates a sequence of polling delay values.
  /// The sequence starts at an initial delay, repeating N times, then backs off 
  /// in an exponential fashion by multiplying by N. Therefore the total
  /// delay at a given backoff step is the same as a single delay at the next-
  /// higher step. After a given number of backoffs, the delay remains fixed
  /// at an upper limit.
  class PollingSchedule
  {
  public:

    PollingSchedule(
      unsigned initialDelayInMilliseconds,
      unsigned backoffMultiplier,
      unsigned numberOfBackoffs);

    unsigned nextDelayInMilliseconds();

    void sleep();

    void reset();

  private:
    std::vector<unsigned> delayValues;
    std::vector<unsigned>::size_type currentIndex;

  };

} // namespace Dakota

#endif

