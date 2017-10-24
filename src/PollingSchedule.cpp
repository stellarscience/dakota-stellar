//==============================================================================
// This software developed by Stellar Science Ltd Co.
// Copyright (C) 2016 Stellar Science.
// See README-Stellar for license
//------------------------------------------------------------------------------
#include <stdexcept>
#include <sys/types.h> // MAY REQUIRE ifndef(HPUX)

// Note: refer to SysCallApplicInterface.cpp for these includes.
#if defined(_WIN32) || defined(_MSC_VER) || defined(__MINGW32__)
#define NOMINMAX
#include <windows.h> // for Sleep()
#elif defined(HAVE_UNISTD_H)
#include <unistd.h> // for usleep()
#endif

#include "PollingSchedule.hpp"

namespace Dakota {

  PollingSchedule::PollingSchedule(
    unsigned initialDelayInMilliseconds,
    unsigned backoffMultiplier,
    unsigned numberOfBackoffs)
    :
    delayValues(),
    currentIndex(0)
  {
    if (initialDelayInMilliseconds <= 0) {
      throw std::invalid_argument("initialDelayInMilliseconds must be positive");
    }
    if (backoffMultiplier <= 0) {
      throw std::invalid_argument("backoffMultiplier must be positive");
    }

    unsigned delay = initialDelayInMilliseconds;
    for (int backoffCount = 0; backoffCount < numberOfBackoffs + 1; backoffCount++) {
      for (int attempt = 0; attempt < backoffMultiplier; attempt++) {
        delayValues.push_back(delay);
      }
      delay *= backoffMultiplier;
    }
  }

  void PollingSchedule::sleep()
  {
    const unsigned delay = nextDelayInMilliseconds();

    // Note: refer to SysCallApplicInterface.cpp for code to sleep.
#if defined(_WIN32) || defined(_MSC_VER) || defined(__MINGW32__)
    Sleep(delay);
#elif defined(HAVE_USLEEP)
    usleep(1000 * delay); // 1000 microseconds = 1 millisec
#endif // SLEEP

    if (currentIndex < delayValues.size()) {
      currentIndex++;
    }
  }

  unsigned PollingSchedule::nextDelayInMilliseconds()
  {
    return (currentIndex >= delayValues.size())
            ? delayValues.back()
            : delayValues.at(currentIndex);
  }

  void PollingSchedule::reset()
  {
    currentIndex = 0;
  }

} // namespace Dakota
