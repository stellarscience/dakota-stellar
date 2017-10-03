#ifndef DDaceObservationError_H
#define DDaceObservationError_H

#include <stdexcept>
using std::runtime_error;

class DDaceObservationError : public runtime_error {
  public:
    DDaceObservationError();
    DDaceObservationError(std::string message);
};

#endif
