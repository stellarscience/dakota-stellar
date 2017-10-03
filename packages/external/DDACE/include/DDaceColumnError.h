#ifndef DDaceColumnError_H
#define DDaceColumnError_H

#include <stdexcept>
using std::runtime_error;

class DDaceColumnError : public runtime_error {
  public:
    DDaceColumnError();
    DDaceColumnError(std::string message);
};

#endif
