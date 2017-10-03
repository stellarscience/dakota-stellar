#ifndef DDaceZeroDivide_H
#define DDaceZeroDivide_H

#include <stdexcept>
using std::runtime_error;

class DDaceZeroDivide : public runtime_error {
  public:
    DDaceZeroDivide();
    DDaceZeroDivide(std::string message);
};

#endif
