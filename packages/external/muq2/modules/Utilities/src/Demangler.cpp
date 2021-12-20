#include "MUQ/Utilities/Demangler.h"

#include <cxxabi.h>
#include <memory>

std::string muq::Utilities::demangle(const char* name) {

    int status;

    std::unique_ptr<char, void(*)(void*)> res {
        abi::__cxa_demangle(name, NULL, NULL, &status),
        std::free
    };

    // If the demangling was successful, return the demangled name as a string
    return (status==0) ? res.get() : name ;
}
