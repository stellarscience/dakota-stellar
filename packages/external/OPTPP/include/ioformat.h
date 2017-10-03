// Uncomment this only is you are REALLY sure templates work on your system
//#define WANT_TEMPLATES  

#ifndef ioformat_h
#define ioformat_h

//Convert to string 10/05/01

#include <iostream>
#include <string>

typedef std::ios_base::fmtflags opt_mode;


namespace OPTPP {

class oformatstate
{
public:
  oformatstate(std::ostream& ut);
  oformatstate (char code, int w=0, int p=0, char c=' ', opt_mode f = std::ios_base::fixed);
  friend std::ostream& operator << (std::ostream& ut, oformatstate const& fmt);
  friend std::ostream& operator >> (std::ostream& ut, oformatstate& fmt);
  
private:
  int owidth;
  int oprecision;
  char ofill;
  opt_mode oflags;
};

std::string format(double val, oformatstate const& fmt);
std::string format(int    val, oformatstate const& fmt);

inline std::string
d(int val, int w=0, int p=0, char c=' ', opt_mode f = std::ios_base::fixed)
{
  oformatstate fmt('d', w, p, c, f);
  return format(val, fmt);
}
inline std::string
e(double val, int w=0, int p=0, char c=' ', opt_mode f = std::ios_base::fixed)
{
  oformatstate fmt('e', w, p, c, f);
  return format(val, fmt);
}
inline std::string
f(double val, int w=0, int p=0, char c=' ', opt_mode f = std::ios_base::fixed)
{
  oformatstate fmt('f', w, p, c, f);
  return format(val, fmt);
}

} // namespace OPTPP

#endif
