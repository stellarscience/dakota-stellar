//
// Test program for cout
//
#include <fstream>
int main ()
{
  
  char *status_file = {"cout.out"};
  ofstream opt_fp(status_file);
  cout = opt_fp;
  cout << "Testing reassignment of cout\n";

}
