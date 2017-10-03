#include <iostream>
#include <fstream>

int main ()
{
  char *opt_input  = {"opt.input"};
  ifstream optin(opt_input);

  int count, keyword_count=0;
  char token[80];
  
  optin >> token;
  count = optin.gcount();
  cout << " count = " << count << "\n";

  while (!optin.eof()) {

    keyword_count++;
    cout << keyword_count << " keyword = " << token << "\n";
    optin >> token;
    count = optin.gcount();
    cout << " count = " << count << "\n";
  }

}
