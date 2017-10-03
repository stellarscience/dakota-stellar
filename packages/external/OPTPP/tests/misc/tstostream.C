#include <fstream>

class OptimizeClass {

private:
  int x_optout_fd;

protected:
  int  dim;
  ofstream optout;
  int optout_fd;

public:

  OptimizeClass(): x_optout_fd(-1),dim(0){
    optout.open("optdefault.out");
    if (!optout) {
      cout << "OptimizeClass::SetOutPut: Can't open default output file\n";
      optout_fd = 0;
    }
    else {
      optout_fd = optout.rdbuf()->fd();
    }
  };

  virtual ~OptimizeClass() {
    cout << "OptimizeClass destructor\n";
    Cleanup();
  }

  int SetOutputFile(const char *filename) { 

    if (x_optout_fd == -1) {  // Change the default output file
      optout.close();
      optout.open(filename);
      if (!optout) {
	cout << "OptimizeClass::SetOutPut: Can't open " << filename << endl;
	optout_fd = 0;
      }
      else {
	optout_fd = optout.rdbuf()->fd();
      }
    }
    else {
      cout << "OptimizeClass::SetOutPut: File already attached\n";
    }
    return optout_fd;
  }
  int SetOutputFile(int FileDescriptor) { 

    if (x_optout_fd == -1) {  // Change the default output file
      optout.close();
      optout.attach(FileDescriptor);
      if (optout.good()) {
	optout_fd = FileDescriptor;
	x_optout_fd = FileDescriptor;
	optout.flush();
      }
      else {
	cout << "OptimizeClass::SetOutputFile: Can't attach optout file\n";
	cout << "OptimizeClass::SetOutputFile: file descriptor " 
	     << optout_fd << endl;
	cout << "OptimizeClass::SetOutputFile: Are you sure the file " 
	     << "is open\n";
      }
    }
    else {
      cout << "OptimizeClass::SetOutputFile: optout file already set\n";
    }
    return optout_fd;
  }
  int SetOutputFile(ofstream& fout) { 

    int FileDescriptor;
    cout << "OptimizeClass::SetOutputFile: File stream version\n";
    if (x_optout_fd == -1) {  // Change the default output file
      optout.close();
      FileDescriptor = fout.rdbuf()->fd();
      optout.attach(FileDescriptor);
      if (optout.good()) {
	optout_fd = FileDescriptor;
	x_optout_fd = FileDescriptor;
	optout.flush();
      }
      else {
	cout << "OptimizeClass::SetOutputFile: Can't attach optout file\n";
	cout << "OptimizeClass::SetOutputFile: file descriptor " 
	     << optout_fd << endl;
	cout << "OptimizeClass::SetOutputFile: Are you sure the file " 
	     << "is open\n";
      }
    }
    else {
      cout << "OptimizeClass::SetOutputFile: optout file already set\n";
    }
    return optout_fd;
  }

  void  PrintStatus(char *s) {
    optout << "\n\n=========  " << s << "  ===========\n\n";
    optout << "Optimize Class\n";
    optout << "\n\n=========  " << s << "  ===========\n\n";
  }

  void  Cleanup() {optout.flush();};
};

void main () 
{
  char filename[80];
  OptimizeClass tstopt1, tstopt2, tstopt3, tstopt4;
  int fd, fd1, fd2;
  int FileDescriptor;

  //
  //  1. Test of default
  //
  tstopt1.PrintStatus("SetOutputFile not called");
  tstopt1.Cleanup();

  //
  //  2. Call with file name
  //
  cout << "Enter file name:";
  cin  >> filename;
  fd1 = tstopt2.SetOutputFile(filename);
  cout << "Test case 2: File descriptor = " << fd1 << endl;

  tstopt2.PrintStatus("Test of SetOutputFile using filename constructor");
  tstopt2.Cleanup();

  //
  // 3a. Call with file descriptor
  //
  
  cout << "Enter second file name:";
  cin  >> filename;
  ofstream fout(filename);
  FileDescriptor = fout.rdbuf()->fd();
  fd2 = tstopt3.SetOutputFile(FileDescriptor);
  cout << "Test case 3: FileDescriptor = " << FileDescriptor << endl;
  cout << "Test case 3: fd2            = " << fd2 << endl;

  tstopt3.PrintStatus("Test of SetOutputFile using file descriptor");

  //
  // 3b. Call again: This should generate an error
  //
  
  fd2 = tstopt3.SetOutputFile(fd1);
  cout << "You should have gotten an error message right above this one\n";
  tstopt3.Cleanup();

  //
  // 4. Call with ofstream
  //
  
  cout << "Enter third (and final) file name:";
  cin  >> filename;
  ofstream fout4(filename);
  FileDescriptor = fout4.rdbuf()->fd();
  fd2 = tstopt4.SetOutputFile(fout4);
  cout << "Test case 4: FileDescriptor = " << FileDescriptor << endl;
  cout << "Test case 4: fd2            = " << fd2 << endl;
  tstopt4.PrintStatus("Test of SetOutputFile using ofstream");
  tstopt4.Cleanup();

  cout << "End of test\n";
    
}
