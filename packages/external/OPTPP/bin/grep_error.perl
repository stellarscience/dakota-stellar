#!/usr/bin/perl
#
#  This is a perl script that is used to determine if an error has
#  occured within a src build. The scripts searches the output
#  for specific keywords.
#
################################################################
#  Usage:  	grep_error.perl
#  Written:	Mario Alleva
#  Date:	Jan, 2003
################################################################

# for all dakota*.in files, create a baseline test file
$input = @ARGV[0]; # base file name to parse, e.g. make_utilib.out

# create necessary filenames
$output = $input;
substr($output, -3, 3) = "err"; # e.g., create make_utilib.err

# remove old error file, if present
if (-e $output) {
  unlink $output;
}

#define strings to search for
$search = " ERROR |\tERROR| Error |\tError| Error:| error:";

#define strings to ignore
$ignore = "ignored";

# read input file until EOF
open (INPUT_FILE, $input) || die "cannot open file $!" ;
print "Grepping for Errors in $input \n";
$found = 0;
while (<INPUT_FILE>) { # read each line of file
  if (/$search/) {
    if (/$ignore/) { # check if error should be ignored
      print ("Ignore error \n");
    }
    else { # actual build error
      if ($found == 0) { # if err file not yet open, open it
	open (ERROR_FILE, ">$output") || die "cannot open file $!" ;
      }
      $found = 1;
      print ERROR_FILE $_;
    }
  }
}
# close both files
close (INPUT_FILE);
if ($found == 1) {
  close (ERROR_FILE);
}
