#!/bin/sh

/bin/rm -f failures.log details.log results.dat
touch failures.log details.log results.dat

echo "-------------------------------------------------------" >> details.log;
echo " " >> details.log;
echo "OPT++ hock test results: " `date` >> details.log ;
echo " " >> details.log;
echo "-------------------------------------------------------" >> details.log;

for reg_test in ./tsthock1 ./tsthock2 ./tsthock5 ./tsthock6 ./tsthock7 ./tsthock10 ./tsthock13 ./tsthock14 ./tsthock26 ./tsthock28 ./tsthock35 ./tsthock65 ./tsthock77 ./tsthock78
do
#
#   First check to see if the program runs at all
#
    if $reg_test; then
	:
    else
	echo -e "$reg_test crashed\n" >> failures.log;
    fi
#
#   Now check to see if the test problems passed/failed
#
    set status = `grep PASSED $reg_test.out >> details.log`
    set status = `grep FAILED $reg_test.out >> details.log`
    set status = `grep PASSED $reg_test.out >> results.dat`
    set status = `grep FAILED $reg_test.out >> results.dat`
done
#
#   If any of the tests resulted in a failure store result in failures log
#   The flag -s checks for existence and nonzero size file
#
set status = `grep -l FAILED *.out >> failures.log`

#
#  Print out messages depending on whether all tests passed or not
#
if [ -s failures.log ]; then
    echo "----------------------------------------------------------------";
    echo " ";
    echo "  OPT++ hock test results: FAILURES DETECTED";
    echo "  ";
    echo "      " `date` ;
    echo " ";
    echo "The following tests failed:";
    echo `cat failures.log`;
    echo " "
    echo "Please check the *.out files for details";
    echo "----------------------------------------------------------------";

    echo "----------------------------------------------------" >> details.log;
    echo "OPT++ hock test results: FAILURES DETECTED" >> details.log;
    echo "----------------------------------------------------" >> details.log;

else
    echo "----------------------------------------------------------------";
    echo " ";
    echo "  OPT++ hock test results: all tests PASSED";
    echo "  ";
    echo "      " `date` ;
    echo "----------------------------------------------------------------";

    echo "----------------------------------------------------" >> details.log;
    echo "OPT++ hock test results: all tests PASSED" >> details.log;
    echo "-----------------------
-----------------------------" >> details.log;

fi

exit 0
