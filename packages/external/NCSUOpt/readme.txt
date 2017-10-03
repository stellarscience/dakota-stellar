DIRECT Version 2.0.4 Readme

Changes from Version 2.0.3

1.) Added another error message / condition:
    Ierror =   -6   Error in DIRDoubleInsert, that is an error occured   
                    DIRECT tried to add all hyperrectangles with the same
                    size and function value at the center. Either        
                    increase maxdiv or use our modification (Jones = 1).

2.) Changed DIRDivide to prevent array bounds error when array bounds checking
    is enabled in compilers.

3.) Changed the test for the volume termination criteria. 

4.) Removed output of cdata(1) in DIRheader since cdata(1) does not need
    to be initialised.

5.) Changed if and while statements with .AND. conditions. There are compilers
    which evaluate the second expression in an .AND. statement even if the
    first statement is already false.

6.) Changed read statements in the example main program. Some computers did not
    read in file names from main.ini correctly without FORMAT statements.

7.) Note: Some computers don't like Maxfunc set to 90000 in DIRect.f. Reduce
          Maxfunc to 40000 (or another smaller numer) should solve this problem.

I want to thank Jian He, Mario Innorta and Manfred Stickel for pointing out
these possibilities for me to improve the code.

Joerg Gablonsky

Raleigh, NC, 7/16/01
