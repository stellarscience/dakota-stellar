% Program to create a table with the number of function evaluations
% needed to find a function value within DIRperc of the global optimal
% function value. The table created shows these numbers for the original
% DIRECT algorithm and our modified algorithm, DIRECT-l.
%
% Joerg Gablonsky, 10/31/00
% Last Changed     04/16/01

numproblems = 12;
startprob = 0;
setdirect;
DIRnumf = 20000;
DIRperc = 0.01;
result = zeros(numproblems,7);

problem = 0;
for problem = startprob:numproblems+startprob
  setproblem;
  disp(sprintf('Problem : %s ',problemname))
  for DIRJones = 0:1
    writeDIRECT(DIReps,DIRnumf,numberT,DIRJones, DIRperc);
    writemain(problem,0);
    ! \rm ../result.out
    ! cd ..; TestDIRect > nul ; \rm nul; cd matlab
    helpi = counting;
    i1 = 2+DIRJones*3;
    i2 = i1 + 2;
    result(problem+1-startprob,i1:i2) = helpi(2:4);
    result(problem+1-startprob,1) = helpi(1);
  end
end
fileoutput
disp(sprintf('The results of the optimizations can be found in the file results.txt'))
