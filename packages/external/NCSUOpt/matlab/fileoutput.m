% Function to write down a LaTex table with the number of function evaluations
% needed for the testproblems.
%
% Joerg Gablonsky, 10/29/00
% Last changed     04/15/01

title1 = '\\begin{tabular}{lccccccc} \n';
title2 = 'Problem        & \\multicolumn{3}{c}{DIRECT} &  \\multicolumn{3}{c}{DIRECT-l} \\\\ \n';
title3 = '               &  f.eval. & Time & $p$  &  f.eval. & Time & $p$  & \\\\ \n';
title4 = '\\end{tabular}\n';
form1  = '%28s &   %5i  &  %8.3f  &  %10.2E & %5i  &  %8.3f  &  %10.2E  \\\\ \n';

n = size(result,1);

FID1 = fopen('results.txt','w');  

fprintf(FID1,title1);
fprintf(FID1,title2);
fprintf(FID1,title3);
for i=1:n
  problem = result(i,1);
  setproblem;
  fprintf(FID1,form1,problemname,result(i,2:7));
end
fprintf(FID1,title4);
fclose(FID1);

