function writeDirect(DIReps,DIRnumf,numberT,DIRJones, DIRperc)
% function to store the parameters of the problem and 
% the wanted number of iterations / functionevaluations
% in the respective files.
%
% Joerg Gablonsky, 10/31/00
% Last Changed     03/14/01


FID1 = fopen('../ini/DIRECT.ini','w');  

fprintf(FID1,'%14.10e    eps \n',DIReps);     % the eps-value
fprintf(FID1,'%8i            maximal function evaluations \n',DIRnumf); % the approx. max. # of functions
fprintf(FID1,'%8i            maximal number of iterations \n',numberT); % the approx. max. # of iterations
fprintf(FID1,'%8i            Use JONES way (0) or our way (1) \n',DIRJones);
fprintf(FID1,'%14.10e    Percentage when to stop. \n',DIRperc);     % the eps-value
fprintf(FID1,'%14.10e    Volume Percentage when to stop. \n',-1);     % the eps-value
fprintf(FID1,'%14.10e    Bound on the measure. Stop when mesure is below this bound \n',-1);     % the eps-value

fclose(FID1); % close the file with problem-definition.
