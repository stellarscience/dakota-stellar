function writemain(problem,alg)
% Function to store the parameters of the problem and 
% the wanted number of iterations / functionevaluations
% in the respective files.
%
% Joerg Gablonsky, 10/29/00
% Last Changed     04/15/01


FID1 = fopen('../ini/main.ini','w');  

fprintf(FID1,'DIRECT.ini           Name of DIRECT-Init file \n');
fprintf(FID1,'%8i             Problem \n',problem); 
fprintf(FID1,'%8i             Which algorithm \n',alg); 
fclose(FID1); % close the file with problem-definition.
