% Set the name of the problem.
%
% Joerg Gablonsky, 10/29/00
% Last Changed     03/14/01

fid = fopen('../ini/problems.ini');
help = fgetl(fid);
numproblems1 = sscanf(help,'%6i');
for q = 1:problem+1
  help1 = fgetl(fid);
  problemname = sscanf(help1,'%s40');
end
fclose(fid);
problemname = strcat('../problem/',problemname);

fid = fopen(problemname);
problemname = fgetl(fid);
help1 = fgetl(fid);
n = str2num(help1(1:2));
fclose(fid);

return

      if (problem == 0)  
         problemname = 'Constant';
      elseif (problem == 1)  
         problemname = 'Linear';
      elseif (problem == 2)  
         problemname = 'Quadratic';
      elseif (problem == 3)  
         problemname = 'Gomez 3';
      elseif (problem == 4)  
         problemname = 'Branin';
      elseif (problem == 5)  
         problemname = 'Shekel-5';
      elseif (problem == 6)  
         problemname = 'Shekel-7';
      elseif (problem == 7)  
         problemname = 'Shekel-10';
      elseif (problem == 8)  
         problemname = 'Hartman-3';
      elseif (problem == 9)  
         problemname = 'Hartman-6';
      elseif (problem == 10)  
         problemname = 'Goldprice';
      elseif (problem == 11)  
         problemname = 'Sixhump';
      elseif (problem == 12)  
         problemname = 'Shubert';
      elseif (problem == 13)  
         problemname = 'Sphere';
      elseif (problem == 14)  
         problemname = 'Saddle';
      elseif (problem == 15)  
         problemname = 'Step';
      end
