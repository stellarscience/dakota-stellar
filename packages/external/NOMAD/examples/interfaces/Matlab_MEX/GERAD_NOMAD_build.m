%% GERAD NOMAD Build for Matlab

% This file will help you compile NOMAD for use with MATLAB. 

% To recompile you will need to get / do the following:

% 1) Get NOMAD
% NOMAD is available from http://www.gerad.ca/NOMAD/PHP_Forms/Download.php.
% Complete the download form then download the latest version. Define the
% $NOMAD_HOME environment variable.

% 2) NOMAD MEX Interface
% The NOMAD MEX Interface is a simple MEX interface written to use NOMAD.

% 3) Compile the MEX File
% The code below will build the NOMAD MEX file. Once you have completed all the
% above steps, simply run this file to compile NOMAD! You MUST BE in the 
% base directory of OPTI!

clear nomad

switch(computer)
case 'PCWIN'
libdir = ' -Lwin32\';
case 'PCWIN64' 
libdir = ' -Lwin64\';        
case 'GLNX86'
libdir = 'glnx86/';
case 'GLNXA64'
libdir = 'glnxa64/';
case 'MACI64'
libdir = 'maci64/';
end

clear nomad_home nomad_src;


fprintf('\n------------------------------------------------\n');
fprintf('NOMAD MEX FILE BUILD --- GERAD VERSION \n\n');

nomad_home = getenv('NOMAD_HOME');
if (length(nomad_home)<1)
    error('opti:nomad','Please set NOMAD_HOME variables properly with the command setenv(''NOMAD_HOME'',ARG1)!');
end
nomad_src=[nomad_home filesep 'src' filesep];


%Get NOMAD Libraries
post = [' -I.  -I' nomad_src ' -lm -lut -output nomad'];

%CD to Source Directory
cdir = cd;
cd 'Source';

%Compile & Move
pre = ['mex -v -g -largeArrayDims nomadmex.cpp ' nomad_src 'Parameters.cpp ' nomad_src 'Barrier.cpp ' nomad_src 'Cache.cpp '...
nomad_src 'Cache_File_Point.cpp ' nomad_src 'Cache_Point.cpp ' nomad_src 'Cache_Search.cpp ' nomad_src 'Clock.cpp '...
nomad_src 'Direction.cpp ' nomad_src 'Directions.cpp ' nomad_src 'Display.cpp '...
nomad_src 'Double.cpp ' nomad_src 'Eval_Point.cpp ' nomad_src 'Evaluator.cpp ' nomad_src 'Evaluator_Control.cpp ' nomad_src 'Exception.cpp '...
nomad_src 'Extended_Poll.cpp ' nomad_src 'L_Curve.cpp ' nomad_src 'LH_Search.cpp ' nomad_src 'OrthogonalMesh.cpp ' nomad_src 'Mads.cpp ' nomad_src 'Model_Sorted_Point.cpp '...
nomad_src 'Model_Stats.cpp ' nomad_src 'Multi_Obj_Evaluator.cpp ' nomad_src 'Parameter_Entries.cpp '...
nomad_src 'Parameter_Entry.cpp ' nomad_src 'Pareto_Front.cpp ' nomad_src 'Pareto_Point.cpp ' nomad_src 'Phase_One_Evaluator.cpp '...
nomad_src 'Phase_One_Search.cpp ' nomad_src 'Point.cpp ' nomad_src 'Priority_Eval_Point.cpp ' nomad_src 'Quad_Model.cpp '...
nomad_src 'Quad_Model_Evaluator.cpp ' nomad_src 'Quad_Model_Search.cpp ' nomad_src 'Random_Pickup.cpp ' nomad_src 'RNG.cpp '...
nomad_src 'Signature.cpp ' nomad_src 'Slave.cpp ' nomad_src 'SMesh.cpp ' nomad_src 'Speculative_Search.cpp ' nomad_src 'Stats.cpp ' nomad_src 'utils.cpp '...
nomad_src 'Variable_Group.cpp ' nomad_src 'VNS_Search.cpp ' nomad_src 'XMesh.cpp'];

try
    eval([pre post])
    movefile(['nomad.' mexext],['..' filesep],'f')
    cd(cdir);
    clear nomad_home nomad_src cdir post pre libdir;
    fprintf('Done!\n');
catch ME
    cd(cdir);
	clear nomad_home nomad_src cdir post pre libdir;
    error('opti:nomad','Error Compiling NOMAD!\n%s',ME.message);
end
fprintf('------------------------------------------------\n');
