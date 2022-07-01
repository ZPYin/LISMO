projectDir = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(projectDir, 'include')));
addpath(genpath(fullfile(projectDir, 'lib')));
addpath(genpath(fullfile(projectDir, 'data')));

disp('Install toolbox successfully!');