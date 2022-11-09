global LISMO_VARS;

LISMO_VARS.author = 'Zhenping Yin';
LISMO_VARS.email = 'zp.yin@whu.edu.cn';
LISMO_VARS.program = 'LISMO';
LISMO_VARS.version = '1.0';
LISMO_VARS.updated_time = datenum(2022, 11, 9);
LISMO_VARS.description = 'LISMO: LIdar for Sea-fog MOnitoring Simulation Toolbox';
LISMO_VARS.projectDir = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(projectDir, 'include')));
addpath(genpath(fullfile(projectDir, 'lib')));
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'analysis')));

disp('Install toolbox successfully!');