global LISMO_VARS;

%% Parameter Initialization
distArr = 7.5:15:20000;   % distance array. (m)
eleAngle = 3;   % elevation angle. (degree)
laserWL = 1030;   % laser wavelength. (nm)
seaFogType = 'sea-fog-none';
saveFile = 'sea-fog-none-nighttime-test.mat';

%% Atmosphere Module
height = distArr * sin(eleAngle / 180 * pi);
[mBsc, mExt] = MolModel(height, laserWL, 'meteor', 'standard_atmosphere');
[aBsc, aExt] = AerModel(height, 'scene', 'marine-moderate');
[fBsc, fExt] = SeaFogModel(distArr, laserWL, 'scene', seaFogType);
tBsc = mBsc + aBsc + fBsc;
tExt = mExt + aExt + fExt;

%% Signal Simulation
dataSim = LISMO_Model(distArr, ...
    'tBsc', tBsc, ...
    'tExt', tExt, ...
    'eleAngle', eleAngle, ...
    'laserWL', laserWL);

%% Save results
save(fullfile(LISMO_VARS.projectDir, 'results', saveFile), 'distArr', 'eleAngle', 'aBsc', 'aExt', 'fBsc', 'fExt', 'mBsc', 'mExt', 'dataSim');