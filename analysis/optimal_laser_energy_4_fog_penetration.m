%------------------------------------------------------------------------------%
% Simulate penetration depth for a given laser energy.
%------------------------------------------------------------------------------%

global LISMO_VARS;

%% Parameter Definition
distArr = 7.5:15:20000;
eleAngle = 0;
laserWL = 1064;
pulseEn = 0.15;
runSimulation = true;
fogType = 'fog-lsw';

%% Signal Simulation
height = distArr * sin(eleAngle / 180 * pi);
[mBsc, mExt] = MolModel(height, laserWL, 'meteor', 'standard_atmosphere');
[aBsc, aExt] = AerModel(height, 'scene', 'marine-moderate');
[fBsc, fExt] = SeaFogModel(distArr, laserWL, 'scene', fogType, 'distLayerFront', 0, 'distLayerBack', 10000);
dataSim = LISMO_Model(distArr, 'tBsc', mBsc + aBsc + fBsc, ...
                               'tExt', mExt + aExt + fExt, ...
                               'eleAngle', eleAngle, ...
                               'laserWL', laserWL, ...
                               'pulseEn', pulseEn, ...
                               'accShots', 25000, ...
                               'PB', 0, ...
                               'visible', 'on', ...
                               'ylim', [0, 4]);