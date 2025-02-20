function [extOut] = extRet_Fernald(range, sig, bg, mBsc, mExt, varargin)
% EXTRET_FERNALD obtain extinction coefficient with Fernald method for 
% horizontal-pointintg lidar.
%
% USAGE:
%    [extOut] = extRet_Fernald(range, sig, bg, mBsc, mExt, varargin)
%
% INPUTS:
%    range: numeric
%        horizontal range bin. (m)
%    sig: numeric
%        signal array. (photon count)
%    bg: numeric
%        background signal. (photon count)
%    mBsc: numeric
%        horizontal molecular backscatter coefficient array. (m-1 sr-1)
%    mExt: numeric
%        horizontal molecular extinction coefficient array. (m-1)
%
% KEYWORDS:
%    minSNR: numeric
%    hFullOL: numeric
%    calibration_constant: numeric
%    lr: numeric
%    maxDecomRange: numeric
%    maxDecomHeight: numeric
%    maxDecomThickness: numeric
%    minRefThickness: numeric
%    snr: numeric
%
% OUTPUTS:
%    extOut: numeric
%        aerosol extinction coefficient. (m-1)
%
% HISTORY:
%    2024-12-23: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'range', @isnumeric);
addRequired(p, 'sig', @isnumeric);
addRequired(p, 'bg', @isnumeric);
addRequired(p, 'mBsc', @isnumeric);
addRequired(p, 'mExt', @isnumeric);
addParameter(p, 'minSNR', 1, @isnumeric);
addParameter(p, 'hFullOL', 500, @isnumeric);
addParameter(p, 'calibration_constant', 1, @isnumeric);
addParameter(p, 'lr', 50, @isnumeric);
addParameter(p, 'maxDecomRange', 20000, @isnumeric);
addParameter(p, 'maxDecomThickness', 1000, @isnumeric);
addParameter(p, 'decomSmWin', 10, @isnumeric);
addParameter(p, 'minRefThickness', 300, @isnumeric);
addParameter(p, 'snr', [], @isnumeric);

parse(p, range, sig, bg, mBsc, mExt, varargin{:});

isInOL = (range < p.Results.hFullOL);
nBinsInRef = 50;

if isempty(p.Results.snr)
    snr = lidarSNR(sig, bg);
else
    snr = p.Results.snr;
end
isLowSNR = (snr <= p.Results.minSNR) & (~ isInOL);
idxLowSNR = find(isLowSNR, 1);

if (idxLowSNR > nBinsInRef)
    aBsc = transpose(fernald(range, sig, bg, p.Results.lr, [range(idxLowSNR - 50), range(idxLowSNR)], 1e-6, mBsc, 16));
else
    aBsc = transpose(fernald(range, sig, bg, p.Results.lr, range(idxLowSNR), 1e-6, mBsc, 16));
end
aBsc(isLowSNR) = NaN;
extOut = aBsc * p.Results.lr;

end