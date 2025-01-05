function [extOut] = extRet_Fernald_plus(range, sig, bg, mBsc, mExt, varargin)
% EXTRET_FERNALD_plus obtain extinction coefficient with Fernald method for 
% horizontal-pointintg lidar.
%
% USAGE:
%    [extOut] = extRet_Fernald_plus(range, sig, bg, mBsc, mExt, varargin)
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

parse(p, range, sig, bg, mBsc, mExt, varargin{:});

isInOL = (range < p.Results.hFullOL);
idxFullOL = find(~ isInOL, 1);
rcs = smooth(sig .* range.^2, 8);
mSig = mBsc .* exp(-2 * nancumsum(mExt .* [range(1); diff(range)]));
dz = (range(2) - range(1));

%% SNR
snr = lidarSNR(sig, bg);
isLowSNR = (snr <= p.Results.minSNR) & (~ isInOL);
idxLowSNR = find(isLowSNR, 1);

%% Scattering ratio
sr = sig .* range.^2 ./ mSig;

isNoisy = false;
if isempty(idxLowSNR) || (idxLowSNR < 2)
    isNoisy = true;
end

%% Douglas-Peucker Decomposition
hRefIdxMat = [];
if ~ isNoisy
    % DPIdx = DouglasPeucker(sr(idxFullOL:idxLowSNR), range(idxFullOL:idxLowSNR), 0.1, ...
    %     p.Results.hFullOL, p.Results.maxDecomRange, ...
    %     p.Results.maxDecomThickness, p.Results.decomSmWin) + idxFullOL - 1;

    sr(sr <= 0) = NaN;
    deriv1 = movingslope(log(sr), 40) / dz * (-0.5);
    deriv2 = movingslope(deriv1, 10) / dz;
    isStable = abs(smooth(deriv2, 10)) <= (0.1e-4 / dz);
    [L, num] = bwlabel(isStable);

    for iL = 1:num
        startIdx = find(L == iL, 1);
        endIdx = find(L == iL, 1, 'last');

        if (range(endIdx) - range(startIdx)) >= p.Results.minRefThickness
            hRefIdxMat = cat(1, hRefIdxMat, [startIdx, endIdx]);
        end
    end

    % isGoodCom = false;
    % for iDP = 1:(length(DPIdx) - 1)
    %     if (range(DPIdx(iDP + 1)) - range(DPIdx(iDP))) > p.Results.minRefThickness
    %         hRefIdxMat = cat(1, hRefIdxMat, [DPIdx(iDP), DPIdx(iDP + 1)]);
    %     end
    % end

    if ~ isempty(hRefIdxMat)
        isGoodCom = true;
    end
end

%% Reference Value
extRef = [];
extRefStd = [];
if isGoodCom
    for iRef = 1:size(hRefIdxMat, 1)
        isGoodFit = true;
        isInRef = hRefIdxMat(iRef, 1):hRefIdxMat(iRef, 2);
        try
            [~, thisExtRef, ~, thisExtRefStd] = chi2fit(range(isInRef), log(rcs(isInRef)), sqrt(sig(isInRef)) .* range(isInRef).^2);
        catch
            isGoodFit = false;
        end

        if isGoodFit
            isPotentialRef = false;

            thisExtRef = thisExtRef * (-0.5);
            thisExtRefStd = thisExtRefStd * 0.5;
            if (thisExtRef > 0) && (thisExtRef <= 2e-3)
                isPotentialRef = true;
            end

            if isPotentialRef
                extRef = cat(2, extRef, thisExtRef);
                extRefStd = cat(2, extRefStd, thisExtRefStd);
            else
                extRef = cat(2, extRef, NaN);
                extRefStd = cat(2, extRefStd, NaN);
            end
        end
    end
end

finalRefH = [];
finalRefVal = [];
isGoodRef = any(~ isnan(extRef));
if isGoodRef
    % [~, bestRefIdx] = min(abs(extRefStd ./ extRef));   % 找消光中位数对应的参考高度区间
    bestRefIdx = find(~ isnan(extRef), 1, 'last');   %找最后的参考高度区间
    finalRefH = [range(hRefIdxMat(bestRefIdx, 1)), range(hRefIdxMat(bestRefIdx, 2))];
    finalRefVal = extRef(bestRefIdx);
end

%% Retrieval
if isGoodRef
    aBsc = transpose(fernald(range, sig, bg, p.Results.lr, finalRefH, finalRefVal / p.Results.lr, mBsc, 8));
else
    aBsc = transpose(fernald(range, sig, bg, p.Results.lr, range(idxLowSNR), 1e-6, mBsc, 16));
end
aBsc(isLowSNR) = NaN;
extOut = aBsc * p.Results.lr;

end