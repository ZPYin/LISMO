function [tExt] = extRet_Xian(range, signal, bg, varargin)
% EXTRET_XIAN calculate extinction coefficient based on Xian's method.
%
% USAGE:
%    [tExt] = extRet_Xian(range, rcs, varargin)
%
% INPUTS:
%    range: numeric
%        distance in range bins. (m)
%    signal: numeric
%        background removed backscatter signal in photon counts.
%    bg: numeric
%        background in photon counts
%
% KEYWORDS:
%    rangeFullOverlap: numeric
%        range of full overlap achieved. (m)
%    minSNR: numeric
%        minimum SNR for signal retrieval.
%    flagConvert550: logical
%        convert to 550 nm extinction coefficient. (bug in Xian's algorithm)
%
% OUTPUTS:
%    tExt: numeric
%        total extinction coefficient. (m^-1)
%
% REFERENCES:
%    Xian, J., Han, Y., Huang, S., Sun, D., Zheng, J., Han, F., Zhou, A., Yang, S., Xu, W., Song, Q., Wei, L., Tan, Q., and Li, X.: Novel Lidar algorithm for horizontal visibility measurement and sea fog monitoring, Opt. Express, 26, 34853-34863, 10.1364/OE.26.034853, 2018.
%
% HISTORY:
%    2023-03-07: first edition by Zhenping
%    2024-06-10: update by Zhenping - bug fix (2 times difference).
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'range', @isnumeric);
addRequired(p, 'signal', @isnumeric);
addRequired(p, 'bg', @isnumeric);
addParameter(p, 'rangeFullOverlap', 500.0, @isnumeric);
addParameter(p, 'minSNR', 0.5, @isnumeric);
addParameter(p, 'flagConvert550', true, @islogical);

parse(p, range, signal, bg, varargin{:});

tExt = NaN(size(signal));   % total extinction (aerosol + molecule)

% search the first range index with full overlap
[~, idxMinOL] = min(abs(range - p.Results.rangeFullOverlap));

% obtaining SNR profile
SNR = lidarSNR(signal, bg, 'bgBins', [length(signal) - 15, length(signal)]);

RCS = signal .* range.^2;   % range corrected signal

RCS_A = RCS(idxMinOL);
ratioAB = RCS / RCS_A;

isFullOverlap = (range >= p.Results.rangeFullOverlap);
isLowSNR = (SNR <= p.Results.minSNR);
detectLimitIdx = find(isFullOverlap & isLowSNR, 1);
isInvalid = (~ isFullOverlap) | isLowSNR;

if (~ isempty(detectLimitIdx))
    % find the index with limited signal
    ratioAB(isInvalid | (ratioAB <= 0)) = NaN;
    [~, idxB] = min(ratioAB(1:detectLimitIdx));

else
    % too strong signal (which is not expected)
    warning('Signal saturation.');
    return;

end

avgLen = 20;
if ((idxB - ceil(avgLen / 2)) > 1) && ((idxB + ceil(avgLen / 2)) < length(range))
    RCS_B = nanmean(RCS((idxB - ceil(avgLen / 2)):(idxB + ceil(avgLen / 2))));
elseif ((idxB - ceil(avgLen / 2)) < 1)
    RCS_B = nanmean(RCS(1:(idxB + ((idxB - ceil(avgLen / 2)) > 1))));
else
    RCS_B = nanmean(RCS((idxB - ((idxB - ceil(avgLen / 2)) > 1)):end));
end

if (idxB ~= idxMinOL)
    I_A_B = trapz(range(idxMinOL:idxB), RCS(idxMinOL:idxB));
else
    I_A_B = 0;
end

% wavelength conversion (to 550 nm based on Angstroem exponent of 1)
if (p.Results.flagConvert550)
    wlConvRatio = (550 / 1030) .^ 1;
else
    wlConvRatio = 1;
end

tExt(idxB) = RCS_B / (I_A_B / (1 - RCS_B / RCS_A) - I_A_B);
I_A_x = I_A_B;
for iBin = (idxB - 1):-1:idxMinOL
    I_A_x = I_A_x - (RCS(iBin) + RCS(iBin + 1)) / 2 * (range(iBin + 1) - range(iBin));
    tExt(iBin) = wlConvRatio * RCS(iBin) / (I_A_B / (1 - RCS_B / RCS_A) - I_A_x);
end

end