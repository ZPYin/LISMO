function [aExt] = extRet_Xian(range, signal, bg, varargin)
% EXTRET_XIAN calculate extinction coefficient based on Xian's method.
%
% USAGE:
%    [aExt] = extRet_Xian(range, rcs, varargin)
%
% INPUTS:
%    range: numeric
%        distance in range bins. (m)
%    signal: numeric
%        backscatter signal.
%    bg: numeric
%        background
%
% KEYWORDS:
%    rangeFullOverlap: numeric
%        range of full overlap achieved. (m)
%    minSNR: numeric
%        minimum SNR for signal retrieval.
%
% OUTPUTS:
%    aExt: numeric
%        extinction coefficient. (m^-1)
%
% HISTORY:
%    2023-03-07: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'range', @isnumeric);
addRequired(p, 'signal', @isnumeric);
addRequired(p, 'bg', @isnumeric);
addParameter(p, 'rangeFullOverlap', 500.0, @isnumeric);
addParameter(p, 'minSNR', 0.5, @isnumeric);

parse(p, range, signal, bg, varargin{:});

[~, idxMinOL] = min(abs(range - p.Results.rangeFullOverlap));
SNR = lidarSNR(signal, bg, 'bgBins', [length(signal) - 100, length(signal) - 50]);
RCS = signal .* range.^2;
RCS_A = RCS(idxMinOL);
ratioAB = RCS / RCS_A;
isInvalid = (range < p.Results.rangeFullOverlap) | (SNR <= p.Results.minSNR);
ratioAB(isInvalid | (ratioAB <= 0)) = NaN;
[~, idxB] = min(ratioAB);
RCS_B = RCS(idxB);

aExt = NaN(size(RCS));
if (idxB ~= idxMinOL)
    I_A_B = trapz(range(idxMinOL:idxB), RCS(idxMinOL:idxB));
else
    I_A_B = 0;
end

for iBin = (idxMinOL + 1):idxB
    aExt(iBin) = RCS(iBin) / (I_A_B / (1 - RCS_B / RCS_A) - trapz(range(idxMinOL:iBin), RCS(idxMinOL:iBin)));
end

end