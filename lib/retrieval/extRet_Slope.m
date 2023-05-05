function [tExt, tExtStd] = extRet_Slope(range, sig, varargin)
% EXTRET_SLOPE Retrieve total extinction coefficient based on Slope method.
%
% USAGE:
%    [tExt, tExtStd] = extRet_Slope(range, sig)
%
% INPUTS:
%    range: numeric
%        range array. (m)
%    sig: numeric
%        background removed signal (photon count)
%
% OUTPUTS:
%    tExt: numeric
%        total extinction coefficient. (m-1)
%    tExtStd: numeric
%        uncertainty of total extinction coefficient. (m-1)
%
% HISTORY:
%    2023-04-26: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'range', @isnumeric);
addRequired(p, 'sig', @isnumeric);

parse(p, range, sig, varargin{:});

%% Calculate Range-Corrected Signal
rcs = sig .* range.^2;
sigNonNegtive = sig;
sigNonNegtive(sigNonNegtive < 0) = NaN;
rcsStd = sqrt(sigNonNegtive) .* range.^2;

%% Remove Non-Positive Signal (Filled with NaN)
rcs(rcs <= 0) = NaN;

[~, slope, ~, slopeStd] = chi2fit(range, log(rcs), rcsStd);

tExt = slope / -2;
tExtStd = slopeStd / 2;

end