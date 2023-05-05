function [snr0] = lidarSNR(sig, bg, varargin)
% LIDARSNR lidar signal SNR
% USAGE:
%    [snr0] = lidarSNR(sig, bg)
% INPUTS:
%    sig: array
%    bg: array | numeric
% KEYWORDS:
%    bgBins: 2-element array
% OUTPUTS:
%    snr0: array
% HISTORY:
%    2021-09-21: first edition by Zhenping
% .. Authors: - zhenping@tropos.de

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'sig', @isnumeric);
addRequired(p, 'bg', @isnumeric);
addParameter(p, 'bgBins', [], @isnumeric);

parse(p, sig, bg, varargin{:});

if isempty(p.Results.bgBins)
    bgBins = [length(sig) - 50, length(sig)];
else
    bgBins = p.Results.bgBins;
end

snr0 = NaN(size(sig));
if bg > 0
    ratio = nanstd(sig(bgBins(1):bgBins(2))) / sqrt(bg);

    flagZero = ((sig + bg) <= 0) | (sig <= 0);
    snr0(~ flagZero) = sig(~ flagZero) ./ (ratio .* sqrt(sig(~ flagZero) + bg));
    snr0(flagZero) = 0;
else
    tot = sig + bg;
    snr0 = sig ./ ADSigStd(tot, 5);
    % snr0 = sig ./ nanstd(sig((end - 50):end));
end

snr0(isnan(snr0)) = 0;

end