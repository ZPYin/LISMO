function [signalCor] = overlapCor(range, signal, ol, varargin)
% OVERLAPCOR overlap correction.
%
% USAGE:
%    [output] = overlapCor(range, signal, ol, varargin)
%
% INPUTS:
%    range: numeric
%    signal: numeric (height x time)
%    ol: numeric
%        overlap function
%
% KEYWORDS:
%    glueRange: numeric
%
% OUTPUTS:
%    signalCor: numeric
%
% HISTORY:
%    2023-11-23: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'range', @isnumeric);
addRequired(p, 'signal', @isnumeric);
addRequired(p, 'ol', @isnumeric);
addParameter(p, 'glueRange', [1000, 1200], @isnumeric);

parse(p, range, signal, ol, varargin{:});

weightH = (range >= p.Results.glueRange(2)) + ...
          ((range >= p.Results.glueRange(1)) & (range < p.Results.glueRange(1))) .* ((range - p.Results.glueRange(1)) ./ (p.Results.glueRange(2) - p.Results.glueRange(1)));
weightL = 1 - weightH;

sigFirstCor = signal ./ repmat(reshape(ol, size(signal, 1), 1), 1, size(signal, 2));

signalCor = repmat(reshape(weightL, size(signal, 1), 1), 1, size(sigFirstCor, 2)) .* (sigFirstCor) + repmat(reshape(weightH, size(signal, 1), 1), 1, size(sigFirstCor, 2)) .* signal;

end