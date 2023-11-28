function [quasi_par_bsc, quasi_par_ext] = quasiRetrieval(height, att_beta, ...
                    molExt, molBsc, LRaer, varargin)
% QUASIRETRIEVAL Retrieve the aerosol optical properties with quasi retrieving method
%
% USAGE:
%    [quasi_par_bsc, quasi_par_ext] = quasiRetrieval(height, att_beta, ...
%                                        molExt, molBsc, LRaer)
%
% INPUTS:
%    height: array
%        height. [m] 
%    att_beta: matrix (height x time)
%        attenuated backscatter. [m^{-1}Sr^{-1}] 
%    molExt: matrix (height x time)
%        molecule extinction coefficient. [m^{-1}] 
%    molBsc: matrix (height x time)
%        molecule backscatter coefficient. [m^{-1}Sr^{-1}]
%    LRaer: float | matrix (height x time)
%        aerosol lidar ratio. [Sr]
%
% KEYWORDS:
%    nIters: numeric
%        iteration times. (default: 2)
%    flagAutoConverge: logical
%        whether to determine iteration numbers based on convergence. (default: false)
%    devHRange: numeric
%        height range for calculating deviation. (m)
%
% OUTPUTS:
%    quasi_par_bsc: matrix
%        quasi particle backscatter coefficient. [m^{-1}Sr^{-1}] 
%    quasi_par_ext: matrix
%        quasi particle extinction coefficient. [m^{-1}]
%
% REFERENCES:
%    Baars, H., Seifert, P., Engelmann, R. & Wandinger, U. Target categorization of aerosol and clouds by continuous multiwavelength-polarization lidar measurements. Atmospheric Measurement Techniques 10, 3175-3201, doi:10.5194/amt-10-3175-2017 (2017).
%
% HISTORY:
%    - 2018-12-25: First Edition by Zhenping
%    - 2019-03-31: Add the keywork of 'nIters' to control the iteration times.
%    - 2023-11-28: Add auto converge.
%
% .. Authors: - zhenping@tropos.de

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'height', @isnumeric);
addRequired(p, 'att_beta', @isnumeric);
addRequired(p, 'molExt', @isnumeric);
addRequired(p, 'molBsc', @isnumeric);
addRequired(p, 'LRaer', @isnumeric);
addParameter(p, 'nIters', 2, @isnumeric);
addParameter(p, 'flagAutoConverge', false, @islogical);
addParameter(p, 'devHRange', [0, 10000], @isnumeric);

parse(p, height, att_beta, molExt, molBsc, LRaer, varargin{:});

height = reshape(height, 1, []);
if length(LRaer) == 1
    LRaer = LRaer * ones(size(att_beta));
end

diffHeight = repmat(transpose([height(1), diff(height)]), 1, size(att_beta, 2));
mol_att = exp(- nancumsum(molExt .* diffHeight, 1));
quasi_par_ext = zeros(size(molBsc));
quasi_par_bsc = zeros(size(molBsc));
quasi_par_att = zeros(size(molBsc));

if p.Results.flagAutoConverge
    MAXITER = 20;

    % iteration untill convergence
    isInDevHeight = (height > p.Results.devHRange(1)) & (height < p.Results.devHRange(2));

    for iProf = 1:size(att_beta, 2)
        isConverge = false;
        isMaxIter = false;
        iterCounter = 1;

        while((~ isConverge) && (~ isMaxIter))

            quasiBsc = quasi_par_bsc(:, iProf);
            quasi_par_att(:, iProf) = exp(-nancumsum(quasi_par_ext(:, iProf) .* diffHeight(:, iProf)));
            quasi_par_bsc(:, iProf) = att_beta(:, iProf) ./ (mol_att(:, iProf) .* quasi_par_att(:, iProf)).^2 - molBsc(:, iProf);
            quasi_par_bsc(quasi_par_bsc(:, iProf) <= 0, iProf) = 0;
            quasi_par_ext(:, iProf) = quasi_par_bsc(:, iProf) .* LRaer(:, iProf);

            sumBsc = nansum(abs(quasi_par_bsc(isInDevHeight, iProf)));
            sumDev = nansum(abs(quasi_par_bsc(isInDevHeight, iProf) - quasiBsc(isInDevHeight)));
            isConverge = (sumDev / sumBsc) < 0.001;
            isMaxIter = (iterCounter > MAXITER);

            iterCounter = iterCounter + 1;

        end
    end
else
    for iLoop = 1:p.Results.nIters
        quasi_par_att = exp(-nancumsum(quasi_par_ext .* diffHeight, 1));
        quasi_par_bsc = att_beta ./ (mol_att .* quasi_par_att).^2 - molBsc;
        quasi_par_bsc(quasi_par_bsc <= 0) = 0;
        quasi_par_ext = quasi_par_bsc .* LRaer;
    end
end

end