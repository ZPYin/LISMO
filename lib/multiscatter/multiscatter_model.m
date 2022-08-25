function data = multiscatter_model(filename, varargin)
%MULTISCATTER_MODEL running the mulple scattering model to calculate the 
%multiple scattering factor. Detailed information you can find in 
%'../lib/multiscatter.m'
% Inputs: 
%    filename: char
%        the absolute path for the lidar configuration file. You can find 
%        an detailed information about the format and the multiple 
%        scattering model [here](http://www.met.reading.ac.uk/clouds/multiscatter/)
% Outputs:
%    data: struct
%        range: array
%            apparent range above the ground. [m]
%        cloudExt: array
%            cloud extintion coefficient. [m^{-1}]
%        cloudRadius: array
%            cloud effective mean radius. [microns]
%        att_total: array
%            total attenuated backscatter. [m^{-1}*Sr^{-1}]
%        att_single: array
%            attenuated backscatter with single backscattering. 
%            [m^{-1}*Sr^{-1}]
% Note:
%    Make sure you've compiled the multiple scattering model and added it to 
%    the searching path before you run the function.
% History:
%    2019-08-14. Add comments by Zhenping Yin
% Contact: 
%    zhenping@tropos.de

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'filename', @ischar);
addParameter(p, 'cygwinPath', '', @ischar);

parse(p, filename, varargin{:});

if exist(filename, 'file') ~= 2
    error('The configuration file does not exist.');
end

if isempty(p.Results.cygwinPath)
    warning('Cygwin does not exist.');
    cygwinPath = 'C:\\cygwin64\\bin\\bash';
else
    cygwinPath = p.Results.cygwinPath;
end

output_file1 = [tempname '.txt'];
output_file2 = [tempname '.txt'];

[s1, ~] = system([cygwinPath, ' --login -c ', '"', strrep([fullfile(fileparts(fileparts(fileparts(which('multiscatter_model')))), ...
'include', 'multiscatter.exe') ' -algorithms original lidar < ' ...
filename ' > ' output_file1], '\', '/'), '"']);
[s2, ~] = system([cygwinPath, ' --login -c ', '"', strrep([fullfile(fileparts(fileparts(fileparts(which('multiscatter_model')))), ...
'include', 'multiscatter.exe') ' -algorithms single lidar < ' ...
filename ' > ' output_file2], '\', '/'), '"']);

if ~ (s1 == 0) || ~ (s2 == 0)
    error('error in excuting multiscatter program.');
end

fid1 = fopen(output_file1, 'r');
data1 = fscanf(fid1, '%f %f %f %f %f', [5, Inf]);
data1 = data1';
fclose(fid1);
fid2 = fopen(output_file2, 'r');
data2 = fscanf(fid2, '%f %f %f %f %f', [5, Inf]);
data2 = data2';
fclose(fid2);

delete(output_file1);
delete(output_file2);

data.range = data1(:, 2);
data.cloudExt = data1(:, 3);
data.cloudRadius = data1(:, 4);
data.att_total = data1(:, 5);
data.att_single = data2(:, 5);

end