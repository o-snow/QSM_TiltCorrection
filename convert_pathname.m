function out_pathname = convert_pathname(in_pathname, flags)
% function to convert the windows pathname into a linux pathname that works
% in the linux subsystem for windows 
%
% INPUT: 
% windows_pathname = the pathname to a file\directory in windows syntax
%
% OUTPUT:
% linux_pathname = the pathname converted to linux with /mnt/c/ 
%
% Author: Oliver Kiersnowski
% Date: November 2019

if nargin < 2
    flags = '-u'; % default from Windows to WSL path
end

% [filepath, name, ext] = fileparts(in_pathname);
% 
% drive = filepath(1);
% filepath = strrep(filepath, [drive ':\'], ['/mnt/' lower(drive) '/']);
% 
% % replace all of the backslashes with forward slashes
% 
% filepath = strrep(filepath, '\', '/');
% 
% out_pathname = [filepath '/' name ext];

convert_command = sprintf('wsl wslpath %s "%s"',flags, in_pathname);
[status, out_pathname] = system(convert_command);

if status~=0
    error('path conversion failed with result:\n %s',out_pathname);
end

out_pathname = deblank(out_pathname); % Fix this neater?

end
