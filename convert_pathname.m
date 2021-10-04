function linux_pathname = convert_pathname(windows_pathname)
% function to convert the windows pathname into a linux pathname that works
% in the linux subsystem for windows 
%
% INPUT: 
% windows_pathname = the pathname to a file/directory in windows syntax
%
% OUTPUT:
% linux_pathname = the pathname converted to linux with /mnt/c/ 
%
% Author: Oliver Kiersnowski
% Date: November 2019

[filepath, name, ext] = fileparts(windows_pathname);

drive = filepath(1);
filepath = strrep(filepath, [drive ':\'], ['/mnt/' lower(drive) '/']);

% replace all of the backslashes with forward slashes

filepath = strrep(filepath, '\', '/');

linux_pathname = [filepath '/' name ext];

end