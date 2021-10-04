function [RotatedImage, OriginalAxis] = RotateNiftiWindows(Parameters)
%
% Description:
%               Function to rotate an image volume in Nifti format into
%               alignment with the scanner (world) axes and therefore align
%               the k-axis of the image frame with the real direction of
%               the main magnetic field of the MRI scanner. 
%
%               Nifti orientation information:
%               (https://brainder.org/2012/09/23/the-nifti-file-format/)
%
%               WINDOWS VERSION. Requires: Windows Subsystem for Linux (WSL)
%
% Inputs:
%   Parameters (struct): Parameters.Image(struct) - Nifti struct containing the image to be rotated and the appropriate header information
%                        Parameters.DesiredAxis(vector) - vector of the desired orientation of the B0 direction in the image frame
%                        Parameters.Reverse(boolean) - whether to rotate towards the desired axis or back the other way (if reorienting back to original)
%                        Parameters.Interpolation(int) - 'trilinear', 'nearestneighbour', 'sinc', 'spline'
%                        Parameters.FSLPath(string) - string to fsl directory (e.g. '/usr/local/fsl')
%                        Parameters.Padding(int or array) - the amount to pad in all directions (if int) or an array to pad (if array)
%                                                           to ensure no cropping of volume of interest (optional)
%       
% Outputs:
%   RotatedImage (double array) - rotated image in alignment with scanner axes
%   OriginalAxis (vector) - vector of the original B0 direction in the
%                           image frame. To be used if user wants to reorient 
%                           image in original orienation
%
% Dependencies:
%   FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSL)
%   Windows Subsystem for Linux (WSL) (https://docs.microsoft.com/en-us/windows/wsl/install-win10)
%   convert_pathname.m (converts Windows paths into Linux readable paths)
%
% Author:
%   Magnetic Resonance Imaging Group, 
%   Department of Medical Physics and Biomedical Engineering, 
%   University College London, UK, 2021

% Read input parameters
if isfield(Parameters, 'Image')
    if ~isstruct(Parameters.Image)
        warndlg('Please specify a Nifti struct for Parameters.Image!', 'Warning')
        return;
    else
        Image = Parameters.Image;
    end
else
    warndlg('Please specify Parameters.Volume!', 'Warning')
    return;
end
if isfield(Parameters, 'DesiredAxis')
    DesiredAxis = Parameters.DesiredAxis;
    if all(size(DesiredAxis) == [1 3]) || all(size(DesiredAxis)==[3 1])
        if all(size(DesiredAxis)==[1 3])
            DesiredAxis = reshape(DesiredAxis, [3 1]);
        end
    else
        warndlg('Please input Parameters.DesiredAxis in the correct format, e.g. [0;0;1]!', 'Warning')
        return;
    end
else
    DesiredAxis = [0;0;1];
end
if isfield(Parameters, 'Reverse')
    Reverse = Parameters.Reverse;
else
    Reverse = 0; % default not to reverse orientation
end
if isfield(Parameters, 'Interpolation')
    Interpolation = Parameters.Interpolation;
else
    Interpolation = 'trilinear'; % default trilinear interpolation
end
if isfield(Parameters, 'FSLPath')
   if ~isstring(Parameters.FSLPath) && ~ischar(Parameters.FSLPath)
       warndlg('Please enter path to FSL bin in string format for Parameters.FSLPath!', 'Warning')
       return;
   else
       FSLPath = char(Parameters.FSLPath);
   end
else
    warndlg('Please enter the path to the FSL bin for Parameters.FSL!', 'Warning')
    return;
end
if isfield(Parameters, 'Padding')
    if ~isa(Parameters.Padding, 'int') && ~isa(Parameters.Padding, 'double')
        warndlg('Please enter an integer value for Parameters.Padding!', 'Warning')
        return;
    else
        Padding = Parameters.Padding;
            if size(Padding) == 1
                PaddingArray = [1,1,1]*Padding; % both sides of each dimension
                if ~Reverse
                    Image.img = padarray(Image.img, PaddingArray, 0, 'both');
                    Image.hdr.dime.dim(2:4) = size(Image.img, 1:3);
                end
            elseif size(Padding, 1) == 3 || size(Padding, 2) == 3
                PaddingArray = Padding;
                if ~Reverse
                    Image.img = padarray(Image.img, PaddingArray, 0, 'both');
                    Image.hdr.dime.dim(2:4) = size(Image.img, 1:3);
                end
        end
    end
end

% Create a temp folder in current directory for the outputs
TempDir = [pwd '\temp_rotate_nifti'];
mkdir(TempDir)

% Read the qfac value
q = Image.hdr.dime.pixdim(1);
if q ~= -1 && q ~= 1
    q = 1;
end

% Method 2:
if Image.hdr.hist.qform_code > 0
    % Define the a,b,c,d quaternion parameters and form rotation matrix
    b = Image.hdr.hist.quatern_b;
    c = Image.hdr.hist.quatern_c;
    d = Image.hdr.hist.quatern_d;
    if abs((b^2 + c^2 + d^2) - 1) < 10*eps('single')
        a = 0;
    else
        a = sqrt(1 - (b^2 + c^2 + d^2));
    end

    R = [a*a+b*b-c*c-d*d, 2*b*c-2*a*d, 2*b*d+2*a*c;...
         2*b*c+2*a*d, a*a+c*c-b*b-d*d, 2*c*d-2*a*b;...
         2*b*d-2*a*c, 2*c*d+2*a*b, a*a+d*d-c*c-b*b];
     
    B0direction = R*DesiredAxis; 
    B0direction(3) = B0direction(3)/q; % the true B0 direction in the image frame
% Method 3:
elseif Image.hdr.hist.sform_code > 0
    R = [Image.hdr.hist.srow_x/Image.hdr.dime.pixdim(2);...
         Image.hdr.hist.srow_y/Image.hdr.dime.pixdim(3);...
         Image.hdr.hist.srow_z/Image.hdr.dime.pixdim(4);...
         0, 0, 0, 1];
    R = R(1:3,1:3);
    R(:,3) = R(:,3)*q;
    
    B0direction = R*DesiredAxis; 
    B0direction(3) = B0direction(3)/q; % the true B0 direction in the image frame
end

% Rotate the image such that the true B0 direction is aligned with the k-axis in the image frame

% Calculate rotation vector and angle
v_c = cross(B0direction, DesiredAxis);
v_c = v_c/norm(v_c);
v_angle = acos(dot(B0direction,DesiredAxis));

% get sine and cosine of rotation angle
s=sin(v_angle);
c=cos(v_angle);
t=1-c;

% build rotation matrix
TransOrigin = [1 0 0 -round(size(Image.img,1)/2)*Image.hdr.dime.pixdim(2);...
                0 1 0 -round(size(Image.img,2)/2)*Image.hdr.dime.pixdim(3);...
                0 0 1 -round(size(Image.img,3)/2)*Image.hdr.dime.pixdim(4);...
                0 0 0 1];
            
TransCentre = [1 0 0 round(size(Image.img,1)/2)*Image.hdr.dime.pixdim(2);...
                0 1 0 round(size(Image.img,2)/2)*Image.hdr.dime.pixdim(3);...
                0 0 1 round(size(Image.img,3)/2)*Image.hdr.dime.pixdim(4);...
                0 0 0 1];
RotMatrix = [t*v_c(1)*v_c(1) + c,        t*v_c(1)*v_c(2) - s*v_c(3), t*v_c(1)*v_c(3) + s*v_c(2), 0.0; ...
             t*v_c(1)*v_c(2) + s*v_c(3), t*v_c(2)*v_c(2) + c,        t*v_c(2)*v_c(3) - s*v_c(1), 0.0;...
             t*v_c(1)*v_c(3) - s*v_c(2), t*v_c(2)*v_c(3) + s*v_c(1), t*v_c(3)*v_c(3) + c,        0.0; ...
             0.0, 0.0, 0.0, 1.0];

if Reverse
    RotMatrix = TransCentre*RotMatrix'*TransOrigin;
else
    RotMatrix = TransCentre*RotMatrix*TransOrigin;
end

% Save everything into the temp folder
dlmwrite([TempDir '\Matrix.txt'], RotMatrix, 'delimiter',' ');

try
    save_untouch_nii(Image, [TempDir '\Image.nii']);
catch
    save_nii(Image, [TempDir '\Image.nii']);  
end

% Convert the pathnames into linux 
Matrix_path = convert_pathname([TempDir '\Matrix.txt']);
Image_path = convert_pathname([TempDir '\Image.nii']);
output = convert_pathname([TempDir '\RotatedImage.nii']);

unix_cmd = sprintf(['wsl . %s' '/etc/fslconf/fsl.sh;'...
                    'FSLOUTPUTTYPE=NIFTI; %s' '/bin/flirt -noresample -nosearch -paddingsize 100 -setbackground 0 -in %s ' ...
                    '-ref %s -applyxfm -init %s -interp %s -out %s'],...
                    FSLPath, FSLPath, Image_path, Image_path,...
                    Matrix_path, Interpolation, output); 

system(unix_cmd);

% Load in the image and delete temp
try
    RotatedImage = load_untouch_nii([TempDir '\RotatedImage.nii']);
catch
    RotatedImage = load_nii([TempDir '\RotatedImage.nii']);
end

pause(5);
rmdir(TempDir, 's');

% Unpad the images if reversed back to original orientation and if padded
if Reverse
    if isfield(Parameters, 'Padding')
        CropWdw = centerCropWindow3d(size(RotatedImage.img), size(RotatedImage.img) - 2*PaddingArray);
        RotatedImage.img = imcrop3(RotatedImage.img, CropWdw);
        RotatedImage.hdr.dime.dim(2:4) = size(RotatedImage.img, 1:3);
    end
end

OriginalAxis = B0direction;


end