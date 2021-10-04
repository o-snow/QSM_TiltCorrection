# QSM_TiltCorrection
## Description
Software to carry out alignment of the image axes with the MRI scanner axes, such that the image is aligned with the main magnetic field. For QSM processing we recommend rotating the total field map to the scanner axes prior to any background field removal to ensure accurate background field removal and susceptibility calculation. 

Folder:      `<Rotation_Software>`

Functions:   `<RotateNiftiWindows.m>`
             `<RotateNiftiLinux.m>`
             `<convert_pathname.m>`

For Windows users please use the Windows version and for Linux/Mac users please use the Linux version. 

## Requirements

Windows version requires Windows Subsystem for Linux installed (https://docs.microsoft.com/en-us/windows/wsl/install-win10
All versions require FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSL)

## Usage

Parameters struct is used for this function. The struct must have the following fields:

```
Parameters.Image(struct) - Nifty struct containing the image to be rotated and the appropriate header information

Parameters.DesiredAxis(vector) - vector of the desired orientation of the B0 direction in the image frame

Parameters.Reverse(boolean) - whether to rotate towards the desired axis or back the other way (if reorienting back to original)

Parameters.Interpolation(int) - 'trilinear', 'nearestneighbour', 'sinc', 'spline'

Parameters.FSLPath(string) - string to fsl directory (e.g. '/usr/local/fsl')

Parameters.Padding(int or array) - the amount to pad in all directions (if int) or an array to pad (if array) to ensure no cropping of volume of interest (optional)

```

*Recommended use:*

Rotate the total field map prior to PDF background field removal. Carry out the susceptibility calculation in the same rotated frame.

If it is desired for the final QSM to be in the same orientation as the acquired data (for ease in segmentation e.g.) then, keeping all other Parameters the same, change the Parameters.Input to be the QSM, and add an additional field Parameters.Reverse = true. Setting this as true will calculate the rotation matrix as before, but apply the inverse of it, rotating the QSM in the opposite direction to the total field map.

For ease in use, we recommend using a Parameters struct that stays constant throughout the pipeline. Therefore name it differently to other Parameter structs in your pipeline. E.g. ParametersForRot.Input = FieldMap; and so on... This allows reversal of the rotation without having to define all of the initial fields again, if you have had to clear the Parameters struct. 

