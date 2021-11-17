# QSM_TiltCorrection

From 2020, this toolbox is being mantained and expanded at the Department of Medical Physics and Biomedical Engineering, University College London, UK, by Oliver Kiersnowski under the supervision of Prof. Karin Shmueli.

## Description
Software to carry out alignment of the image axes with the MRI scanner axes, such that the image is aligned with the main magnetic field. For QSM processing we recommend rotating the total field map to the scanner axes prior to any background field removal to ensure accurate background field removal and susceptibility calculation [1].

Folder:      `<Rotation_Software>`

Functions:   `<RotateNiftiWindows.m>`
             `<RotateNiftiLinux.m>`
             `<convert_pathname.m>`

For Windows users please use the Windows version and for Linux/Mac users please use the Linux version. 

## Requirements

Windows version requires Windows Subsystem for Linux installed (https://docs.microsoft.com/en-us/windows/wsl/install-win10
All versions require FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSL) [2,3,4]

Tools for NIfTI and ANALYZE image (https://uk.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image) for inputting the correct nifti structure into the `Parameters` struct.

## Usage

Parameters struct is used for this function. The struct must have the following fields:

```
Parameters.Image(struct) - Nifti struct containing the image to be rotated and the appropriate header information

Parameters.DesiredAxis(vector) - vector of the desired orientation of the B0 direction in the image frame

Parameters.Reverse(boolean) - whether to rotate towards the desired axis or back the other way (if reorienting back to original)

Parameters.Interpolation(int) - 'trilinear', 'nearestneighbour', 'sinc', 'spline'

Parameters.FSLPath(string) - string to fsl directory (e.g. '/usr/local/fsl')

Parameters.Padding(int or array) - the amount to pad in all directions (if int) or an array to pad (if array) to ensure no cropping of volume of interest (optional)

```
The inputted Parameters.Image is a Nifti struct containing the fields: Image.hdr that contains the header information and Image.img that contains the image volume. We recommend reading the image volume from file using  `load_untouch_nii.m` from the Tools for NIfTI and ANALYZE image (https://uk.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image).

*Recommended use:*

Rotate the total field map prior to PDF background field removal. Carry out the susceptibility calculation in the same rotated frame.

If it is desired for the final QSM to be in the same orientation as the acquired data (for ease in segmentation e.g.) then, keeping all other Parameters the same, change the Parameters.Input to be the QSM, and add an additional field Parameters.Reverse = true. Setting this as true will calculate the rotation matrix as before, but apply the inverse of it, rotating the QSM in the opposite direction to the total field map.

For ease in use, we recommend using a Parameters struct that stays constant throughout the pipeline. Therefore name it differently to other Parameter structs in your pipeline. E.g. ParametersForRot.Input = FieldMap; and so on... This allows reversal of the rotation without having to define all of the initial fields again, if you have had to clear the Parameters struct. 

# References
[1]: Kiersnowski OC, Karsa A, Thornton JS, Shmueli K. The Effect of Oblique Image Slices on the Accuracy of Quantitative Susceptibility Mapping and a Robust Tilt Correction Method #0794. Proc. Int. Soc. Magn. Reson. Med. 2021 doi: 10.1002/mrm.22135.

[2]: M.W. Woolrich, S. Jbabdi, B. Patenaude, M. Chappell, S. Makni, T. Behrens, C. Beckmann, M. Jenkinson, S.M. Smith. Bayesian analysis of neuroimaging data in FSL. NeuroImage, 45:S173-86, 2009

[3]: S.M. Smith, M. Jenkinson, M.W. Woolrich, C.F. Beckmann, T.E.J. Behrens, H. Johansen-Berg, P.R. Bannister, M. De Luca, I. Drobnjak, D.E. Flitney, R. Niazy, J. Saunders, J. Vickers, Y. Zhang, N. De Stefano, J.M. Brady, and P.M. Matthews. Advances in functional and structural MR image analysis and implementation as FSL. NeuroImage, 23(S1):208-19, 2004

[4]: M. Jenkinson, C.F. Beckmann, T.E. Behrens, M.W. Woolrich, S.M. Smith. FSL. NeuroImage, 62:782-90, 2012

# Disclaimer
THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL OLIVER KIERSNOWSKI OR HIS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, BUSINESS INTERRUPTION; PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; AND LOSS OF USE, DATA OR PROFITS) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
