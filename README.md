# OpenCapsule - Pendant Capsule Elastometry
C/C++ Software / Command Line Interface / Linux  
Developed at TU Dortmund University  
Author: Jonas Hegemann (jonas.hegemann@tu-dortmund.de)  
Framework: Sebastian Knoche (sebastian.knoche@tu-dortmund.de)  
Supervisor: Jan Kierfeld (jan.kierfeld@tu-dortmund.de)  

## Description:
We provide a C/C++ software for the shape analysis of 
deflated elastic capsules in a pendant capsule geometry, 
which is based on the proper elastic description of the 
capsule material as a quasi two-dimensional elastic layer 
using shell theory. This pendant capsule elastometry 
provides a new tool for interfacial rheology of capsules 
that goes beyond determination of the Gibbs- or dilatational 
modulus from area-dependent measurements of the surface 
tension using pendant drop tensiometry, which can only 
give a rough estimate of the elastic capsule properties 
as they are based on a purely liquid interface model. 
The software allows to determine the following quantities 
for each digitized image of a deflated capsule shape 
given another image of its undeformed reference shape: 
the surface tension characterizing the reference shape, 
Young's surface modulus (or alternatively area compression 
modulus) and Poisson's ratio. If series of images are 
available, these moduli can be determined as a function 
of the capsule volume to analyze viscoelastic effects 
for different volume change rates and aging effects 
depending on the deformation history. An additional 
wrinkling wavelength measurement allows to determine 
the bending modulus, from which the layer thickness 
can be derived. We verified the method by analyzing 
several materials and comparing the results to available 
rheological measurements. We make the software 
available under the GPL license. 

## Publications

* Pendant capsule elastometry  
  J. Hegemann, S. Knoche, S. Egger, H. Rehage, and J. Kierfeld  
  Journal of Colloid and Interface Science, ... to be published

* Elastometry of deflated capsules: Elastic moduli from shape and wrinkle analysis  
  S. Knoche, D. Vella, E. Aumaitre, P. Degen, H. Rehage, P. Cicuta, and J. Kierfeld  
  Langmuir, Vol. 29, No. 40, Pages 12463--12471, 2013, ACS Publications

## Installation Guide

OpenCapsule was tested (successful compilation and execution) on 
Ubuntu 16.04 LTS and Scientific Linux 6.6 with GNU/GCC compiler. 
To install OpenCapsule first download the tarball and extract 
it to some arbitrary direction or simply clone the repository. 
Before compiling install the required dependencies
Boost, OpenCV and Gnu Scientific Library (GSL) e.g. 
on Ubuntu typing something like

```
sudo apt-get install gcc g++ make libboost-all-dev libopencv-dev libgsl-dev libgsl2
```

in a login-shell will be sufficient. On Scientific Linux, 
i.e. Fedora system, use Yum package manager or install 
the libraries manually. OpenCapsule uses OpenCV for image 
handling and therefore only supports the image formats 
supported by your specific installation. Check if all PATH 
variables are set appropriately or adapt the Makefile by 
informing g++ explicitly where to find the libraries by 
using the -L option. Finally, simply type

```
make
```

to compile the project. If you like to 
set up a more usable working environment also type

```
make install
```

to make OpenCapsule accessible in all working directories.

## Image requirements

To ensure images are compatible with OpenCapsule 
and can be analyzed correctly, they have to 
meet several requirements, which are given in the following. 

* The needle should be visible in the upper part of the image.
* As far as possible, the capsule should to be aligned along the vertical axis.
* Gravitational forces should act along the vertical axis.
* The background should be equally colored with no other objects or artefacts in there troubling contour detection.
* The camera should not be moved during recording of the sequence.
* Inner and outer capillary diameter should differ significantly.

If you are not sure what your images should look like,
take a look at the publications listed above.

## Getting Started

To give an impression of the typical command line 
workflow some examples are given in the following. 
A typical task would be to perform a Laplace-Young 
analysis, which is nowadays the standard algorithm 
for analyzing droplets and capsules. Based on this, 
OpenCapsule provides an additional algorithm based 
on Hookean elasticity, which is capable to describe 
shapes of elastic shells more accurately than the 
Laplace-Young analysis. This model exhibits much 
more information, such as the elastic moduli, stress 
distributions and details about the wrinkling. 
OpenCapsule differs from other implementations 
by a high degree of automation and numerical 
robustness. Change to your prefered working 
directory and call OpenCapsule by simply typing
```
OpenCapsule
```
in the shell. If your working environment is set 
up correctly, OpenCapsule gives you some general 
information about itself. Additionally, all 
necessary directories and a standard configuration 
file are created. Please edit the configuration 
file `./config/config.cfg` before you proceed. 
List your input files, which should placed 
(by default) in the `./input/` folder just created, 
under `REFERENCE_SHAPE` and `ELASTIC_SHAPE`. 
OpenCapsule averages over reference images, 
which should all show the same, undeformed 
state of the capsule. Make sure you have 
specified the fields `EXPERIMENT_CAPDIAMETER`
and `EXPERIMENT_DENSITY` in SI units. If you 
are not sure how to handle the configuration 
file, see the corresponding section for 
more details about this. Next, typing
```
OpenCapsule -h
```
will give you an overview over the command line 
options as well as a list of common examples. 
The results of the following standard commands 
are placed in the `./global_out/` folder, where 
you will find the results listed in plain 
`*.txt` files, images in `*.png` format, and 
wrapped in a comprehensive html-report. 

### Edge Detection

```
OpenCapsule -i [path to image]
```
This command performs an isolated image analysis 
of the image specified by the relative 
path `./path-to-image/image.png`. 
If your image processing fails or results are 
obviously bad, you can adapt OpenCapsule to your 
specific images. Try to tune the related values 
in the configuration file to get an optimal 
result. If everything works fine, perform the shape analysis. 

### Laplace-Young Analysis
```
OpenCapsule -r
```
This command triggers the standard Laplace-Young analysis. 
No initial values have to be specified. 

### Elastic/Hookean Analysis
```
OpenCapsule -s
```
This command triggers the Hookean shape analysis 
and will automatically use the specified reference 
shapes. _This sets OpenCapsule apart from
available commercial software packages._
You can also provide initial values by using the call
```
OpenCapsule -s --poisson 0.5 --compression 10.0 --pressure 1.0
```
These three parameters are optional. The first 
one specifies an initial guess of Poisson's 
ratio. The second one specifies an initial 
guess of the area compression modulus in units 
of the surface tension obtained by the reference 
regression. The third one specifies an initial 
guess for the dimensionless pressure at the apex 
of the deflated capsule in units of the surface 
tension and the inverse capillary diameter. Only 
if the simple call without initial values gives 
no adequate result, you should provide intitial 
values for the fit parameters. 

## Generated Output
The essential results are saved to the 
`GLOBAL_OUT_FOLDER`. Most interesting are the 
summarized results written to the files 
`reference.dat` and `sequence.dat`. Moreover 
a lot of miscellaneous output is generated 
including the determined shapes, extracted 
contours, binary images and general image 
information, i.e. quantities measured from 
the binary image. The most instructive output 
is the input image with the determined shape 
and wrinkle domain overlayed. From these 
images you can immediately judge if the 
analysis was successful. With `GNUPLOT_SUPPORT` 
activated you also get all relevant plots 
as `*.eps` files. Don't worry about the large 
number of output files. They are named 
intuitively and self explanatory. The file 
`report.html` merges a lot of information 
into a single file, including all determined 
values in SI units and images with shapes overlayed.

## License
OpenCapsule is GPLv3 licensed and thereby 100% 
Open Source Software. Feel free to use it 
in-house in your research department or company. 
Feel free to improve or adapt OpenCapsule, 
but notice that the resulting code has to be 
GPLv3 licensed, too. You are NOT allowed to 
include the OpenCapsule code or a derivate in 
proprietary or closed source applications. 
You are also NOT allowed to compile OpenCapsule 
into a library and link your proprietary or 
closed source software against this library. 
However, it is legitimate to access the OpenCapsule 
binary (as a black box) from your proprietary 
software. For details see the original GPLv3 license. 

## Frequent Issues

### OpenCapsule does not detect a closed contour...

Most probably it helps to increase the variable 
`GAUSSIAN_SIGMA` in 1.0 steps to smooth the image. 
Edge detection strongly depends on the contrast 
and brightness of your image. Sometimes it helps 
to adjust the thresholds `T_HIGH` and `T_LOW` for 
the hysteresis alogrithm performed as the last 
step of edge detection. Try some other (probably higher) 
values most preferable at ratios 3:1, 2:1 
or something similar. You can switch off the option 
`CHECK_IF_CLOSED`, but this is not recommended. Note 
that both thresholds are in the range between zero 
and one. Very small values force the edge detection 
to be very sensitive whereas values close to one 
will find only very few contour points.

### Some of the numerical algorithms don't converge...

First, turn on the corresponding `WATCH_XXX` variable 
to get an impression of what happens. Then, carefully 
tune the corresponding numerical thresholds, e.g. 
lower the required accuracies `EPS_XXX`. Also think 
about tuning your initial values and make sure 
`NELDERMEAD_PREFITTING` is enabled.

### OpenCapsule simply aborts before doing anything...

This is probably an issue of your selected paths 
or input files. Make sure all filenames are valid 
and backed by a real file. Check if your 
configuration file contains all necessary values. 
Consider a rollback to the standard configuration file.

## Upcoming Changes
In addition to Hookean elasticity, we plan to implement 
more elastic material laws, such as Mooney-Rivlin 
elasticity. This could help to decide which material 
law fits best for a specific capsule type. We also 
aim to implement extensions in order to include 
viscoelastic effects. Having analyzed plenty of materials 
in the future, it would be useful to have a large
database collecting all the different materials/capsules 
similar to already existing protein- or gen-databases.
This would help to classify materials based on their elastic properties. 

## Configuration File
All settings are controlled by the configuration 
file `./config/config.cfg` relative to the current 
working directory. Settings can be divided in 
several major classes, which are explained 
in the following. NOTE: only the sections 
**1)**+**2)** and maybe **4)**+**5)** are important 
for the average user. These sections are marked 
as "user-space". Sections, which must not be left 
unspecified are marked as "mandatory". Values 
given behind the properties are standard values. 
Reasonable ranges are given in bracktes.

### 1) Input files (user-space, mandatory)

Put the file names of your input images in a list 
seperated by colons. For the undeformed reference 
shapes give as many (similar) images as possible, 
because OpenCapsule will average over all these images. 
No specific ordering is required for the reference 
shapes. For the deformed (deflated) shapes make sure 
they are ordered chronologically, such that parameters 
can be traced and convergence will be more stable 
and faster. If some deflated shapes are very close 
or similar to the reference shape, please discard them.
```
REFERENCE_SHAPE reference1.png;reference2.png
```
Unordered list of files (in input-folder) for reference images separated by colons.
```
ELASTIC_SHAPE elastic1.png;elastic2.png;elastic3.png;elastic4.png
```
(Chronlogically ordered) list of files (in input-folder) for elastic images separated by colons.

### 2) Experimental quantities (user-space, mandatory)

OpenCapsule won't run or give any result if you don't 
give the outer capillary diameter in m and the density 
difference between inner and outer phase in kg/m^3.

```
EXPERIMENT_DENSITY (100.0 .. 900.0)
```
Density difference between inner and outer phase in kg/m^3.
```
EXPERIMENT_CAPDIAMETER (0.0005 .. 0.002)
```
Outer diameter of the used needle in m.

### 3) Numerical thresholds/constants
Normally, OpenCapsule should run smoothly with its default 
values. However, to achieve convergence of the fitting 
or shooting algorithms, it can be necessary 
(in some rare cases) to adapt the numerics to your 
specific input data. If you dare/need to change 
some of the numerical constants, you should exactly 
know what you are doing and be very careful, 
because values are tested/optimized intensively 
and small changes can dramatically change the behaviour of the algorithms.
```
EPS_RMS 1e-16 (1.0e-16 .. 1.0e-12)
```
Maximum arc-length deviation in bisection used 
to determine shortest distance between contour points and shape.
```
EPS_NEWTON_LAPLACE 1e-6, EPS_NEWTON_HOOKE 1e-6 (1.0e-8 .. 1.0e-3)
```
Newton minimization used for regression reference/elastic 
equations stops if euclidian norm of parameter update falls below this threshold.
```
EPS_SINGLE_SHOOTING 1e-6 (1.0e-8 .. 1.0e-4)
```
Shooting method used to solve elastic equations 
stops if deviation to boundary falls below 
this threshold. Deviation of 0.01 corresponds to 1%.
```
EPS_MULTI_SHOOTING 1e-6 (1.0e-8 .. 1.0e-4)
```
Shooting method used to improve single shooting 
stops if euclidian norm of residual falls 
below this threshold. Residual aggregates 
all distance vectors between individual segments of solution.
```
DX_FIT 1e-6 (1.0e-8 .. 1.0e-4)
```
Interval width used to determine numerical 
derivatives during regression, i.e. 
change of shape with parameters (moduli, pressure, ..).
```
DX_SHOOTING 1e-6 (1.0e-8 .. 1.0e-4)
```
Interval width used to determine numerical 
derivatives in multiple shooting method, 
i.e., change of segment end-points with initial values.
```
INTEGRATION_STEP_LAPLACE 1e-4 (1.0e-4 .. 1.0e-3)
```
Step width of RK4 method to solve reference equations.
```
INTEGRATION_STEP_HOOKE 1e-4 (1.0e-4 .. 1.0e-3)
```
Step width of RK4 method to solve elastic equations.
```
EPS_IMPLICIT_RK 1e-10
```
Implicit RK4 integration used to solve elastic 
equations stops if euclidian norm of 
slope update falls below this threshold.
```
IMPLICIT_INTEGRATION 0
```
Activate/deactivate implicit RK4 integration
(obsolete, slow, introduced for testing of stability).
```
EPS_NELDER_MEAD 1e-4
```
Nelder-Mead (downhill simplex) method stops if error 
difference between worst and best vertex falls below this threshold.

### 4) Image processing (user-space)
If your image processing fails, you might want to tune some of the following values.
```
SCALE_IMAGE 1e-2
```
Extracted contour points from image are simply 
multiplied by this value to do some prescaling 
before fitting reference shapes.
```
T_LOW 0.3, T_HIGH 0.6 (0.0 .. 1.0)
```
Lower and higher threshold of hysteresis 
algorithm performed as last step of edge 
detection. you can safely use higher values, 
OpenCapsule will automatically lower 
the thresholds if no closed contour is found.
```
GAUSSIAN_SIGMA 1.0 (1.0 .. 10.0)
```
Width of gaussian smoothing image before edge detection.
```
R_MIN 1, R_MAX 5 (1 .. 10)
```
Minimum and maximum radius of applied 
gaussian filters (OpenCapsule averages 
over binary images resulting from 
smoothing with differently sized gaussians).
```
THRESHOLD_CAPILLARY 5 (1 .. 20)
```
Used for tracing outer capillary from upper 
image border to inner capillary. searching 
terminates if contour pixel deviates more 
than this threshold from average horizontal position.
```
TOP_BUFFER 1 (0 .. 20)
```
Vertical distance from detected capillary 
height to nearest contour point (necessary 
because of possible surfactant/polymer 
sediments at the inner capillary).
```
CHECK_IF_CLOSED 1
```
Activates/deactivate check for closed 
contour in binary image (keep activated if possible, 
deactivation will skip abortion if 
contour not closed, behaviour will be undefined).
```
FIX_CHARACTERISTICS 1
```
Capillary position will be detected only from 
reference state, elastic states then use the 
same position. this option is useful if capsule 
shapes occur, which are parallel to the needle 
at the capsule mounting point. You can keep 
`FIX_CHARACTERISTICS` activated if you are sure 
that the camera or the needle does not move 
during the whole video.

### 5) Functional switches/features (user-space)

These switches can be set to 1 (active) and 0 (inactive) 
and are used to control program flow. For example you can 
activate `NELDERMEAD_PREFITTING` (recommended if you can 
only roughly imagine your initial values for elastic 
moduli and pressure). If you have measured for example 
the pressure or know poisson's ratio you can switch of 
the regression of these specific parameters. Note that 
this is a constraint and will probably lead to higher 
fit errors. In case of only weakly deflated capsules 
you can deactivate `EXTENDED_SHOOTING` to improve speed of the algorithm.

```
FIT_PRESSURE 1
```
Sets pressure as free parameter of the regression.
```
FIT_POISSON 1
```
Sets poisson ratio as free parameter of the regression.
```
FIT_COMPRESSION 1
```
Sets area compression modulus as free parameter of the regression.
```
EXTENDED_SHOOTING 1
```
Activates/deactivates multiple shooting method used to improve solution of single shooting method.
```
FORCE_SYMMETRY 0 (experimental)
```
Activates/deactivates rotation of contour points to align them to symmetry axis.
```
FORCE_BOUNDARY_SYMMETRY 1
```
Activates/deactivates symmetric aligning of two the contour points associated with capillary.
```
NELDERMEAD_PREFITTING 0
```
Activates/deactivates nelder-mead regression before newton regression.
```
PARAMETER_TRACING 0
```
Actiates/deactivates use of fit result in next image as initial values.
```
WRINKLING_WAVELENGTH 0.0
```
Specifies a manually measured wrinkle wavelength in units of pixels.

### 6) General directories
```
INPUT_FOLDER ./input/
```
Specifies folder, from which images are read.
```
CONFIG_FOLDER ./config/
```
Specifies folder, in which configuration file is placed.
```
OUT_FOLDER ./out/
```
Specifies folder for miscellaneous output files.
```
GLOBAL_OUT_FOLDER ./global_out/
```
Specifies folder for final output provided to user.
```
TMP_FOLDER ./tmp/
```
Specifies folder, in which automatically generated gnuplot scripts are placed.

### 7) Debugging:

These switches are very useful for debugging purposes. If something goes seriously wrong, you should activate some of the following switches to see what happens.
```
WATCH_SINGLE_SHOOTING 0
```
Activates/deactivates stdout information of single shooting method.
```
WATCH_MULTI_SHOOTING 0
```
Activates/deactivates stdout information of multiple shooting method.
```
WATCH_LAPLACE_FITTING 1
```
Activates/deactivates stdout information of reference Newton regression.
```
WATCH_HOOKE_FITTING 1
```
Activates/deactivates stdout information of elastic Newton regression.
```
WATCH_QR_DECOMPOSITION 0
```
Activates/deactivates stdout information of QR-decomposition.
```
WATCH_NELDER_MEAD 1
```
Activates/deactivates stdout information of elastic Nelder-Mead regression.
