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

OpenCapsule was tested (successful compiling and running) on 
Ubuntu 16.04 LTS and Scientific Linux 6.6 with GNU/GCC compiler. 
To install OpenCapsule first download the tarball and extract 
it to some arbitrary direction or simply clone the repository. 
Before compiling install the required dependencies
Boost, OpenCV and Gnu Scientific Library (GSL) e.g. 
on Ubuntu typing something like

```
sudo apt-get install gcc g++ make libboost-all-dev libopencv-dev libgsl0-dev libgsl0ldbl
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
file "config/config.cfg" before you proceed. 
List your input files, which should placed 
(by default) in the "input" folder just created, 
under REFERENCE_SHAPE and ELASTIC_SHAPE. 
OpenCapsule averages over reference images, 
which should all show the same, undeformed 
state of the capsule. Make sure you have 
specified the fields EXPERIMENT_CAPDIAMETER 
and EXPERIMENT_DENSITY in SI units. If you 
are not sure how to handle the configuration 
file, see the corresponding section for 
more details about this. Next, typing
```
OpenCapsule -h
```
will give you an overview over the command line 
options as well as a list of common examples. 
The results of the following standard commands 
are placed in the "global_out" folder, where 
you will find the results listed in plain 
*.txt files, images in *.png format, and 
wrapped in a comprehensive html-report. 

### Edge Detection

```
OpenCapsule -i [path to image]
```
This command performs an isolated image analysis 
of the image specified by the option --image. 
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

### Hookean Analysis
```
OpenCapsule -s
```
This command triggers the Hookean shape analysis 
and will automatically use the specified reference 
shapes. You can also provide initial values by using the call
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
