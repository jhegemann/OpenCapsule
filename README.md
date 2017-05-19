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
Boost, OpenCV, Gnu Scientific Library (GSL) e.g. 
on Ubuntu typing something like

```
sudo apt-get install gcc g++ make libboost-all-dev libopencv-dev libgsl0-dev libgsl0ldbl
```

in a login-shell will be sufficient. On Scientific Linux, i.e. Fedora system, use Yum package manager or install the libraries manually. OpenCapsule uses OpenCV for image handling and therefore only supports the image formats supported by your specific installation. Check if all PATH variables are set appropriately or adapt the Makefile by informing g++ explicitly where to find the libraries by using the -L option. Finally, simply type

```
make
```

to compile the project. If you like to set up a more usable working environment also type

```
make install
```

to make OpenCapsule accessible in all working directories.
