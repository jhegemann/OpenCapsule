/* OpenCapsule
 * Copyright (C) 2017 Jonas Hegemann, TU Dortmund
 * mail: jonas.hegemann@hotmail.de 
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#ifndef IMAGE_H 
#define IMAGE_H

#define _USE_MATH_DEFINES

#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#include "objects.h"
#include "config.h"
#include "capsule.h"

using namespace std;
using namespace cv;

typedef vector< pair<double,double> > PointSet;

class ImageInfo
{
public:
        ImageInfo() : a(0.0), b(0.0), h(0.0), y_top(0.0), y_bottom(0.0), x_left(0.0), x_right(0.0) {}
        ~ImageInfo() {}
        double a;
        double b;
        double h;
        double y_top;
        double y_bottom;
        double x_left;
        double x_right;
        string path;
};

class Image
{
public:
        
        ImageData foto;
        bool fix_capheight;
        int capheight;
        
        Image();
        ~Image();
        
        /* clean up */
        void FreeImage();
        
        /* used for Graham Scan algorithm */
        static bool ComparePoints(const pair<double,double> p1, const pair<double,double> p2);
        
        /* specified to width & height */
        double **AllocateMatrix();
        
        /* specified to width & height */
        void FreeMatrix(double **matrix);
        
        /* measure height of capsule from binary image */
        double MeasureH(double **image);
        
        /* measure inner capillary diameter from binary image */
        double MeasureA(double **image);
        
        /* delete duplicates from list */
        void DeleteDuplicates(PointSet *v);
        
        /* align principal axis along vertical y-axis */
        void ForceSymmetry(PointSet *v);
        
        /* extract contour points from binary image */
        PointSet *PointCloud(double **image);
        
        /* determine average wrinkle wavelength from edge detection in capsule inner */
        double WrinkleWavelength(HookeData *shape);
        
        /* use with ClosedContour
         * recursive coloring algorithm
         * coloring stored in c
         * output parameter : c
         * c should be initialized as c[x][y] = 0 for all x,y
         * input parameter : binary-image 
         * output parameter : c */
        void Color(double **image, double **c, int x, int y);
        
        /* check if closed contour exists in binary image */
        bool ClosedContour(double **image);
        
        /* measure distance from left image border to capillary */
        double MeasureXLeft(double **image);
        
        /* measure distance from right image border to capillary */
        double MeasureXRight(double **image);
        
        /* measure inner capillary diameter from binary image */
        double MeasureB(double **image);
        
        /* measure capsule distance from lower image border from binary image */
        double MeasureYBottom(double **image);
        
        /* measure capillary height from binary image */
        double MeasureYTop(double **image);
        
        /* measure capillary height on left side from binary image */
        int MeasureYTopLeft(double **image);
        
        /* measure capillary height on right side from binary image */
        int MeasureYTopRight(double **image);
        
        /* convolution
         * apply arbitrary filter with specified radius to image */
        double **ApplyFilter(double **image, double **filter, int radius);
        
        /* combines images r = r_min, ...,r_max
         * trace edges with hysteresis 
         * final part of Canny edge detection */
        double **CombinedHysteresis(string filename, const bool wrinkles = false);
        
        /* compute gradient image 
         * output parameters : grad, gradphi, gradphi0 */
        void Gradient(double **grad, double **gradphi, double **gradphi0, double **dx, double **dy);
        
        /* r : radius of smoothing operator
         * read image
         * smooth image
         * calculate derivatives
         * output parameter : grad, gradphi0 */
        void ConvertImage(double **grad, double **gradphi0, int r);
        
        /* iterative edge-tracing
         * applies threshold t_high and t_low
         * output parameter : c */
        void TraceEdge(double **grad, double **gradphi0, double **binary, double max_brightness);
        
        /* check if colored neighbor exists */
        bool NeighborColored(double **c, int x, int y);
        
        /* hysteris 
         * output parameter : binary */
        void ApplyHysteresis(double **grad, double **gradphi0, double **binary);
        
        /* determine maximum brightness of grey scale image */
        double MaxBrightness(double **image);
        
        /* check if gradient image has maximum perpendicular to edge */
        bool CheckMaxInDir(double **grad, double **gradphi0,  int x, int y);
        
        /* check if gradient image has maximum at pixel (x,y) */
        bool CheckMax(double **grad, int x, int y);
        
        /* check for maximum in 135 degree direction */
        bool CheckMax3(double** gradient, int x, int y);
        
        /* check for maximum in 90 degree direction */
        bool CheckMax2(double **gradient, int x, int y);
        
        /* check for maximum in 45 degree direction */
        bool CheckMax1(double **gradient, int x, int y);
        
        /* check for maximum in 0 degree direction */
        bool CheckMax0(double **gradient, int x, int y);
        
        /* write image to file */
        int Write(const char* filename, double **i);
        
        /* read image from file */
        int Read(const char* filename_in, const char* filename_out);
        
        /* write laplace capsule image together with shape */
        static void WriteLaplace(string filename_out, ImageInfo *img_data, LaplaceData *shape, double conversion);	
        
        /* write hooke capsule image together with shape */
        static void WriteHooke(string filename_out, ImageInfo *img_data, HookeData *shape, double conversion);
        
};



#endif
