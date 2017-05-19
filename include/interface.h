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

#ifndef INTERFACE_H
#define INTERFACE_H

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <omp.h>

#include "image.h"
#include "capsule.h"
#include "config.h"
#include "script.h"
#include "objects.h"

using namespace std;

template <typename T>
string
to_string(T value)
{
  std::ostringstream buffer;
  buffer << value;
  return buffer.str();
}

extern vector<PointSet*> _reference_shape; 
extern vector<PointSet*> _elastic_shape;
extern vector<ImageInfo*> _reference_image;
extern vector<ImageInfo*> _elastic_image;

vector<ImageInfo*> LoadImages(string preamble, vector<string> &file, vector<PointSet*> &contour, bool elastic);
void FinalizeContours(string preamble, double scaling, vector<PointSet*> &contour);
void InitializeSolver();
void FinalizeSolver();

void SimpleLaplace(string filename, double pressure, double density, double scaling);
void SimpleHooke(string filename, double pressure_laplace, double density_laplace, double scaling_laplace, double pressure, double poisson, double compression);
LaplaceData *LaplaceMean(double pressure, double density, double scaling);
void Sequence(LaplaceData *reference, double pressure, double poisson, double compression);
void EdgeDetection(string filename);

#endif

