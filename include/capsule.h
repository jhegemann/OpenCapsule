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

#ifndef CAPSULE_H
#define CAPSULE_H

#define _USE_MATH_DEFINES

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <vector>
#include <string>
#include <sstream>
#include <omp.h>
#include <math.h>
#include <algorithm>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include "objects.h"
#include "config.h"

using namespace std;

#define FMT setw(20) << setprecision(10)

typedef vector< pair<double,double> > PointSet;

/* ODE evaluation and stability testing */
/* right side evaluation of Laplace-Young system 
 * s : arc-length
 * y : solution
 * f : right side
 * data contains parameter set, solution and splines */
void LaplaceEquation(double s, double *y, double *f, LaplaceData *data);
/* right side evaluation of Hooke system
 * s : arc-length
 * y : solution
 * y_ext : curvatures and strains
 * f : right side
 * data contains parameter set, solution and splines */
void HookeEquation(double s, double *y, double *y_ext, double *f, HookeData *data);
/* jacobian calculation of Hooke system
 * jacobian stored in J
 * used for implicit integration and stability analysis */
void HookeJacobi(double s, double *y, double *y_ext, gsl_matrix *J, HookeData *data);
/* stability estimation of Hooke system */
void GerschGorin(double s, double *y, double *y_ext, HookeData *data);
/* calculates eigenvalues of matrix A and stores result in eval */
void EigenValues(gsl_matrix *A, double *eva, int N);

/* Integration */
/* runge kutta 4 integration step of Laplace-Young system
 * s : arc-length
 * h : step-size
 * solution stored in data */
void LaplaceRungeKutta(double s, double h, LaplaceData *data);
/* runge kutta 4 integration step of Hooke system
 * s : arc-length
 * h : step-size
 * solution stored in data */
void HookeRungeKutta(double s, double h, HookeData *data);
/* implicit runge kutta 4 integration 
 * solution stored in data
 * finds slopes with minimization */
void ImplicitHookeRungeKutta(double s, double h, HookeData *data);
/* implicit runge kutta 4 integration
 * solution stored in data
 * finds slopes with fixed point iteration */
void ImplicitHookeRungeKuttaFixpoint(double s, double h, HookeData *data);

/* Solve-methods */
/* find Laplace solution to parameters specified in data and store solution back to data */
void SolveLaplace(LaplaceData *data);
/* measure boundary deviation of Hooke shape in dependence of initial values tau_s(0) */
void ErrorFunctionSingleShooting(HookeData *data);
/* find Hooke solution with single shooting method and store in data */
bool SingleShooting(HookeData *data, double t_max, double tau);
bool SingleShooting(HookeData *data); /* if accuracy not reached, improve solution with parallel shooting */
/* find Hooke solution without any shooting by tracing from the reference form (does not perform well, experimental) */
bool ParallelShootingTracing(HookeData *data); 
/* find Hooke solution based on solution found by single shooting
* imporve solution found by single shooting by applying multiple shooting method */
bool ParallelShooting(HookeData *data, int p_interval); 
/* find solution of Hooke equation 
* boundary condition not checked here */
void SolveHooke(HookeData *data, double t_min, double t_max);
void SolveHooke(HookeData *data);

/* Plotting, interpolation, scaling, translation */
/* write data file from Laplace shape */
void PlotLaplace(string filename, LaplaceData *data);
/* write data file from Hooke shape */
void PlotHooke(string filename, HookeData *data);
/* scale Laplace shape to dimensionless units */
void ScaleLaplace(LaplaceData *data);
/* interpolate Laplace shape */
void InterpolateLaplace(LaplaceData *data);
/* interpolate Hooke shape */
void InterpolateHooke(HookeData *data);
/* shifting, moving, align capsules in coordinate system */
void MoveLaplaceUp(LaplaceData *data);
void MoveLaplaceDown(LaplaceData *data);
void MoveHookeUp(HookeData *data);
void MoveHookeDown(HookeData *data);
/* determine wrinkling region */
bool WrinklingRegion(HookeData *data, double &s1, double &s2);
/* determine average meridional tension in wrinkling region */
double AverageWrinklingTension(HookeData *data, double s1, double s2);
/* maximum radius */
double MaximumRadius(HookeData *data);

/* Fit-methods */
/* fitting Laplace shape to contour points (fast) */
void LaplaceFit(LaplaceData *data, PointSet *points);
/* fitting Hooke shape to contour points (fast) */
void HookeFit(HookeData *data, PointSet *points, const bool fit_p, const bool fit_nu, const bool fit_k);
bool HookeFitIteration(HookeData *data, PointSet *points, const bool *fit, vector<int> pind, int nparams, int maxparams, 
                double *parameters, double *params, LaplaceData *l, HookeData *c, gsl_matrix *J, gsl_vector *F, gsl_vector *D, gsl_vector *PI);
/* Simulated annealing fit of Laplace shape (slow) */
void Metropolis(LaplaceData *data, PointSet *points);
/* Simulated annealing fit of Hooke shape (slow) */
void Metropolis(HookeData *data, PointSet *points);
/* Nelder Mead method (derivation free) for prefitting (medium fast) 
* recommended for bad initial guess */
void NelderMead(HookeData *data, PointSet *points);

/* Nelder Mead methods */
bool CompareVertices(const Vertex *v1, const Vertex *v2);
void ReplaceWorst(vector<Vertex*> *v, Vertex *newone);
void CalcReflection(vector<Vertex*> *v, Vertex *r, Vertex *c, HookeData *dummy, PointSet *points, double alpha);
void CalcCentroid(Vertex *c, vector<Vertex*> *v);
void SortVertices(vector<Vertex*> *v, HookeData *dummy, PointSet *points);
double Error(Vertex *v, HookeData *dummy, PointSet *points);
double Error(Data *shape, PointSet *points);

/* Linear algebra */
void QR(gsl_matrix *A, gsl_vector *x, gsl_vector *b, int M, int N);
void PrintVector(gsl_vector *x, const char *name);
void PrintMatrix(gsl_matrix *x, const char *name);

/* Errors and redsiduals */
double Volume(Data *shape);
double Area(Data *shape);
double DistDeriv(Data *shape, double s, double z_i, double r_i);
double Error(Data *shape, PointSet *points);
double HookeError(HookeData *shape, PointSet *points);
Residual Dist(Data *shape, double s, double z_i, double r_i);
Residual SingleError(Data *shape, pair<double,double> point);

/* Runge-Kutta Implicit Table */
static const double _A[4] = {0.06943184420297371, 0.3300094782075719, 0.6699905217924281, 0.9305681557970262};
static const double _B[4][4] = {{0.08696371128436346, -0.0266041800849988, 0.01262746268940472, -0.003555149685795692},
{0.1881181174998681, 0.1630362887156365, -0.02788042860247092, 0.00673550059453816},
{0.1671919219741888, 0.353953006033744, 0.1630362887156365, -0.01419069493114115},
{0.1774825722545226, 0.3134451147418683, 0.3526767575162718, 0.08696371128436346}};
static const double _C[4] = {0.1739274225687269, 0.326072577431273, 0.326072577431273, 0.1739274225687269};
 
 
 
 #endif
 