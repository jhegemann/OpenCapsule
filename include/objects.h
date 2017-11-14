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

#ifndef OBJECTS_H
#define OBJECTS_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "config.h"

using namespace std;

/* classes */
class Vertex;
class Residual;
class InterpolationObject;
class Data;
class LaplaceData;
class HookeData;
class ImageData;


class Vertex
{
	public:
		
		Vertex();
		Vertex(double p, double nu, double k);
		~Vertex();
		Vertex *Copy();
		void Print();
		
		double p, nu, k;
		double error;
		bool changed;
};

class Residual
{
	public:
		
		Residual();
		~Residual();
		double Norm();
		
		double delta_r;
		double delta_z;
		double s;
};

class InterpolationObject
{
	public:
		
		InterpolationObject() ;
		~InterpolationObject();
		void Clean();
		
		gsl_interp_accel *acc;
		gsl_spline *r_s, *z_s, *psi_s;
		gsl_spline *tau_s;
		gsl_spline *r_z, *r2_z;
		
		/* external quantities */
		gsl_spline *l_s, *l_phi, *k_phi, *t_phi;
		
};

class Data
{
	public:
		
		Data() ;
		virtual ~Data() ;
		void AllocateData(int dim);
		
		/* general data */
		double *parameters;
		InterpolationObject splines;
		vector<double> s, r, r2, z, psi;
		double L0;
		int size;
		int dim;
		bool invalid;
		
		/* explicit runge kutta */
		double *y;
		double *yrk4;
		double *k1, *k2, *k3, *k4;
		double *y2, *y3, *y4;
};

class LaplaceData : public Data
{
	public:
		
		LaplaceData();
		~LaplaceData();
		void Clear();
		void Push(double t, double y0, double y02, double y1, double y2);
		double GetPressure();
		double GetDensity();
		double GetScaling();
		void SetParameters(double p, double rho, double a);
		double GetConversion();
		void SetConversion(double x);
		
		double V0;
		double A0;
		double surface_tension;
		double inner_diameter;
		double conversion;
		
};

class HookeData : public Data
{
	public:
		
		HookeData();
		~HookeData();
		void SetParameters(double p, double nu, double EH0);
		void Clear();
		void Push(double t, double y0, double y02, double y1, double y2, double y3, double yext0, double yext1, double yext2, double yext3);
		void SetInitialConditions(double r, double z, double psi, double tau_s);
		void SetExternConditions(double ls, double lp, double kp, double tp);
		double GetPressure();
		double GetPoisson();
		double GetCompression();
		double GetEH0();
		void SetLaplace(LaplaceData *data);
		
		LaplaceData *undeformed;
		vector<double> tau_s;
		
		double *y_ext; /* lambda_s, lambda_phi, kappa_phi, tau_phi */
		double mu;
		
		vector<double> lambda_s;
		vector<double> lambda_phi;
		vector<double> kappa_phi;
		vector<double> tau_phi;
		
		/* pointer for implicit rk */
		gsl_matrix *JP;
		gsl_matrix *J;
		gsl_vector *X;
		gsl_vector *Y;
		gsl_vector *K;
		gsl_vector *F;
		gsl_vector *BF;
		gsl_vector *R;
		
		gsl_vector *min_vec;
		double lambda_s_trace;
		double lambda_phi_trace;
		double tau_phi_trace_1;
		double tau_phi_trace_2;
		bool integration_mark_0;
		bool integration_mark_1;
		bool integration_mark_2;
		bool integration_mark_3;
		
};

class ImageData
{
	public:
		ImageData();
		~ImageData();
		double t_high, t_low;
		int width, height, step, channels;
		double **matrix;
};



#endif
