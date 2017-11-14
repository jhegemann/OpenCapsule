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

#include "objects.h"


///////////////////////////////// VERTEX CLASS ////////////////////////////////////
Vertex::Vertex() 
{ 
        changed = false;
}

Vertex::Vertex(double p, double nu, double k) 
{
        changed = false;
        this->p = p;
        this->nu = nu;
        this->k = k;
}

Vertex::~Vertex()
{
        
}

Vertex *Vertex::Copy()
{
        Vertex *cpy = new Vertex(this->p, this->nu, this->k);
        cpy->error = this->error;
        return cpy;
}

void Vertex::Print()
{
        cout << setw(20) << "Set" << setw(20) << p << setw(20) << nu << setw(20) << k << setw(20) << "Error" << setw(20) << error << endl;
}



///////////////////////////////// RESIDUAL CLASS ////////////////////////////////////

Residual::Residual()
{
        
}

Residual::~Residual()
{
        
}

double Residual::Norm()
{
        return (delta_r * delta_r + delta_z * delta_z);
}


///////////////////////////////// INTERPOLATIONOBJECT CLASS ////////////////////////////////////


InterpolationObject::InterpolationObject() 
{ 
        acc = NULL;
        r_s = NULL;
        z_s = NULL;
        psi_s = NULL;
        tau_s = NULL;
        r_z = NULL;
        r2_z = NULL;
        l_s = NULL; 
        l_phi = NULL;
        k_phi = NULL; 
        t_phi = NULL;
}

InterpolationObject::~InterpolationObject()
{
        Clean();
}

void InterpolationObject::Clean()
{
        if (r_s) gsl_spline_free(r_s);
        if (z_s) gsl_spline_free(z_s);
        if (psi_s) gsl_spline_free(psi_s);
        if (tau_s) gsl_spline_free(tau_s);
        if (r_z) gsl_spline_free(r_z);
        if (r2_z) gsl_spline_free(r2_z);
        if (l_s) gsl_spline_free(l_s);
        if (l_phi) gsl_spline_free(l_phi);
        if (k_phi) gsl_spline_free(k_phi);
        if (t_phi) gsl_spline_free(t_phi);
        if (acc) gsl_interp_accel_free(acc);
}


///////////////////////////////// DATA CLASS ////////////////////////////////////

Data::Data() 
{ 
        parameters = (double*)calloc(3, sizeof(double));
        size = 0;
        invalid = true;
}

Data::~Data() 
{ 
        free(parameters);
        
        free(y);
        
        free(yrk4);
        
        free(k1);
        free(k2);
        free(k3);
        free(k4);
        
        free(y2);
        free(y3);
        free(y4);
}

void Data::AllocateData(int dim)
{
        y = (double*)calloc(dim, sizeof(double));
        
        yrk4 = (double*)calloc(dim, sizeof(double));
        
        k1 = (double*)calloc(dim, sizeof(double));
        k2 = (double*)calloc(dim, sizeof(double));
        k3 = (double*)calloc(dim, sizeof(double));
        k4 = (double*)calloc(dim, sizeof(double));
        
        y2 = (double*)calloc(dim, sizeof(double));
        y3 = (double*)calloc(dim, sizeof(double));
        y4 = (double*)calloc(dim, sizeof(double));
}


///////////////////////////////// LAPLACEDATA CLASS ////////////////////////////////////


LaplaceData::LaplaceData() : Data()
{ 
        conversion = 1.0/SCALE_IMAGE;
        dim = 3;
        AllocateData(dim);
}

LaplaceData::~LaplaceData()
{
        Clear();
}

void LaplaceData::Clear()
{
        s.clear();
        r.clear();
        r2.clear();
        z.clear();
        psi.clear();
        size = 0;
}

void LaplaceData::Push(double t, double y0, double y02, double y1, double y2)
{
        s.push_back(t);
        r.push_back(y0);
        r2.push_back(y02);
        z.push_back(y1);
        psi.push_back(y2);
        size++;
}

double LaplaceData::GetPressure()
{
        return parameters[0];
}

double LaplaceData::GetDensity()
{
        return parameters[1];
}

double LaplaceData::GetScaling()
{
        return parameters[2];
}

void LaplaceData::SetParameters(double p, double rho, double a)
{
        parameters[0] = p;
        parameters[1] = rho;
        parameters[2] = a;
} 	

double LaplaceData::GetConversion()
{
        return conversion;
}

void LaplaceData::SetConversion(double x)
{
        conversion = x;
}


///////////////////////////////// HOOKEDATA CLASS ////////////////////////////////////


HookeData::HookeData() : Data()
{ 
        dim = 4;
        y_ext = (double*)calloc(4, sizeof(double));
        /* implicit rk containers */
        if (IMPLICIT_INTEGRATION)
        {
                JP = gsl_matrix_calloc(dim, dim);
                J = gsl_matrix_calloc(dim * 4, dim);
                X = gsl_vector_calloc(dim * 4);
                Y = gsl_vector_calloc(dim * 4);
                K = gsl_vector_calloc(dim * 4);
                F = gsl_vector_calloc(dim * 4);
                BF = gsl_vector_calloc(dim * 4);
                R = gsl_vector_calloc(dim * 4);
        }
        /* explicit rk containers */
        AllocateData(dim);
	/* allocate minimizer */
	min_vec = gsl_vector_alloc(2);
}

HookeData::~HookeData()
{
        Clear();
        free(y_ext);
        if (IMPLICIT_INTEGRATION)
        {
                gsl_matrix_free(JP);
                gsl_matrix_free(J);
                gsl_vector_free(X);
                gsl_vector_free(Y);
                gsl_vector_free(K);
                gsl_vector_free(F);
                gsl_vector_free(BF);
                gsl_vector_free(R);
        }
        undeformed = NULL;
}

void HookeData::SetParameters(double p, double nu, double k)
{
        /* pressure, poisson, compression-modulus */
        parameters[0] = p;
        parameters[1] = nu;
        parameters[2] = k;
}

void HookeData::Clear()
{
        s.clear();
        r.clear();
        r2.clear();
        z.clear();
        psi.clear();
        tau_s.clear();
        lambda_s.clear();
        lambda_phi.clear();
        kappa_phi.clear();
        tau_phi.clear();
        size = 0;
}

void HookeData::Push(double t, double y0, double y02, double y1, double y2, double y3, double yext0, double yext1, double yext2, double yext3)
{
        s.push_back(t);	
        r.push_back(y0);
        r2.push_back(y02);
        z.push_back(y1);
        psi.push_back(y2);
        tau_s.push_back(y3);
        lambda_s.push_back(yext0);
        lambda_phi.push_back(yext1);
        kappa_phi.push_back(yext2);
        tau_phi.push_back(yext3);
        size++;
}

void HookeData::SetInitialConditions(double r, double z, double psi, double tau_s)
{
        y[0] = r;
        y[1] = z;
        y[2] = psi;
        y[3] = tau_s;
        SetExternConditions(0.0, 0.0, 0.0, 0.0);
}

void HookeData::SetExternConditions(double ls, double lp, double kp, double tp)
{
        y_ext[0] = ls;
        y_ext[1] = lp;
        y_ext[2] = kp;
        y_ext[3] = tp;
}

double HookeData::GetPressure()
{
        return parameters[0]; 
}

double HookeData::GetPoisson()
{
        return parameters[1];
}

double HookeData::GetCompression()
{
        return parameters[2]; 
}

double HookeData::GetEH0()
{
        return GetCompression() * (2.0 * (1.0 - GetPoisson()));
}

void HookeData::SetLaplace(LaplaceData *data)
{
        undeformed = data; 
}


///////////////////////////////// IMAGEDATA CLASS ////////////////////////////////////

ImageData::ImageData()
{
        
}

ImageData::~ImageData()
{
        
}

