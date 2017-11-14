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

#include "capsule.h"

void (*ShapeEquation)(double, double*, double*, double*, HookeData*);

void SetConstitutiveLaw() 
{
	switch (CONSTITUTIVE_LAW) 
	{
		case 0:
			ShapeEquation = &HookeEquationLinear;
			break;
		case 1:
			ShapeEquation = &HookeEquation;
			break;
		case 2:
			ShapeEquation = &HookeEquationMooneyRivlin;
			break;
		default:
			printf("invalid constitutive law!\n"); fflush(stdout);
			exit(1);
			break;
	}	
}

struct rivlin_params { double ts; double tp; double ls; double lp; double em; double psi; };
struct hooke_params { double ts; double tp; double ls; double lp; double eh0; double nu; };

bool CompareVertices(const Vertex *v1, const Vertex *v2)
{
        return (v1->error < v2->error);
}

double TauS(double ls, double lp, double EM, double PSI)
{
	return (EM/(3.0*ls*lp))*(ls*ls-1.0/(ls*ls*lp*lp))*(PSI+(1.0-PSI)*lp*lp) + GAMMA_SCALE;
}

double TauP(double ls, double lp, double EM, double PSI)
{
	return (EM/(3.0*ls*lp))*(lp*lp-1.0/(ls*ls*lp*lp))*(PSI+(1.0-PSI)*ls*ls) + GAMMA_SCALE;
}

double TauSMinWrinkle(double ls, double lp, double tsbar, double lpbar, double EM, double PSI) 
{
	return tsbar - (EM/(3.0*ls*lpbar))*(ls*ls-1.0/(ls*ls*lp*lp))*(PSI+(1.0-PSI)*lp*lp) - (lp/lpbar)*GAMMA_SCALE;
}

double TauPMinWrinkle(double ls, double lp, double EM, double PSI)
{
	return (EM/(3.0*ls*lp))*(lp*lp-1.0/(ls*ls*lp*lp))*(PSI+(1.0-PSI)*ls*ls) + GAMMA_SCALE;
}

double TauSMin(double ts, double tp, double ls, double lp, double EM, double PSI) 
{
	return ts - (EM/(3.0*ls*lp))*(ls*ls-1.0/(ls*ls*lp*lp))*(PSI+(1.0-PSI)*lp*lp) - GAMMA_SCALE;
}

double TauPMin(double ts, double tp, double ls, double lp, double EM, double PSI)
{
	return tp - (EM/(3.0*ls*lp))*(lp*lp-1.0/(ls*ls*lp*lp))*(PSI+(1.0-PSI)*ls*ls) - GAMMA_SCALE;
}

double TauMinSing(double t, double l, double EM, double PSI) 
{
	return t - (EM/3.0)*(1.0-1.0/pow(l, 6))*(PSI+(1.0-PSI)*l*l) - GAMMA_SCALE;
}

int StrainMooneyRivlin(const gsl_vector *x, void *params, gsl_vector *f) 
{
	struct rivlin_params *p = (struct rivlin_params *) params;
	
	const double x0 = gsl_vector_get(x, 0);
	const double x1 = gsl_vector_get(x, 1);

	const double y0 = TauSMin(p->ts, x0, x1, p->lp, p->em, p->psi);
	const double y1 = TauPMin(p->ts, x0, x1, p->lp, p->em, p->psi);

	gsl_vector_set(f, 0, y0);
	gsl_vector_set(f, 1, y1);

	return GSL_SUCCESS;
}

double StrainMooneyRivlinSing(double x, void *params) 
{
	struct rivlin_params *p = (struct rivlin_params *) params;
	return TauMinSing(p->ts, x, p->em, p->psi);
}

int StrainMooneyRivlinWrinkle(const gsl_vector *x, void *params, gsl_vector *f) 
{
	struct rivlin_params *p = (struct rivlin_params *) params;
	
	const double x0 = gsl_vector_get(x, 0);
	const double x1 = gsl_vector_get(x, 1);

	const double y0 = TauSMinWrinkle(x0, x1, p->ts, p->lp, p->em, p->psi);
	const double y1 = TauPMinWrinkle(x0, x1, p->em, p->psi);

	gsl_vector_set(f, 0, y0);
	gsl_vector_set(f, 1, y1);

	return GSL_SUCCESS;
}

double TauSMinHookeWrinkle(double ls, double lp, double tsbar, double lpbar, double EH0, double NU) 
{
	return tsbar - (lp/lpbar)*(EH0/(1.0-NU*NU))*((ls-1.0)+NU*(lp-1.0)) - (lp/lpbar)*GAMMA_SCALE;
}

double TauPMinHookeWrinkle(double ls, double lp, double EH0, double NU) 
{
	return (EH0/(1.0-NU*NU))*((lp-1.0)+NU*(ls-1.0)) + GAMMA_SCALE;
}

int StrainHookeLinearWrinkle(const gsl_vector *x, void *params, gsl_vector *f) 
{
	struct hooke_params *p = (struct hooke_params *)params;
	
	const double x0 = gsl_vector_get(x, 0);
	const double x1 = gsl_vector_get(x, 1);
	
	const double y0 = TauSMinHookeWrinkle(x0, x1, p->ts, p->lp, p->eh0, p->nu);
	const double y1 = TauPMinHookeWrinkle(x0, x1, p->eh0, p->nu);
	
	gsl_vector_set(f, 0, y0);
	gsl_vector_set(f, 1, y1);
	
	return GSL_SUCCESS;
}

/* output parameter: right side f */
void LaplaceEquation(double s, double *y, double *f, LaplaceData *data)
{
        /* params[0] : reference pressure - p0
         * params[1] : density - rho0
         * params[2] : scaling factor - a */
        int dim = data->dim;
        double *pi = data->parameters;
        double p = pi[0] / pi[2];
        double r = pi[1] / pi[2] / pi[2];
        if (fabs(s) < 1.0e-8)
        {
                f[0] = cos(y[2]);
                f[1] = sin(y[2]);
                f[2] = 0.5*(p-r*y[1]);
        }
        else
        {
                f[0] = cos(y[2]);
                f[1] = sin(y[2]);
                f[2] = (p-r*y[1])-sin(y[2])/(y[0]);
        }
}

const gsl_multiroot_fsolver_type *T2d = gsl_multiroot_fsolver_hybrids;
gsl_multiroot_fsolver *s2d = gsl_multiroot_fsolver_alloc(T2d, 2);  
const gsl_root_fsolver_type *T1d = gsl_root_fsolver_brent;
gsl_root_fsolver *s1d = gsl_root_fsolver_alloc(T1d);

/* output parameter: right side f */
void HookeEquationMooneyRivlin(double s, double *y, double *y_ext, double *f, HookeData *data)
{
        /* params[0] : deformed pressure - p
         * params[1] : poisson ratio - nu
         * params[2] : area compression modulus - k */
	struct rivlin_params params;
        int dim = data->dim;
        double *pi = data->parameters;
        double a = data->undeformed->parameters[2];
        double p = pi[0] / a;
        double nu = pi[1];
	double EH0 = pi[2];
	double rho = data->undeformed->parameters[1] / a / a;
        double r0 = gsl_spline_eval(data->undeformed->splines.r_s, s, data->undeformed->splines.acc);
        if (fabs(s) < 1.0e-8)
        {
                /* singularity evaluation in apex/on symmetry axis */
		gsl_function strain_function;
		strain_function.function = &StrainMooneyRivlinSing;
		params.ts = y[3];
		params.em = EH0;
		params.psi = nu;
		strain_function.params = &params;
		double x_lo = 0.5;
		double x_hi = 1.5;
		double r = 1;
		gsl_root_fsolver_set(s1d, &strain_function, x_lo, x_hi);
		int status;
		int iter = 0;
		do
		{
			iter++;
			status = gsl_root_fsolver_iterate(s1d);
			r = gsl_root_fsolver_root(s1d);
			x_lo = gsl_root_fsolver_x_lower(s1d);
			x_hi = gsl_root_fsolver_x_upper(s1d);
			status = gsl_root_test_interval(x_lo, x_hi, 1.0e-4, 1.0e-4);		 
		}
		while (status == GSL_CONTINUE && iter < 1000);
		y_ext[0] = r;
		y_ext[1] = y_ext[0];
                f[0] = y_ext[0] * cos(y[2]);
                f[1] = y_ext[0] * sin(y[2]);
                f[2] = (y_ext[0] * p) / (2.0 * y[3]);
                f[3] = 0.0;
		
		data->lambda_s_trace = y_ext[0];
		data->tau_phi_trace_1 = y[3];
		data->integration_mark_0 = true;
        }
        else
        {
                y_ext[2] = sin(y[2]) / y[0];
                y_ext[1] = y[0] / r0;
		params.ts = y[3];
		params.lp = y_ext[1];
		params.em = EH0;
		params.psi = nu;
		gsl_multiroot_function strain_function = {&StrainMooneyRivlin, 2, &params};
		if (!data->integration_mark_1) {
			/* before wrinkling */
			gsl_vector_set(data->min_vec, 0, data->tau_phi_trace_1); /* TRACER */
		} else {
			if (data->integration_mark_2) {
				/* after wrinkling */
				gsl_vector_set(data->min_vec, 0, data->tau_phi_trace_2); /* TRACER */
			}
		}
		gsl_vector_set(data->min_vec, 1, data->lambda_s_trace); /* TRACER */
		gsl_multiroot_fsolver_set(s2d, &strain_function, data->min_vec);
		int status;
		int iter = 0;
		do {
			iter++;
			status = gsl_multiroot_fsolver_iterate(s2d);
			if(status) {
				break;
			}
			status = gsl_multiroot_test_residual(s2d->f, 1e-16);
		} while(status == GSL_CONTINUE && iter < 1000);
		y_ext[0] = gsl_vector_get(s2d->x, 1); /* lambda_s */
		y_ext[3] = gsl_vector_get(s2d->x, 0); /* tau_phi */
                if (y_ext[3] < 0.0)
                {
			gsl_multiroot_function strain_function = {&StrainMooneyRivlinWrinkle, 2, &params};
			gsl_vector_set(data->min_vec, 0, data->lambda_s_trace); /* TRACER */ /* lambda_s */
			gsl_vector_set(data->min_vec, 1, data->lambda_phi_trace); /* TRACER */ /* lambda_phi */
			gsl_multiroot_fsolver_set(s2d, &strain_function, data->min_vec);
			int status;
			int iter = 0;
			do {
				iter++;
				status = gsl_multiroot_fsolver_iterate(s2d);
				if(status) {
					break;
				}
				status = gsl_multiroot_test_residual(s2d->f, 1e-16);
			} while(status == GSL_CONTINUE && iter < 1000);
			y_ext[0] = gsl_vector_get(s2d->x, 0);
			data->lambda_phi_trace = gsl_vector_get(s2d->x, 1);
		} 
		data->tau_phi_trace_1 = y_ext[3]; /* TRACER */
		data->tau_phi_trace_2 = y_ext[3]; /* TRACER */
 		data->lambda_s_trace = y_ext[0];
                f[0] = y_ext[0] * cos(y[2]);
                f[1] = y_ext[0] * sin(y[2]);
                if (y_ext[3] < 0.0)
                {
                        /* wrinkling domain */
			if(!data->integration_mark_1) {
				data->lambda_phi_trace = y_ext[1];
				data->integration_mark_1 = true;
			}
                        f[2] = (y_ext[0] / y[3]) * (p - rho * y[1]);
                        f[3] = (-1.0) * y_ext[0] * cos(y[2]) * y[3] / y[0];
                }
                else
                {
                        /* non wrinkling domain */
			if(data->integration_mark_1) {
				data->tau_phi_trace_2 = y_ext[3];
				data->integration_mark_2 = true;
			}
                        f[2] = (y_ext[0] / y[3]) * (p - rho * y[1] - y_ext[2] * y_ext[3]);
                        f[3] = (-1.0) * y_ext[0] * cos(y[2]) * (y[3] - y_ext[3]) / y[0];
                }
        }
}

/* output parameter: right side f */
void HookeEquationLinear(double s, double *y, double *y_ext, double *f, HookeData *data)
{
        /* params[0] : deformed pressure - p
         * params[1] : poisson ratio - nu
         * params[2] : area compression modulus - k */
	struct hooke_params params;
        int dim = data->dim;
        double *pi = data->parameters;
        double a = data->undeformed->parameters[2];
        double p = pi[0] / a;
        double nu = pi[1];
        double EH0 = pi[2] * (2.0 * (1.0 - nu));
	double rho = data->undeformed->parameters[1] / a / a;
        double r0 = gsl_spline_eval(data->undeformed->splines.r_s, s, data->undeformed->splines.acc);
        if (fabs(s) < 1.0e-8)
        {
                y_ext[0] = 1.0 + (y[3]-GAMMA_SCALE)*(1.0-nu)/EH0;
		y_ext[1] = y_ext[0];
                f[0] = y_ext[0] * cos(y[2]);
                f[1] = y_ext[0] * sin(y[2]);
                f[2] = (y_ext[0] * p) / (2.0 * y[3]);
                f[3] = 0.0;
        }
        else
        {
                y_ext[2] = sin(y[2]) / y[0];
                y_ext[1] = y[0] / r0;
                y_ext[0] = ((1.0 - nu * nu) * (y[3] - GAMMA_SCALE)) / EH0 + 1.0 - nu * (y_ext[1] - 1.0);
                y_ext[3] = (EH0 / (1.0 - nu * nu)) * ((y_ext[1] - 1.0) + nu * (y_ext[0] - 1.0)) + GAMMA_SCALE;
                if (y_ext[3] < 0.0) 
                {
			double coeff1 = EH0*(EH0+GAMMA_SCALE*(-1.0+nu))*(1.0+2.0*nu);
			double coeff2 = EH0*EH0*(EH0*EH0 + GAMMA_SCALE*GAMMA_SCALE*(-1.0 + nu)*(-1.0 + nu) + 2.0*EH0*(GAMMA_SCALE*(-1.0+nu) - 2.0*y_ext[1]*y[3]*nu));
			double ls2 = (coeff1-sqrt(coeff2))/(2.0*EH0*EH0*nu);
			y_ext[0] = ls2;
                }		
                f[0] = y_ext[0] * cos(y[2]);
                f[1] = y_ext[0] * sin(y[2]);
                if (y_ext[3] < 0.0)
                {
                        f[2] = (y_ext[0] / y[3]) * (p - rho * y[1]);
                        f[3] = (-1.0) * y_ext[0] * cos(y[2]) * y[3] / y[0];
                }
                else
                {
                        /* non wrinkling domain */
                        f[2] = (y_ext[0] / y[3]) * (p - rho * y[1] - y_ext[2] * y_ext[3]);
                        f[3] = (-1.0) * y_ext[0] * cos(y[2]) * (y[3] - y_ext[3]) / y[0];
                }
        }
}

/* output parameter: right side f */
void HookeEquation(double s, double *y, double *y_ext, double *f, HookeData *data)
{
        /* params[0] : deformed pressure - p
         * params[1] : poisson ratio - nu
         * params[2] : area compression modulus - k */
        int dim = data->dim;
        double *pi = data->parameters;
        double a = data->undeformed->parameters[2];
        double p = pi[0] / a;
        double nu = pi[1];
        double EH0 = pi[2] * (2.0 * (1.0 - nu));
	double rho = data->undeformed->parameters[1] / a / a;
        double r0 = gsl_spline_eval(data->undeformed->splines.r_s, s, data->undeformed->splines.acc);
        if (fabs(s) < 1.0e-8)
        {
                /* singularity evaluation in apex/on symmetry axis */
                y_ext[0] = (1.0 + nu) / ((1.0 + nu) - (y[3] - GAMMA_SCALE) * (1.0 - nu * nu) / EH0); 
		y_ext[1] = y_ext[0];
                f[0] = y_ext[0] * cos(y[2]);
                f[1] = y_ext[0] * sin(y[2]);
                f[2] = (y_ext[0] * p) / (2.0 * y[3]);
                f[3] = 0.0;
        }
        else
        {
                y_ext[2] = sin(y[2]) / y[0];
                y_ext[1] = y[0] / r0;
                y_ext[0] = ((1.0 - nu * nu) * y_ext[1] * (y[3] - GAMMA_SCALE)) / EH0 + 1.0 - nu * (y_ext[1] - 1.0);
                y_ext[3] = (EH0 / (1.0 - nu * nu)) * (1.0 / y_ext[0]) * ((y_ext[1] - 1.0) + nu * (y_ext[0] - 1.0)) + GAMMA_SCALE;
                if (y_ext[3] < 0.0) 
                {
                        y_ext[0] = (y[3]*y_ext[1]+EH0-GAMMA_SCALE*(1.0+nu))/(EH0-2.0*nu*GAMMA_SCALE-(1.0-nu*nu)*GAMMA_SCALE*GAMMA_SCALE/EH0);
                }		
                f[0] = y_ext[0] * cos(y[2]);
                f[1] = y_ext[0] * sin(y[2]);
                if (y_ext[3] < 0.0)
                {
                        /* wrinkling domain */
                        f[2] = (y_ext[0] / y[3]) * (p - rho * y[1]);
                        f[3] = (-1.0) * y_ext[0] * cos(y[2]) * y[3] / y[0];
                }
                else
                {
                        /* non wrinkling domain */
                        f[2] = (y_ext[0] / y[3]) * (p - rho * y[1] - y_ext[2] * y_ext[3]);
                        f[3] = (-1.0) * y_ext[0] * cos(y[2]) * (y[3] - y_ext[3]) / y[0];
                }
        }
}


/* calculate eigenvalues of matrix A and store in eva */
void EigenValues(gsl_matrix *A, double *eva, int N)
{
        gsl_vector *eval = gsl_vector_calloc(N);
        gsl_matrix *evec = gsl_matrix_calloc(N,N);
        
        gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(N);
        gsl_eigen_symmv(A, eval, evec, w);
        gsl_eigen_symmv_free(w);
        
        for (int i = 0; i < N; i++) eva[i] = gsl_vector_get(eval, i);
        
        gsl_vector_free(eval);
        gsl_matrix_free(evec);
}

/* calculate gersch-gorin circles to estimate stability */
void GerschGorin(double s, double *y, double *y_ext, HookeData *data)
{
        bool print = false;
        gsl_matrix *J = gsl_matrix_calloc(4,4);
        HookeJacobi(s, y, y_ext, J, data);
        double c1, c2, c3, c4;
        double ival[4][2];
        for (int i = 0; i < 4; i++)
        {
                c1 = gsl_matrix_get(J, i, i);
                /* row and column sum */
                ival[i][0] = 0.0;
                ival[i][1] = 0.0;
                for (int k = 0; k < 4; k++)
                {
                        if (k != i)
                        {
                                ival[i][0] += fabs(gsl_matrix_get(J, i, k));
                                ival[i][1] += fabs(gsl_matrix_get(J, k, i));
                        }
                }
                /* swap if neccessary */
                if (ival[i][0] > ival[i][1])
                {
                        double tmp = ival[i][0];
                        ival[i][0] = ival[i][1];
                        ival[i][1] = tmp;
                }
                /* condition of matrix bad? */
                if ((ival[i][0] > 1e4) || (ival[i][0] < 1e-4)) print = true;
        }
        if (print) /* debug purposes */
        {
                cout << FMT << "+++++" << FMT << "+++++" << flush << endl;
                for (int i = 0; i < 4; i++)
                {
                        cout << FMT << ival[i][0] << FMT << ival[i][1] << flush << endl;
                }
        }
        gsl_matrix_free(J);
}

/* output parameter: jacobian J 4 x 4 */
void HookeJacobi(double s, double *y, double *y_ext, gsl_matrix *J, HookeData *data)
{
        int dim = data->dim;
        double *pi = data->parameters;
        double a = data->undeformed->parameters[2];
        double p = pi[0] / a;
        double nu = pi[1];
        double EH0 = pi[2] * (2.0 * (1.0 - nu));
        double rho = data->undeformed->parameters[1] / a / a;
        double r0 = gsl_spline_eval(data->undeformed->splines.r_s, s, data->undeformed->splines.acc);
        
        if (s <= 0.0)
        {
                /* singularity evaluation in apex/on symmetry axis */
                y_ext[0] = (1.0 + nu) / ((1.0 + nu) - (y[3] - 1.0) * (1.0 - nu * nu) / EH0); 
        }
        else
        {
                y_ext[2] = sin(y[2]) / y[0];
                y_ext[1] = y[0] / r0;
                y_ext[0] = ((1.0 - nu * nu) * y_ext[1] * (y[3] - 1.0)) / EH0 + 1.0 - nu * (y_ext[1] - 1.0);
                y_ext[3] = (EH0 / (1.0 - nu * nu)) * (1.0 / y_ext[0]) * ((y_ext[1] - 1.0) + nu * (y_ext[0] - 1.0)) + 1.0;
                if (y_ext[3] < 0.0) 
                {
                        y_ext[0] = (y[3]*y_ext[1]+EH0-1.0*(1.0+nu))/(EH0-2.0*nu*1.0-(1.0-nu*nu)*1.0*1.0/EH0);
                }
        }
        
        double ls = y_ext[0];
        double lphi = y_ext[1];
        double kphi = y_ext[2];
        double tphi = y_ext[3];
        double r = y[0];
        double z = y[1];
        double psi = y[2];
        double ts = y[3];
        double kphi_dr = (-1.0) * sin(psi) / r / r;
        double lphi_dr = 1.0 / r0;
        double ls_dr = -nu/r0 + (ts-1.0)*(1-nu*nu)/r0/EH0;
        double tphi_dr = EH0*(nu*ls_dr + 1.0/r0)/ls/(1-nu*nu) + ls_dr*tphi/ls;
        double kphi_dpsi = cos(psi)/r;
        double ls_dts = lphi * (1-nu*nu)/EH0;
        double tphi_dts = lphi*nu/ls - lphi*(lphi+nu*(ls-1.0)-1.0)/ls/ls;
        double theta = 1.0;
        theta /= EH0 - 1.0*1.0*(1.0-nu*nu)/EH0 - 2.0*nu*1.0;
        double sigma = 1.0 + nu + (1.0-nu*nu)*(ts-1.0)/EH0;
        double tmp;
        
        if (s <= 0.0)
        {
                /* singularity evaluation in apex/on symmetry axis */
                
                /* r - column */
                gsl_matrix_set(J, 0, 0, 0.0);
                gsl_matrix_set(J, 1, 0, 0.0);
                gsl_matrix_set(J, 2, 0, 0.0);
                gsl_matrix_set(J, 3, 0, 0.0);
                
                /* z - column */
                gsl_matrix_set(J, 0, 1, 0.0);
                gsl_matrix_set(J, 1, 1, 0.0);
                gsl_matrix_set(J, 2, 1, 0.0);
                gsl_matrix_set(J, 3, 1, 0.0);
                
                /* psi - column */
                tmp = (-1.0)*ls*sin(psi);
                gsl_matrix_set(J, 0, 2, tmp);
                tmp = ls*cos(psi);
                gsl_matrix_set(J, 1, 2, tmp);
                gsl_matrix_set(J, 2, 2, 0.0);
                gsl_matrix_set(J, 3, 2, 0.0);
                
                /* tau_s - column */
                tmp = (1.0+nu)*(1-nu*nu)*cos(psi)/EH0/sigma/sigma;
                gsl_matrix_set(J, 0, 3, tmp);
                tmp = (1.0+nu)*(1-nu*nu)*sin(psi)/EH0/sigma/sigma;
                gsl_matrix_set(J, 1, 3, tmp);
                tmp = (-1.0)*p*ls/2.0/ts/ts + p*(1.0+nu)*(1-nu*nu)/2.0/EH0/sigma/sigma/ts;
                gsl_matrix_set(J, 2, 3, tmp);
                gsl_matrix_set(J, 3, 3, 0.0);
        }
        else
        {
                if (tphi < 0.0)
                {
                        /* wrinkling domain */
                        
                        /* r - column */
                        gsl_matrix_set(J, 0, 0, 0.0);
                        gsl_matrix_set(J, 1, 0, 0.0);
                        gsl_matrix_set(J, 2, 0, 0.0);
                        tmp = ls*cos(psi)*ts/r/r;
                        gsl_matrix_set(J, 3, 0, tmp);
                        
                        /* z - column */
                        gsl_matrix_set(J, 0, 1, 0.0);
                        gsl_matrix_set(J, 1, 1, 0.0);
                        tmp = (-1.0)*ls*rho/ts;
                        gsl_matrix_set(J, 2, 1, tmp);
                        gsl_matrix_set(J, 3, 1, 0.0);
                        
                        /* psi - column */
                        tmp = (-1.0)*ls*sin(psi);
                        gsl_matrix_set(J, 0, 2, tmp);
                        tmp = ls*cos(psi);
                        gsl_matrix_set(J, 1, 2, tmp);
                        gsl_matrix_set(J, 2, 2, 0.0);
                        tmp = ls*sin(psi)*ts/r;
                        gsl_matrix_set(J, 3, 2, tmp);
                        
                        /* tau_s - column */
                        tmp = theta*lphi*cos(psi);
                        gsl_matrix_set(J, 0, 3, tmp);
                        tmp = theta*lphi*sin(psi);
                        gsl_matrix_set(J, 1, 3, tmp);
                        tmp = (-1.0)*ls*(p-rho*z)/ts/ts + theta*lphi*(p-rho*z)/ts;
                        gsl_matrix_set(J, 2, 3, tmp);
                        tmp = (-1.0)*ls*cos(psi)/r-theta*lphi*cos(psi)*ts/r;
                        gsl_matrix_set(J, 3, 3, tmp);
                }
                else
                {
                        /* non wrinkling domain */
                        
                        /* r - column */
                        tmp = cos(psi) * ls_dr;
                        gsl_matrix_set(J, 0, 0, tmp);
                        tmp = sin(psi) * ls_dr;
                        gsl_matrix_set(J, 1, 0, tmp);
                        tmp = (1.0/ts)*(ls_dr*(p-rho*z) - ls_dr*kphi*tphi - kphi_dr*ls*tphi - tphi_dr*ls*kphi);
                        gsl_matrix_set(J, 2, 0, tmp);
                        tmp = (1.0/r)*(ls_dr*cos(psi)*(tphi-ts) + ls*cos(psi)*tphi_dr) + (1.0/r/r)*ls*cos(psi)*(ts-tphi);
                        gsl_matrix_set(J, 3, 0, tmp);
                        
                        /* z - column */
                        gsl_matrix_set(J, 0, 1, 0.0);
                        gsl_matrix_set(J, 1, 1, 0.0);
                        tmp = ls*rho/ts;
                        gsl_matrix_set(J, 2, 1, tmp);
                        gsl_matrix_set(J, 3, 1, 0.0);
                        
                        /* psi - column */
                        tmp = (-1.0)*ls*sin(psi);
                        gsl_matrix_set(J, 0, 2, tmp);
                        tmp = ls*cos(psi);
                        gsl_matrix_set(J, 1, 2, tmp);
                        tmp = (-1.0)*ls*tphi*cos(psi)/ts/r;
                        gsl_matrix_set(J, 2, 2, tmp);
                        tmp = ls*sin(psi)*(ts-tphi)/r;
                        gsl_matrix_set(J, 3, 2, tmp);
                        
                        /* tau_s - column */
                        tmp = lphi*(1-nu*nu)*cos(psi)/EH0;
                        gsl_matrix_set(J, 0, 3, tmp);
                        tmp = lphi*(1-nu*nu)*sin(psi)/EH0;
                        gsl_matrix_set(J, 1, 3, tmp);
                        tmp = ((1.0/ts)*ls_dts-ls/ts/ts)*(p-rho*z-kphi*tphi) - ls*kphi*tphi_dts/ts;
                        gsl_matrix_set(J, 2, 3, tmp);
                        tmp = (cos(psi)/r)*(tphi*ls_dts+ls*tphi_dts-ts*ls_dts-ls);
                        gsl_matrix_set(J, 3, 3, tmp);
                }
        }
}

/* runge kutta 4 integration step of Laplace-Young equation 
 * output parameter: data */
void LaplaceRungeKutta(double s, double h, LaplaceData *data)
{
        int dim = data->dim;
        
        double *y = data->y;
        double *yrk4 = data->yrk4;
        
        double *k1, *k2, *k3, *k4;
        double *y2, *y3, *y4;
        k1 = data->k1;
        k2 = data->k2;
        k3 = data->k3;
        k4 = data->k4;
        
        y2 = data->y2;
        y3 = data->y3;
        y4 = data->y4;
        
        LaplaceEquation(s, y, k1, data);
        for (int i = 0; i < dim; i++) y2[i] = y[i] + 0.5*h*k1[i];
        LaplaceEquation(s + 0.5*h, y2, k2, data);
        for (int i = 0; i < dim; i++) y3[i] = y[i] + 0.5*h*k2[i];
        LaplaceEquation(s + 0.5*h, y3, k3, data);
        for (int i = 0; i < dim; i++) y4[i] = y[i] + h*k3[i];
        LaplaceEquation(s + h, y4, k4, data);
        
        for (int i = 0; i < dim; i++) yrk4[i] = y[i] + h*(1.0/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
        for (int i = 0; i < dim; i++) y[i] = yrk4[i];
} 

/* runge kutta 4 integration step of Hooke equation 
 * output parameter: data */
void HookeRungeKutta(double s, double h, HookeData *data)
{
        int dim = data->dim;
        
        double *y = data->y;
        double *y_ext = data->y_ext;
        double *yrk4 = data->yrk4;
        
        double *k1, *k2, *k3, *k4;
        double *y2, *y3, *y4;
        k1 = data->k1;
        k2 = data->k2;
        k3 = data->k3;
        k4 = data->k4;
        
        y2 = data->y2;
        y3 = data->y3;
        y4 = data->y4;       
        ShapeEquation(s, y, y_ext, k1, data);
        for (int i = 0; i < dim; i++) y2[i] = y[i] + 0.5*h*k1[i];
        ShapeEquation(s + 0.5*h, y2, y_ext, k2, data);
        for (int i = 0; i < dim; i++) y3[i] = y[i] + 0.5*h*k2[i];
        ShapeEquation(s + 0.5*h, y3, y_ext, k3, data);
        for (int i = 0; i < dim; i++) y4[i] = y[i] + h*k3[i];
        ShapeEquation(s + h, y4, y_ext, k4, data);
        
        for (int i = 0; i < dim; i++) yrk4[i] = y[i] + h*(1.0/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
        for (int i = 0; i < dim; i++) y[i] = yrk4[i];
}

/* implicit runge kutta step using fixed point iteration to determine slopes
 * output parameter : data */
void ImplicitHookeRungeKuttaFixpoint(double s, double h, HookeData *data)
{
        int dim = data->dim;
        int ord = 4;
        double *y = data->y;
        double *y_ext = data->y_ext;
        gsl_matrix *JP = data->JP;
        gsl_matrix *J = data->J;
        gsl_vector *X = data->X;
        gsl_vector *Y = data->Y;
        gsl_vector *K = data->K;
        gsl_vector *F = data->F;
        gsl_vector *BF = data->BF;
        gsl_vector *R = data->R;
        bool converged = false;
        
        /* initialize */
        for (int i = 0; i < ord; i++)
        {
                for (int j = 0; j < dim; j++)
                {
                        gsl_vector_set(Y, i*dim + j, y[j]);
                        gsl_vector_set(K, i*dim + j, y[j] + h);
                }
        }
        
        int it = 0;
        /* iteration */
        while (!converged)
        {
                /* set F */
                for (int i = 0; i < ord; i++)
                {
                        double f[4];
                        double k[4];
                        for (int j = 0; j < dim; j++)
                        {
                                k[j] = gsl_vector_get(K, i*dim + j);
                        }
                        HookeEquation(s + _A[i] * h, k, y_ext, f, data);
                        for (int j = 0; j < dim; j++)
                        {
                                gsl_vector_set(F, i*dim + j, f[j]);
                        }
                        HookeJacobi(s + _A[i] * h, k, y_ext, JP, data);
                        for (int j = 0; j < dim; j++)
                        {
                                for (int l = 0; l < dim; l++)
                                {
                                        gsl_matrix_set(J, i*dim + j, l, h * gsl_matrix_get(JP, j, l));
                                }
                        }
                }
                /* set BF */
                for (int i = 0; i < ord; i++)
                {
                        for (int j = 0; j < dim; j++)
                        {
                                double tmp = 0.0;
                                for (int l = 0; l < dim; l++)
                                {
                                        tmp += h * _B[i][l] * gsl_vector_get(F, l*dim + j);
                                }
                                gsl_vector_set(BF, i*dim + j, tmp);
                        }
                }
                /* set R */
                for (int i = 0; i < ord; i++)
                {
                        for (int j = 0; j < dim; j++)
                        {
                                gsl_vector_set(R, i*dim + j, gsl_vector_get(Y, i*dim + j) + gsl_vector_get(BF, i*dim + j));
                                gsl_vector_set(X, i*dim + j, gsl_vector_get(R, i*dim + j) - gsl_vector_get(K, i*dim + j));
                        }
                }
                /* set K */
                for (int i = 0; i < ord; i++)
                {
                        for (int j = 0; j < dim; j++)
                        {
                                gsl_vector_set(K, i*dim + j, gsl_vector_get(R, i*dim + j));
                        }
                }
                if (gsl_blas_dnrm2(X) < EPS_IMPLICIT_RK) converged = true;
        }
        
        /* converged, update */
        double tmp[4];
        for (int i = 0; i < dim; i++)
        {
                tmp[i] = 0.0;
                for (int j = 0; j < dim; j++)
                {
                        tmp[i] += h * _C[j] * gsl_vector_get(F, j*dim + i);
                }
        }
        for (int i = 0; i < dim; i++)
        {
                y[i] += tmp[i];
        }
}

/* implicit runge kutta step using minimization to determine slopes
 * output parameter : data */
void ImplicitHookeRungeKutta(double s, double h, HookeData *data)
{
        int dim = data->dim;
        int ord = 4;
        double *y = data->y;
        double *y_ext = data->y_ext;
        gsl_matrix *JP = data->JP;
        gsl_matrix *J = data->J;
        gsl_vector *X = data->X;
        gsl_vector *Y = data->Y;
        gsl_vector *K = data->K;
        gsl_vector *F = data->F;
        gsl_vector *BF = data->BF;
        gsl_vector *R = data->R;
        
        bool converged = false;
        
        /* set Y, K, F */
        for (int i = 0; i < ord; i++)
        {
                for (int j = 0; j < dim; j++)
                {
                        gsl_vector_set(Y, i*dim + j, y[j]);
                        gsl_vector_set(K, i*dim + j, y[j] + h);
                }
        }
        int it = 0;
        
        while (!converged)
        {
                it++;
                /* set F, J */
                for (int i = 0; i < ord; i++)
                {
                        double f[4];
                        double k[4];
                        for (int j = 0; j < dim; j++)
                        {
                                k[j] = gsl_vector_get(K, i*dim + j);
                        }
                        HookeEquation(s + _A[i] * h, k, y_ext, f, data);
                        for (int j = 0; j < dim; j++)
                        {
                                gsl_vector_set(F, i*dim + j, f[j]);
                        }
                        HookeJacobi(s + _A[i] * h, k, y_ext, JP, data);
                        for (int j = 0; j < dim; j++)
                        {
                                for (int l = 0; l < dim; l++)
                                {
                                        gsl_matrix_set(J, i*dim + j, l, h * gsl_matrix_get(JP, j, l));
                                }
                        }
                }
                
                /* set BF */
                for (int i = 0; i < ord; i++)
                {
                        for (int j = 0; j < dim; j++)
                        {
                                double tmp = 0.0;
                                for (int l = 0; l < dim; l++)
                                {
                                        tmp += h * _B[i][l] * gsl_vector_get(F, l*dim + j);
                                }
                                gsl_vector_set(BF, i*dim + j, tmp);
                        }
                }
                
                /* set R */
                for (int i = 0; i < ord; i++)
                {
                        for (int j = 0; j < dim; j++)
                        {
                                gsl_vector_set(R, i*dim + j, gsl_vector_get(K, i*dim + j) - gsl_vector_get(Y, i*dim + j) - gsl_vector_get(BF, i*dim + j));
                        }
                }
                
                gsl_vector_scale(R, -(1.0));
                
                QR(J, X, R, ord * dim, dim);
                
                bool continue_linesearch = true;
                double error_prev, error_next;
                error_prev = gsl_blas_dnrm2(R);
                
                while (continue_linesearch)
                {
                        if (gsl_blas_dnrm2(X) < EPS_IMPLICIT_RK)
                        {
                                continue_linesearch = false;
                                converged = true;
                                break;
                        }
                        for (int i = 0; i < ord; i++)
                        {
                                for (int j = 0; j < dim; j++)
                                {
                                        gsl_vector_set(R, i*dim + j, gsl_vector_get(K, i*dim + j) + gsl_vector_get(X, i*dim + j) - gsl_vector_get(Y, i*dim + j) - gsl_vector_get(BF, i*dim + j));
                                }
                        }
                        error_next = gsl_blas_dnrm2(R);
                        if (error_next < error_prev)
                        {
                                continue_linesearch = false;
                        }
                        else
                        {
                                gsl_vector_scale(X, 0.5);
                        }
                }
                
                if (!converged)
                {
                        gsl_vector_add(K, X);
                }
                else
                {
                        double tmp[4];
                        for (int i = 0; i < dim; i++)
                        {
                                tmp[i] = 0.0;
                                for (int j = 0; j < dim; j++)
                                {
                                        tmp[i] += h * _C[j] * gsl_vector_get(F, j*dim + i);
                                }
                        }
                        for (int i = 0; i < dim; i++)
                        {
                                y[i] += tmp[i];
                        }
                }
        }
}

/* find complete solution of Laplace-Young equation and store in data */
void SolveLaplace(LaplaceData *data)
{
        /* no shooting method applied here 
         * searching instead for convex shape with width larger than capillary diameter 
         * thus needs to be adapted for smaller shapes */
        int dim = data->dim;
        for (int i = 0; i < dim; i++) data->y[i] = 0.0;
        double *y = data->y;
        double t = 0.0;
        double tMax = 10.0;
        double h = INTEGRATION_STEP_LAPLACE;
        int count_bc_crossing_from_left = 0;
        int count_bc_crossing_from_right = 0;
        double r0_before_step = 0.0;
        double r0_after_step = 0.0;
        double z0_before_step = 0.0;
        double z0_after_step = 0.0;
        double psi0_before_step = 0.0;
        double psi0_after_step = 0.0;
        double psi_dev = 0.0;
        double r_dev = 0.0;
        data->invalid = false;
        
        /* store initial conditions */
        data->Clear();
        data->Push(t, y[0], y[0]*y[0], y[1], y[2]);
        
        double a = data->parameters[2];
        
        /* scaling */
        tMax *= a;
        h *= a;
        
        while (t < tMax)
        {
                r0_before_step = y[0];
                z0_before_step = y[1];
                psi0_before_step = y[2];
                
                LaplaceRungeKutta(t, h, data);
                /* step applied -> update s0 */
                t += h;
                
                r0_after_step = y[0];
                z0_after_step = y[1];
                psi0_after_step = y[2];
                
                /* derivatives */
                r_dev = (r0_after_step - r0_before_step) / h;
                psi_dev = (psi0_after_step - psi0_before_step) / h;
                
                /* check for invalid solution */
                if (!( ((y[2] < M_PI/2) || (y[0] > 0.1)) && ((y[2] > -M_PI/2) || (y[0] > 0.1)) && (fabs(psi_dev) < 50) && (y[2] < 2*M_PI) ) )
                {
                        data->invalid = true;
                        data->Clear();
                        return;
                }
                
                /* store current solution */
                data->Push(t, y[0], y[0]*y[0], y[1], y[2]);
                
                /* check boundary conditions */
                if ((r0_before_step < a * 0.5) && (r0_after_step > a * 0.5))
                {
                        count_bc_crossing_from_left++;
                }
                if ((r0_before_step > a * 0.5) && (r0_after_step < a * 0.5))
                {
                        count_bc_crossing_from_right++;
                }
                
                /* matching boundary condition? */
                if (count_bc_crossing_from_right == 1)
                {
                        /* cubic interpolation */
                        InterpolateLaplace(data);
                        /* calculate L0 via splines */
                        double h0 = 1e-08;
                        /* backtrace through spline with higher precision till boundary condition is reached exactly */
                        while (gsl_spline_eval(data->splines.r_s, t, data->splines.acc) < a * 0.5) t -= h0;
                        data->L0 = t;
                        data->invalid = false;
                        
                        /* translate theory -> z(L0) == z(capillary) */
                        MoveLaplaceDown(data);
                        return;
                }
                
        }
        /* no boundary match -> fail */
        data->invalid = true;
        data->Clear();
        return;
}

/* plot error function of single shooting and determine zero crossings
 * debug purposes only */
void ErrorFunctionSingleShooting(HookeData *data)
{
        double p = data->parameters[0];
        double nu = data->parameters[1];
        double k = data->parameters[2];
        double pref1, pref2, pref3;
        pref1 = data->undeformed->parameters[0];
        pref2 = data->undeformed->parameters[1];
        pref3 = data->undeformed->parameters[2];
        
        HookeData f;
        LaplaceData l;
        l.SetParameters(pref1, pref2, pref3);
        SolveLaplace(&l);
        f.SetLaplace(&l);
        f.SetParameters(p, nu, k);
        
        PointSet func;
        
        fstream errfunc("errfuncsingle.dat", ios::out);
        
        for (double t = 1e-2; t <= 2.0; t += 1e-4)
        {
                f.SetInitialConditions(0.0, 0.0, 0.0, t);
                SolveHooke(&f);
                double d = f.y[0] - 0.5;
                errfunc << FMT << t << FMT << d << endl;
                func.push_back(make_pair(t, d));
        }
        
        /* determine zero crossings */
        for (int i = 1; i < func.size() - 1; i++)
        {
                if (((func[i-1].second > 0) && (func[i+1].second < 0)) || ((func[i-1].second < 0) && (func[i+1].second > 0)))
                {
                        cout << "zero-crossing at " << func[i].first << flush << endl;
                }
        }
        
        errfunc.close();
}

/* find complete solution of Hooke equation and store in data */
bool SingleShooting(HookeData *data)
{
        return SingleShooting(data, data->undeformed->L0, 0.5);
}

/* applies single shooting method
 * adapt initial value tau_s(0) until boundary deviation sufficient small 
 * performs improved bisection (subdivision in multiple intervals) to find valid solution */
bool SingleShooting(HookeData *data, double t_max, double target)
{
        bool done = false;
        bool success = false;
        if (WATCH_SINGLE_SHOOTING) cout << "Single shooting.." << endl;
        
        double p = data->parameters[0];
        double nu = data->parameters[1];
        double k = data->parameters[2];
        double pref1, pref2, pref3;
        pref1 = data->undeformed->parameters[0];
        pref2 = data->undeformed->parameters[1];
        pref3 = data->undeformed->parameters[2];
        
        if (WATCH_SINGLE_SHOOTING) 
        {
                cout << "Single shooting laplace-parameter p = " << pref1 << " rho = " << pref2 << " alpha = " << pref3 << endl;
                cout << "Single shooting hooke-parameter: p = " << p << " nu = " << nu << " k = " << k << flush << endl;
        }
        
        int m = 3;
        int oc = 0;
        
        HookeData f[m];
        LaplaceData l[m];
        
        double s[m];
        double d[m];
        double d_min;
        double s_min;
	int i_min;
        double a, b;
        a = 1e-2;
        b = 2.0;
        double h;
        h = (b - a) / ((double)m - 1.0);
        while (h > 1e-15)
        {
                h = (b - a) / ((double)m - 1.0);
                if (WATCH_SINGLE_SHOOTING) cout << h << flush << endl;
                for (int i = 0; i < m; i++)
                {
                        s[i] = a + i * h;
                }
                
                for (int i = 0; i < m; i++)
                {
                        l[i].SetParameters(pref1, pref2, pref3);
                        SolveLaplace(&l[i]);
                        f[i].SetLaplace(&l[i]);
                        f[i].SetParameters(p, nu, k);
                        f[i].SetInitialConditions(0.0, 0.0, 0.0, s[i]);
                        SolveHooke(&f[i], 0.0, t_max);
                        if (f[i].invalid)
                        {
                                d[i] = -1.0;
                        }
                        else
                        {			
                                d[i] = f[i].y[0] - target;
                        }
                }
                
                /* find minimal boundary deviation */
                d_min = fabs(d[0]);
                s_min = s[0];
                for (int i = 1; i < m; i++)
                {
                        if (fabs(d[i]) < d_min)
                        {
                                d_min = fabs(d[i]);
                                s_min = s[i];
				i_min = i;
                        }
                }
                
                if (WATCH_SINGLE_SHOOTING)
                {
                        cout << "Single shooting Current Interval: " << setw(15) << a << setw(15) << b << endl;
                        cout << "Single shooting Evaluation at     ";
                        for (int i = 0; i < m; i++) cout << setw(15) << s[i];
                        cout << endl;
                        cout << "Single shooting Errorfunction     ";
                        for (int i = 0; i < m; i++) cout << setw(15) << d[i];
                        cout << endl;
                        cout << "Single shooting Minimum           " << setw(15) << d_min << endl;
                        cout << endl;
                }
                
                /*
                if (d_min < EPS_SINGLE_SHOOTING) {
			if (WATCH_SINGLE_SHOOTING) cout << "accuracy reached... single shooting converged." << endl;
			break;
		}*/
                
                bool zero_crossing = false;
                bool interval_changed = true;
                /* search interval where deviation crosses zero */
                for (int i = 0; i < m - 1; i++)
                {
                        if ((d[i] <= 0 && d[i+1] >= 0) || (d[i] >= 0 && d[i+1] <= 0)) 
                        {
                                zero_crossing = true;
                                // if ((a == s[i]) && (b == s[i+1]))
                                if (fabs(a - s[i]) < 1.0e-15 && fabs(b - s[i+1]) < 1.0e-15)
                                {
                                        interval_changed = false;
                                }
                                a = s[i];
                                b = s[i+1];
                                break;
                        }
                }
                
                if (!zero_crossing)
                {
                        if (WATCH_SINGLE_SHOOTING) cout << "no zero crossing found in single shooting.. boundary deviation " << d_min << flush << endl; 
                        s_min = (a+b)/2.0;
                        break;
                }
                
                if (!interval_changed)
                {
                        if (WATCH_SINGLE_SHOOTING) cout << "interval did not change in single shooting.. boundary deviation " << d_min << flush << endl;
                        s_min = (a+b)/2.0;
                        break;
                }
                
        }
        data->SetInitialConditions(0.0, 0.0, 0.0, s_min);
        data->mu = s_min;
        SolveHooke(data, 0.0, t_max);
        if (data->invalid)
        {
                if (WATCH_SINGLE_SHOOTING) cout << "invalid.." << flush << endl;
                return false;
        }
        double boundary_deviation = fabs(data->y[0] - target);
        if (boundary_deviation > 0.5) 
        {
                if (WATCH_SINGLE_SHOOTING) cout << "bd not reached.." << flush << endl;
                return false;
        }
        
        /* check if required accuracy reached */
        if (boundary_deviation < EPS_SINGLE_SHOOTING) 
        {
                if (WATCH_SINGLE_SHOOTING) cout << "Single shooting converged.." << flush << endl;
                success = true;
        }
        else
        {
                /* improve accuracy with multiple shooting method if required */
                if (EXTENDED_SHOOTING)
                {
                        if (WATCH_SINGLE_SHOOTING) cout << "Single shooting (Boundary Deviation = " << boundary_deviation << ") -> Extended shooting.." << flush << endl;
                        return ParallelShooting(data, 4);
                }
        }
        if (!EXTENDED_SHOOTING)
        {
                if (!success)
                {
                        if (WATCH_SINGLE_SHOOTING) cout << "Single shooting incomplete.. Boundary Deviation = " << boundary_deviation << flush << endl;
                        return false;
                }
        }
        MoveHookeDown(data);
        if (WATCH_SINGLE_SHOOTING) cout << "Single shooting successful.." << flush << endl;
        return success;
}

/* calculate numerical differential quotient of Hooke equation with respect to parameter param for use in jacobian 
 * output parameter : dy */
void FiniteDifference(HookeData *dummy, double s1, double s2, int param, double *y, double *dy)
{
        double dplus[4], dminus[4];
        switch (param)
        {
                case 1:
                        dummy->SetInitialConditions(y[0] + DX_SHOOTING, y[1], y[2], y[3]);
                        break;
                case 2:
                        dummy->SetInitialConditions(y[0], y[1] + DX_SHOOTING, y[2], y[3]);
                        break;
                case 3:
                        dummy->SetInitialConditions(y[0], y[1], y[2] + DX_SHOOTING, y[3]);
                        break;
                case 4:
                        dummy->SetInitialConditions(y[0], y[1], y[2], y[3] + DX_SHOOTING);
                        break;
                default:
                        break;
        }
        SolveHooke(dummy, s1, s2);
        for (int k = 0; k < 4; k++) dplus[k] = dummy->y[k];
        switch (param)
        {
                case 1:
                        dummy->SetInitialConditions(y[0] - DX_SHOOTING, y[1], y[2], y[3]);
                        break;
                case 2:
                        dummy->SetInitialConditions(y[0], y[1] - DX_SHOOTING, y[2], y[3]);
                        break;
                case 3:
                        dummy->SetInitialConditions(y[0], y[1], y[2] - DX_SHOOTING, y[3]);
                        break;
                case 4:
                        dummy->SetInitialConditions(y[0], y[1], y[2], y[3] - DX_SHOOTING);
                        break;
                default:
                        break;
        }
        SolveHooke(dummy, s1, s2);
        for (int k = 0; k < 4; k++) dminus[k] = dummy->y[k];
        for (int k = 0; k < 4; k++) dy[k] = (1.0/(2.0*DX_SHOOTING))*(dplus[k] - dminus[k]);
}

/* method to solve hooke-membrane without single shooting
 * by tracing from reference shape 
 * use only recommended if any shooting method fails to converge 
 * slow and no convergence guaranteed (experimental) 
 * output parameter : data */
bool ParallelShootingTracing(HookeData *data)
{
        double *pi = data->parameters;
        LaplaceData *udef = data->undeformed;
        double root[3] = {udef->parameters[0], 0.5, 5.0};
        double target[3] = {pi[0], pi[1], pi[2]};
        data->SetParameters(root[0], target[1], target[2]);
        data->SetInitialConditions(0.0, 0.0, 0.0, 1.0);
        int N = 100;
        double h[3][N];
        for (int i = 0; i < 3; i++)
        {
                for (int k = 0; k < N; k++)
                {
                        h[i][k] = (target[i] - root[i])/(double)N;
                        h[i][k] *= (2.0)/(double)N;
                        h[i][k] *= ((double)N-((double)k+1));
                }
        }
        
        /* solve undeformed configuration */
        bool success = ParallelShooting(data, 4);
        PlotHooke("0.dat", data);
        if (!success) return false;
        data->SetParameters(root[0] + h[0][0], root[1], root[2]);
        success = ParallelShooting(data, 4);
        PlotHooke("1.dat", data);
        if (!success) return false;
        double h_ts = data->tau_s[0] - 1.0;
        
        for (int i = 2; i <= N; i++)
        {
                stringstream ss;
                ss << i;
                string filename = ss.str();
                filename.append(".dat");
                cout << i << ": " << data->tau_s[0] << flush << endl;
                data->SetInitialConditions(0.0, 0.0, 0.0, data->tau_s[0]);
                data->SetParameters(root[0] + i * h[0][i-1], target[1], target[2]);
                success = ParallelShooting(data, 4);
                PlotHooke(filename.c_str(), data);
                if (!success) return false;
        }
}

/* multiple shooting method
 * output parameter : data */
bool ParallelShooting(HookeData *data, int p_interval)
{
        double *pi = data->parameters;
        LaplaceData *udef = data->undeformed;
        double pref1 = udef->parameters[0];
        double pref2 = udef->parameters[1];
        double pref3 = udef->parameters[2];
        double L0 = udef->L0;
        
        HookeData dummy;
        dummy.SetLaplace(udef);
        dummy.SetParameters(pi[0], pi[1], pi[2]);
        
        dummy.SetInitialConditions(0.0, 0.0, 0.0, data->mu);
        SolveHooke(&dummy);
        
        if (dummy.invalid) 
        {
                cout << "invalid dummy in ParallelShooting.." << flush << endl;
                return false;
        }
        
        if (WATCH_MULTI_SHOOTING) 
        {
                cout << "Parallel shooting with " << p_interval << " intervals" << endl;
                cout << "Parallel shooting laplace-parameter p = " << udef->parameters[0] << " rho = " << udef->parameters[1] << " alpha = " << udef->parameters[2] << endl;
                cout << "Parallel shooting hooke-parameter: p = " << pi[0] << " nu = " << pi[1] << " k = " << pi[2] << flush << endl;
        }
        
        /* number intervals will be increased if multiple shooting not successful */
        int p = p_interval;
        /* number gridpoints */
        int m = p + 1;
        int n = 4;
        
        HookeData f[m-3];
        LaplaceData l[m-3];
        
        double *s = (double*)calloc(m, sizeof(double));
        double h = udef->L0 / (double)p;
        for (int i = 0; i < m; i++) s[i] = (double)i * h;
        double **y = (double**)calloc(m, sizeof(double*));
        for (int i = 0; i < m; i++) y[i] = (double*)calloc(4, sizeof(double));
        double **yy = (double**)calloc(m, sizeof(double*));
        for (int i = 0; i < m; i++) yy[i] = (double*)calloc(4, sizeof(double));
        
        /* initialize gridpoints */
        y[0][0] = 0.0;
        y[0][1] = 0.0;
        y[0][2] = 0.0;
        y[0][3] = gsl_spline_eval(dummy.splines.tau_s, s[0], dummy.splines.acc);
        for (int i = 1; i < m-1; i++)
        {
                y[i][0] = gsl_spline_eval(dummy.splines.r_s, s[i], dummy.splines.acc);
                y[i][1] = gsl_spline_eval(dummy.splines.z_s, s[i], dummy.splines.acc);
                y[i][2] = gsl_spline_eval(dummy.splines.psi_s, s[i], dummy.splines.acc);
                y[i][3] = gsl_spline_eval(dummy.splines.tau_s, s[i], dummy.splines.acc);
        }
        y[m-1][0] = 0.5;
        y[m-1][1] = gsl_spline_eval(dummy.splines.z_s, s[m-1], dummy.splines.acc);
        y[m-1][2] = gsl_spline_eval(dummy.splines.psi_s, s[m-1], dummy.splines.acc);
        y[m-1][3] = gsl_spline_eval(dummy.splines.tau_s, s[m-1], dummy.splines.acc);
        
        /* prepare method */
        double dy[4];
        double error_prev, error_next;
        bool continue_linesearch;
        bool converged = false;
        
        int it = 0;
        int M = p * n - 3;
        int N = p * n - 3;
        gsl_matrix *J = gsl_matrix_calloc(M, N); /* jacobian */
        gsl_vector *F = gsl_vector_calloc(M); /* residual */
        gsl_vector *D = gsl_vector_calloc(M); /* parameter shift */
        
        
        while (!converged)
        {
                it++;
                /* set up jacobian */
                dummy.SetInitialConditions(y[0][0], y[0][1], y[0][2], y[0][3]);
                SolveHooke(&dummy, s[0], s[1]);
                
                for (int i = 0; i < n; i++) gsl_vector_set(F, i, dummy.y[i] - y[1][i]);
                
                FiniteDifference(&dummy, s[0], s[1], 4, y[0], dy);
                
                for (int i = 0; i < n; i++) gsl_matrix_set(J, i, 0, dy[i]);
                for (int i = 0; i < n; i++) gsl_matrix_set(J, i, i + 1, -1.0);
                
                /* compute segements in parallel */
                //#pragma omp parallel for private(i) shared(f, l) num_threads(4)
                for (int i = 1; i < m-2; i++)
                {
                        double loc_dy[4];
                        int ii = i-1;
                        l[ii].SetParameters(pref1, pref2, pref3);
                        SolveLaplace(&l[ii]);
                        f[ii].SetLaplace(&l[ii]);
                        f[ii].SetParameters(pi[0], pi[1], pi[2]);
                        f[ii].SetInitialConditions(y[i][0], y[i][1], y[i][2], y[i][3]);
                        SolveHooke(&f[ii], s[i], s[i+1]);
                        
                        for (int k = 0; k < n; k++) gsl_vector_set(F, i*n + k, f[ii].y[k] - y[i+1][k]);
                        
                        for (int pindex = 0; pindex < 4; pindex++)
                        {
                                FiniteDifference(&f[ii], s[i], s[i+1], pindex + 1, y[i], loc_dy);
                                for (int k = i*n; k < i*n+4; k++) gsl_matrix_set(J, k, 1 + (i-1)*n + pindex, loc_dy[k-i*n]);
                        }    
                        
                        for (int k = 0; k < n; k++) 
                        {
                                gsl_matrix_set(J, i*n + k, 1 + i*n + k, -1.0);
                        }
                }
                
                dummy.SetInitialConditions(y[m-2][0], y[m-2][1], y[m-2][2], y[m-2][3]);
                SolveHooke(&dummy, s[m-2], s[m-1]);
                gsl_vector_set(F, M-1, dummy.y[0] - y[m-1][0]);
                
                for (int k = 1; k <= n; k++)
                {
                        FiniteDifference(&dummy, s[m-2], s[m-1], k, y[m-1], dy);
                        gsl_matrix_set(J, M-1, N-(5-k), dy[0]);
                }
                
                gsl_vector_scale(F, -1.0);
                QR(J, D, F, M, N);
                
                /* linesearch */
                continue_linesearch = true;
                error_prev = gsl_blas_dnrm2(F);
                int cit = 0;
                while (continue_linesearch)	/* begin linesearch */
                {
                        yy[0][3] = y[0][3] + gsl_vector_get(D, 0);
                        for (int i = 1; i < m-1; i++)
                        {
                                for (int k = 0; k < n; k++) yy[i][k] = y[i][k] + gsl_vector_get(D, (i-1)*n + (k+1));
                        }
                        
                        dummy.SetInitialConditions(yy[0][0], yy[0][1], yy[0][2], yy[0][3]);
                        SolveHooke(&dummy, s[0], s[1]);
                        
                        for (int i = 0; i < n; i++) gsl_vector_set(F, i, dummy.y[i] - y[1][i]);
                        
                        for (int i = 1; i < m-2; i++)
                        {
                                dummy.SetInitialConditions(yy[i][0], yy[i][1], yy[i][2], yy[i][3]);
                                SolveHooke(&dummy, s[i], s[i+1]);
                                for (int k = 0; k < n; k++) gsl_vector_set(F, i*n + k, dummy.y[k] - y[i+1][k]);
                        }
                        
                        dummy.SetInitialConditions(yy[m-2][0], yy[m-2][1], yy[m-2][2], yy[m-2][3]);
                        
                        SolveHooke(&dummy, s[m-2], s[m-1]);
                        gsl_vector_set(F, M-1, dummy.y[0] - y[m-1][0]);
                        
                        error_next = gsl_blas_dnrm2(F);
                        
                        if (error_next < error_prev)
                        {
                                continue_linesearch = false;
                        }
                        else
                        {
                                gsl_vector_scale(D, 0.1);
                        }
                        
                        if (gsl_blas_dnrm2(D) <= 1e-16) 
                        {
                                continue_linesearch = false;
                                /* ITERATION RUNS IN CIRCLE if still converged == false */
                                /* converged = true; */
                                /* ??????????? */
                                if (WATCH_MULTI_SHOOTING) cout << "Linesearch not successful.." << flush << endl;
                                if (p_interval < 16) 
                                {
                                        return ParallelShooting(data, p + 4);
                                }
                                else 
                                {
                                        return false;
                                }
                        }
                        
                } /* end line search */
                
                /* update grid points */
                y[0][3] += gsl_vector_get(D, 0);
                for (int i = 1; i < m-1; i++)
                {
                        for (int k = 0; k < n; k++) y[i][k] += gsl_vector_get(D, (i-1)*n + (k+1));
                }
                
                if (WATCH_MULTI_SHOOTING) cout << "Parallel shooting Error" << setw(20) << gsl_blas_dnrm2(F) << endl;
                
                if (gsl_blas_dnrm2(F) <= EPS_MULTI_SHOOTING) converged = true;
                
                if (WATCH_MULTI_SHOOTING) 
                {
                        if (converged) 
                        {
                                if (WATCH_MULTI_SHOOTING) cout << "Parallel shooting converged.." << endl;
                                break;
                        }
                }
                
                /* check if error has increased
                 * perform 100 iterations on max */
                if ((error_next >= error_prev) || (it > 100))
                {
                        if (WATCH_MULTI_SHOOTING) cout << "Error increased.." << flush << endl;
                        if (p_interval < 16) 
                                return ParallelShooting(data, p + 4);
                        else 
                                return false;
                }
                
        }/* end iteration */
        
        gsl_matrix_free(J);
        gsl_vector_free(D);
        gsl_vector_free(F);
        
        dummy.SetInitialConditions(y[m-2][0], y[m-2][1], y[m-2][2], y[m-2][3]);
        SolveHooke(&dummy, s[m-2], s[m-1]);
        for (int i = 0; i < n; i++)
        {
                y[m-1][i] = dummy.y[i];
        }
        
        /* updating data */
        data->Clear();
        data->SetInitialConditions(y[0][0], y[0][1], y[0][2], y[0][3]);
        data->L0 = L0;
        /* merge segments to complete solution 
         * should be continuous here */
        for (int i = 0; i < p; i++)
        {
                dummy.SetInitialConditions(y[i][0], y[i][1], y[i][2], y[i][3]);
                SolveHooke(&dummy, s[i], s[i+1]);
                int max, min;
                if (i == 0)
                {
                        min = 0;
                        max = dummy.size;
                }
                else
                {
                        min = 1;
                        max = dummy.size;
                }
                for (int i = min; i < max; i++)
                {
                        data->Push(dummy.s[i], dummy.r[i], dummy.r2[i], dummy.z[i], dummy.psi[i], dummy.tau_s[i], dummy.lambda_s[i], dummy.lambda_phi[i], dummy.kappa_phi[i], dummy.tau_phi[i]);
                }
        }
        data->invalid = false;
        
        /* interpolate solution and align to horizontal axis */
        InterpolateHooke(data);
        MoveHookeDown(data);
        
        if (WATCH_MULTI_SHOOTING)
        {
                cout << "Boundary Deviation: " << fabs(data->r[data->size-1] - 0.5) << flush << endl;
        }
        
        /* clean up */
        for (int i = 0; i < m; i++) 
        {
                free(y[i]);
                free(yy[i]);
        }
        free(y);
        free(yy);
        free(s);
        if (WATCH_MULTI_SHOOTING) cout << "Parallel shooting successful.." << flush << endl;
        return true;
}

/* find solution to Hooke equation and store in data
 * solution is not checked for boundary condition here 
 * function used in shooting method */
void SolveHooke(HookeData *data)
{
// 	gsl_vector_set(data->min_trace, 0, 1.0);
// 	gsl_vector_set(data->min_trace, 1, 1.0);
        SolveHooke(data, 0.0, data->undeformed->L0); 
}

/* find solution to Hooke equation and store in data
 * solution is not checked for boundary condition here 
 * function used in shooting method 
 * also used for integration of single segments */
void SolveHooke(HookeData *data, double t_min, double t_max)
{
        int dim = data->dim;
        double *y = data->y;
        double *y_ext = data->y_ext;
        
        data->Clear();
        double a = data->undeformed->parameters[2];
        
        double t = t_min; 
        double h = INTEGRATION_STEP_HOOKE;
        data->L0 = data->undeformed->L0;
        double z_before_step = 0.0;
        double z_after_step = 0.0;
        double psi_before_step = 0.0; 
        double psi_after_step = 0.0;
        double z_dev = 0.0;
        double psi_dev = 0.0;
        int count_steps = (int)((t_max - t_min) / h + 0.5);
        
        double y0, y1, y2, y3;
        y0 = y[0];
        y1 = y[1];
        y2 = y[2];
        y3 = y[3];
        
        data->Push(t, y0, y0*y0, y1, y2, y3, y_ext[0], y_ext[1], y_ext[2], y_ext[3]);
        
        double ac = data->L0; /* threshold for valid solution */
        double tmp1, tmp2, tmp3, tmp4;
        
        int it = 0;
	
	data->lambda_s_trace = 1.0;
	data->lambda_phi_trace = 1.0;
	data->tau_phi_trace_1 = 1.0;
	data->tau_phi_trace_2 = 1.0;
	data->integration_mark_0 = false;
	data->integration_mark_1 = false;
	data->integration_mark_2 = false;
	data->integration_mark_3 = false;
	
        while (t < t_max)
        {
                tmp1 = y[0];
                tmp2 = y[1];
                tmp3 = y[2];
                tmp4 = y[3];
                
                z_before_step = y[1];
                psi_before_step = y[2];
                
                if (t_max - t < h) 
                {
                        h = t_max - t;
                }
                else 
                {
                        h = INTEGRATION_STEP_HOOKE;
                }
                
                if (t <= 0.01*data->L0)
                {
                        HookeRungeKutta(t, h, data);
		}
                else
                {
                        if (IMPLICIT_INTEGRATION)
                        {
                                ImplicitHookeRungeKutta(t, h, data);
                        }
                        else
                        {
                                HookeRungeKutta(t, h, data);
                        }
                }
                
                t += h;
                
                z_after_step = y[1];
                psi_after_step = y[2];
                psi_dev = (psi_after_step - psi_before_step) / h;
               
                /* check for valid solution */
                if ((y[0] > 10.0*ac) || (y[1] > 10.0*ac) || (y[0] <= 0.0) || (y[1] < 0.0) || (y[2] < 0.0) || (y[2] > M_PI) || (y[3] < 0.0))
                {
                        data->invalid = true;
                        data->Clear();
                        return;
                } 
                
                if (isnan(y[0])) 
                {
                        cout << "prev\t" << tmp1 << "\t" << tmp2 << "\t" << tmp3 << "\t" << tmp4 << flush << endl;
                        cout << "curr\t" << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << y[3] << flush << endl;
                        cout << "r is nan" << flush << endl;
                        return;
                }
                
                if (isnan(y[1]))
                {
                        cout << "prev\t" << tmp1 << "\t" << tmp2 << "\t" << tmp3 << "\t" << tmp4 << flush << endl;
                        cout << "curr\t" << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << y[3] << flush << endl;
                        cout << "z is nan" << flush << endl;
                        return;
                }
                
                if (isnan(y[2]))
                {
                        cout << "prev\t" << tmp1 << "\t" << tmp2 << "\t" << tmp3 << "\t" << tmp4 << flush << endl;
                        cout << "curr\t" << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << y[3] << flush << endl;
                        cout << "psi is nan" << flush << endl;
                        return;
                }
                
                if (isnan(y[3]))
                {
                        cout << "prev\t" << tmp1 << "\t" << tmp2 << "\t" << tmp3 << "\t" << tmp4 << flush << endl;
                        cout << "curr\t" << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << y[3] << flush << endl;
                        cout << "tau_s is nan" << flush << endl;
                        return;
                }
                
                it++;
                /* store current solution vector */
                if ((it % (count_steps / count_steps) == 0) || (t_max - t <= h))
                {
                        y0 = y[0];
                        y1 = y[1];
                        y2 = y[2];
                        y3 = y[3];
                        
                        data->Push(t, y0, y0*y0, y1, y2, y3, y_ext[0], y_ext[1], y_ext[2], y_ext[3]);
                }
        } 
        /* valid solution found, interpolate and return */
        data->invalid = false;
        InterpolateHooke(data);
}

/* translation methods
 * used for alignment of finished solutions
 * capillary on horizontal axis
 * capsule in lower half space */
void MoveLaplaceUp(LaplaceData *data)
{
        double dz = (-1.0) * gsl_spline_eval(data->splines.z_s, 0.0, data->splines.acc);
        for (int i = 0; i < data->size; i++) data->z[i] += dz;
        InterpolateLaplace(data);
}

void MoveLaplaceDown(LaplaceData *data)
{
        double dz = gsl_spline_eval(data->splines.z_s, data->L0, data->splines.acc);
        for (int i = 0; i < data->size; i++) data->z[i] -= dz;
        InterpolateLaplace(data);
}

void MoveHookeUp(HookeData *data)
{
        double dz = (-1.0) * gsl_spline_eval(data->splines.z_s, 0.0, data->splines.acc);
        for (int i = 0; i < data->size; i++) data->z[i] += dz;
        InterpolateHooke(data);
}

void MoveHookeDown(HookeData *data)
{
        double dz = gsl_spline_eval(data->splines.z_s, data->L0, data->splines.acc);
        for (int i = 0; i < data->size; i++) data->z[i] -= dz;
        InterpolateHooke(data);
}

bool WrinklingRegion(HookeData *data, double &s1, double &s2)
{
        double zero_begin, zero_end;
        double h = 1e-2;
        bool abort = false;
        bool wrinkle_begin = false, wrinkle_end = false;
        double s0 = 0.0;
        while (s0 <= data->L0)
        {
                if ((gsl_spline_eval(data->splines.t_phi, s0, data->splines.acc) > 0.0) && (gsl_spline_eval(data->splines.t_phi, s0 + h, data->splines.acc) < 0.0))
                {
                        wrinkle_begin = true;
                        double s0_backup = s0;
                        double h_loc = 1e-8;
                        bool done = false;
                        while (!done)
                        {
                                s0 += h_loc;
                                if (gsl_spline_eval(data->splines.t_phi, s0, data->splines.acc) < 0.0)
                                {
                                        zero_begin = (2.0*s0 + h_loc)/2.0;
                                        done = true;
                                }
                        }
                        s0 = s0_backup;
                }
                if ((gsl_spline_eval(data->splines.t_phi, s0, data->splines.acc) < 0.0) && (gsl_spline_eval(data->splines.t_phi, s0 + h, data->splines.acc) > 0.0))
                {
                        wrinkle_end = true;
                        double s0_backup = s0;
                        double h_loc = 1e-8;
                        bool done = false;
                        while (!done)
                        {
                                s0 += h_loc;
                                if (gsl_spline_eval(data->splines.t_phi, s0, data->splines.acc) > 0.0)
                                {
                                        zero_end = (2.0*s0 + h_loc)/2.0;
                                        done = true;
                                }
                        }
                        s0 = s0_backup;
                }	
                s0 += h;
                if (abort) break;
        }
        if (wrinkle_begin || wrinkle_end) {
                if (!zero_begin) s1 = 0.0;
                else s1 = zero_begin;
                if (!zero_end) s2 = data->L0;
                else s2 = zero_end;
                return true;
        } else return false;
}

double AverageWrinklingTension(HookeData *data, double s1, double s2)
{
        double avg_tau = 0.0;
        int counter = 0;
        for (double s = s1; s <= s2; s += 1.0e-2) {
                avg_tau += gsl_spline_eval(data->splines.tau_s, s, data->splines.acc);
                counter++;
        }
        avg_tau /= (double)counter;
        return avg_tau;
}

double MaximumRadius(HookeData *data)
{
        double max = 0.0;
        for (int i = 0; i < data->r.size(); i++) 
        {
                if (data->r[i] > max) max = data->r[i];
        }
        return max;
}

/* transform length units of Laplace shape */
void ScaleLaplace(LaplaceData *data)
{
        double scal = data->parameters[2];
        for (int i = 0; i < data->size; i++)
        {
                data->s[i] /= scal;
                data->r[i] /= scal;
                data->r2[i] /= scal * scal;
                data->z[i] /= scal;
        }
        data->L0 /= scal;
        data->parameters[2] = 1.0;
        data->conversion *= scal;
        InterpolateLaplace(data);
}

/* plot Hooke shape, tension, strain and curvature profiles */
void PlotHooke(string filename, HookeData *data)
{	
        ofstream outfile(filename.c_str(), ios::out);
        // header
        outfile << "#" << "L0=" << data->L0 << ",delta_z=" << gsl_spline_eval(data->splines.z_s, 0.0, data->splines.acc) << endl;
        outfile
        << "#"
        << setw(19) << "s"
        << setw(20) << "r"
        << setw(20) << "z"
        << setw(20) << "psi"
        << setw(20) << "tau_s"
        << setw(20) << "lambda_s"
        << setw(20) << "lambda_phi"
        << setw(20) << "kappa_phi"
        << setw(20) << "tau_phi" << endl;
        double s;
        for (s = 0.0; s <= data->L0; s += 1e-4)
        {
                outfile 
                << FMT << s
                << FMT << gsl_spline_eval(data->splines.r_s, s, data->splines.acc)
                << FMT << gsl_spline_eval(data->splines.z_s, s, data->splines.acc)
                << FMT << gsl_spline_eval(data->splines.psi_s, s, data->splines.acc)
                << FMT << gsl_spline_eval(data->splines.tau_s, s, data->splines.acc) 
                << FMT << gsl_spline_eval(data->splines.l_s, s, data->splines.acc)
                << FMT << gsl_spline_eval(data->splines.l_phi, s, data->splines.acc)
                << FMT << gsl_spline_eval(data->splines.k_phi, s, data->splines.acc)
                << FMT << gsl_spline_eval(data->splines.t_phi, s, data->splines.acc) << endl;
        }
        s = data->L0;
        outfile 
        << FMT << s
        << FMT << gsl_spline_eval(data->splines.r_s, s, data->splines.acc)
        << FMT << gsl_spline_eval(data->splines.z_s, s, data->splines.acc)
        << FMT << gsl_spline_eval(data->splines.psi_s, s, data->splines.acc)
        << FMT << gsl_spline_eval(data->splines.tau_s, s, data->splines.acc)
        << FMT << gsl_spline_eval(data->splines.l_s, s, data->splines.acc)
        << FMT << gsl_spline_eval(data->splines.l_phi, s, data->splines.acc)
        << FMT << gsl_spline_eval(data->splines.k_phi, s, data->splines.acc)
        << FMT << gsl_spline_eval(data->splines.t_phi, s, data->splines.acc) << endl;
        outfile.close();
}

/* plot Laplace shape */
void PlotLaplace(string filename, LaplaceData *data)
{
        ofstream outfile(filename.c_str(), ios::out);
        // header
        outfile << "#" << "L0=" << data->L0 << ",delta_z=" << gsl_spline_eval(data->splines.z_s, 0.0, data->splines.acc) << endl;
        outfile
        << "#"
        << setw(19) << "s"
        << setw(20) << "r"
        << setw(20) << "z"
        << setw(20) << "psi" << endl;
        double s;
        for (s = 0.0; s <= data->L0; s += 1e-2)
        {
                outfile 
                << FMT << s 
                << FMT << gsl_spline_eval(data->splines.r_s, s, data->splines.acc) 
                << FMT << gsl_spline_eval(data->splines.z_s, s, data->splines.acc) 
                << FMT << setprecision(15) << gsl_spline_eval(data->splines.psi_s, s, data->splines.acc) << endl;
                s += 1e-2;
        }
        s = data->L0;
        outfile 
        << FMT << s 
        << FMT << gsl_spline_eval(data->splines.r_s, s, data->splines.acc) 
        << FMT << gsl_spline_eval(data->splines.z_s, s, data->splines.acc) 
        << FMT << gsl_spline_eval(data->splines.psi_s, s, data->splines.acc) << endl;
        outfile.close();
}

/* interpolate Laplace shape using gsl methods */
void InterpolateLaplace(LaplaceData *data)
{
        data->splines.Clean();
        data->splines.acc = gsl_interp_accel_alloc();
        /* intialize gsl spline objects */
        data->splines.r_s = gsl_spline_alloc(gsl_interp_cspline, data->size);
        data->splines.z_s = gsl_spline_alloc(gsl_interp_cspline, data->size);
        data->splines.psi_s = gsl_spline_alloc(gsl_interp_cspline, data->size);
        data->splines.r_z = gsl_spline_alloc(gsl_interp_cspline, data->size);
        data->splines.r2_z = gsl_spline_alloc(gsl_interp_cspline, data->size);
        /* initialize splines */
        gsl_spline_init(data->splines.r_s, &(data->s[0]), &(data->r[0]), data->size);
        gsl_spline_init(data->splines.z_s, &(data->s[0]), &(data->z[0]), data->size);
        gsl_spline_init(data->splines.psi_s, &(data->s[0]), &(data->psi[0]), data->size);
        gsl_spline_init(data->splines.r_z, &(data->z[0]), &(data->r[0]), data->size);
        gsl_spline_init(data->splines.r2_z, &(data->z[0]), &(data->r2[0]), data->size);
}

/* interpolate Hooke shape using gsl methods */
void InterpolateHooke(HookeData *data)
{
        data->splines.Clean();
        data->splines.acc = gsl_interp_accel_alloc();
        /* initialize gsl spline objects */
        data->splines.r_s = gsl_spline_alloc(gsl_interp_cspline, data->size);
        data->splines.z_s = gsl_spline_alloc(gsl_interp_cspline, data->size);
        data->splines.psi_s = gsl_spline_alloc(gsl_interp_cspline, data->size);
        data->splines.tau_s = gsl_spline_alloc(gsl_interp_cspline, data->size);
        data->splines.r_z = gsl_spline_alloc(gsl_interp_cspline, data->size);
        data->splines.r2_z = gsl_spline_alloc(gsl_interp_cspline, data->size);
        data->splines.l_s = gsl_spline_alloc(gsl_interp_cspline, data->size);
        data->splines.l_phi = gsl_spline_alloc(gsl_interp_cspline, data->size);
        data->splines.k_phi = gsl_spline_alloc(gsl_interp_cspline, data->size);
        data->splines.t_phi = gsl_spline_alloc(gsl_interp_cspline, data->size);
        /* initialize splines */
        gsl_spline_init(data->splines.r_s, &(data->s[0]), &(data->r[0]), data->size);
        gsl_spline_init(data->splines.z_s, &(data->s[0]), &(data->z[0]), data->size);
        gsl_spline_init(data->splines.psi_s, &(data->s[0]), &(data->psi[0]), data->size);
        gsl_spline_init(data->splines.tau_s, &(data->s[0]), &(data->tau_s[0]), data->size);
        gsl_spline_init(data->splines.r_z, &(data->z[0]), &(data->r[0]), data->size);
        gsl_spline_init(data->splines.r2_z, &(data->z[0]), &(data->r2[0]), data->size);
        gsl_spline_init(data->splines.l_s, &(data->s[0]), &(data->lambda_s[0]), data->size);
        gsl_spline_init(data->splines.l_phi, &(data->s[0]), &(data->lambda_phi[0]), data->size);
        gsl_spline_init(data->splines.k_phi, &(data->s[0]), &(data->kappa_phi[0]), data->size);
        gsl_spline_init(data->splines.t_phi, &(data->s[0]), &(data->tau_phi[0]), data->size);
}

/* fit laplace equation given in data to measurement given in points */
void LaplaceFit(LaplaceData *data, PointSet *points)
{
        if (data->invalid) 
        {
                cout << "invalid input shape in laplace fit!" << endl;
                exit(0);
        }
        
        int N = points->size();
        double *parameters = data->parameters;
        int dim = data->dim;
        
        LaplaceData c[2*dim];
        
        gsl_matrix *J = gsl_matrix_calloc(2*N, dim); /* jacobian */
        gsl_vector *F = gsl_vector_calloc(2*N); /* residuals */
        gsl_vector *D = gsl_vector_calloc(dim); /* parameter shift */
        gsl_vector *PI = gsl_vector_calloc(dim); /* current parameters */
        
        double xi, err_prev, err_next;
        
        pair<double,double> *p;
        Residual error, error_y1, error_x1, error_y2, error_x2;
        int s;
        bool continue_linesearch;
        bool converged = false;
        
        for (int i = 0; i < dim; i++) gsl_vector_set(PI, i, parameters[i]);
        
        while (!converged)
        { 
                /* finite differences */
                c[0].SetParameters(gsl_vector_get(PI, 0) + DX_FIT, gsl_vector_get(PI, 1), gsl_vector_get(PI, 2));
                SolveLaplace(&c[0]);
                c[1].SetParameters(gsl_vector_get(PI, 0), gsl_vector_get(PI, 1) + DX_FIT, gsl_vector_get(PI, 2));
                SolveLaplace(&c[1]);
                c[2].SetParameters(gsl_vector_get(PI, 0), gsl_vector_get(PI, 1), gsl_vector_get(PI, 2) + DX_FIT);
                SolveLaplace(&c[2]);
                c[3].SetParameters(gsl_vector_get(PI, 0) - DX_FIT, gsl_vector_get(PI, 1), gsl_vector_get(PI, 2));
                SolveLaplace(&c[3]);
                c[4].SetParameters(gsl_vector_get(PI, 0), gsl_vector_get(PI, 1) - DX_FIT, gsl_vector_get(PI, 2));
                SolveLaplace(&c[4]);
                c[5].SetParameters(gsl_vector_get(PI, 0), gsl_vector_get(PI, 1), gsl_vector_get(PI, 2) - DX_FIT);
                SolveLaplace(&c[5]);
                
                /* set up jacobian */
                for (int y = 0; y < N; y++)
                {
                        p = NULL;
                        p = &(*points)[y];
                        error = SingleError(data, *p);
                        gsl_vector_set(F, y, - error.delta_r);
                        gsl_vector_set(F, y + N, - error.delta_z);
                        for (int x = 0; x < dim; x++)
                        {
                                error_x1 = SingleError(&c[x], *p);
                                error_x2 = SingleError(&c[x + dim], *p);
                                gsl_matrix_set(J, y, x, (error_x1.delta_r - error_x2.delta_r) / (2*DX_FIT));
                                gsl_matrix_set(J, y + N, x, (error_x1.delta_z - error_x2.delta_z) / (2*DX_FIT));
                        }
                }
                
                QR(J, D, F, 2*N, 3);
                
                /* perform linesearch on update */
                xi = 0.5;
                continue_linesearch = true;
                err_prev = Error(data, points);
                while (continue_linesearch)
                {
                        data->SetParameters(gsl_vector_get(PI, 0) + gsl_vector_get(D, 0), gsl_vector_get(PI, 1) + gsl_vector_get(D, 1), gsl_vector_get(PI, 2) + gsl_vector_get(D, 2));								
                        
                        SolveLaplace(data);
                        if (!data->invalid)
                        {
                                err_next = Error(data, points);
                                if (err_next < err_prev)
                                {
                                        continue_linesearch = false;
                                        if (WATCH_LAPLACE_FITTING) cout << setw(20) << "Error" << setw(20) << data->GetConversion() * err_next << endl;
                                        break;
                                }
                        }
                        
                        gsl_vector_scale(D, xi);
                        
                        if (gsl_blas_dnrm2(D) < EPS_NEWTON_LAPLACE)
                        {
                                continue_linesearch = false;
                                converged = true;
                        }
                        
                        if (WATCH_LAPLACE_FITTING)
                        {
                                cout << setw(20) << "shift" << endl;
                                cout << setw(20) << gsl_vector_get(D, 0) << setw(20) << gsl_vector_get(D, 1) << setw(20) << gsl_vector_get(D, 2) << endl;
                        }
                }
                
                if (!converged)
                {
                        gsl_vector_add(PI, D);
                        data->SetParameters(gsl_vector_get(PI, 0), gsl_vector_get(PI, 1), gsl_vector_get(PI, 2));
                        SolveLaplace(data);
                }
                
                if (WATCH_LAPLACE_FITTING)
                {
                        cout << setw(20) << "pressure" << setw(20) << "density" << setw(20) << "scaling" << endl;
                        cout << setw(20) << gsl_vector_get(PI, 0) << setw(20) << gsl_vector_get(PI, 1) << setw(20) << gsl_vector_get(PI, 2) << endl;
                        cout << endl;
                }
                
        }
        
        /* clean up */
        parameters = NULL;
        p = NULL;
        gsl_vector_free(F);
        gsl_vector_free(D);
        gsl_vector_free(PI);
        gsl_matrix_free(J);
}

/* perform single iteration in fitting Hooke equation stored in data to set of points stored in points
 * booleans can be chosen freely to fit any desired subset of possible parameters or the full set */
bool HookeFitIteration(HookeData *data, PointSet *points, bool fit_p, bool fit_nu, bool fit_k)
{
        int nparams = 0;
        int maxparams = 3;
        bool fit[3] = {fit_p, fit_nu, fit_k};
        vector<int> pind;
        
        for (int i = 0; i < maxparams; i++) if (fit[i]) nparams++;
        for (int i = 0; i < maxparams; i++) if (fit[i]) pind.push_back(i);
        
        if (nparams == 0)
        {
                cout << "fit aborted, all parameters discarded due to invalid derivatives.." << flush << endl;
                return true;
        }
        
        int N = points->size();
        double *parameters = data->parameters;
        double params[3];
        
        HookeData c[2*nparams];
        LaplaceData l[nparams];
        
        gsl_matrix *J = gsl_matrix_calloc(2*N, nparams); /* jacobian */
        gsl_vector *F = gsl_vector_calloc(2*N); /* residual */
        gsl_vector *D = gsl_vector_calloc(nparams); /* parameter shift */
        gsl_vector *PI = gsl_vector_calloc(nparams); /* current parameter set */
        
        for (int i = 0; i < pind.size(); i++) gsl_vector_set(PI, i, parameters[pind[i]]);
        
        bool converged = false;
        int it;
        double err = Error(data, points);
        
        if (WATCH_HOOKE_FITTING)
        {
                cout << setw(20) << "PI(p,nu,k) ->";
                for (int i = 0; i < nparams; i++) cout << setw(20) << setprecision(10) << gsl_vector_get(PI, i);
                cout << endl;
                cout << setw(20) << "ER(px) ->" << setw(20) << data->undeformed->GetConversion() * err << endl;
                cout << setw(20) << "V(unit) ->" << setw(20) << Volume(data) << endl;
        }
        
        /* get current parameters and check if shall be fitted */
        it = 0;
        for (int i = 0; i < maxparams; i++)
        {
                if (fit[i])
                {
                        params[i] = gsl_vector_get(PI, it);
                        it++;
                }
                else
                {
                        params[i] = parameters[i];
                }
        }
        
        /* finite differences */
        bool para_shoot_success[nparams];
        for (int i = 0; i < nparams; i++) para_shoot_success[i] = true;
        
        for (int i = 0; i < nparams; i++)
        {
                if (WATCH_HOOKE_FITTING) cout << "." << flush;
                
                double para = params[pind[i]];
                
                double local[3];
                local[0] = params[0];
                local[1] = params[1];
                local[2] = params[2];
                
                l[i].SetParameters(data->undeformed->GetPressure(), data->undeformed->GetDensity(), data->undeformed->GetScaling());
                SolveLaplace(&l[i]);
                
                local[pind[i]] = para + DX_FIT;
                c[i].SetParameters(local[0], local[1], local[2]);
                c[i].SetLaplace(&l[i]);
                
                bool success1, success2;
                
                success1 = SingleShooting(&c[i]);
                
                local[pind[i]] = para - DX_FIT;
                c[i+nparams].SetParameters(local[0], local[1], local[2]);
                c[i+nparams].SetLaplace(&l[i]);
                
                success2 = SingleShooting(&c[i+nparams]);
                
                if ((!success1) || (!success2)) 
                {
                        para_shoot_success[i] = false;
                }
                else 
                {
                        if (WATCH_HOOKE_FITTING) cout << "." << flush;
                }
        }
        if (WATCH_HOOKE_FITTING) cout << endl << flush;
        
        /* discard fit parameter temporarily, if shooting method has failed */  
        for (int i = 0; i < nparams; i++) 
        {
                if (!para_shoot_success[i])
                {
                        fit[pind[i]] = false;
                        cout << "discard problematic fit-parameter" << flush << endl;
                        cout << "active parameters:\t" << flush;
                        cout << fit[0] << "\t" << fit[1] << "\t" << fit[2] << flush << endl;
                        return HookeFitIteration(data, points, fit[0], fit[1], fit[2]);
                }
        }
        
        pair<double,double> *p;    
        Residual error, error_x1, error_x2;
        /* set up jacobian and residual */
        for (int y = 0; y < N; y++)
        {
                p = &(*points)[y];
                error = SingleError(data, *p);
                gsl_vector_set(F, y, - error.delta_r);
                gsl_vector_set(F, y + N, - error.delta_z);
                for (int x = 0; x < nparams; x++)
                {
                        error_x1 = SingleError(&c[x], *p);
                        error_x2 = SingleError(&c[x + nparams], *p);
                        gsl_matrix_set(J, y, x, (error_x1.delta_r - error_x2.delta_r) / (2*DX_FIT));
                        gsl_matrix_set(J, y + N, x, (error_x1.delta_z - error_x2.delta_z) / (2*DX_FIT));
                }
        }
        
        /* overdetermined system */
        QR(J, D, F, 2*N, nparams);
        
        /* perform linsearch on update */
        double xi = 0.5;
        double err_prev, err_next;
        bool continue_linesearch = true;
        err_prev = Error(data, points);
        gsl_vector_scale(D, 10.0/gsl_blas_dnrm2(D));
        
        while (continue_linesearch)
        {
                if (WATCH_HOOKE_FITTING)
                {
                        cout << setw(20) << "D ->" << flush;
                        for (int i = 0; i < nparams; i++)
                        {
                                cout << setw(20) << setprecision(10) << gsl_vector_get(D, i) << flush;
                        }
                        cout << endl << flush;
                }    
                
                if (gsl_blas_dnrm2(D) < EPS_NEWTON_HOOKE)
                {
                        continue_linesearch = false;
                        converged = true;
                        break;
                }
                
                it = 0;
                for (int i = 0; i < maxparams; i++)
                {
                        if (fit[i])
                        {
                                params[i] = gsl_vector_get(PI, it) + gsl_vector_get(D, it);
				if((i == 1) && (params[i] > 1.0)) {
					params[i] = 1.0;
				}
                                it++;
                        }
                        else 
                        {
                                params[i] = parameters[i];
                        }
                }
                /* check if parameters in still in between valid boundaries */
                if ((params[0] <= 0.0) || (params[2] <= 0.0) || (params[1] < ((CONSTITUTIVE_LAW == 2) ? 0.0 : -1.0)) || (params[1] > 1.0))
                {
                        gsl_vector_scale(D, xi);
                }
                else
                {
                        data->SetParameters(params[0], params[1], params[2]);
                        
                        bool success;
                        success = SingleShooting(data);
                        if ((!data->invalid) && success)
                        {
                                err_next = Error(data, points);
                                if (err_next < err_prev)
                                {
                                        continue_linesearch = false;
                                        break;
                                }
                        }
                        
                        if ((data->invalid) || (!success) || (continue_linesearch))
                        {
                                gsl_vector_scale(D, xi);
                        }
                }
                
        } /* end of linesearch */
        
        if (!converged)
        {
                gsl_vector_add(PI, D);
                it = 0;
                for (int i = 0; i < maxparams; i++)
                {
                        if (fit[i])
                        {
                                params[i] = gsl_vector_get(PI, it);
                                it++;
                        }
                        else
                        {
                                params[i] = parameters[i];
                        }
                }
                data->SetParameters(params[0], params[1], params[2]);
                SingleShooting(data);
        } 
        else
        {
                it = 0;
                for (int i = 0; i < maxparams; i++)
                {
                        if (fit[i])
                        {
                                params[i] = gsl_vector_get(PI, it);
                                it++;
                        }
                        else
                        {
                                params[i] = parameters[i];
                        }
                }
                
                data->SetParameters(params[0], params[1], params[2]);
                SingleShooting(data);
        }
        
        /* clean up */
        p = NULL;
        parameters = NULL;
        gsl_vector_free(F);
        gsl_vector_free(D);
        gsl_vector_free(PI);
        gsl_matrix_free(J);
        
        return converged;
}

/* perform complete fit of elastic shape equations */
void HookeFit(HookeData *data, PointSet *points, const bool fit_p, const bool fit_nu, const bool fit_k)
{
        bool converged = false;
        
        while (!converged)
        {
                converged = HookeFitIteration(data, points, fit_p, fit_nu, fit_k);
        }
        
        if (WATCH_HOOKE_FITTING)
        {
                cout << setw(20) << "final set";
                cout << FMT << data->GetPressure();
                cout << FMT << data->GetPoisson();
                cout << FMT << data->GetCompression();
                cout << FMT << data->undeformed->GetConversion() * Error(data, points) << endl;
                cout << endl;
        }
}


/* nelder mead method to fit elastic shape equations
 * recommended to find initial guess for further improvement with newton method 
 * needs no derivatives -> much more stable far from minimium */
void NelderMead(HookeData *data, PointSet *points)
{
        cout << endl << endl;
        cout << setw(20) << "NELDER MEAD FIT" << endl;
        
        vector<Vertex*> v;
        Vertex redu1, redu2, redu3;
        Vertex refl, expa, cont;
        Vertex cent;
        double edge = 1.0e-1;
        
        double r, e, c;
        r = 1.0;
        e = 3.0;
        c = 1.0/3.0;
        
        /* setup initial simplex */
        /* step 0 */
        v.push_back(new Vertex(data->parameters[0] - 0.1*data->parameters[0], data->parameters[1] - 0.1, data->parameters[2] - 0.1*data->parameters[2]));
        v.push_back(new Vertex(data->parameters[0] + 0.1*data->parameters[0], data->parameters[1], data->parameters[2]));
        v.push_back(new Vertex(data->parameters[0], data->parameters[1] + 0.1, data->parameters[2]));
        v.push_back(new Vertex(data->parameters[0], data->parameters[1], data->parameters[2] + 0.1*data->parameters[2]));
        
        HookeData dummy;
        dummy.SetLaplace(data->undeformed);
        
        /* iteration */
        /* implementation with GOTOs probably bad style, but works */
        while (true)
        {
                STEP1:
                SortVertices(&v, &dummy, points);
                if (v[3]->error - v[0]->error < EPS_NELDERMEAD)
                {
                        break;
                }
                if (WATCH_NELDERMEAD)
                {
                        cout << endl;
                        for (int i = 0; i < v.size(); i++)
                        {
                                cout << setw(20) << "parameter" << setw(20) 
                                << v[i]->p << setw(20) << v[i]->nu << setw(20) 
                                << v[i]->k << setw(20) << "rms" << setw(20) << v[i]->error << endl; 
                        }
                }
                STEP2:
                CalcCentroid(&cent, &v);
                STEP3:
                CalcReflection(&v, &refl, &cent, &dummy, points, r); 
                if ((v[0]->error < refl.error) && (refl.error < v[2]->error))
                {
                        ReplaceWorst(&v, &refl);
                        goto STEP1;
                }
                STEP4:
                if (refl.error < v[0]->error)
                {
                        CalcReflection(&v, &expa, &cent, &dummy, points, e);
                        if (expa.error < refl.error)
                        {
                                ReplaceWorst(&v, &expa);
                                goto STEP1;
                        }
                        else
                        {
                                ReplaceWorst(&v, &refl);
                                goto STEP1;
                        }
                }
                else if (v[2]->error < refl.error)
                {
                        goto STEP5;
                }
                STEP5:
                CalcReflection(&v, &cont, &cent, &dummy, points, -c);
                if (cont.error < v[3]->error)
                {
                        ReplaceWorst(&v, &cont);
                        goto STEP1;
                }
                else 
                {
                        goto STEP6;
                }
                STEP6:
                CalcReflection(&v, &redu1, v[0], &dummy, points, c);
                CalcReflection(&v, &redu2, v[0], &dummy, points, c);
                CalcReflection(&v, &redu3, v[0], &dummy, points, c);
                if (redu1.error > v[3]->error || redu2.error > v[3]->error || redu3.error > v[3]->error)
                {
                        break;
                }
                for (int i = 0; i < 3; i++) 
                {
                        delete v[v.size()-1];
                        v.pop_back();
                }
                v.push_back(redu1.Copy());
                v.push_back(redu2.Copy());
                v.push_back(redu3.Copy());
        }
        
        cout << setw(20) << "FINAL PARAMETERS: " << setw(20) 
        << v[0]->p << setw(20) << v[0]->nu << setw(20) 
        << v[0]->k << setw(20) << "ERROR" << setw(20) << v[0]->error << endl;
        
        data->SetParameters(v[0]->p, v[0]->nu, v[0]->k); 
        SingleShooting(data);
        
        for (int i = 0; i < v.size(); i++) 
        {
                delete v[v.size()-1];
                v.pop_back();
        }
        v.clear();
}

/* helper method for nelder mead */
void ReplaceWorst(vector<Vertex*> *v, Vertex *newone)
{
        delete (*v)[v->size()-1];
        v->pop_back();
        v->push_back(newone->Copy());
} 

/* helper method for nelder mead */
void CalcReflection(vector<Vertex*> *v, Vertex *r, Vertex *c, HookeData *dummy, PointSet *points, double alpha)
{
        // calculate reflection
        double p, nu, k;
        Vertex *worst = (*v)[v->size() - 1];
        r->p = c->p + alpha * (c->p - worst->p);
        r->nu = c->nu + alpha * (c->nu - worst->nu);
        r->k = c->k + alpha * (c->k - worst->k);
        r->error = Error(r, dummy, points);
}

/* helper method for nelder mead */
void CalcCentroid(Vertex *c, vector<Vertex*> *v)
{
        // calculate centroid
        int size = v->size();
        c->p = 0.0; c->nu = 0.0; c->k = 0.0;
        for (int j = 0; j < size - 1; j++)
        {
                c->p += (*v)[j]->p; 
                c->nu += (*v)[j]->nu;
                c->k += (*v)[j]->k;
        }
        c->p /= (double)size - 1.0; c->nu /= (double)size - 1.0; c->k /= (double)size - 1.0;
}

/* helper method for nelder mead */
void SortVertices(vector<Vertex*> *v, HookeData *dummy, PointSet *points)
{
        /* calculate errors */
        double err;
        for (int i = 0; i < v->size(); i++)
        {
                (*v)[i]->error = Error((*v)[i], dummy, points);
        }
        sort(v->begin(), v->end(), CompareVertices);
}

/* overloaded for vertex */
double Error(Vertex *v, HookeData *dummy, PointSet *points)
{
        dummy->SetParameters(v->p, v->nu, v->k);
        bool success = SingleShooting(dummy);
        if (success)
        {
                if ((dummy->GetPoisson() > 1.0) || (dummy->GetPoisson() < ((CONSTITUTIVE_LAW == 2) ? 0.0 : -1.0)) || (dummy->GetCompression() < 0.0) || (dummy->GetPressure() < 0.0))
                {
                        /* ivalid parameters -> assign very large error */
                        return 1000.0;
                }
                else 
                {
                        return dummy->undeformed->GetConversion() * Error(dummy, points);
                }
        }
        else
        {
                /* no shooting success return large error */
                return 1000.0;
        }
}

/* calculate dimensionless volume for capsule scaled width a != 1 */
double UnitVolume(LaplaceData *data)
{
        double scal = data->parameters[2];
        for (int i = 0; i < data->size; i++)
        {
                data->r[i] /= scal;
                data->z[i] /= scal;
                data->r2[i] /= scal * scal;
        }
        
        gsl_interp_accel *accel = gsl_interp_accel_alloc();
        gsl_spline *r2_z = gsl_spline_alloc(gsl_interp_cspline, data->size);
        gsl_spline_init(r2_z, &(data->z[0]), &(data->r2[0]), data->size);
        
        double z_max = data->z[0];
        double volume = gsl_spline_eval_integ(r2_z, z_max, 0.0, accel) * M_PI;
        
        gsl_interp_accel_free(accel);
        gsl_spline_free(r2_z);
        
        return volume;
}

/* calculate capsule volume */
double Volume(Data *shape)
{
        double z_max = gsl_spline_eval(shape->splines.z_s, 0.0, shape->splines.acc);
        return gsl_spline_eval_integ(shape->splines.r2_z, z_max, 0.0, shape->splines.acc) * M_PI;
}

/* calculate capsule area */
double Area(Data *shape)
{
        double z_max = gsl_spline_eval(shape->splines.z_s, 0.0, shape->splines.acc);
        double A = 2.0 * M_PI * gsl_spline_eval_integ(shape->splines.r_z, z_max, 0.0, shape->splines.acc);
        return A;
}

/* calculate distance between shape and contour point */
Residual Dist(Data *shape, double s, double r_i, double z_i)
{
        double z = gsl_spline_eval(shape->splines.z_s, s, shape->splines.acc);
        double r = gsl_spline_eval(shape->splines.r_s, s, shape->splines.acc);
        if (r_i < 0.0)
        {
                r *= (-1);
        }
        Residual result;
        result.s = s;
        result.delta_z = z_i - z;
        result.delta_r = r_i - r;
        return result;
}

/* calculate derivative of distance between shape and contour */
double DistDeriv(Data *shape, double s, double r_i, double z_i)
{
        double r_deriv = gsl_spline_eval_deriv(shape->splines.r_s, s, shape->splines.acc);
        double z_deriv = gsl_spline_eval_deriv(shape->splines.z_s, s, shape->splines.acc);
        double r = gsl_spline_eval(shape->splines.r_s, s, shape->splines.acc);
        double z = gsl_spline_eval(shape->splines.z_s, s, shape->splines.acc);
        return 2*(z*z_deriv - z_deriv*z_i + r*r_deriv - r_deriv*r_i);
}

/* calculate fit error */
double Error(Data *shape, PointSet *points)
{
        int size = points->size();
        double rms = 0.0;
        for (int i = 0; i < size; i++)
        {
                rms += SingleError(shape, (*points)[i]).Norm();
        }
        return (1.0/sqrt((double)size + 1.0))*sqrt(rms);
}

/* calculate residual vs s 
 * use for debug purposes only */
double HookeError(HookeData *shape, PointSet *points)
{
        int size = points->size();
        double rms = 0.0;
        
        ofstream res("residuum.dat", ios::out);
        for (int i = 0; i < size; i++)
        {
                rms += SingleError(shape, (*points)[i]).Norm();
                double signum = -1.0;
                if ((*points)[i].first == 0.0)
                {
                        if (fabs((*points)[i].second) < gsl_spline_eval(shape->splines.z_s, SingleError(shape, (*points)[i]).s, shape->splines.acc))
                        {
                                signum *= (-1.0);
                        }
                }
                else if (fabs((*points)[i].first) < gsl_spline_eval(shape->splines.r_s, SingleError(shape, (*points)[i]).s, shape->splines.acc))
                {
                        signum *= (-1.0);
                }
                res << setw(20) << i + 1 << setw(20) << shape->undeformed->GetConversion() * signum * sqrt(SingleError(shape, (*points)[i]).Norm()) << endl;
        }
        res.close();
        return rms;
}

/* calculate minimum distance between shape and single contour point */
Residual SingleError(Data *shape, pair<double,double> point)
{
        double a = 0.0, b = shape->L0;
        double pos = 0.5 * b;
        double step = 0.5 * b;
        bool done = false;
        double d1, d2;
        double d;
        bool rightend, leftend;
        double delta = 1e-4;
        double left, right;
        while (!done)
        {
                rightend = false;
                leftend = false;
                step /= 2.0;
                if (step < EPS_RMS)
                {
                        done = true;
                        break;
                }
                if (pos + delta > b)
                { 
                        d1 = Dist(shape, b, point.first, point.second).Norm();
                        rightend = true;
                }
                else
                {
                        d1 = Dist(shape, pos + delta, point.first, point.second).Norm();  
                }
                if (pos - delta < a)  
                {
                        d2 = Dist(shape, a, point.first, point.second).Norm();
                        leftend = true;
                }
                else
                {
                        d2 = Dist(shape, pos - delta, point.first, point.second).Norm(); 
                }
                if (d1 < d2)
                {
                        if (rightend)
                                pos = b;
                        else
                                pos += step;
                }
                else if (d1 > d2)
                {
                        if (leftend)
                                pos = a;
                        else
                                pos -= step;
                }
                else
                {
                        d = Dist(shape, pos, point.first, point.second).Norm();
                        if ((d1 >= d) && (d2 >= d))
                        {
                                done = true;
                        }
                        else if ((d1 < d) && (d2 < d))
                        {
                                left = Dist(shape, a, point.first, point.second).Norm();
                                right = Dist(shape, b, point.first, point.second).Norm();
                                if (left < right)
                                {
                                        pos = left;
                                        done = true;
                                }
                                else
                                {
                                        pos = right;
                                        done = true;
                                }
                        }
                        else cout << "error in rms calculation" << flush << endl;
                }
        }
        return Dist(shape, pos, point.first, point.second);
}

/* M X N MATRIX A
 * N VECTOR X
 * M VECTOR B
 * QR decomposition via householder transforms
 * solves quadratic as well as over-determined systems */
void QR(gsl_matrix *A, gsl_vector *x, gsl_vector *b, int M, int N)
{
        bool square = (M == N) ? true : false;
        
        gsl_matrix *J = gsl_matrix_calloc(M, N);
        gsl_matrix *Qi = gsl_matrix_calloc(M, M);
        gsl_matrix *vvT = gsl_matrix_calloc(M, M);
        
        gsl_matrix *Q = gsl_matrix_calloc(M, M);
        gsl_matrix *R = gsl_matrix_calloc(M, N);
        
        // intialize Q
        gsl_matrix_set_identity(Q);
        
        double alpha, abs_v;
        int cltrafo = N;
        if (square) cltrafo--;
        gsl_matrix_memcpy(J, A);
        
        for (int i = 0; i < cltrafo; i++)
        {
                gsl_vector *v = gsl_vector_calloc(M - i);
                for (int k = 0; k < M - i; k++)
                {
                        gsl_vector_set(v, k, gsl_matrix_get(J, i + k, i));
                }
                
                if (WATCH_QR_DECOMPOSITION)
                {
                        PrintVector(v, "v");
                        PrintMatrix(J, "J");
                        PrintVector(b, "b");
                }
                
                alpha = gsl_blas_dnrm2(v);
                if (gsl_vector_get(v, 0) > 0) alpha *= (-1);
                
                gsl_vector *e = gsl_vector_calloc(M - i);
                gsl_vector_set(e, 0, alpha);  
                gsl_vector_add(v, e);
                
                if (WATCH_QR_DECOMPOSITION)
                {
                        PrintVector(v, "v + e");
                        PrintVector(e, "e");
                }
                
                gsl_vector_free(e);
                
                abs_v = gsl_blas_dnrm2(v);
                gsl_vector_scale(v, 1/abs_v);
                
                if (WATCH_QR_DECOMPOSITION) PrintVector(v, "v + e unit");
                
                for (int m = 0; m < M; m++) 
                {
                        for (int n = 0; n < M; n++)
                        {
                                if (m == n)
                                {
                                        gsl_matrix_set(Qi, m, n, 1.0);
                                }
                                else
                                {
                                        gsl_matrix_set(Qi, m, n, 0.0);
                                }
                        }
                }
                
                if (WATCH_QR_DECOMPOSITION) PrintMatrix(Qi, "Qi");
                
                for (int m = 0; m < M; m++) 
                {
                        for (int n = 0; n < M; n++)
                        {
                                if ((m >= i) && (n >= i))
                                {
                                        gsl_matrix_set(vvT, m, n, gsl_vector_get(v, m - i) * gsl_vector_get(v, n - i));
                                }
                                else
                                {
                                        gsl_matrix_set(vvT, m, n, 0.0);
                                }
                        }
                }
                
                gsl_vector_free(v);
                
                gsl_matrix_scale(vvT, 2.0);
                
                if (WATCH_QR_DECOMPOSITION) PrintMatrix(vvT, "2*v*v^T");
                
                gsl_matrix_sub(Qi, vvT);
                
                if (WATCH_QR_DECOMPOSITION) PrintMatrix(Qi, "Qi");
                
                gsl_matrix *tmp = gsl_matrix_calloc(M, M);
                gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, Q, Qi, 0.0, tmp);
                gsl_matrix_memcpy(Q, tmp);
                gsl_matrix_free(tmp);
                
                gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Q, A, 0.0, R);
                
                if (WATCH_QR_DECOMPOSITION) PrintMatrix(R, "R");
                
                gsl_matrix_memcpy(J, R);
        }
        
        if (WATCH_QR_DECOMPOSITION)
        {
                PrintMatrix(Q, "Q");
                PrintMatrix(R, "R");
                gsl_matrix *T = gsl_matrix_calloc(M, N);
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Q, R, 0.0, T);
                gsl_matrix_sub(T, A);
                PrintMatrix(T, "Test");
                gsl_matrix_free(T);
        }
        
        gsl_matrix_free(J);
        gsl_matrix_free(Qi);
        gsl_matrix_free(vvT);
        
        gsl_vector *u = gsl_vector_calloc(M);
        gsl_blas_dgemv(CblasTrans, 1.0, Q, b, 0.0, u);
        if (WATCH_QR_DECOMPOSITION) PrintVector(u, "QTF");
        gsl_matrix *U = gsl_matrix_calloc(N, N);
        for (int m = 0; m < N; m++)
        {
                gsl_vector_set(x, m, gsl_vector_get(u, m));
        }
        for (int m = 0; m < N; m++)
        {
                for (int n = m; n < N; n++)
                {
                        gsl_matrix_set(U, m, n, gsl_matrix_get(R, m, n));
                }
        }
        
        /* solve by backward substitution */
        gsl_linalg_R_svx(U, x);
        
        if (WATCH_QR_DECOMPOSITION)
        {
                PrintMatrix(U, "U");
                PrintVector(x, "x");
        }
        
        /* clean up */
        gsl_vector_free(u);
        gsl_matrix_free(U);
        gsl_matrix_free(R);
        gsl_matrix_free(Q);
}

/* helper methods for debug purposes only */
void PrintVector(gsl_vector *x, const char *name)
{
        cout << "VECTOR: " << name << endl;
        for (int m = 0; m < x->size; m++)
        {
                cout << setw(20) << gsl_vector_get(x, m) << endl;
        }
        cout << endl << endl;
}

/* helper methods for debug purposes only */
void PrintMatrix(gsl_matrix *x, const char *name)
{
        cout << "MATRIX: " << name << endl;
        for (int m = 0; m < x->size1; m++)
        {
                for (int n = 0; n < x->size2; n++)
                {
                        cout << setw(20) << gsl_matrix_get(x, m, n);
                }
                cout << endl;
        }
        cout << endl << endl;
}


