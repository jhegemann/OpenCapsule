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

#ifndef CONFIG_H
#define CONFIG_H

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <string>

using namespace std;

#define STD_EPS_RMS 1e-16
#define STD_EPS_NEWTON_LAPLACE 1e-6
#define STD_EPS_NEWTON_HOOKE 1e-6
#define STD_EPS_SINGLE_SHOOTING 1e-6
#define STD_EPS_MULTI_SHOOTING 1e-6
#define STD_DX_FIT 1e-6
#define STD_DX_SHOOTING 1e-6
#define STD_INTEGRATION_STEP_LAPLACE 1e-4
#define STD_INTEGRATION_STEP_HOOKE 1e-4
#define STD_EPS_NELDERMEAD 1e-4
#define STD_SCALE_IMAGE 1e-2
#define STD_T_HIGH 0.6
#define STD_T_LOW 0.3
#define STD_GAUSSIAN_SIGMA 1.0
#define STD_R_MIN 1
#define STD_R_MAX 5
#define STD_THRESHOLD_CAPILLARY 5
#define STD_TOP_BUFFER 0
#define STD_CHECK_IF_CLOSED 1
#define STD_IMPLICIT_INTEGRATION 0
#define STD_EPS_IMPLICIT_RK 1e-10
#define STD_INPUT_FOLDER "./input/"
#define STD_CONFIG_FOLDER "./config/"
#define STD_OUT_FOLDER "./out/"
#define STD_GLOBAL_OUT_FOLDER "./global_out/"
#define STD_TMP_FOLDER "./tmp/"
#define STD_WATCH_SINGLE_SHOOTING 0
#define STD_WATCH_MULTI_SHOOTING 0
#define STD_WATCH_LAPLACE_FITTING 1
#define STD_WATCH_HOOKE_FITTING 1
#define STD_WATCH_QR_DECOMPOSITION 0
#define STD_WATCH_METROPOLIS_LAPLACE 0
#define STD_WATCH_METROPOLIS_HOOKE 0
#define STD_WATCH_NELDERMEAD 1
#define STD_FIT_PRESSURE 1
#define STD_FIT_POISSON 1
#define STD_FIT_COMPRESSION 1
#define STD_FIX_CHARACTERISTICS 1
#define STD_EXTENDED_SHOOTING 1
#define STD_FORCE_SYMMETRY 0
#define STD_FORCE_BOUNDARY_SYMMETRY 1
#define STD_NELDERMEAD_PREFITTING 0
#define STD_GNUPLOT_SUPPORT 0
#define STD_PARAMETER_TRACING 0
#define STD_EXPERIMENT_DENSITY 140
#define STD_EXPERIMENT_CAPDIAMETER 0.00165
#define STD_GAMMA_SCALE 1.0
#define STD_WRINKLING_WAVELENGTH 0.0

extern vector<string> REFERENCE_SHAPE;
extern vector<string> ELASTIC_SHAPE;
extern double EPS_RMS;
extern double EPS_NEWTON_LAPLACE;
extern double EPS_NEWTON_HOOKE;
extern double EPS_SINGLE_SHOOTING;
extern double EPS_MULTI_SHOOTING;
extern double DX_FIT;
extern double DX_SHOOTING;
extern double INTEGRATION_STEP_LAPLACE;
extern double INTEGRATION_STEP_HOOKE;
extern double EPS_NELDERMEAD;
extern double SCALE_IMAGE;
extern double T_HIGH;
extern double T_LOW;
extern double GAUSSIAN_SIGMA;
extern int R_MIN;
extern int R_MAX;
extern double THRESHOLD_CAPILLARY;
extern int TOP_BUFFER;
extern double EPS_IMPLICIT_RK;
extern bool IMPLICIT_INTEGRATION;
extern bool CHECK_IF_CLOSED;
extern string INPUT_FOLDER;
extern string CONFIG_FOLDER;
extern string OUT_FOLDER;
extern string GLOBAL_OUT_FOLDER;
extern string TMP_FOLDER;
extern bool WATCH_SINGLE_SHOOTING;
extern bool WATCH_MULTI_SHOOTING;
extern bool WATCH_LAPLACE_FITTING;
extern bool WATCH_HOOKE_FITTING;
extern bool WATCH_QR_DECOMPOSITION;
extern bool WATCH_METROPOLIS_LAPLACE;
extern bool WATCH_METROPOLIS_HOOKE;
extern bool WATCH_NELDERMEAD;
extern bool FIT_PRESSURE;
extern bool FIT_POISSON;
extern bool FIT_COMPRESSION;
extern bool FIX_CHARACTERISTICS;
extern bool EXTENDED_SHOOTING;
extern bool FORCE_SYMMETRY;
extern bool FORCE_BOUNDARY_SYMMETRY;
extern bool NELDERMEAD_PREFITTING;
extern bool GNUPLOT_SUPPORT;
extern bool PARAMETER_TRACING;
extern double EXPERIMENT_DENSITY;
extern double EXPERIMENT_CAPDIAMETER;
extern double GAMMA_SCALE;
extern double WRINKLING_WAVELENGTH;

size_t replace_all(string& str, const string from, const string to);
void trim_wss(string& str);
vector<string> explode(string text, string delimiter);
bool contains(map<string, string> &cfg, string key);
void parse_configuration();


#endif 
