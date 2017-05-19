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

#include "config.h"

vector<string> REFERENCE_SHAPE;
vector<string> ELASTIC_SHAPE;
double EPS_RMS = STD_EPS_RMS;
double EPS_NEWTON_LAPLACE = STD_EPS_NEWTON_LAPLACE;
double EPS_NEWTON_HOOKE = STD_EPS_NEWTON_HOOKE;
double EPS_SINGLE_SHOOTING = STD_EPS_SINGLE_SHOOTING;
double EPS_MULTI_SHOOTING = STD_EPS_MULTI_SHOOTING;
double DX_FIT = STD_DX_FIT;
double DX_SHOOTING = STD_DX_SHOOTING;
double INTEGRATION_STEP_LAPLACE = STD_INTEGRATION_STEP_LAPLACE;
double INTEGRATION_STEP_HOOKE = STD_INTEGRATION_STEP_HOOKE;
double EPS_NELDERMEAD = STD_EPS_NELDERMEAD;
double SCALE_IMAGE = STD_SCALE_IMAGE;
double T_HIGH = STD_T_HIGH;
double T_LOW = STD_T_LOW;
double GAUSSIAN_SIGMA = STD_GAUSSIAN_SIGMA;
int R_MIN = STD_R_MIN;
int R_MAX = STD_R_MAX;
double THRESHOLD_CAPILLARY = STD_THRESHOLD_CAPILLARY;
int TOP_BUFFER = STD_TOP_BUFFER;
bool CHECK_IF_CLOSED = STD_CHECK_IF_CLOSED;
double EPS_IMPLICIT_RK = STD_EPS_IMPLICIT_RK;
bool IMPLICIT_INTEGRATION = STD_IMPLICIT_INTEGRATION;
string INPUT_FOLDER = STD_INPUT_FOLDER;
string CONFIG_FOLDER = STD_CONFIG_FOLDER;
string OUT_FOLDER = STD_OUT_FOLDER;
string GLOBAL_OUT_FOLDER = STD_GLOBAL_OUT_FOLDER;
string TMP_FOLDER = STD_TMP_FOLDER;
bool WATCH_SINGLE_SHOOTING = STD_WATCH_SINGLE_SHOOTING;
bool WATCH_MULTI_SHOOTING = STD_WATCH_MULTI_SHOOTING;
bool WATCH_LAPLACE_FITTING = STD_WATCH_LAPLACE_FITTING;
bool WATCH_HOOKE_FITTING = STD_WATCH_HOOKE_FITTING;
bool WATCH_QR_DECOMPOSITION = STD_WATCH_QR_DECOMPOSITION;
bool WATCH_METROPOLIS_LAPLACE = STD_WATCH_METROPOLIS_LAPLACE;
bool WATCH_METROPOLIS_HOOKE = STD_WATCH_METROPOLIS_HOOKE;
bool WATCH_NELDERMEAD = STD_WATCH_NELDERMEAD;
bool FIT_PRESSURE = STD_FIT_PRESSURE;
bool FIT_POISSON = STD_FIT_POISSON;
bool FIT_COMPRESSION = STD_FIT_COMPRESSION;
bool FIX_CHARACTERISTICS = STD_FIX_CHARACTERISTICS;
bool EXTENDED_SHOOTING = STD_EXTENDED_SHOOTING;
bool FORCE_SYMMETRY = STD_FORCE_SYMMETRY;
bool FORCE_BOUNDARY_SYMMETRY = STD_FORCE_BOUNDARY_SYMMETRY;
bool NELDERMEAD_PREFITTING = STD_NELDERMEAD_PREFITTING;
bool GNUPLOT_SUPPORT = STD_GNUPLOT_SUPPORT;
bool PARAMETER_TRACING = STD_PARAMETER_TRACING;
double EXPERIMENT_DENSITY = STD_EXPERIMENT_DENSITY;
double EXPERIMENT_CAPDIAMETER = STD_EXPERIMENT_CAPDIAMETER;
double GAMMA_SCALE = STD_GAMMA_SCALE;
double WRINKLING_WAVELENGTH = STD_WRINKLING_WAVELENGTH;

size_t replace_all(string& str, const string from, const string to)
{
        size_t start_pos;
        size_t number = 0;
        while((start_pos = str.find(from)) != string::npos) {
                str.replace(start_pos, from.length(), to);
                number++;
        }
        return number;
}

void trim_wss(string& str)
{
        str.erase(0, str.find_first_not_of(" "));
        str.erase(str.find_last_not_of(" ") + 1);
        replace_all(str, "  ", " ");
}

vector<string> explode(string text, string delimiter)
{
        vector<string> result;
        replace_all(text, delimiter+delimiter, delimiter);
        text.erase(0, text.find_first_not_of(delimiter));
        text.erase(text.find_last_not_of(delimiter) + 1);
        size_t found;
        string segment;
        while((found = text.find(delimiter)) != string::npos) {
                segment = text.substr(0, found);
                trim_wss(segment);
                result.push_back(segment);
                text = text.substr(found + delimiter.length());
        }
        trim_wss(text);
        result.push_back(text);
        return result;
}

bool contains(map<string, string> &cfg, string key)
{
        map<string, string>::iterator it = cfg.find(key);
        if (it != cfg.end()) 
        {
                return true;
        } else 
        {
                return false;
        }
}

void parse_configuration()
{
        string path = CONFIG_FOLDER;
        path.append("config.cfg");
        ifstream cfg(path.c_str(), ios::in);
        
        map<string, string> fs;
        
        string line;
        int line_counter = 0;
        while (getline(cfg, line)) {
                line_counter++;
                vector<string> line_content = explode(line, " ");
                if (line_content.size() != 2) 
                {
                        cout << "configuration problem: invalid key value pair in line " << line_counter << "!" << endl;
                        exit(0);
                }
                fs.insert(make_pair(line_content[0], line_content[1]));
        }
        
        REFERENCE_SHAPE = explode(fs["REFERENCE_SHAPE"], ";");
        ELASTIC_SHAPE = explode(fs["ELASTIC_SHAPE"], ";");
        
        if (contains(fs, "EPS_RMS")) EPS_RMS = atof(fs["EPS_RMS"].c_str());
        if (contains(fs, "EPS_NEWTON_LAPLACE")) EPS_NEWTON_LAPLACE = atof(fs["EPS_NEWTON_LAPLACE"].c_str());
        if (contains(fs, "EPS_NEWTON_HOOKE")) EPS_NEWTON_HOOKE = atof(fs["EPS_NEWTON_HOOKE"].c_str());
        if (contains(fs, "EPS_SINGLE_SHOOTING")) EPS_SINGLE_SHOOTING = atof(fs["EPS_SINGLE_SHOOTING"].c_str());
        if (contains(fs, "EPS_MULTI_SHOOTING")) EPS_MULTI_SHOOTING = atof(fs["EPS_MULTI_SHOOTING"].c_str());
        if (contains(fs, "DX_FIT")) DX_FIT = atof(fs["DX_FIT"].c_str());
        if (contains(fs, "DX_SHOOTING")) DX_SHOOTING = atof(fs["DX_SHOOTING"].c_str());
        if (contains(fs, "INTEGRATION_STEP_LAPLACE")) INTEGRATION_STEP_LAPLACE = atof(fs["INTEGRATION_STEP_LAPLACE"].c_str()); 
        if (contains(fs, "INTEGRATION_STEP_HOOKE")) INTEGRATION_STEP_HOOKE = atof(fs["INTEGRATION_STEP_HOOKE"].c_str());
        if (contains(fs, "EPS_NELDERMEAD")) EPS_NELDERMEAD = atof(fs["EPS_NELDERMEAD"].c_str()); 
        if (contains(fs, "SCALE_IMAGE")) SCALE_IMAGE = atof(fs["SCALE_IMAGE"].c_str()); 
        if (contains(fs, "T_HIGH")) T_HIGH = atof(fs["T_HIGH"].c_str());
        if (contains(fs, "T_LOW")) T_LOW = atof(fs["T_LOW"].c_str()); 
        if (contains(fs, "GAUSSIAN_SIGMA")) GAUSSIAN_SIGMA = atof(fs["GAUSSIAN_SIGMA"].c_str());
        if (contains(fs, "R_MIN")) R_MIN = atoi(fs["R_MIN"].c_str()); 
        if (contains(fs, "R_MAX")) R_MAX = atoi(fs["R_MAX"].c_str()); 
        if (contains(fs, "THRESHOLD_CAPILLARY")) THRESHOLD_CAPILLARY = atof(fs["THRESHOLD_CAPILLARY"].c_str());
        if (contains(fs, "TOP_BUFFER")) TOP_BUFFER = atoi(fs["TOP_BUFFER"].c_str());
        if (contains(fs, "CHECK_IF_CLOSED")) CHECK_IF_CLOSED = (bool)atoi(fs["CHECK_IF_CLOSED"].c_str());
        if (contains(fs, "EPS_IMPLICIT_RK")) EPS_IMPLICIT_RK = atof(fs["EPS_IMPLICIT_RK"].c_str());
        if (contains(fs, "IMPLICIT_INTGERATION")) IMPLICIT_INTEGRATION = (bool)atoi(fs["IMPLICIT_INTEGRATION"].c_str());
        if (contains(fs, "INPUT_FOLDER")) INPUT_FOLDER = fs["INPUT_FOLDER"];
        if (contains(fs, "CONFIG_FOLDER")) CONFIG_FOLDER = fs["CONFIG_FOLDER"];
        if (contains(fs, "OUT_FOLDER")) OUT_FOLDER = fs["OUT_FOLDER"]; 
        if (contains(fs, "GLOBAL_OUT_FOLDER")) GLOBAL_OUT_FOLDER = fs["GLOBAL_OUT_FOLDER"];
        if (contains(fs, "TMP_FOLDER")) TMP_FOLDER = fs["TMP_FOLDER"]; 
        if (contains(fs, "WATCH_SINGLE_SHOOTING")) WATCH_SINGLE_SHOOTING = (bool)atoi(fs["WATCH_SINGLE_SHOOTING"].c_str());
        if (contains(fs, "WATCH_MULTI_SHOOTING")) WATCH_MULTI_SHOOTING = (bool)atoi(fs["WATCH_MULTI_SHOOTING"].c_str());
        if (contains(fs, "WATCH_LAPLACE_FITTING")) WATCH_LAPLACE_FITTING = (bool)atoi(fs["WATCH_LAPLACE_FITTING"].c_str());
        if (contains(fs, "WATCH_HOOKE_FITTING")) WATCH_HOOKE_FITTING = (bool)atoi(fs["WATCH_HOOKE_FITTING"].c_str());
        if (contains(fs, "WATCH_QR_DECOMPOSITION")) WATCH_QR_DECOMPOSITION = (bool)atoi(fs["WATCH_QR_DECOMPOSITION"].c_str());
        if (contains(fs, "WATCH_METROPOLIS_LAPLACE")) WATCH_METROPOLIS_LAPLACE = (bool)atoi(fs["WATCH_METROPOLIS_LAPLACE"].c_str()); 
        if (contains(fs, "WATCH_METROPOLIS HOOKE")) WATCH_METROPOLIS_HOOKE = (bool)atoi(fs["WATCH_METROPOLIS_HOOKE"].c_str());
        if (contains(fs, "WATCH_NELDERMEAD")) WATCH_NELDERMEAD = (bool)atoi(fs["WATCH_NELDERMEAD"].c_str());
        if (contains(fs, "FIT_PRESSURE")) FIT_PRESSURE = (bool)atoi(fs["FIT_PRESSURE"].c_str()); 
        if (contains(fs, "FIT_POISSON")) FIT_POISSON = (bool)atoi(fs["FIT_POISSON"].c_str()); 
        if (contains(fs, "FIT_COMPRESSION")) FIT_COMPRESSION = (bool)atoi(fs["FIT_COMPRESSION"].c_str()); 
        if (contains(fs, "FIX_CHARACTERISTICS")) FIX_CHARACTERISTICS = (bool)atoi(fs["FIX_CHARACTERISTICS"].c_str());
        if (contains(fs, "EXTENDED_SHOOTING")) EXTENDED_SHOOTING = (bool)atoi(fs["EXTENDED_SHOOTING"].c_str()); 
        if (contains(fs, "FORCE_SYMMETRY")) FORCE_SYMMETRY = (bool)atoi(fs["FORCE_SYMMETRY"].c_str()); 
        if (contains(fs, "FORCE_BOUNDARY_SYMMETRY")) FORCE_BOUNDARY_SYMMETRY = (bool)atoi(fs["FORCE_BOUNDARY_SYMMETRY"].c_str());
        if (contains(fs, "NELDERMEAD_PREFITTING")) NELDERMEAD_PREFITTING = (bool)atoi(fs["NELDERMEAD_PREFITTING"].c_str());
        if (contains(fs, "GNUPLOT_SUPPORT")) GNUPLOT_SUPPORT = (bool)atoi(fs["GNUPLOT_SUPPORT"].c_str());
        if (contains(fs, "PARAMETER_TRACING")) PARAMETER_TRACING = (bool)atoi(fs["PARAMETER_TRACING"].c_str());
        if (contains(fs, "GAMMA_SCALE")) GAMMA_SCALE = atof(fs["GAMMA_SCALE"].c_str());
        if (contains(fs, "WRINKLING_WAVELENGTH")) WRINKLING_WAVELENGTH = atof(fs["WRINKLING_WAVELENGTH"].c_str());
        
        if (!contains(fs, "EXPERIMENT_DENSITY"))
        {
                cout << "please specify the density difference between inner and outer phase in your configuration file!" << endl;
                exit(0);
        }
        EXPERIMENT_DENSITY = atof(fs["EXPERIMENT_DENSITY"].c_str()); 
        if (!contains(fs, "EXPERIMENT_CAPDIAMETER")) 
        {
                cout << "please specify the outer capillary diameter in cour configuration file!" << endl;
                exit(0);
        }
        EXPERIMENT_CAPDIAMETER = atof(fs["EXPERIMENT_CAPDIAMETER"].c_str()); 
        
        for(int i=0; i<REFERENCE_SHAPE.size(); i++) 
        {
                REFERENCE_SHAPE[i] = INPUT_FOLDER + REFERENCE_SHAPE[i];
        }
        for(int i=0; i<ELASTIC_SHAPE.size(); i++) 
        {
                ELASTIC_SHAPE[i] = INPUT_FOLDER + ELASTIC_SHAPE[i];
        }
        
        if (REFERENCE_SHAPE.size() == 0) 
        {
                cout << "please specify your reference files in your configuration file!" << endl;
                exit(0);
        }
        
}
