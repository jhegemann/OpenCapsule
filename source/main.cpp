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

#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <fstream>
#include <time.h>
#include <string>
#include <sstream>
#include <cmath>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include "interface.h"
#include "config.h"

using namespace std;
namespace flsys = boost::filesystem;
namespace propt = boost::program_options;

int main(int argc, char **argv)
{
        cout << "OpenCapsule" << endl;
        cout << "Version 1.0" << endl;
        cout << "Copyright (C) 2017 Jonas Hegemann, TU Dortmund" << endl;
        cout << "C++ Software / Command Line Interface" << endl;
        cout << "Author: Jonas Hegemann (jonas.hegemann@tu-dortmund.de)" << endl;
        cout << "Framework: Sebastian Knoche (sebastian.knoche@tu-dortmund.de)" << endl;
        cout << "Supervisor: Jan Kierfeld (jan.kierfeld@tu-dortmund.de)" << endl;
        cout << "This program comes with ABSOLUTELY NO WARRANTY" << endl;
        cout << "This is free software, and you are welcome to redistribute it under certain conditions given in license.txt" << endl;
        
        gsl_set_error_handler_off();
        
        flsys::path path_global_out(GLOBAL_OUT_FOLDER.c_str()); 
        if (flsys::create_directory(path_global_out)) cout << "directory ./global_out created..." << endl;
        else cout << "directory ./global_out exists..." << endl;
        
        flsys::path path_input(INPUT_FOLDER.c_str());
        if (flsys::create_directory(path_input)) cout << "directory ./input created..." << endl;
        else cout << "directory ./input exists..." << endl;
        
        flsys::path path_out(OUT_FOLDER.c_str());
        if (flsys::create_directory(path_out)) cout << "directory ./out created..." << endl;
        else cout << "directory ./out exists..." << endl;
        
        flsys::path path_tmp(TMP_FOLDER.c_str());
        if (flsys::create_directory(path_tmp)) cout << "directory ./tmp created..." << endl;
        else cout << "directory ./tmp exists..." << endl;
        
        flsys::path path_config(CONFIG_FOLDER.c_str());
        if (flsys::create_directory(path_config)) 
        {
                cout << "directory ./config created..." << endl;
                fstream cfg;
                string configuration_name = CONFIG_FOLDER + "config.cfg";
                cfg.open(configuration_name.c_str(), ios::out);
                cfg << "REFERENCE_SHAPE                                         1.png" << endl;
                cfg << "ELASTIC_SHAPE                                           2.png" << endl;
                cfg << "EXPERIMENT_DENSITY                                      140.0" << endl;
                cfg << "EXPERIMENT_CAPDIAMETER                                  0.001" << endl;
                cfg.close();
                cout << "smallest possible configuration written..." << endl;
                cout << "please adapt configuration file " << configuration_name << " to your input data!" << endl;
                cout << "list reference and deformed shapes, set density difference and needle diameter!" << endl;
                cout << "terminating..." << endl;
                exit(0);
        }
        else cout << "directory ./config exists..." << endl;
        
        cout << "parsing config..." << endl << endl;
        parse_configuration();
        
        string filename_image;
        
        double pressure = 1.0;
        double poisson = 0.5;
        double compression = 10.0;
        
        propt::options_description desc("Allowed Options");
        desc.add_options()
        ("help,h", "produce help message")
        
        ("reference,r", "Reference Regression")
        ("sequence,s", "Sequence Regression")
        ("image,i", propt::value<string>(&filename_image), "Edge Detection")
        ("poisson", propt::value<double>(&poisson), "Initial Guess Poisson's Ratio")
        ("compression", propt::value<double>(&compression), "Initial Guess Area Compression Modulus")
        ("pressure", propt::value<double>(&pressure), "Initial Guess Pressure Deflated Capsule");
        
        propt::variables_map vm;
        propt::store(propt::parse_command_line(argc, argv, desc), vm);
        propt::notify(vm);
        
        if (vm.count("help")) 
        {
                cout << desc << endl;
                cout << "Usage Examples:" << endl;
                cout << "  OpenCapsule --image ./input/test.png" << endl;
                cout << "  OpenCapsule --reference" << endl;
                cout << "  OpenCapsule --sequence" << endl;
                cout << "  OpenCapsule --sequence --pressure 0.5 --poisson 0.5 --compression 10.0" << endl;
                return 0;
        }
        
        if (vm.count("reference") && vm.count("sequence")) 
        {
                cout << "invalid command line options...\n" << endl;
                return 0;
        }
        
        if (vm.count("reference")  || vm.count("sequence")) InitializeSolver();
        LaplaceData *reference = NULL;
        
        if (vm.count("reference")  || vm.count("sequence"))
        {
                reference = LaplaceMean(2.5, 0.25, 2.0);
        }
        
        if (vm.count("sequence"))
        {
                Sequence(reference, pressure, poisson, compression);
        }
        
        if (vm.count("image")) 
        {
                EdgeDetection(filename_image);
        }
        
        if (reference) delete reference;
        
        if (vm.count("reference")) FinalizeSolver();
        
        return 0;
}

