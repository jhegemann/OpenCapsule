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

#include "interface.h"

/* global variables */
vector<PointSet*> _reference_shape;
vector<PointSet*> _elastic_shape;
vector<ImageInfo*> _reference_image;
vector<ImageInfo*> _elastic_image;

/* methods */
vector<ImageInfo*> LoadImages(string preamble, vector<string> &file, vector<PointSet*> &contour, bool elastic)
{
        vector<ImageInfo*> result;
        cout << "loading images...\n" << endl;
        for(int i=0; i<file.size(); i++) {
                cout << file[i] << endl;
                string filename2, filename3, filename4;
                /* output files */
                filename2 = GLOBAL_OUT_FOLDER + preamble + "_" + to_string<int>(i) + "_image.dat";
                filename3 = GLOBAL_OUT_FOLDER + preamble + "_" + to_string<int>(i) + "_image_info.dat";
                filename4 = GLOBAL_OUT_FOLDER + preamble + "_" + to_string<int>(i) + "_image_binary.png";
                /* load image */
                Image image;
                if (elastic && FIX_CHARACTERISTICS) {
                        double avg_y_top = 0.0;
                        for (int i = 0; i < _reference_image.size(); i++) avg_y_top += _reference_image[i]->y_top;
                        avg_y_top /= (double)_reference_image.size();
                        image.fix_capheight = true;
                        image.capheight = avg_y_top;
                }
                image.Read(file[i].c_str(), filename2.c_str());
                double **binary = image.CombinedHysteresis(filename4.c_str());
                if (!binary)
                {
                        cout << "error during image processing, probably no closed contour detected..." << endl;
                        exit(0);
                }
                /* writing info data */
                ofstream imagedata(filename3.c_str(), ios::out);
                ImageInfo *info = new ImageInfo();
                info->a = image.MeasureA(binary);
                info->b = image.MeasureB(binary);
                info->h = image.MeasureH(binary);
                info->y_top = image.MeasureYTop(binary);
                info->y_bottom = image.MeasureYBottom(binary);
                info->x_left = image.MeasureXLeft(binary);
                info->x_right = image.MeasureXRight(binary);
                info->path = file[i];
                result.push_back(info);
                imagedata << "a[px]" << "\t" << info->a << endl;
                imagedata << "b[px]" << "\t" << info->b << endl;
                imagedata << "h[px]" << "\t" << info->h << endl;
                imagedata << "y_top[px]" << "\t" << info->y_top << endl;
                imagedata << "y_bottom[px]" << "\t" << info->y_bottom << endl;
                imagedata << "x_left[px]" << "\t" << info->x_left << endl;
                imagedata << "x_right[px]" << "\t" << info->x_right << endl;
                imagedata.close();
                /* save point set for regression */
                contour.push_back(image.PointCloud(binary));
                /* clean up */
                image.FreeMatrix(binary);
                image.FreeImage();
        }
        cout << endl;
        return result;
}

void FinalizeContours(string preamble, double scaling, vector<PointSet*> &contour)
{
        /* scale all images with determined scaling-factor */
        for (int i = 0; i < contour.size(); i++)
        {
                /* scaling */
                for (int k = 0; k < contour[i]->size(); k++)
                {
                        (*contour[i])[k].first /= scaling;
                        (*contour[i])[k].second /= scaling;
                }
                /* pipe all contours in unit-scaling */
                string filename = GLOBAL_OUT_FOLDER + preamble + "_" + to_string<int>(i) + "_contour.dat";
                ofstream output(filename.c_str(), ios::out);
                for (int k = 0; k < contour[i]->size(); k++)
                {
                        output 
                        << setw(20) << (*contour[i])[k].first 
                        << setw(20) << (*contour[i])[k].second << endl;
                }
                output.close();
        }
}

void InitializeSolver()
{
        _reference_image = LoadImages("reference", REFERENCE_SHAPE, _reference_shape, false);
        _elastic_image = LoadImages("elastic", ELASTIC_SHAPE, _elastic_shape, true);
}

void FinalizeSolver()
{
        for(int i = 0; i < _reference_shape.size(); i++) delete _reference_shape[i];
        _reference_shape.clear();
        for(int i = 0; i < _elastic_shape.size(); i++) delete _elastic_shape[i];
        _elastic_shape.clear();
        for(int i = 0; i < _reference_image.size(); i++) delete _reference_image[i];
        _reference_image.clear();
        for(int i = 0; i < _elastic_image.size(); i++) delete _elastic_image[i];
        _elastic_image.clear();
}

void SimpleLaplace(string filename, double pressure, double density, double scaling)
{
        string extended_filename = GLOBAL_OUT_FOLDER;
        extended_filename.append(filename);
        extended_filename.append(".dat");
        
        LaplaceData shape;
        shape.SetParameters(pressure, density, scaling);
        SolveLaplace(&shape);
        PlotLaplace(extended_filename.c_str(), &shape);
        
        if (GNUPLOT_SUPPORT) vislap(filename.c_str(), 1.0);
}

void SimpleHooke(string filename, double pressure_laplace, double density_laplace, double scaling_laplace, double pressure, double poisson, double compression)
{
        string extended_filename1 = GLOBAL_OUT_FOLDER;
        extended_filename1.append(filename);
        extended_filename1.append(".dat");
        string extended_filename2 = OUT_FOLDER;
        extended_filename2.append("reference_");
        extended_filename2.append(filename);
        extended_filename2.append(".dat");
        
        LaplaceData ref_shape;
        ref_shape.SetParameters(pressure_laplace, density_laplace, scaling_laplace);
        SolveLaplace(&ref_shape);
        ScaleLaplace(&ref_shape);
        PlotLaplace(extended_filename2.c_str(), &ref_shape);
        
        HookeData shape;
        shape.SetLaplace(&ref_shape);
        shape.SetParameters(pressure, poisson, compression);
        shape.SetInitialConditions(0.0, 0.0, 0.0, 1.0);
        
        bool success = SingleShooting(&shape);
        
        if (!success)
        {
                cout << "No solution found.." << flush << endl;
                return;
        }
        PlotHooke(extended_filename1.c_str(), &shape);
        
        string ref_filename = "reference_";
        ref_filename.append(filename);
        
        if (GNUPLOT_SUPPORT)
        {
                vislap(ref_filename.c_str(), 1.0);
                vishok(filename.c_str(), 1.0, 1.0);
        }
}

LaplaceData *LaplaceMean(double pressure, double density, double scaling)
{
        /* average over all reference contours */
        double mean_pressure = 0.0;
        double mean_density = 0.0;
        double mean_scaling = 0.0;
        
        LaplaceData shape[_reference_shape.size()];
        
        string tmp = GLOBAL_OUT_FOLDER;
        tmp.append("reference.dat");
        ofstream global_output(tmp.c_str(), ios::out);
        global_output 
        << setw(20) << "index" 
        << setw(20) << "area [m^2]" 
        << setw(20) << "volume [m^3]" 
        << setw(20) << "pressure [N/m^2]" 
        << setw(20) << "density [kg/m^3]" 
        << setw(20) << "scaling [m]" 
        << setw(20) << "error [px]" 
        << setw(20) << "tension [N/m]" 
        << endl;
        
        /* calculate average scaling factor */
        double cap_convert = 0.0;
        for (int i = 0; i < _reference_shape.size(); i++) cap_convert += _reference_image[i]->a/_reference_image[i]->b;
        cap_convert /= (double)_reference_shape.size();
        double inner_capillary_m = EXPERIMENT_CAPDIAMETER * cap_convert;
        
        cout << setw(40) << "Reference..." << endl;
        #pragma omp parallel for shared(shape) reduction(+: mean_pressure, mean_density, mean_scaling)
        for (int i = 0; i < _reference_shape.size(); i++)
        {
                shape[i].SetParameters(pressure, density, scaling);
                SolveLaplace(&shape[i]);
                LaplaceFit(&shape[i], _reference_shape[i]);
                mean_pressure += shape[i].GetPressure();
                mean_density += shape[i].GetDensity();
                mean_scaling += shape[i].GetScaling();
        }
        
        mean_pressure /= (double)_reference_shape.size();
        mean_density /= (double)_reference_shape.size();
        mean_scaling /= (double)_reference_shape.size();
        
        /* surface tension */
        double surface_tension = pow(inner_capillary_m, 2)*EXPERIMENT_DENSITY*9.81/mean_density;
        
        LaplaceData *mean = new LaplaceData();
        string filename = GLOBAL_OUT_FOLDER;
        filename.append("basic_shape.dat");
        mean->SetParameters(mean_pressure, mean_density, mean_scaling);
        SolveLaplace(mean);
        ScaleLaplace(mean);
        mean->V0 = Volume(mean);
        mean->A0 = Area(mean);
        mean->surface_tension = surface_tension;
        mean->inner_diameter = inner_capillary_m;
        PlotLaplace(filename.c_str(), mean);
        
        FinalizeContours("reference", mean_scaling, _reference_shape);
        FinalizeContours("elastic", mean_scaling, _elastic_shape);
        
        /* plot laplace shapes */
        for (int k = 0; k < _reference_shape.size(); k++)
        {
                string filename = GLOBAL_OUT_FOLDER + "reference_" + to_string<int>(k) + "_shape.dat";
                ScaleLaplace(&shape[k]);
                PlotLaplace(filename.c_str(), &shape[k]);
                double local_inner_capillary_m = EXPERIMENT_CAPDIAMETER * _reference_image[k]->a/_reference_image[k]->b;
                double tension = pow(local_inner_capillary_m, 2)*EXPERIMENT_DENSITY*9.81/shape[k].GetDensity();
                global_output 
                << setw(20) << k 
                << setw(20) << Area(&shape[k])*pow(local_inner_capillary_m, 2)
                << setw(20) << Volume(&shape[k])*pow(local_inner_capillary_m, 3)
                << setw(20) << shape[k].GetPressure()*tension/local_inner_capillary_m
                << setw(20) << shape[k].GetDensity()*tension/pow(local_inner_capillary_m, 2)/9.81
                << setw(20) << local_inner_capillary_m /* from dimensionless to m */
                << setw(20) << mean->GetConversion()*Error(&shape[k], _reference_shape[k]) 
                << setw(20) << tension
                << endl;
                string out_image = GLOBAL_OUT_FOLDER + "final_reference_" + to_string<int>(k) + ".png";
                Image::WriteLaplace(out_image, _reference_image[k], &shape[k], mean->GetConversion());
        }
        global_output << endl << endl;
        cout << setw(40) << "done!" << flush << endl << endl;
        
        global_output << endl << endl;
        global_output << setw(20) << "pressure [N/m^2]" << setw(20) << "density [kg/m^3]" << setw(20) << "scaling-factor" << endl;
        global_output << setw(20) << mean_pressure*surface_tension/inner_capillary_m << setw(20) << mean_density*surface_tension/pow(inner_capillary_m, 2)/9.81 << setw(20) << mean_scaling/SCALE_IMAGE << endl << endl;
        global_output << setw(20) << "surface tension" << setw(20) << surface_tension << "N/m" << endl << endl;
        global_output.close();
        
        cout << setw(20) << "pressure [N/m^2]" << setw(20) << "density [kg/m^3]" << setw(20) << "scaling-factor" << endl;
        cout << setw(20) << mean_pressure*surface_tension/inner_capillary_m << setw(20) << mean_density*surface_tension/pow(inner_capillary_m, 2)/9.81 << setw(20) << mean_scaling/SCALE_IMAGE << endl << endl;
        cout << setw(20) << "surface tension" << setw(20) << surface_tension << "N/m" << endl << endl;
        
        if (GNUPLOT_SUPPORT)
        {
                mergecnt();
                mergelap(mean->inner_diameter);
        }
        
        return mean;
}

void Sequence(LaplaceData *reference, double pressure, double poisson, double compression)
{
        string tmp = GLOBAL_OUT_FOLDER;
        tmp.append("sequence.dat");
        ofstream sequence(tmp.c_str(), ios::out);
        tmp = GLOBAL_OUT_FOLDER;
        tmp.append("wrinkling.dat");
        ofstream wrinkling(tmp.c_str(), ios::out);
        
        /* prepare full html report */
        string report_name = GLOBAL_OUT_FOLDER + "report.html";
        ofstream report(report_name.c_str(), ios::out);
        report << "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" >" << endl;
        report << "<html><head><title>Capsule Analysis</title><style>body { width: 75%; margin: 0 auto; font-size:14; } .bottom_space{margin-bottom:2cm;} h1 { padding-left: 0px; color: #004488; } h2 { padding-left: 20px; color: #004488; }</style></head><body>" << endl;
        report << "<h1>Laplace-Young Analysis</h1>" << endl;
        report << "<table class=\"bottom_space\">" << endl;
        report << "<tr><td>surface tension</td><td>" << reference->surface_tension << "N/m</td></tr>" << endl;
        report << "<tr><td>pressure at the apex</td><td>" << reference->GetPressure()*reference->surface_tension/reference->inner_diameter << "N/m^2</td></tr>" << endl;
        report << "<tr><td>density difference</td><td>" << reference->GetDensity()*reference->surface_tension/pow(reference->inner_diameter, 2)/9.81 << "kg/m^3</td></tr>" << endl;
        report << "<tr><td>scaling</td><td>" << reference->GetConversion() << "</td></tr>" << endl;
        report << "<tr><td>surface</td><td>" << reference->A0*pow(reference->inner_diameter, 2) << "</td></tr>" << endl;
        report << "<tr><td>volume</td><td>" << reference->V0*pow(reference->inner_diameter, 3) << "</td></tr>" << endl;
        report << "</table>" << endl;
        vector<string> final_reference_image = get_files_containing("final_reference");
        for (int i = 0; i < final_reference_image.size(); i++) 
        {
                vector<string> decomposed = explode(final_reference_image[i], "/");
                report << "<img src=\"" << decomposed[decomposed.size()-1] << "\"/>" << endl;
        }
        
        double V0 = reference->V0;
        double A0 = reference->A0;
        double surface_tension = reference->surface_tension;
        double length_scale = reference->inner_diameter;
        double pref1, pref2, pref3;
        pref1 = reference->GetPressure();
        pref2 = reference->GetDensity();
        pref3 = reference->GetScaling();
        double pdef1, pdef2, pdef3;
        pdef1 = pressure;
        pdef2 = poisson;
        pdef3 = compression;
        double conversion = reference->GetConversion();
        
        LaplaceData udef;
        HookeData shape;
        
        cout << setw(40) << "Sequence..." << flush << endl;
        sequence << "#" << setw(11) << "image" << setw(12) << "p[N/m^2]" << setw(12) << "nu" << setw(12) << "k2d[N/m]" << setw(12) << "rms[px]" << setw(12) << "rel.vol." << setw(12) << "rel.area" << setw(12) << "y2d[N/m]" << setw(12) << "Lw[m]" << setw(12) << "tau_s[N/m]" << setw(12) << "lambda[m]" << setw(12) << "eb[Nm]" << setw(12) << "fppl" << setw(12) << "H0[m]"<< endl;
        wrinkling << "#" << setw(19) << "image" << setw(20) << "rel. volume" << setw(20) << "L0" << setw(20) << "s1/L0" << setw(20) << "s2/L0" << endl;
        report << "<h1>Elastic Analysis</h1>" << endl;
        
        for (int i = 0; i < _elastic_shape.size(); i++)
        {
                string filename = GLOBAL_OUT_FOLDER + "elastic_" + to_string<int>(i) + "_shape.dat";
                udef.SetParameters(pref1, pref2, pref3);
                SolveLaplace(&udef); 
                shape.SetLaplace(&udef);
                shape.SetParameters(pdef1, pdef2, pdef3);
                shape.SetInitialConditions(0.0, 0.0, 0.0, 1.0);
                bool success = true;
                success = SingleShooting(&shape);
                
                if (!success) 
                {
                        cout << "vary initial conditions.. no shooting success.." << flush << endl;
                        return;
                }
                
                if (NELDERMEAD_PREFITTING) NelderMead(&shape, _elastic_shape[i]);
                HookeFit(&shape, _elastic_shape[i], FIT_PRESSURE, FIT_POISSON, FIT_COMPRESSION);		
                
                /* determine wrinkling domain by finding zero crossings of circumferential tension */
                double s1, s2;
                double wrinkle_length = 0.0;
                double avg_tau = 0.0;
                double lambda = 0.0;
                double bending_modulus = 0.0;
                if (WrinklingRegion(&shape, s1, s2))
                {
                        cout << "Wrinkling analysis...\n" << endl;
                        wrinkling << setw(20) << i << setw(20) << Volume(&shape) / V0 << setw(20) << shape.L0 << setw(20) << s1 << setw(20) << s2 << flush << endl;
                        wrinkle_length = s2 - s1;
                        avg_tau = AverageWrinklingTension(&shape, s1, s2);
                        Image image;
                        string wrinkle_image_name = GLOBAL_OUT_FOLDER + "wrinkle_" + to_string<int>(i) + "_image.dat";
                        image.Read(ELASTIC_SHAPE[i].c_str(), wrinkle_image_name.c_str());
                        if(fabs(WRINKLING_WAVELENGTH) < 1.0e-14) lambda = image.WrinkleWavelength(&shape);
                        else lambda = WRINKLING_WAVELENGTH;
                        image.FreeImage();
                        bending_modulus = pow(lambda*length_scale/reference->GetConversion(), 4)*avg_tau*surface_tension/(16.0*M_PI*M_PI*pow(wrinkle_length*length_scale, 2));
                }
                
                double error = reference->GetConversion()*Error(&shape, _elastic_shape[i]);
                double volume = Volume(&shape)/V0;
                double area = Area(&shape)/A0;
                double eh0 = shape.GetEH0();
                
                if (PARAMETER_TRACING) {
                        pdef1 = shape.GetPressure();
                        pdef2 = shape.GetPoisson();
                        pdef3 = shape.GetCompression();
                }
                
                sequence 
                << setw(12) << i 
                << setw(12) << shape.GetPressure()*surface_tension/length_scale
                << setw(12) << shape.GetPoisson() 
                << setw(12) << shape.GetCompression()*surface_tension
                << setw(12) << error 
                << setw(12) << volume 
                << setw(12) << area
                << setw(12) << shape.GetEH0()*surface_tension
                << setw(12) << wrinkle_length*length_scale
                << setw(12) << avg_tau*surface_tension
                << setw(12) << lambda*length_scale/reference->GetConversion()
                << setw(12) << bending_modulus
                << setw(12) << shape.GetEH0()*surface_tension*pow(MaximumRadius(&shape)*length_scale, 2)/bending_modulus
                << setw(12) << sqrt(12.0*bending_modulus*(1.0-shape.GetPoisson())/(shape.GetEH0()*surface_tension))
                << endl;
                PlotHooke(filename.c_str(), &shape);
                string out_image = GLOBAL_OUT_FOLDER + "final_elastic_" + to_string<int>(i) + ".png";
                Image::WriteHooke(out_image, _elastic_image[i], &shape, reference->GetConversion());
                
                report << "<h2>file: " << out_image << "</h2>" << endl;
                report << "<table><tr><td>" << endl;
                vector<string> decomposed = explode(out_image, "/");
                report << "<img src=\"" << decomposed[decomposed.size()-1] << "\"/></td><td>" << endl;
                report << "<table>" << endl;
                report << "<tr><td>pressure at the apex</td><td>" << shape.GetPressure()*surface_tension/length_scale << "N/m^2</td></tr>" << endl;
                report << "<tr><td>poisson's ratio</td><td>" << shape.GetPoisson() << "</td></tr>" << endl;
                report << "<tr><td>area compression modulus</td><td>" << shape.GetCompression()*surface_tension << "N/m</td></tr>" << endl;
                report << "<tr><td>fit error</td><td>" << error << "px</td></tr>" << endl;
                report << "<tr><td>relative surface</td><td>" << area << "</td></tr>" << endl;
                report << "<tr><td>relative volume</td><td>" << volume << "</td></tr>" << endl;
                report << "<tr><td>young's modulus</td><td>" << shape.GetEH0()*surface_tension << "N/m</td></tr>" << endl;
                report << "<tr><td>wrinkle length</td><td>" << wrinkle_length*length_scale << "m</td></tr>" << endl;
                report << "<tr><td>average meridional tension in wrinkled domain</td><td>" << avg_tau*surface_tension << "N/m</td></tr>" << endl;
                report << "<tr><td>average wrinkle wavelength</td><td>" << lambda*length_scale/reference->GetConversion() << "m</td></tr>" << endl;
                report << "<tr><td>bending modulus</td><td>" << bending_modulus << "Nm</td></tr>" << endl;
                report << "<tr><td>f&oumlppl von karman number</td><td>" << shape.GetEH0()*surface_tension*pow(MaximumRadius(&shape)*length_scale, 2)/bending_modulus << "</td></tr>" << endl;
                report << "<tr><td>thickness of layer</td><td>" << sqrt(12.0*bending_modulus*(1.0-shape.GetPoisson())/(shape.GetEH0()*surface_tension)) << "m</td></tr>" << endl;
                report << "</table></td></table>" << endl;
        }
        
        wrinkling.close();
        sequence.close();
        
        /* finalize report */
        report << "</body></html>" << endl;
        report.close();
        
        cout << setw(40) << "Sequence done!" << flush << endl << endl;
        
        if (GNUPLOT_SUPPORT)
        {
                mergehok(reference->inner_diameter, reference->surface_tension);
                visseq(reference->inner_diameter, reference->surface_tension);
        }
}

void EdgeDetection(string filename)
{
        Image image;
        string out_filename = GLOBAL_OUT_FOLDER + "image.dat";
        string bin_filename = GLOBAL_OUT_FOLDER + "binary.png";
        image.Read(filename.c_str(), out_filename.c_str());
        CHECK_IF_CLOSED = false;
        double **binary = image.CombinedHysteresis(bin_filename.c_str());
        if (!binary)
        {
                cout << "error during image processing, probably no closed contour detected..." << endl;
                exit(0);
        }
}


