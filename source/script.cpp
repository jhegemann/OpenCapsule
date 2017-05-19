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

#include "script.h"

/* compare different filenames */
bool compare_filenames(const string file1, const string file2)
{
        int num1 = file_num(file1);
        int num2 = file_num(file2);
        return (num1 < num2);
}

/* plot elastic sequence in dependence of relative volume V/V0 */
void visseq(double length_scale, double surface_tension)
{
        flsys::remove("./tmp/seq.gp");
        ofstream seq("./tmp/seq.gp", ios::out);
        
        seq << "reset" << endl;
        seq << "set terminal postscript enhanced solid color" << endl;
        
        seq << "set output './global_out/sequence_pressure.eps'" << endl;
        seq << "set xlabel 'relative volume'" << endl;
        seq << "set ylabel 'pressure [N/m^2]'" << endl;
        seq << "plot './global_out/sequence.dat' u 6:2 t '' w p pt 54 lc 1" << endl;
        
        seq << "set output './global_out/sequence_poisson.eps'" << endl;
        seq << "set xlabel 'relative volume'" << endl;
        seq << "set ylabel 'poisson'" << endl;
        seq << "plot './global_out/sequence.dat' u 6:3 t '' w p pt 54 lc 1" << endl;
        
        seq << "set output './global_out/sequence_compression.eps'" << endl;
        seq << "set xlabel 'relative volume'" << endl;
        seq << "set ylabel 'area compression modulus [N/m]'" << endl;
        seq << "plot './global_out/sequence.dat' u 6:4 t '' w p pt 54 lc 1" << endl;
        
        seq << "set output './global_out/sequence_young.eps'" << endl;
        seq << "set xlabel 'relative volume'" << endl;
        seq << "set ylabel 'young modulus [N/m]'" << endl;
        seq << "plot './global_out/sequence.dat' u 6:8 t '' w p pt 54 lc 1" << endl;
        
        seq << "set output './global_out/sequence_wrinkling_length.eps'" << endl;
        seq << "set xlabel 'relative volume'" << endl;
        seq << "set ylabel 'wrinkling length [m]'" << endl;
        seq << "plot './global_out/sequence.dat' u 6:9 t '' w p pt 54 lc 1" << endl;
        
        seq << "set output './global_out/sequence_wrinkling_tension.eps'" << endl;
        seq << "set xlabel 'relative volume'" << endl;
        seq << "set ylabel 'wrinkling tension [N/m]'" << endl;
        seq << "plot './global_out/sequence.dat' u 6:10 t '' w p pt 54 lc 1" << endl;
        
        seq << "set output './global_out/sequence_wrinkling_wavelength.eps'" << endl;
        seq << "set xlabel 'relative volume'" << endl;
        seq << "set ylabel 'wrinkling wavelength [m]'" << endl;
        seq << "plot './global_out/sequence.dat' u 6:11 t '' w p pt 54 lc 1" << endl;
        
        seq << "set output './global_out/sequence_bending.eps'" << endl;
        seq << "set xlabel 'relative volume'" << endl;
        seq << "set ylabel 'bending modulus [Nm]'" << endl;
        seq << "plot './global_out/sequence.dat' u 6:12 t '' w p pt 54 lc 1" << endl;
        
        seq << "set output './global_out/sequence_foeppl.eps'" << endl;
        seq << "set xlabel 'relative volume'" << endl;
        seq << "set ylabel 'foeppl von karman number'" << endl;
        seq << "plot './global_out/sequence.dat' u 6:13 t '' w p pt 54 lc 1" << endl;
        
        seq << "set output './global_out/sequence_thickness.eps'" << endl;
        seq << "set xlabel 'relative volume'" << endl;
        seq << "set ylabel 'layer thickness [m]'" << endl;
        seq << "plot './global_out/sequence.dat' u 6:14 t '' w p pt 54 lc 1" << endl;
        
        seq.close();
        
        int ret = system("gnuplot ./tmp/seq.gp &> /dev/null");
}

/* plot single laplace-form in SI units */
void vislap(const char *file, double length_scale)
{
        flsys::remove("./tmp/lap.gp");
        ofstream lap("./tmp/lap.gp", ios::out);
        
        lap << "reset" << endl;
        lap << "set terminal postscript enhanced solid color" << endl;
        lap << "set output './global_out/" << file << ".eps'" << endl;
        lap << "set key bottom right" << endl;
        lap << "set xlabel 'x [mm]'" << endl;
        lap << "set ylabel 'y [mm]'" << endl;
        lap << "set title '" << file << " z(r)'" << endl;
        lap << "plot './out/" << file << ".dat' u ($2*" << length_scale << "):($3*" << length_scale << ") t 'spline' w l lw 2 lt 7" << endl;
        lap << "set title '" << file << " r(s)'" << endl;
        lap << "plot './out/" << file << ".dat' u ($1*" << length_scale << "):($2*" << length_scale << ") t 'spline' w l lw 2 lt 7" << endl;
        lap << "set title '" << file << " z(s)'" << endl;
        lap << "plot './out/" << file << ".dat' u ($1*" << length_scale << "):($3*" << length_scale << ") t 'spline' w l lw 2 lt 7" << endl;
        lap << "set title '" << file << " {/Symbol \\131}(s)'" << endl;
        lap << "plot './out/" << file << ".dat' u ($1*" << length_scale << "):4 t 'spline' w l lw 2 lt 7" << endl;
        
        lap.close();
        
        int ret = system("gnuplot ./tmp/lap.gp &> /dev/null");
}

// plot single hooke-form
void vishok(const char *file, double length_scale, double surface_tension)
{
        flsys::remove("./tmp/hok.gp");
        ofstream hok("./tmp/hok.gp", ios::out);
        
        hok << "reset" << endl;
        hok <<  "set terminal postscript enhanced solid color" << endl;
        hok << "set output './global_out/" << file << ".eps'" << endl;
        hok << "set key bottom right" << endl;
        
        hok << "set title '" << file << " z(r)'" << endl;
        hok << "plot './out/" << file << ".dat' u ($2*" << length_scale << "):($3*" << length_scale << ") t 'spline' w l lw 2 lt 7" << endl;
        hok << "set title '" << file << " r(s)'" << endl;
        hok << "plot './out/" << file << ".dat' u ($1*" << length_scale << "):($2*" << length_scale << ") t 'spline' w l lw 2 lt 7" << endl;
        hok << "set title '" << file << " z(s)'" << endl;
        hok << "plot './out/" << file << ".dat' u ($1*" << length_scale << "):($3*" << length_scale << ") t 'spline' w l lw 2 lt 7" << endl;
        hok << "set title '" << file << " {/Symbol \\131}(s)'" << endl;
        hok << "plot './out/" << file << ".dat' u ($1*" << length_scale << "):4 t 'spline' w l lw 2 lt 7" << endl;
        hok << "set title '" << file << " {/Symbol t}_s(s)'" << endl;
        hok << "plot './out/" << file << ".dat' u ($1*" << length_scale << "):($5*" << surface_tension << ") t 'spline' w l lw 2 lt 7" << endl;
        hok << "set title '" << file << " {/Symbol l}_s(s)'" << endl;
        hok << "plot './out/" << file << ".dat' u ($1*" << length_scale << "):6 t 'spline' w l lw 2 lt 7" << endl;
        hok << "set title '" << file << " {/Symbol l}_{/Symbol f}(s)'" << endl;
        hok << "plot './out/" << file << ".dat' u ($1*" << length_scale << "):7 t 'spline' w l lw 2 lt 7" << endl;
        hok << "set title '" << file << " {/Symbol k}_{/Symbol f}(s)'" << endl;
        hok << "plot './out/" << file << ".dat' u ($1*" << length_scale << "):($8/" << length_scale << ") t 'spline' w l lw 2 lt 7" << endl;
        hok << "set title '" << file << " {/Symbol t}_{/Symbol f}(s)'" << endl;
        hok << "plot './out/" << file << ".dat' u ($1*" << length_scale << "):($9*" << surface_tension << ") t 'spline' w l lw 2 lt 7" << endl;
        
        hok.close();
        
        int ret = system("gnuplot ./tmp/hok.gp &> /dev/null");
}

/* merge laplace data into one file and plot solutions together with contour */
void mergelap(double length_scale)
{
        vector<string> red = get_files_containing("reference", "shape");
        merge_files("reference_shapes.dat", red);
        vector<string> cnt = get_files_containing("reference", "contour");
        
        flsys::remove("/tmp/ref.gp");
        // gnu script
        ofstream ref("./tmp/ref.gp", ios::out);
        ref << "reset" << endl;
        ref << "set terminal postscript enhanced solid color" << endl;
        ref << "set output './global_out/reference_shapes.eps'" << endl;
        ref << "set key right bottom" << endl;
        
        for (int i = 0; i < red.size(); i++)
        {
                string title = red[i];
                replace_all(title, "_", " ");
                replace_all(title, " ", "\\_");
                ref << "set title '" << title << " z(r)'" << endl;
                ref << "plot '" << cnt[i] << "' u ($1*" << length_scale << "):($2*" << length_scale << ") t 'contour' w p pt 7 lc 1 ps 1.0, '" << red[i] << "' u ($2*" << length_scale << "):($3*" << length_scale << ") t 'theory' w l lw 3 lt 7" << endl;
                ref << "set title '" << title << " r(s)'" << endl;
                ref << "plot '" << red[i] << "' u ($1*" << length_scale << "):($2*" << length_scale << ") t '' w l lw 2 lt 7" << endl;
                ref << "set title '" << title << " z(s)'" << endl;
                ref << "plot '" << red[i] << "' u ($1*" << length_scale << "):($3*" << length_scale << ") t '' w l lw 2 lt 7" << endl;
                ref << "set title '" << title << " {/Symbol \\131}(s)'" << endl;
                ref << "plot '" << red[i] << "' u ($1*" << length_scale << "):4 t '' w l lw 2 lt 7" << endl;
        }
        
        ref << "set output './global_out/basic_shape.eps'" << endl;
        ref << "set title 'reference z(r)'" << endl;
        ref << "plot './global_out/basic_shape.dat' u ($2*" << length_scale << "):($3*" << length_scale << ") t '' w l lw 2 lt 7" << endl;
        ref << "set title 'reference r(s)'" << endl;
        ref << "plot './global_out/basic_shape.dat' u ($1*" << length_scale << "):($2*" << length_scale << ") t '' w l lw 2 lt 7" << endl;
        ref << "set title 'reference z(s)'" << endl;
        ref << "plot './global_out/basic_shape.dat' u ($1*" << length_scale << "):($3*" << length_scale << ") t '' w l lw 2 lt 7" << endl;
        ref << "set title 'reference {/Symbol \\131}(s)'" << endl;
        ref << "plot './global_out/basic_shape.dat' u ($1*" << length_scale << "):4 t '' w l lw 2 lt 7" << endl;
        ref.close();
        
        int ret = system("gnuplot ./tmp/ref.gp &> /dev/null");
}

/* merge hooke data into one file and plot solutions together with contour */
void mergehok(double length_scale, double surface_tension)
{
        vector<string> red = get_files_containing("elastic", "shape");
        vector<string> cnt = get_files_containing("elastic", "contour");
        merge_files("elastic_shapes.dat", red);
        
        // gnu script
        flsys::remove("./tmp/def.gp");
        ofstream def("./tmp/def.gp", ios::out);
        def << "reset" << endl;
        def << "set terminal postscript enhanced solid color" << endl;
        def << "set output './global_out/elastic_shapes.eps'" << endl;
        def << "set key right bottom" << endl;
        
        for (int i = 0; i < red.size(); i++)
        {
                string title = red[i];
                replace_all(title, "_", " ");
                replace_all(title, " ", "\\_");
                def << "set title '" << title << " z(r)'" << endl;
                def << "plot '" << cnt[file_num(red[i])] << "' u ($1*" << length_scale << "):($2*" << length_scale << ") t 'contour' w p pt 7 lc 1 ps 1.0, '" << red[i] << "' u ($2*" << length_scale << "):($3*" << length_scale << ") t 'theory' w l lw 3 lt 7" << endl;
                def << "set title '" << title << " r(s)'" << endl;
                def << "plot '" << red[i] << "' u ($1*" << length_scale << "):($2*" << length_scale << ") t '' w l lw 2 lt 7" << endl;
                def << "set title '" << title << " z(s)'" << endl;
                def << "plot '" << red[i] << "' u ($1*" << length_scale << "):($3*" << length_scale << ") t '' w l lw 2 lt 7" << endl;
                def << "set title '" << title << " {/Symbol \\131}(s)'" << endl;
                def << "plot '" << red[i] << "' u ($1*" << length_scale << "):4 t '' w l lw 2 lt 7" << endl;
                def << "set title '" << title << " {/Symbol t}_s(s)'" << endl;
                def << "plot '" << red[i] << "' u ($1*" << length_scale << "):($5*" << surface_tension << ") t '' w l lw 2 lt 7" << endl;
                def << "set title '" << title << " {/Symbol l}_s(s)'" << endl;
                def << "plot '" << red[i] << "' u ($1*" << length_scale << "):6 t '' w l lw 2 lt 7" << endl;
                def << "set title '" << title << " {/Symbol l}_{/Symbol f}(s)'" << endl;
                def << "plot '" << red[i] << "' u ($1*" << length_scale << "):7 t '' w l lw 2 lt 7" << endl;
                def << "set title '" << title << " {/Symbol k}_{/Symbol f}(s)'" << endl;
                def << "plot '" << red[i] << "' u ($1*" << length_scale << "):($8/" << length_scale << ") t '' w l lw 2 lt 7" << endl;
                def << "set title '" << title << " {/Symbol t}_{/Symbol f}(s)'" << endl;
                def << "plot '" << red[i] << "' u ($1*" << length_scale << "):($9*" << surface_tension << ") t '' w l lw 2 lt 7" << endl;
        }
        def.close();
        
        int ret = system("gnuplot ./tmp/def.gp &> /dev/null");
}

/* merge contour data into one file */
void mergecnt()
{
        vector<string> red = get_files_containing("contour");
        merge_files("contours.dat", red);
}

/* merge image data into one file */
void mergeimg()
{
        vector<string> red = get_files_containing("imageinfo");
        merge_files("imginfo.dat", red);
}

/* return files from GLOBAL_OUT_FOLDER containg some certain string */
vector<string> get_files_containing(string str)
{
        flsys::path dir(GLOBAL_OUT_FOLDER);
        vector<string> all;
        vector<string> red;
        flsys::directory_iterator end_itr;
        for (flsys::directory_iterator itr(dir); itr != end_itr; ++itr)
        {
                all.push_back(itr->path().string());
        }
        for (int i = 0; i < all.size(); i++)
        {
                size_t pos = all[i].find(str);
                if (pos != string::npos)
                {
                        red.push_back(all[i]);
                }
        }
        all.clear();
        sort(red.begin(), red.end(), &compare_filenames);
        return red;
}

vector<string> get_files_containing(string str1, string str2)
{
        vector<string> query = get_files_containing(str1);
        vector<string> result;
        for(int i = 0; i < query.size(); i++) 
        {
                size_t pos = query[i].find(str2);
                if (pos != string::npos) 
                {
                        result.push_back(query[i]);
                }
        }
        return result;
}

/* merge content of files listed in vector files into one file */
void merge_files(string str, vector<string> files)
{
        string path = "./global_out/";
        path.append(str);
        ofstream merge(path.c_str(), ios::out);
        for (int i = 0; i < files.size(); i++)
        {
                ifstream incoming(files[i].c_str(), ios::in);
                string line;
                merge << "#" << files[i] << endl;
                while (getline(incoming, line) && (line != ""))
                {
                        merge << line << endl;
                }
                merge << endl << endl << endl << endl;
                incoming.close();
        }
        merge.close();
}

/* returns n from filename {reference,elastic}_{n}_name */
int file_num(string str)
{
        string strnum = str.substr(str.find("/") + 1);
        strnum = strnum.substr(strnum.find("/") + 1);
        strnum = strnum.substr(strnum.find("_") + 1);
        strnum = strnum.substr(0, strnum.find("_"));
        int num = atoi(strnum.c_str());
        return num;
}

