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

#ifndef SCRIPT_H
#define SCRIPT_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <boost/filesystem.hpp>

#include "interface.h"

using namespace std;
namespace flsys = boost::filesystem;

bool compare_filenames(const string file1, const string file2);

void visseq(double length_scale, double surface_tension);
void vislap(const char *file, double length_scale);  
void vishok(const char *file, double length_scale, double surface_tension);
void mergelap(double length_scale);
void mergehok(double length_scale, double surface_tension);
void mergecnt();
void mergeimg();

vector<string> get_files_containing(string str); 
vector<string> get_files_containing(string str1, string str2);
void merge_files(string str, vector<string> files);  
int file_num(string str);








#endif
