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

#include "image.h"

Image::Image()
{
        fix_capheight = false;
        capheight = 0;
}

Image::~Image()
{
        
}

void Image::FreeImage()
{
        FreeMatrix(foto.matrix);
}

/* used for Graham Scan sorting algorithm */
bool Image::ComparePoints(const pair<double,double> p1, const pair<double,double> p2)
{
        double w1, w2;
        double d1, d2;
        d1 = sqrt(p1.first * p1.first + p1.second * p1.second);
        d2 = sqrt(p2.first * p2.first + p2.second * p2.second);
        w1 = atan2(p1.second, p1.first);
        w2 = atan2(p2.second, p2.first);
        if (fabs(w1 - w2) > 1.0e-15)
        {
                return (w1 < w2);
        }
        else
        {
                return (d1 < d2);
        }
}

/* allocate memory for image on heap */
double **Image::AllocateMatrix()
{
        double **f;
        f = (double**)calloc(foto.width, sizeof(double*));
        for (int i = 0; i < foto.width; i++)
        {
                f[i] = (double*)calloc(foto.height, sizeof(double));
        }
        return f;
}

/* free image memory */
void Image::FreeMatrix(double **matrix)
{
        for (int i = 0; i < foto.width; i++)
        {
                free(matrix[i]);
        }
        free(matrix);
}

double Image::MeasureH(double **image)
{
        double y_top = MeasureYTop(image);
        double y_bottom = MeasureYBottom(image);
        return (foto.height - y_top - y_bottom);
}

double Image::MeasureA(double **image)
{
        double y_top = MeasureYTop(image);
        double **shape = AllocateMatrix();
        /* delete capillary
         * shape will contain capsule without capillary */
        for (int x = 0; x < foto.width; x++)
        {
                for (int y = 0; y < foto.height; y++)
                {
                        /* cut at y = foto.height - 1 - y_top - top_buffer */
                        if (y > foto.height - 1 - (int)y_top - TOP_BUFFER)
                        {
                                shape[x][y] = 0.0;
                        }
                        else
                        {
                                shape[x][y] = image[x][y];
                        }
                }
        }
        int x_left = 0.0, x_right = 0.0;
        /* measure inner capillary diameter */
        /*for (int x = 0; x <= foto.width/2; x++)*/
        for (int x = 0; x <= foto.width; x++)
        {
                if ((int)(shape[x][foto.height - 1 - (int)y_top - TOP_BUFFER] + 0.5) == 1)
                {
                        x_left = x;
                        break;
                }
        } 
        /*for (int x = foto.width - 1; x > foto.width/2; x--)*/
        for (int x = foto.width - 1; x > 0; x--)
        {
                if ((int)(shape[x][foto.height - 1 - (int)y_top - TOP_BUFFER] + 0.5) == 1)
                {
                        x_right = x;
                        break;
                }
        }
        FreeMatrix(shape);
        return fabs((double)x_right - (double)x_left);
}

void Image::DeleteDuplicates(PointSet *v)
{
        PointSet temp;
        bool in_temp;
        for (int i = 0; i < v->size(); i++)
        {
                in_temp = false;
                for (int j = 0; j < temp.size(); j++)
                {
                        if ((fabs(v->at(i).first - temp[j].first) < 1.0e-15) && (fabs(v->at(i).second - temp[j].second) < 1.0e-15)) in_temp = true;
                }
                if (!in_temp) 
                {
                        temp.push_back(v->at(i));
                }
        }
}

/* calculate principal axis and rotate image 
 * y-axis will be aligned with symmetry axis of contour point set */
void Image::ForceSymmetry(PointSet *v)
{
        gsl_matrix *I = gsl_matrix_calloc(2, 2);
        gsl_matrix *evec = gsl_matrix_calloc(2, 2);
        gsl_vector *eval = gsl_vector_calloc(2);
        gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(2);
        
        for (int k = 0; k < 2; k++)
        {
                for (int l = 0; l < 2; l++)
                {
                        double ICom = 0.0;
                        for (int i = 0; i < v->size(); i++)
                        {
                                double x[2] = {(*v)[i].first, (*v)[i].second};
                                if (k == l) ICom += x[0]*x[0] + x[1]*x[1];
                                ICom -= x[k]*x[l];
                        }
                        gsl_matrix_set(I, k, l, ICom);
                }
        }
        
        gsl_eigen_symmv(I, eval, evec, w);
        gsl_eigen_symmv_free(w);
        
        double eval_1 = gsl_vector_get(eval, 0);
        double eval_2 = gsl_vector_get(eval, 1);
        int index = (eval_1 < eval_2) ? 1 : 0;
        double axis[2] = {gsl_matrix_get(evec, 0, index), gsl_matrix_get(evec, 1, index)};
        double ax_nrm = sqrt(axis[0]*axis[0] + axis[1]*axis[1]);
        double arc = acos(axis[1]/ax_nrm);
        if (axis[0] < 0) arc *= (-1.0);
        
        gsl_matrix_free(I);
        gsl_matrix_free(evec);
        gsl_vector_free(eval);
        
        double x_c = 0.0;
        double y_c = 0.0;
        for (int i = 0; i < v->size(); i++) 
        {
                x_c += (*v)[i].first;
                y_c += (*v)[i].second;
        }
        x_c /= (double)v->size();
        y_c /= (double)v->size();
        
        /* translate to zero */
        for (int i = 0; i < v->size(); i++)
        {
                (*v)[i].first -= x_c;
                (*v)[i].second -= y_c;
        }
        
        printf("rotate image by %e degree\n", (arc/(2.0*M_PI))*360.0);
        
        /* rotate image */
        for (int i = 0; i < v->size(); i++)
        {
                double x = (*v)[i].first;
                double y = (*v)[i].second;
                (*v)[i].first = x*cos(arc) - y*sin(arc);
                (*v)[i].second = x*sin(arc) + y*cos(arc);
        }
        
        /* translate to center of mass */
        for (int i = 0; i < v->size(); i++)
        {
                (*v)[i].first += x_c;
                (*v)[i].second += y_c;
        }
        
}

/* extract set of contour points from binary image */
PointSet *Image::PointCloud(double **image)
{
        PointSet cloud;
        double y_top = MeasureYTop(image);
        
        double **shape = AllocateMatrix();
        
        /* x-positions of connection between capillary and cloud */
        double x_cap_left, x_cap_right;
        bool x_cap_left_found = false, x_cap_right_found = false;
        
        /* delete capillary
         * shape will contain capsule without capillary */
        for (int x = 0; x < foto.width; x++)
        {
                for (int y = 0; y < foto.height; y++)
                {
                        /* cut at y = height - 1 - y_top - top_buffer */
                        if (y > foto.height - 1 - (int)y_top - TOP_BUFFER)
                        {
                                shape[x][y] = 0;
                        }
                        else
                        {
                                shape[x][y] = image[x][y];
                        }
                }
        }
        
        /* scan image for contour points */
        for (int y = foto.height - 1; y >= 0; y--)
        {
                /* scan from left side */
                for (int x = 0; x < foto.width; x++)
                {
                        if (shape[x][y] == 1)
                        {
                                if (!x_cap_left_found)
                                {
                                        x_cap_left = (double)x;
                                        x_cap_left_found = true;
                                }
                                cloud.push_back(make_pair((double)x, (double)y));
                                break;
                        }
                }
                /* scan from right side */
                for (int x = foto.width - 1; x >= 0; x--)
                {
                        if (shape[x][y] == 1)
                        {
                                if (!x_cap_right_found)
                                {
                                        x_cap_right = (double)x;
                                        x_cap_right_found = true;
                                }
                                cloud.push_back(make_pair((double)x, (double)y));
                                break;
                        }
                }
        }
        
        /* scan from bottom */
        for (int x = 0; x < foto.width; x++)
        {
                for (int y = 0; y < foto.height; y++)
                {
                        if (shape[x][y] == 1)
                        {
                                cloud.push_back(make_pair((double)x, (double)y));
                                break;
                        }
                }
        }
        
        /* scan from top left */
        for (int x = 0; x <= x_cap_left; x++)
        {
                for (int y = foto.height - 1; y >= 0; y--)
                {
                        if (shape[x][y] == 1)
                        {
                                cloud.push_back(make_pair((double)x, (double)y));
                                break;
                        }
                }
        }
        
        /* scan from top right */
        for (int x = x_cap_right; x < foto.width; x++)
        {
                for (int y = foto.height - 1; y >= 0; y--)
                {
                        if (shape[x][y] == 1)
                        {
                                cloud.push_back(make_pair((double)x, (double)y));
                                break;
                        }
                }
        }
        
        /* optimize symmetry */
        if (FORCE_SYMMETRY) 
        {
                ForceSymmetry(&cloud);       
        }
        
        
        /* delete doubly inserted points from set */
        DeleteDuplicates(&cloud);
        
        double x_s = 0.0, y_s = 0.0;
        int N = 0;
        for (PointSet::iterator it = cloud.begin(); it != cloud.end(); it++)
        {
                x_s += it->first;
                y_s += it->second;
                N++;
        }
        x_s /= (double)N;
        y_s /= (double)N;
        
        for (int i = 0; i < cloud.size(); i++)
        {
                cloud[i].first -= x_s;
                cloud[i].second -= y_s;
        }
        
        /* mirror axis at 45 degrees to prepare for Graham Scan sorting */
        double tmp;
        for (int i = 0; i < cloud.size(); i++)
        {
                tmp = cloud[i].first;
                cloud[i].first = cloud[i].second;
                cloud[i].second = tmp;
        }
        
        /* Graham Scan sorting */
        sort(cloud.begin(), cloud.end(), ComparePoints);
        
        /* scaling, translation and rotation */
        double y_bottom = MeasureYBottom(image);
        for (int i = 0; i < cloud.size(); i++)
        {
                cloud[i].first *= SCALE_IMAGE;
                cloud[i].second *= SCALE_IMAGE;
                cloud[i].first += y_s * SCALE_IMAGE;
                cloud[i].first -= y_bottom * SCALE_IMAGE;
                cloud[i].first -= MeasureH(image) * SCALE_IMAGE;
        }
        
        /* mirror back to normal state */
        for (int i = 0; i < cloud.size(); i++)
        {
                tmp = cloud[i].first;
                cloud[i].first = cloud[i].second;
                cloud[i].second = tmp;
        }
        
        
        /* sparse point set to improve efficiency */
        double max_dist = fabs(MeasureA(image) * SCALE_IMAGE / 2.0);
        PointSet *sparse_region = new PointSet();
        sparse_region->push_back(cloud[0]);
        double sym_trans = 0.0;
        bool keep_checking = true;
        for (int i = 1; i < cloud.size() - 1; i++)
        {
                double dist1 = sqrt((cloud[i-1].first - cloud[i].first) * (cloud[i-1].first - cloud[i].first) + (cloud[i-1].second - cloud[i].second) * (cloud[i-1].second - cloud[i].second)); 
                double dist2 = sqrt((cloud[i+1].first - cloud[i].first) * (cloud[i+1].first - cloud[i].first) + (cloud[i+1].second - cloud[i].second) * (cloud[i+1].second - cloud[i].second)); 
                if (dist1 > max_dist && keep_checking)
                {
                        sym_trans += cloud[i].first;
                        sym_trans += cloud[(i-1+cloud.size())%cloud.size()].first;
                        keep_checking = false;
                }
                if(dist2 > max_dist && keep_checking)
                {
                        sym_trans += cloud[i].first;
                        sym_trans += cloud[(i+1)%cloud.size()].first;
                        keep_checking = false;
                }
                /* take only 1/EVERY of the found points AND points directly at the capillary */
                if(i % EVERY == 0 || dist1 > max_dist || dist2 > max_dist) {
                        sparse_region->push_back(cloud[i]);
                }
        }
        sparse_region->push_back(cloud[cloud.size() - 1]);
        cloud.clear();
        
        /* make points at the capillary symmetric */
        if (FORCE_BOUNDARY_SYMMETRY)
        {
                sym_trans /= 2.0;
                for (int i = 0; i < sparse_region->size(); i++)
                {
                        (*sparse_region)[i].first -= sym_trans;
                }
        }
        
        FreeMatrix(shape);
        
        return sparse_region;
}

double Image::WrinkleWavelength(HookeData *shape)
{
        double **binary = CombinedHysteresis("");
        if (!binary)
        {
                cout << "error during image processing!" << endl;
                exit(0);
        }
        double **wrinkle_binary = CombinedHysteresis("", true);
        Write("./out/wrinkling.png", wrinkle_binary);
        if (!wrinkle_binary)
        {
                cout << "error during image processing!" << endl;
                exit(0);
        }
        double y_top = MeasureYTop(wrinkle_binary);
        double s1, s2;
        int z_upper = 0;
        int z_lower = 0;
        int z_offset = foto.height - (int)(MeasureYBottom(binary));
        double scaling = shape->undeformed->GetConversion();
        double avg_width = 0.0;
        int counter_width = 0;
        if (WrinklingRegion(shape, s1, s2))
        {
                if (s2 > shape->L0) s2 = shape->L0;
                z_upper = (int)(gsl_spline_eval(shape->splines.z_s, s2, shape->splines.acc)*scaling - gsl_spline_eval(shape->splines.z_s, 0.0, shape->splines.acc)*scaling);
                z_lower = (int)(gsl_spline_eval(shape->splines.z_s, s1, shape->splines.acc)*scaling - gsl_spline_eval(shape->splines.z_s, 0.0, shape->splines.acc)*scaling);
        }
        /* delete capillary from binary */
        for (int x = 0; x < foto.width; x++)
        {
                for (int y = 0; y < foto.height; y++)
                {
                        /* cut at y = height - 1 - y_top - top_buffer */
                        if (y > foto.height - 1 - (int)(y_top + 0.5) - TOP_BUFFER)
                        {
                                binary[x][y] = 0;
                        }
                }
        }
        /* scan image and select inner pixels of the wrinkling domain */
        int N = 0;
        int Ne = 0;
        
        int z_0 = max(z_offset - z_upper, 0);
        int z_1 = min(z_offset - z_lower, foto.height);
        
        for (int y = z_0; y <= z_1; y++)
        {
                int xleft = 0;
                int xright = 0;
                /* scan from left side */
                for (int x = 0; x < foto.width; x++)
                {
                        if (binary[x][y] == 1)
                        {
                                xleft = (double)x;
                                break;
                        }
                }
                /* scan from right side */
                for (int x = foto.width - 1; x >= 0; x--)
                {
                        if (binary[x][y] == 1) 
                        {
                                xright = (double)x;
                                break;
                        }
                }
                avg_width += xright - xleft;
                counter_width++;
                for (int x = xleft + 1; x < xright; x++) 
                {
                        N++;
                        if (wrinkle_binary[x][y] == 1) 
                        {
                                Ne++;
                        }
                }
        }
        FreeMatrix(binary);
        FreeMatrix(wrinkle_binary);
        avg_width /= (double)counter_width;
        double rad = avg_width/2.0;
        return M_PI*rad/((double)N/Ne);
}

/* use with ClosedContour
 * recursive coloring algorithm
 * coloring stored in c
 * output parameter : c
 * c should be initialized as c[x][y] = 0 for all x,y
 * input parameter : binary-image 
 * output parameter : c */
void Image::Color(double **image, double **c, int x, int y)
{
        c[x][y] = 1;
        /* proceed if pixel, not visited and not reached the border */
        if (x+1 < foto.width)
        {
                if ((image[x+1][y] == 1) && (c[x+1][y] == 0))
                {
                        Color(image, c, x + 1, y);
                }
        }
        if (x-1 >= 0)
        {
                if ((image[x-1][y] == 1) && (c[x-1][y] == 0))
                {
                        Color(image, c, x - 1, y);
                }
        }
        if (y+1 < foto.height)
        {
                if ((image[x][y+1] == 1) && (c[x][y+1] == 0))
                {
                        Color(image, c, x, y + 1);
                }
        }
        if (y-1 >= 0)
        {
                if ((image[x][y-1] == 1) && (c[x][y-1] == 0))
                {
                        Color(image, c, x, y - 1);
                }
        }
        if ((x+1 < foto.width) && (y+1 < foto.height))
        {
                if ((image[x+1][y+1] == 1) && (c[x+1][y+1] == 0))
                {
                        Color(image, c, x + 1, y + 1);
                }
        }
        if ((x-1 >= 0) && (y+1 < foto.height))
        {
                if ((image[x-1][y+1] == 1) && (c[x-1][y+1] == 0))
                {
                        Color(image, c, x - 1, y + 1);
                }
        }
        if ((x+1 < foto.width) && (y-1 >= 0))
        {
                if ((image[x+1][y-1] == 1) && (c[x+1][y-1] == 0))
                {
                        Color(image, c, x + 1, y - 1);
                }
        }
        if ((x-1 >= 0) && (y-1 >= 0))
        {
                if ((image[x-1][y-1] == 1) && (c[x-1][y-1] == 0))
                {
                        Color(image, c, x - 1, y - 1);
                }
        }
}

bool Image::ClosedContour(double **image)
{
        /* visitor field */
        double **c = AllocateMatrix();
        
        bool result = false;
        
        /* search start and color contour */
        /*for (int x = 0; x < foto.width/2; x++)*/
        for (int x = 0; x < foto.width; x++)
        {
                if (image[x][foto.height - 1] == 1)
                {
                        Color(image, c, x, foto.height - 1);
                        break;
                }
        }
        
        /* search end and check if colored */
        /*for (int x = foto.width - 1; x >= foto.width/2; x--)*/
        for (int x = foto.width - 1; x > 0; x--)
        {
                if (image[x][foto.height - 1] == 1)
                {
                        if (c[x][foto.height - 1] == 1)
                                result = true;
                        break;
                }
        }
        
        FreeMatrix(c);
        
        return result;
}

double Image::MeasureXLeft(double **image)
{
        double y_top = MeasureYTop(image);
        y_top -= 5.0;
        int N = 5;
        int step = y_top / N;
        double x_l;
        double x_l_sum = 0.0;
        for (int i = 0; i <= N; i++)
        {
                int y = foto.height - 1 - i * step;
                /* measure from left side */
                x_l = 0.0;
                /*for (int x = 0; x < foto.width/2; x++)*/
                for (int x = 0; x < foto.width; x++)
                {
                        if (image[x][y] == 0)
                        {
                                x_l += 1.0;
                        } 
                        else
                        {
                                break;
                        }
                }
                x_l_sum += x_l;
        }
        x_l_sum /= (double)N + 1.0;
        return x_l_sum;
}

double Image::MeasureXRight(double **image)
{
        double y_top = MeasureYTop(image);
        y_top -= 5.0;
        int N = 5;
        int step = y_top / N;
        double x_r;
        double x_r_sum = 0.0;
        for (int i = 0; i <= N; i++)
        {
                int y = foto.height - 1 - i * step;
                /* measure from right side */
                x_r = 0.0;
                /*for (int x = foto.width - 1; x > foto.width/2; x--)*/
                for (int x = foto.width - 1; x > 0; x--)
                {
                        if (image[x][y] == 0)
                        {
                                x_r += 1.0;
                        }
                        else
                        {
                                break;
                        }
                }
                x_r_sum += x_r;
        } 
        x_r_sum /= (double)N + 1.0;
        return x_r_sum;
}

double Image::MeasureB(double **image)
{
        double x_l = MeasureXLeft(image);
        double x_r = MeasureXRight(image);
        return (foto.width - 1 - x_r - x_l);
}

double Image::MeasureYBottom(double **image)
{
        int y_bottom = foto.height;
        int it;
        for (int x = 0; x < foto.width; x++)
        {
                it = 0;
                for (int y = 0; y < foto.height; y++)
                {
                        if (image[x][y] == 0)
                        {
                                it++;
                        }
                        else
                        {
                                break;
                        }
                }
                if (it < y_bottom)
                {
                        y_bottom = it;
                }
        }
        return (double)y_bottom;
}

/* return mean from left and right side measurement */
double Image::MeasureYTop(double **image)
{
        if (fix_capheight) return (double)capheight;
        else return (1.0/2.0) * ((double)MeasureYTopLeft(image) + (double)MeasureYTopRight(image));
}

int Image::MeasureYTopLeft(double **image)
{
        /* store x-positions */
        vector<int> xx;
        double mean;
        int d1, d2;
        for (int y = foto.height - 1; y >= 0; y--)
        {
                for (int x = 0; x < foto.width; x++)
                {
                        if (image[x][y] == 1)
                        {
                                xx.push_back(x);
                                break;
                        }
                }
                mean = 0.0;
                for (int i = 0; i < xx.size(); i++)
                {
                        mean += (double)xx[i];
                }
                mean /= xx.size();
                if (fabs((double)xx[xx.size()-1] - mean) > THRESHOLD_CAPILLARY) return foto.height - 1 - y;
        }
        xx.clear();
        return 0;
}

int Image::MeasureYTopRight(double **image)
{
        vector<int> xx;
        double mean;
        int d1, d2;
        for (int y = foto.height - 1; y >= 0; y--)
        {
                for (int x = foto.width - 1; x >= 0; x--)
                {
                        if (image[x][y] == 1)
                        {
                                /* x:absolute -> width - 1 - x: distance from right border */
                                xx.push_back(foto.width - 1 - x);
                                break;
                        }
                }
                mean = 0.0;
                for (int i = 0; i < xx.size(); i++)
                {
                        mean += (double)xx[i];
                }
                mean /= xx.size();
                if (fabs((double)xx[xx.size()-1] - mean) > THRESHOLD_CAPILLARY) return foto.height - 1 - y;
        }
        xx.clear();
        return 0;
}

double **Image::ApplyFilter(double **image, double **filter, int radius)
{
        /* allocate memory for the filtered image */
        double **f = AllocateMatrix();
        int xxm, yym;
        
        for (int y = foto.height - 1; y >= 0; y--)
        {
                for (int x = 0; x < foto.width; x++)
                {
                        f[x][y] = 0.0;
                        for (int xx = -radius; xx <= radius; xx++)
                        {
                                for (int yy = -radius; yy <= radius; yy++)
                                {
                                        /* boundary condition for x + xx, y + yy */
                                        xxm = x + xx;
                                        yym = y + yy;
                                        if (xxm < 0)
                                        {
                                                xxm = 0;
                                        }
                                        if (xxm > foto.width - 1)
                                        {
                                                xxm = foto.width - 1;
                                        }
                                        if (yym < 0)
                                        {
                                                yym = 0;
                                        }
                                        if (yym > foto.height - 1)
                                        {
                                                yym = foto.height - 1;
                                        }
                                        // apply filter
                                        f[x][y] += image[xxm][yym] * filter[radius + xx][radius + yy];
                                }
                        }
                }
        }
        return f;
}

/* combines images r = r_min, ...,r_max
 * trace edges with hysteresis */
double **Image::CombinedHysteresis(string filename, const bool wrinkles)
{
        if (!wrinkles) 
        {
                foto.t_high = T_HIGH;
                foto.t_low = T_LOW;
        } else 
        {
                foto.t_high = 0.10;
                foto.t_low = 0.05;
                CHECK_IF_CLOSED = FALSE;
        }
        
        bool closed = false;
        vector<double**> gradients;
        vector<double**> directions;
        
        /* declare combined image */
        double **com = AllocateMatrix();
        double **bin = AllocateMatrix();
        
        /* calculate and store gradients for r = r_min, ...,r_max */
        for (int i = R_MIN; i <= R_MAX; i++)
        {
                double **grad = AllocateMatrix();
                double **gradphi0 = AllocateMatrix();
                ConvertImage(grad, gradphi0, i);
                gradients.push_back(grad);
                if (wrinkles) 
                {
                        for (int m = 0; m < foto.width; m++) 
                        {
                                for (int n = 0; n < foto.height; n++) 
                                {
                                        gradphi0[m][n] = 0.0;
                                }
                        }
                }
                directions.push_back(gradphi0);
        }
        
        while (!closed)
        {
                /* superposition of images */
                for (int i = R_MIN; i <= R_MAX; i++)
                {
                        ApplyHysteresis(gradients[i - R_MIN], directions[i - R_MIN], bin);
                        for (int y = 0; y < foto.height; y++)
                        {
                                for (int x = 0; x < foto.width; x++)
                                {
                                        /* adding */
                                        com[x][y] += bin[x][y];
                                        /* scaling */
                                        com[x][y] = (int)(com[x][y] > 0);
                                }
                        }
                }
                
                /* check if contour is closed */
                if (CHECK_IF_CLOSED)
                {
                        if (ClosedContour(com))
                        {
                                closed = true;
                        }
                        else
                        {
                                /* adjust thresholds in order to get closed structure */
                                if (foto.t_low <= 0.05) return NULL;
                                foto.t_high = foto.t_high - 0.01;
                                foto.t_low = foto.t_low - 0.01;
                        }
                }
                else
                {
                        closed = true;
                }
        }
        if(!filename.empty()) Write(filename.c_str(), com);
        FreeMatrix(bin);
        for (int i = 0; i < gradients.size(); i++) FreeMatrix(gradients[i]);
        for (int i = 0; i < directions.size(); i++) FreeMatrix(directions[i]);
        gradients.clear();
        directions.clear();
        if (wrinkles) 
        {
                foto.t_high = T_HIGH;
                foto.t_low = T_LOW;
                CHECK_IF_CLOSED = TRUE;
        }
        return com;
}

void Image::Gradient(double **grad, double **gradphi, double **gradphi0, double **dx, double **dy)
{
        /* gradient direction */
        double dir;
        
        for (int y = 0; y < foto.height; y++)
        {
                for (int x = 0; x < foto.width; x++)
                {
                        /* gradient-matrix */
                        grad[x][y] = sqrt(dx[x][y] * dx[x][y] + dy[x][y] * dy[x][y]);
                        
                        /* cases for atan */
                        if ((dx[x][y] == 0) && (dy[x][y] == 0))
                        {
                                gradphi[x][y] = 0.0;
                        }
                        else if ((dx[x][y] == 0) && (dy[x][y] > 0))
                        {
                                gradphi[x][y] = M_PI/2.0;
                        }
                        else 
                        {
                                gradphi[x][y] = atan(dy[x][y] / dx[x][y]);
                        }
                        
                        /* round on pi/4, i.e. discretize direction
                         * 4 directions
                         * x-axis shows in 0 direction 0°
                         * y-axis shows in 2 direction 90° */
                        
                        /* current direction */
                        dir = gradphi[x][y];
                        
                        /* 4 intervals */
                        if ((dir > 3*M_PI/8) || (dir < (-1)*3*M_PI/8))
                        {
                                gradphi0[x][y] = 0;
                        }
                        if ((dir < 3*M_PI/8) && (dir > M_PI/8))
                        {
                                gradphi0[x][y] = 1;
                        }
                        if ((dir < M_PI/8) && (dir > (-1)*M_PI/8))
                        {
                                gradphi0[x][y] = 2;
                        }
                        if ((dir > (-1)*3*M_PI/8) && (dir < -(1)*M_PI/8))
                        {
                                gradphi0[x][y] = 3;
                        }
                }
        }
}

/* r : radius of smoothing operator
 * read image
 * smooth image
 * calculate derivatives
 * output parameter : grad, gradphi0 */
void Image::ConvertImage(double **grad, double **gradphi0, int r)
{
        double **smooth;
        double **dx;
        double **dy;
        double **gradphi;
        double **filter_gauss;
        double **filter_gradx;
        double **filter_grady;
        
        filter_gauss = (double**)calloc(2*r + 1, sizeof(double*));
        for (int i = 0; i < 2*r + 1; i++)
        {
                filter_gauss[i] = (double*)calloc(2*r + 1, sizeof(double));
        }
        
        filter_gradx = (double**)calloc(3, sizeof(double*));
        filter_grady = (double**)calloc(3, sizeof(double*));
        for (int i = 0; i < 3; i++)
        {
                filter_gradx[i] = (double*)calloc(3, sizeof(double));
                filter_grady[i] = (double*)calloc(3, sizeof(double));
        }
        
        /* allocate gaussian matrix */
        for (int xx = -r; xx <= r; xx++)
        {
                for (int yy = -r; yy <= r; yy++)
                {
                        filter_gauss[r + xx][r + yy] = exp(-(xx*xx + yy*yy)/GAUSSIAN_SIGMA/GAUSSIAN_SIGMA);
                }
        }
        
        /* allocate sobel matrices */
        filter_gradx[0][0] = -3.0;  filter_gradx[0][1] = 0.0;    filter_gradx[0][2] = 3.0;
        filter_gradx[1][0] = -10.0; filter_gradx[1][1] = 0.0;    filter_gradx[1][2] = 10.0;
        filter_gradx[2][0] = -3.0;  filter_gradx[2][1] = 0.0;    filter_gradx[2][2] = 3.0;
        
        filter_grady[0][0] = -3.0;  filter_grady[0][1] = -10.0;  filter_grady[0][2] = -3.0;
        filter_grady[1][0] = 0.0;   filter_grady[1][1] = 0.0;    filter_grady[1][2] = 0.0;
        filter_grady[2][0] = 3.0;   filter_grady[2][1] = 10.0;   filter_grady[2][2] = 3.0;
        
        /* apply filter */
        smooth = ApplyFilter(foto.matrix, filter_gauss, r);
        /* calculate derivatives */
        dx = ApplyFilter(smooth, filter_gradx, 1);
        dy = ApplyFilter(smooth, filter_grady, 1);
        
        /* allocate field for gradient direction */
        gradphi = AllocateMatrix();
        
        /* calculate gradient */
        Gradient(grad, gradphi, gradphi0, dx, dy);  
        /* output grad and gradphi0 calculated */
        
        /* write pictures to files (debugging purposes only) */
        // Write("./out/2-smooth.png", smooth);
        // Write("./out/3a-dx.png", dx);
        // Write("./out/3b-dy.png", dy);
        // Write("./out/4a-gradient.png", grad);
        // Write("./out/4b-gradientdirection.png", gradphi);
        // Write("./out/4c-gradientdirectiondiscrete.png", gradphi0);
        
        FreeMatrix(smooth);
        FreeMatrix(dx);
        FreeMatrix(dy);
        FreeMatrix(gradphi);
        
        for (int i = 0; i < 3; i++)
        {
                free(filter_gradx[i]);
                free(filter_grady[i]);
        }
        for (int i = 0; i < 2*r + 1; i++)
        {
                free(filter_gauss[i]);
        }
        
        free(filter_gauss);
        free(filter_gradx);
        free(filter_grady);
}

void Image::TraceEdge(double **grad, double **gradphi0, double **binary, double max_brightness)
{
        /* check for maxima above t_high */
        for (int x = 0; x < foto.width; x++)
        {
                for (int y = 0; y < foto.height; y++)
                {
                        if ((grad[x][y] > foto.t_high * max_brightness) && (CheckMaxInDir(grad, gradphi0, x, y)) && (binary[x][y] == 0))
                        {
                                binary[x][y] = 1;
                        }
                }
        }
        
        /* check for maxima between t_low and t_high with neighbored edges */
        for (int x = 0; x < foto.width; x++)
        {
                for (int y = 0; y < foto.height; y++)
                {
                        if ((grad[x][y] < foto.t_high * max_brightness) && (grad[x][y] > foto.t_low * max_brightness) && (CheckMaxInDir(grad, gradphi0, x, y)) && (binary[x][y] == 0) && (NeighborColored(binary, x, y)))
                        {
                                binary[x][y] = 1;
                        }
                }
        }
}

bool Image::NeighborColored(double **c, int x, int y)
{
        bool result = false;
        
        if (x+1 < foto.width)
        {
                if (c[x+1][y] == 1)
                {
                        result = true;
                }
        }
        if (x-1 >= 0)
        {
                if (c[x-1][y] == 1)
                {
                        result = true;
                }
        }
        if (y+1 < foto.height)
        {
                if (c[x][y+1] == 1)
                {
                        result = true;
                }
        }
        if (y-1 >= 0)
        {
                if (c[x][y-1] == 1)
                {
                        result = true;
                }
        }
        if ((x+1 < foto.width) && (y+1 < foto.height))
        {
                if (c[x+1][y+1] == 1)
                {
                        result = true;
                }
        }
        if ((x-1 >= 0) && (y+1 < foto.height))
        {
                if (c[x-1][y+1] == 1)
                {
                        result = true;
                }
        }
        if ((x+1 < foto.width) && (y-1 >= 0))
        {
                if (c[x+1][y-1] == 1)
                {
                        result = true;
                }
        }
        if ((x-1 >= 0) && (y-1 >= 0))
        {
                if (c[x-1][y-1] == 1)
                {
                        result = true;
                }
        }
        
        return result;
}

void Image::ApplyHysteresis(double **grad, double **gradphi0, double **binary)
{
        double max_brightness = MaxBrightness(grad);
        TraceEdge(grad, gradphi0, binary, max_brightness);
}

double Image::MaxBrightness(double **image)
{
        double max = 0.0;
        for (int y = 0; y < foto.height; y++)
        {
                for (int x = 0; x < foto.width; x++)
                {
                        if (image[x][y] > max) max = image[x][y];
                }
        }
        return max;
}

bool Image::CheckMaxInDir(double **grad, double **gradphi0,  int x, int y)
{
        return ((gradphi0[x][y] == 0 && CheckMax0(grad, x, y)) 
        || (gradphi0[x][y] == 1 && CheckMax1(grad, x, y)) 
        || (gradphi0[x][y] == 2 && CheckMax2(grad, x, y)) 
        || (gradphi0[x][y] == 3 && CheckMax3(grad, x, y)));
}

bool Image::CheckMax(double **grad, int x, int y)
{
        return (CheckMax0(grad, x, y) || CheckMax1(grad, x, y) || CheckMax2(grad, x, y) || CheckMax3(grad, x, y));
}

bool Image::CheckMax3(double** gradient, int x, int y)
{
        bool max = true;
        if ((x > 0) && (y < foto.height - 1))
        {
                if (gradient[x-1][y+1] > gradient[x][y]) max = false;
        }
        if ((x < foto.width - 1) && (y > 0))
        {
                if (gradient[x+1][y-1] > gradient[x][y]) max = false;
        }
        return max;
}

bool Image::CheckMax2(double **gradient, int x, int y)
{
        bool max = true;
        if (y < foto.height - 1)
        {
                if (gradient[x][y+1] > gradient[x][y]) max = false;
        }
        if (y > 0)
        {
                if (gradient[x][y-1] > gradient[x][y]) max = false;
        }
        return max;
}

bool Image::CheckMax1(double **gradient, int x, int y)
{
        bool max = true;
        if ((x > 0) && (y > 0))
        {
                if (gradient[x-1][y-1] > gradient[x][y]) max = false;
        }
        if ((x < foto.width - 1) && (y < foto.height - 1))
        {
                if (gradient[x+1][y+1] > gradient[x][y]) max = false;
        }
        return max;
}

bool Image::CheckMax0(double **gradient, int x, int y)
{
        bool max = true;
        if (x > 0)
        {
                if (gradient[x-1][y] > gradient[x][y]) max = false;
        }
        if (x < foto.width - 1)
        {
                if (gradient[x+1][y] > gradient[x][y]) max = false;
        }
        return max;
}

int Image::Write(const char* filename, double **i)
{
        double g = MaxBrightness(i);
        Mat output(foto.height, foto.width, CV_64FC1);
        for (int x = 0; x < foto.width; x++)
        {
                for (int y = 0; y < foto.height; y++)
                {
                        output.at<double>(foto.height - 1 - y,x) = 255.0 * (i[x][y] / g);
                }
        }
        imwrite(filename, output);
        output.release();
        return 0;
}

int Image::Read(const char* filename_in, const char* filename_out)
{
        IplImage *im = cvLoadImage(filename_in);
        uchar *data;
        
        foto.step = im->widthStep;
        foto.width = im->width;
        foto.height = im->height;
        foto.channels = im->nChannels;
        
        data = (uchar*)im->imageData;
        
        foto.matrix = AllocateMatrix();
        
        ofstream output(filename_out, ios::out);
        
        for (int y = 0; y < foto.height; y++)
        {
                for (int x = 0; x < foto.width; x++)
                {
                        foto.matrix[x][foto.height - 1 - y] = 0.0;
                        for (int k = 0; k < foto.channels; k++)
                        {
                                foto.matrix[x][foto.height - 1 - y] += (double)((int)data[y*foto.step+x*foto.channels+k]);
                        }
                        foto.matrix[x][foto.height - 1 - y] /= foto.channels;
                        output << "\t" << foto.matrix[x][foto.height - 1 - y];
                }
                output << "\n";
        }
        
        output.close();
        cvReleaseImage(&im);
        return 0;
} 

void Image::WriteLaplace(string filename_out, ImageInfo *img_data, LaplaceData *shape, double conversion)
{
        IplImage *im = cvLoadImage(img_data->path.c_str());
        double scaling = conversion;
        int width = im->width;
        int height = im->height;
        
        int offset_r = (int)((img_data->x_left + width - img_data->x_right)/2.0);
        int offset_z = height - (int)(img_data->y_bottom);
        
        int r0 = 0.0;
        int z0 = 0.0;
        int r;
        int z;
        
        cvDrawLine(im, cvPoint(img_data->x_left, 0), cvPoint(img_data->x_left, img_data->y_top), CV_RGB(0, 255, 0), 2);
        cvDrawLine(im, cvPoint(width - 1 - img_data->x_right, 0), cvPoint(width - 1 - img_data->x_right, img_data->y_top), CV_RGB(0, 255, 0), 2);
        cvDrawLine(im, cvPoint(img_data->x_left, img_data->y_top), cvPoint(width - 1 - img_data->x_right, img_data->y_top), CV_RGB(0, 255, 0), 2);
        
        for (double s = 1.0e-2; s <= shape->L0; s += 1.0e-2)
        {
                r = (int)(gsl_spline_eval(shape->splines.r_s, s, shape->splines.acc)*scaling);
                z = (int)(gsl_spline_eval(shape->splines.z_s, s, shape->splines.acc)*scaling - gsl_spline_eval(shape->splines.z_s, 0.0, shape->splines.acc)*scaling);
                cvDrawLine(im, cvPoint(offset_r + r0, offset_z - z0), cvPoint(offset_r + r, offset_z - z), CV_RGB(255, 0, 0), 2);
                r0 = r;
                z0 = z;
        }
        
        r0 = 0.0;
        z0 = 0.0;
        for (double s = 1.0e-2; s <= shape->L0; s += 1.0e-2)
        {
                r = -(int)(gsl_spline_eval(shape->splines.r_s, s, shape->splines.acc)*scaling);
                z = (int)(gsl_spline_eval(shape->splines.z_s, s, shape->splines.acc)*scaling - gsl_spline_eval(shape->splines.z_s, 0.0, shape->splines.acc)*scaling);
                cvDrawLine(im, cvPoint(offset_r + r0, offset_z - z0), cvPoint(offset_r + r, offset_z - z), CV_RGB(255, 0, 0), 2);
                r0 = r;
                z0 = z;
        }
        
        double xi = (int)img_data->b * (0.001 / EXPERIMENT_CAPDIAMETER);
        double xi2 = xi * 0.1;
        double offset_bar = height - (int)(img_data->y_bottom)/2;
        
        cvDrawLine(im, cvPoint(offset_r - xi/2.0, offset_bar), cvPoint(offset_r + xi/2.0, offset_bar), CV_RGB(0, 0, 0), 2);
        cvDrawLine(im, cvPoint(offset_r - xi/2.0, offset_bar - xi2/2.0), cvPoint(offset_r - xi/2.0, offset_bar + xi2/2.0), CV_RGB(0, 0, 0), 2);
        cvDrawLine(im, cvPoint(offset_r + xi/2.0, offset_bar - xi2/2.0), cvPoint(offset_r + xi/2.0, offset_bar + xi2/2.0), CV_RGB(0, 0, 0), 2);
        
        cvSaveImage(filename_out.c_str(), im);
        cvReleaseImage(&im);
}

void Image::WriteHooke(string filename_out, ImageInfo *img_data, HookeData *shape, double conversion)
{
        IplImage *im = cvLoadImage(img_data->path.c_str());
        double scaling = conversion;
        int width = im->width;
        int height = im->height;
        
        int offset_r = (int)((img_data->x_left + width - img_data->x_right - 1)/2.0);
        int offset_z = height - (int)(img_data->y_bottom);
        
        int r0 = 0.0;
        int z0 = 0.0;
        int r;
        int z;
        
        double s1, s2;
        if (WrinklingRegion(shape, s1, s2))
        {
                int z_upper = (int)(gsl_spline_eval(shape->splines.z_s, s2, shape->splines.acc)*scaling - gsl_spline_eval(shape->splines.z_s, 0.0, shape->splines.acc)*scaling);
                int z_lower = (int)(gsl_spline_eval(shape->splines.z_s, s1, shape->splines.acc)*scaling - gsl_spline_eval(shape->splines.z_s, 0.0, shape->splines.acc)*scaling);
                cvDrawLine(im, cvPoint(0, offset_z - z_upper), cvPoint(width, offset_z - z_upper), CV_RGB(0, 0, 255), 2);
                cvDrawLine(im, cvPoint(0, offset_z - z_lower), cvPoint(width, offset_z - z_lower), CV_RGB(0, 0, 255), 2);
        }
        
        cvDrawLine(im, cvPoint(img_data->x_left, 0), cvPoint(img_data->x_left, img_data->y_top), CV_RGB(0, 255, 0), 2);
        cvDrawLine(im, cvPoint(width - 1 - img_data->x_right, 0), cvPoint(width - 1 - img_data->x_right, img_data->y_top), CV_RGB(0, 255, 0), 2);
        cvDrawLine(im, cvPoint(img_data->x_left, img_data->y_top), cvPoint(width - 1 - img_data->x_right, img_data->y_top), CV_RGB(0, 255, 0), 2);
        
        for (double s = 1.0e-2; s <= shape->L0; s += 1.0e-2)
        {
                r = (int)(gsl_spline_eval(shape->splines.r_s, s, shape->splines.acc)*scaling);
                z = (int)(gsl_spline_eval(shape->splines.z_s, s, shape->splines.acc)*scaling - gsl_spline_eval(shape->splines.z_s, 0.0, shape->splines.acc)*scaling);
                if(gsl_spline_eval(shape->splines.t_phi, s, shape->splines.acc)*scaling < 0.0) 
                {
                        cvDrawLine(im, cvPoint(offset_r + r0, offset_z - z0), cvPoint(offset_r + r, offset_z - z), CV_RGB(0, 0, 255), 2);
                } else 
                {
                        cvDrawLine(im, cvPoint(offset_r + r0, offset_z - z0), cvPoint(offset_r + r, offset_z - z), CV_RGB(255, 0, 0), 2);
                }
                r0 = r;
                z0 = z;
        }
        
        r0 = 0.0;
        z0 = 0.0;
        for (double s = 1.0e-2; s <= shape->L0; s += 1.0e-2)
        {
                r = -(int)(gsl_spline_eval(shape->splines.r_s, s, shape->splines.acc)*scaling);
                z = (int)(gsl_spline_eval(shape->splines.z_s, s, shape->splines.acc)*scaling - gsl_spline_eval(shape->splines.z_s, 0.0, shape->splines.acc)*scaling);
                if(gsl_spline_eval(shape->splines.t_phi, s, shape->splines.acc)*scaling < 0.0) 
                {
                        cvDrawLine(im, cvPoint(offset_r + r0, offset_z - z0), cvPoint(offset_r + r, offset_z - z), CV_RGB(0, 0, 255), 2);
                } else 
                {
                        cvDrawLine(im, cvPoint(offset_r + r0, offset_z - z0), cvPoint(offset_r + r, offset_z - z), CV_RGB(255, 0, 0), 2);
                }
                r0 = r;
                z0 = z;
        }
        
        double xi = (int)img_data->b * (0.001 / EXPERIMENT_CAPDIAMETER);
        double xi2 = xi * 0.1;
        double offset_bar = height - (int)(img_data->y_bottom)/2;
        
        cvDrawLine(im, cvPoint(offset_r - xi/2.0, offset_bar), cvPoint(offset_r + xi/2.0, offset_bar), CV_RGB(0, 0, 0), 2);
        cvDrawLine(im, cvPoint(offset_r - xi/2.0, offset_bar - xi2/2.0), cvPoint(offset_r - xi/2.0, offset_bar + xi2/2.0), CV_RGB(0, 0, 0), 2);
        cvDrawLine(im, cvPoint(offset_r + xi/2.0, offset_bar - xi2/2.0), cvPoint(offset_r + xi/2.0, offset_bar + xi2/2.0), CV_RGB(0, 0, 0), 2);
        
        cvSaveImage(filename_out.c_str(), im);
        cvReleaseImage(&im);
}


