#ifndef VPSC_UTIL_H
#define VPSC_UTIL_H

#include <string>
#include <vector>
#include <cmath>

#include <gsl/gsl_matrix.h>

#include "vpsc_core.h"

class Grain;
class SlipSystemsManager;

using namespace std;

void print_matrix(gsl_matrix const * matrix);

double rad2deg(double rad);
double deg2rad(double deg);

double sign(double d);

void matmul(gsl_matrix const *lm, gsl_matrix const *rm, gsl_matrix *res);

vector<Grain*> load_grains(string pathname, const size_t num_of_slipsystems);
void export_grains(const vector<Grain*> grains, const string pathname);

void load_slip_systems(const string& filepath, SlipSystemsManager & manager);

#endif //VPSC_UTIL_H
