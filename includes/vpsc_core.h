#ifndef VPSC_CORE_H
#define VPSC_CORE_H

#include <vector>

#include <gsl/gsl_vector.h>

#include "vpsc_classes.h"
#include "vpsc_util.h"

#define BUNGE_NOTATION 1

//Rotation matrix
Matrix calculate_rotation_matrix(double phi1, double phi, double phi2, int notation);
Matrix calculate_rotation_matrix_bunge(double phi1, double phi, double phi2);

//Compatibile notation x||x y||y* z||z
Matrix get_transformation_matrix(double ctoa); 

Vector hkil2uvw_plane(double h, double k, double i, double l, const Matrix& inv_trans_matrix);
Vector hkil2uvw_dir(double h, double k, double i, double l, const Matrix& inv_trans_matrix);

//Calculates deviatoric stress/strain form the given tensor represented by the matrix
Matrix get_deviatoric(const Matrix& s);

//Calculates von Mises stress, stress must be deviatoric
double evmises_stress(const Matrix& stress);

//Calculates von Mises strain, strain must be deviatoric
double evmises_strain(const Matrix& strain);

//double dot product of two 2nd order tensors 3x3 in matrix notation
double ddot(const Matrix& r, const Matrix& l);

double calculate_shear_strain(const Matrix& stress, const SlipSystem& ssys, double strain_rate);

//Calcualtes rotation matrix after deformation in active slip systems in crystal refernce system
Matrix calc_rm_after_def(const vector<double>& shear_strains, const SlipSystemsManager& manager);

//Calculates the hardening of a grain according to voce law
void voce_grain_hardening(Grain& g, const vector<double>& sstrains, const SlipSystemsManager& manager);

//Calculate shear strains
vector<double> calculate_shear_strains(const Matrix& stress, const Grain& g, const SlipSystemsManager& manager, double strainrate);
vector<double> calculate_shear_strains_helper(const gsl_vector* stress, const Grain& g, const SlipSystemsManager& manager, double strain_rate);
double sum_all_shear_strains(vector<double>& shear_strains);

double calc_shear_strains_sum(const SlipSystemsManager& manager, const vector<double>& shear_strains, size_t m, size_t n);

//Minimization function for finding shear strains during deformation based on stress
int strain_minimization(const gsl_vector* e, void* params, gsl_vector* f);

int simulations(Grain& grain, const Matrix& def, SlipSystemsManager& manager, const double strain_rate, size_t maxiter, double residual);

#endif //VPSC_CORE_H
