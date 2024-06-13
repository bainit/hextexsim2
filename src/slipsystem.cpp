#include <iostream>
#include <string>
#include <cmath>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "../includes/vpsc_classes.h"

using namespace std;

SlipSystem::SlipSystem(initializer_list<double> plane, initializer_list<double> direction, double c, double mm, string l, double t0, double t1, double th0, double th1, initializer_list<double> lh) : slip_plane{plane}, slip_direction{direction,Vector::Orientation::vertical}, crss{c}, m{mm}, label{l}, tau_0{t0}, tau_1{t1}, theta_0{th0}, theta_1{th1}, latent_hardening{lh} {

	calculate();

}

SlipSystem::SlipSystem(Vector& plane, Vector& dir, double c, double mm, string l, double t0, double t1, double th0, double th1, initializer_list<double> lh) : slip_plane{plane}, slip_direction{dir}, crss{c}, m{mm}, label{l}, tau_0{t0}, tau_1{t1}, theta_0{th0}, theta_1{th1}, latent_hardening{lh} {

	calculate();

}

void SlipSystem::calculate(){

	slip_plane.normalize();
	slip_direction.normalize();

	deformation_matrix = slip_direction*slip_plane;

	Matrix def_matrix_transpose = deformation_matrix.give_transpose();

	symetry_def_mat = 0.5*(deformation_matrix+def_matrix_transpose);
	asymetry_def_mat = 0.5*(deformation_matrix-def_matrix_transpose);	

}

size_t SlipSystem::get_lh_size() const {

	return latent_hardening.size();

}

double SlipSystem::get_lh_parameter(size_t p) const {

	try {
		return latent_hardening.at(p);

	}catch (out_of_range const& exp) {

		cout<<exp.what()<<endl;
	}

}

SlipSystem::~SlipSystem() {

}
