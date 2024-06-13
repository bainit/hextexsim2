#include <iostream>
#include <utility>
#include <cmath>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "../includes/vpsc_classes.h"
#include "../includes/vpsc_core.h"
#include "../includes/vpsc_util.h"

using namespace std;

Grain::Grain(double p1, double p, double p2, double w, size_t numssystems) : phi1(p1), phi(p), phi2(p2), weight(w), crsses(numssystems) {

	rotation_matrix = calculate_rotation_matrix_bunge(p1,p,p2);

	T = gsl_multiroot_fsolver_hybrids;
  	multiroot_fsolver = gsl_multiroot_fsolver_alloc (T, 9);

	stress = gsl_vector_alloc(9);
	
}

const Matrix& Grain::get_rotation_matrix() const{

	return rotation_matrix;

}

double Grain::get_crss_at(size_t p) const {

	try{
		return crsses.at(p);

	}catch (out_of_range const& exp){
		cout<<exp.what()<<endl;
	}

}

void Grain::set_crss_at(size_t p, double val) {

	try {
		crsses[p] = val;
	}catch(out_of_range const& exp) {
		cout<<exp.what()<<endl;
	}
}

void Grain::update_rotation(const Matrix& rot) {

	rotation_matrix = rot*rotation_matrix;

}

void Grain::set_stress(gsl_vector* vs) {

	gsl_vector_memcpy(stress,vs);

}

void Grain::print() const {

	cout<<"phi1: "<<get_phi1()<<", phi: "<<get_phi()<<", phi2: "<<get_phi2()<<", weight: "<<get_weight()<<endl;
	cout<<endl;
	cout<<"Rotation matrix: "<<endl;
	rotation_matrix.print();

}

/*
double Grain::get_phi1() const {

	return phi1;
}

double Grain::get_phi() const {

	return phi;

}

double Grain::get_phi2() const {

	return phi2;

}

double Grain::get_weight() const {

	return weight;

}
*/
Grain::~Grain(){

	gsl_multiroot_fsolver_free (multiroot_fsolver);
	gsl_vector_free (stress);

}
