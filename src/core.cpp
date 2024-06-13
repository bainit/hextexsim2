#include <iostream>
#include <cmath>
#include <vector>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "../includes/vpsc_util.h"
#include "../includes/vpsc_classes.h"
#include "../includes/vpsc_core.h"

using namespace std;

Matrix calculate_rotation_matrix(double phi1, double phi, double phi2, int notation) {

	if(notation == BUNGE_NOTATION) {

        	Matrix Z1{cos(phi1),sin(phi1),0,-sin(phi1),cos(phi1),0,0,0,1};
        	Matrix X{1,0,0,0,cos(phi),sin(phi),0,-sin(phi),cos(phi)};
		Matrix Z2{cos(phi2), sin(phi2),0,-sin(phi2), cos(phi2),0,0,0,1};
		
		Matrix operation1 = Z2*X;
		Matrix operation2 = operation1*Z1;

		return operation2;
	}

	throw runtime_error("calculate_rotation_matrix::Notation not implemented");

}

Matrix calculate_rotation_matrix_bunge(double phi1, double phi, double phi2) {

	return calculate_rotation_matrix(phi1,phi,phi2,BUNGE_NOTATION);
}

//Compatibile notation: x||x y||y* z||z
Matrix get_transformation_matrix(double ctoa) {
	
	constexpr static double a = sqrt(3.0)/2.0;
	Matrix m{1,-0.5,0, 0,a,0, 0,0,ctoa};

	return m;

}

Vector hkil2uvw_plane(double h, double k, double i, double l, const Matrix& inv_trans_matrix) {

	Vector p{h,k,l};

	Vector u = p*inv_trans_matrix;
	u.normalize();
	
	return u;

}

Vector hkil2uvw_dir(double h, double k, double i, double l, const Matrix& trans_matrix) {

	Vector d{{h-i,k-i,l},Vector::Orientation::vertical};
	
	Vector v = trans_matrix*d;
	v.normalize();
	
	return v;

}

Matrix get_deviatoric(const Matrix& s) {

	constexpr static double PARAM = 1.0/3.0;

	double p = PARAM*(s.m11()+s.m22()+s.m33());

	Matrix dev{s.m11()-p,s.m12(),s.m13(),s.m21(),s.m22()-p,s.m23(),s.m31(),s.m32(),s.m33()-p};

	return dev;

}

double evmises_stress(const Matrix& stress) {

	constexpr static double PARAM = 3.0/2.0;

	return sqrt(PARAM*ddot(stress,stress));

}

double evmises_strain(const Matrix& strain) {

	constexpr static double PARAM = 2.0/3.0;
	
	return sqrt(PARAM*ddot(strain,strain));

}

double ddot(const Matrix& l, const Matrix& r) {

        double result = 0;
        for(int i=0;i<3;i++) {
                for(int j=0;j<3;j++) {
                        result += l.get(i,j)*r.get(i,j);
                }
        }

        return result;

}

double calculate_shear_strain(const Matrix& stress, const SlipSystem& system, double strain_rate) {

	double tau = ddot(system.get_def_mat(),stress);
	double m = 1.0/system.get_m();
	double gamma = strain_rate*pow(abs(tau/system.get_tau_0()),m)*sign(tau);

	return gamma;


}

//Calculates rotation matrix after deformation in active slip systems in crystal reference system
Matrix calc_rm_after_def(const vector<double>& shear_strains, const SlipSystemsManager& manager){

	if(shear_strains.size() != manager.get_size()) {

		throw runtime_error("calc_rm_after_def(), number of shear strains does not match number of slip systems");

	}

	Matrix rm;
		
	for(int i=0;i<shear_strains.size();i++) {
		const SlipSystem& system = manager.get_slipsystem(i);
	
		rm  = rm + shear_strains.at(i)*system.get_def_mat();

	}
	//Matrix I = {1,0,0,0,1,0,0,0,1};
	rm = Settings::get_instance().IDENTITY_MATRIX + 0.5*(rm-rm.give_transpose());
	//rm = I + 0.5*(rm-rm.give_transpose());
	return rm;

}

void voce_grain_hardening(Grain& g, const vector<double>& sstrains, const SlipSystemsManager& manager) {

	for(size_t i=0;i<g.get_crsses_size();i++) {
		double lh_val = 0;
		const SlipSystem& system = manager.get_slipsystem(i);
		for(size_t j=0;j<system.get_lh_size();j++) {
			lh_val += g.get_crss_at(j) * sstrains.at(j) * system.get_lh_parameter(j); 
		}
			
		double crss_val = system.get_tau_0() + (system.get_tau_1()+system.get_theta_1()*lh_val)*(1-exp(-system.get_theta_0()/system.get_tau_0()*lh_val));
		g.set_crss_at(i,crss_val);
	}

}

vector<double> calculate_shear_strains(const Matrix& stress, const Grain& g, const SlipSystemsManager& manager, double strain_rate) {

	vector<double> rsstrains;

	for(size_t i=0;i<manager.get_size();i++) {

		const SlipSystem& system = manager.get_slipsystem(i);
		
		double tau = ddot(stress,system.get_sym());
		double im = 1/system.get_m();
		double tau_0 = g.get_crss_at(i);		

		double p = pow(abs(tau/tau_0),im);
		double t = tau/abs(tau);

		double gamma = strain_rate*pow(abs(tau/tau_0),im)*(tau/abs(tau));		
		
		rsstrains.push_back(gamma);

	}	

	return rsstrains;

}

vector<double> calculate_shear_strains_helper(const gsl_vector* stress, const Grain& grain, const SlipSystemsManager& manager, double strain_rate) {

	double s11 = gsl_vector_get(stress,0);
        double s12 = gsl_vector_get(stress,1);
        double s13 = gsl_vector_get(stress,2);
        double s21 = gsl_vector_get(stress,3);
        double s22 = gsl_vector_get(stress,4);
        double s23 = gsl_vector_get(stress,5);
        double s31 = gsl_vector_get(stress,6);
        double s32 = gsl_vector_get(stress,7);
        double s33 = gsl_vector_get(stress,8);

        Matrix rstress{s11,s12,s13,s21,s22,s23,s31,s32,s33};
        Matrix dev_stress = get_deviatoric(rstress);
        return calculate_shear_strains(dev_stress, grain, manager, strain_rate);

}

double sum_all_shear_strains(vector<double>& shear_strains) {
	
	double taylor_factor = 0;
	
	for(size_t i=0;i<shear_strains.size();i++) {
		
		taylor_factor += abs(shear_strains.at(i));
	
	}

	return taylor_factor;

}

double calc_shear_strains_sum(const SlipSystemsManager& manager, const vector<double>& shear_strains, size_t m, size_t n) {

	double result = 0;
	for(size_t i=0;i<shear_strains.size();i++) {
		SlipSystem& system = manager.get_slipsystem(i);
		
		result += system.get_sym().get(m,n)*shear_strains.at(i);

	}

	return result;

}

int strain_minimization(const gsl_vector* x, void* params, gsl_vector* f) {

	Matrix& strain = *(((struct smparams*)params)->strain);
	Grain& grain = *(((struct smparams*)params)->grain);
	double srate = ((struct smparams*)params)->strain_rate;
	SlipSystemsManager& manager = *(((struct smparams*)params)->manager);

 	vector<double> shear_strains = calculate_shear_strains_helper(x,grain,manager,srate);

	double result11 = calc_shear_strains_sum(manager,shear_strains,0,0) - strain.m11();
	double result12 = calc_shear_strains_sum(manager,shear_strains,0,1) - strain.m12();
	double result13 = calc_shear_strains_sum(manager,shear_strains,0,2) - strain.m13();
	double result22 = calc_shear_strains_sum(manager,shear_strains,1,1) - strain.m22();
	double result23 = calc_shear_strains_sum(manager,shear_strains,1,2) - strain.m23();

	gsl_vector_set(f,0,result11);
	gsl_vector_set(f,1,result12);

	gsl_vector_set(f,2,result13);
	gsl_vector_set(f,3,result12);

	gsl_vector_set(f,4,result22);
	gsl_vector_set(f,5,result23);
	
	gsl_vector_set(f,6,result13);
	gsl_vector_set(f,7,result23);

	gsl_vector_set(f,8,result11-result22);

	return GSL_SUCCESS;

}

int simulations(Grain& grain, const Matrix& def, SlipSystemsManager& manager, const double strain_rate, size_t maxiter, double residual) {

        const Matrix& a = grain.get_rotation_matrix();
        Matrix at = a.give_transpose();

        Matrix def_l = a*def*at;

        int status = 0;
        size_t iter = 0;

        struct smparams sp{&def_l,&grain,strain_rate,&manager};

        gsl_multiroot_function f = {&strain_minimization, 9, &sp};

        gsl_multiroot_fsolver* s = grain.get_multiroot_fsolver();
        gsl_multiroot_fsolver_set (s, &f, grain.get_stress());

        gsl_multiroot_fsolver_set (s, &f, grain.get_stress());

        gsl_vector* stress = gsl_vector_alloc(9);

        while(1) {
                do {
                        iter++;
                        status = gsl_multiroot_fsolver_iterate (s);

                        if (status)   /* check if solver is stuck */
                                break;

                        status = gsl_multiroot_test_residual (s->f, residual);

                }while (status == GSL_CONTINUE && iter < maxiter);

                if(status==GSL_SUCCESS) {
                        break;
                }

                for(size_t i=0;i<9;i++) {

                        gsl_vector_set(stress, i, rand()%100-50);
                }

                gsl_multiroot_fsolver_set (s, &f, stress);
                iter = 0;

        }

        grain.set_stress(s->x);

        return status;

}
