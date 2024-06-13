#ifndef VPSC_CLASSES_H
#define VPSC_CLASSES_H

#include <string>
#include <memory>
#include <vector>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

using namespace std;

#define FCC_SLIP_SYSTEMS "FCC"
#define HCP_SLIP_SYSTEMS "HCP" 

//gsl_matrix wraper 3x3
class Vector;

class Matrix {

	public:
		Matrix();
		Matrix(initializer_list<double> d);
		Matrix(const Matrix& matrix);
		Matrix(gsl_matrix* m);
		Matrix(Matrix&& other);

		Matrix operator*(const Matrix& m) const;
		Vector operator*(const Vector& v) const;
		Matrix operator+(const Matrix& mr) const;
		Matrix operator-(const Matrix& mr) const;

		Matrix& operator=(Matrix&& other);

		double get(size_t i, size_t j) const;
		
		gsl_matrix* get_matrix();
		gsl_matrix* get_matrix() const;

		double m11() const;
		double m12() const;
		double m13() const;
		double m21() const;
		double m22() const;
		double m23() const;
		double m31() const;
		double m32() const;
		double m33() const;

		Matrix give_inverse() const;
		Matrix give_transpose() const;

		void print() const;

		~Matrix();

		friend class Vector;

	private:
		gsl_matrix* matrix;

		
};

Matrix operator*(double d, const Matrix& m);

//vector wraper over gsl_matrix 1x3
class Vector {

	public:

		enum Orientation {horizontal=0,vertical};
		
		Vector(Orientation o=Orientation::horizontal);
		Vector(initializer_list<double> d, Orientation o=Orientation::horizontal);
		Vector(gsl_matrix* v, Orientation o=Orientation::horizontal);
		Vector(Vector&& other);
		Vector(const Vector& v);

		Vector operator*(const Matrix& m) const;
		Matrix operator*(const Vector& rv) const;
		Vector& operator=(Vector&& other);
	
		double v1() const;
		double v2() const;
		double v3() const;
		
		double get_length();
		void normalize();

		Vector give_transpose() const;
		Vector give_unit_vector(Orientation o=Orientation::horizontal) const;

		void print() const;
		bool ishorizontal() const;
		bool isvertical() const;

		~Vector();

		friend class Matrix;

	private:

		gsl_matrix* vector;
		Orientation orientation;

		double length;

		gsl_matrix* get_vector() const;
		gsl_matrix* create_vector();

};

class Settings {

	public:
		static void create_instance(const string& pathname);
		static const Settings& get_instance();
 
		//identity matrix
		Matrix IDENTITY_MATRIX;

		double VONMISES_STRAIN;
		
		double STRAIN_RATE;
		double DEF_STEP;
		size_t MAX_DEF;
		double TOLERANCE;
		size_t MAXITER;		
		Matrix DEFORMATION;

		~Settings();

	private:
		Settings(double sr, double ds, size_t md, double t, size_t mi, double evm, Matrix& d);
		
		static Settings* instance;	
	
};

class Results {

	public:
		static void create_instance();
		static Results& get_instance();

		void add_to_taylor_factor(double val);
		vector<double> const & get_taylor_factor() const;

	private:
		Results();
		static Results* instance;

		vector<double> taylor_factor;
	
};

class Grain {

	public:
		Grain(double p1, double p, double p2, double weight, size_t numssystems);
		
		const Matrix& get_rotation_matrix() const;
		void print() const;
		double get_phi1() const {return phi1;};
		double get_phi() const {return phi;};
		double get_phi2() const {return phi2;};
		double get_weight() const {return weight;};

		size_t get_crsses_size() const {return crsses.size();};
		double get_crss_at(size_t p) const;
		void set_crss_at(size_t p, double val);
		void update_rotation(const Matrix& rot);
		void set_stress(gsl_vector* vs);
		gsl_vector* get_stress() {return stress;}
		gsl_multiroot_fsolver* get_multiroot_fsolver() {return multiroot_fsolver;}

		~Grain();

	private:
		//Euler angles in radians
		double phi1, phi, phi2, weight;
		vector<double> crsses;
		Matrix rotation_matrix;
		const gsl_multiroot_fsolver_type *T;
  		gsl_multiroot_fsolver *multiroot_fsolver;
		gsl_vector* stress;

	
};

class SlipSystem {

	public:
		//{hkl}<uvw>, {slip plane}<slip direction>
		//tau_0, tau_1, theta_0, theta_1 and latent_hardenign are hardenign parameters for voce hardenig law
		SlipSystem(initializer_list<double> plane, initializer_list<double> direction, double crss, double mm, string label, double tau_0, double tah_1, double theta_0, double theta_1, initializer_list<double> latent_hardening);
		
		SlipSystem(Vector& plane, Vector& dir, double crss, double mm, string label, double tau_0, double tau_1, double theta_0, double theta_1, initializer_list<double> latent_hardening);	
	
		string get_label() const {return label;};

		const Vector& get_slip_plane() const {return slip_plane;};
		const Vector& get_slip_dir() const {return slip_direction;};

		//Schmid matrix
		const Matrix& get_def_mat() const {return deformation_matrix;};
		const Matrix& get_sym() const {return symetry_def_mat;};
		const Matrix& get_asym() const {return asymetry_def_mat;};

		double get_crss() const {return crss;};
		double get_m() const {return m;};

		double get_tau_0() const {return tau_0;};
		double get_tau_1() const {return tau_1;};
		double get_theta_0() const {return theta_0;};
		double get_theta_1() const {return theta_1;};

		double get_lh_parameter(size_t p) const;
		size_t get_lh_size() const;

		~SlipSystem();

	private:
		//normalized 
		//spip_plane (1,3)
		Vector slip_plane;
		//slip_dir (3,1)
		Vector slip_direction;

		string label;

		//Schmid matrix
		Matrix deformation_matrix;

		//Symetry part of deformation matrix
		Matrix symetry_def_mat;

		//Antysymetry part of deformation matrix
		Matrix asymetry_def_mat;

		//critical resolved shear stress
		const double crss;

		//viscoplastic parameter
		double m;

		const double tau_0;
		const double tau_1;
		const double theta_0;
		const double theta_1;
		const vector<double> latent_hardening;

		void calculate();
 
};

class SlipSystemsManager {

	public:
		SlipSystemsManager();
		
		void add_slipsystem(SlipSystem* system);
		size_t get_size() const;
		SlipSystem& get_slipsystem(size_t p) const;		

		~SlipSystemsManager();

	private: 
		vector<SlipSystem*> slipsystems;	
	
}; 

struct smparams {

	smparams(Matrix* e, Grain* g, double sr, SlipSystemsManager* m) : strain{e}, grain{g}, strain_rate{sr}, manager{m} {};
	Matrix* strain;
	Grain* grain;
	double strain_rate;
	SlipSystemsManager* manager;

};
#endif //VPSC_CLASSES_H
