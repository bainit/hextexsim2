#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <utility>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "../includes/vpsc_classes.h"
#include "../includes/vpsc_util.h"
#include "../includes/vpsc_core.h"

using namespace std;

void print_matrix(gsl_matrix const * matrix) {

	for(int i=0;i<matrix->size1;i++) {
		for(int j=0;j<matrix->size2;j++) {
			printf("%.2f ",gsl_matrix_get(matrix,i,j));
		}
		printf("\n");
	}

}

double rad2deg(double rad) {

	return rad*180/M_PI;

}

double deg2rad(double deg) {

	return deg*M_PI/180;	

}

double sign(double d) {

	if(signbit(d)==1) {
		return 1;
	}

	return -1;

}

void matmul(gsl_matrix const *lm, gsl_matrix const *rm, gsl_matrix *res) {
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, lm,rm, 0.0, res);
}

vector<Grain*> load_grains(string pathname, const size_t num_of_slipsystems) {

	printf("Loading grains from: %s\n", pathname.c_str());

   	ifstream input_file;

        input_file.open(pathname);

        if(!input_file.is_open()) {
		printf("Error during opening file: %s\n", pathname.c_str());
  		exit(-1);
	}
	
	vector<Grain*> grains;
	string line;

	while(getline(input_file,line)) {
		string val;
		istringstream stream(line);
		double vals[4];
		int i=0;
		while(getline(stream,val,';') && i<4) {
			vals[i] = stod(val);
			i++;
		}

		grains.push_back(new Grain(vals[0],vals[1],vals[2],vals[3],num_of_slipsystems));
	}

	return grains;

}

void export_grains(const vector<Grain*> grains, const string pathname) {

	printf("Exporting grains to file: %s\n", pathname.c_str());

	stringstream stream(pathname);

	for(int i=0;i<grains.size();i++) {

		Grain& g = *(grains[i]);
		
		const Matrix& rm = g.get_rotation_matrix();

		double phi1=0;
		double PHI=0;
		double phi2=0;

		if(abs(rm.m33())<0.999) {
			PHI = acos(rm.m33());
			double sinphi = sin(PHI);
			phi1 = atan2(rm.m31()/sinphi,-rm.m32()/sinphi);
			phi2 = atan2(rm.m13()/sinphi,rm.m23()/sinphi);
		}else {
			phi1 = atan2(rm.m12(),rm.m11())/2;
			phi2 = phi1;

		}

		stream<<to_string(phi1)<<";"<<to_string(PHI)<<";"<<to_string(phi2)<<";"<<to_string(1)<<endl;

	}

	ofstream outfile(pathname,ofstream::out);
	outfile<<stream.rdbuf();
	outfile.close();

}

void load_slip_systems(const string& pathname, SlipSystemsManager & manager) {

	printf("Loading slip systems from: %s\n", pathname.c_str());

   	ifstream input_file;

        input_file.open(pathname);

        if(!input_file.is_open()) {
		throw runtime_error("load_slip_systems(): cant open slip systems file");
	}

	string line;
	getline(input_file,line);

	if(line.compare(FCC_SLIP_SYSTEMS)==0) {
		printf("Loading fcc slip system family\n");

		while(getline(input_file,line)) {
     			double u,v,w,h,k,l,tau,m;
			string name;

			istringstream input_data(line);
			input_data>>u;
			input_data>>v;
			input_data>>w;
	
			input_data>>h;
			input_data>>k;
			input_data>>l;

			input_data>>tau;
			input_data>>m;

			input_data>>name;

			manager.add_slipsystem(new SlipSystem{{u,v,w},{h,k,l},tau,m,name,1,1,1,1,{1,1,1,1,1,1,1,1,1,1,1,1}});
		
		}

	}else if(line.compare(HCP_SLIP_SYSTEMS)==0) {
		
		printf("Loading hcp slip system family\n");

		double ctoa = 0;
		getline(input_file,line);
		istringstream input_data(line);
		input_data>>ctoa;
		
		Matrix transmat = get_transformation_matrix(ctoa);
		Matrix invtransmat = transmat.give_inverse();

		while(getline(input_file,line)) {

     			double u,v,t,w,h,k,i,l,tau,m;
			string name;

			istringstream input_data(line);
			input_data>>u;
			input_data>>v;
			input_data>>t;
			input_data>>w;
	
			input_data>>h;
			input_data>>k;
			input_data>>i;
			input_data>>l;

			input_data>>tau;
			input_data>>m;

			input_data>>name;

			Vector plane = hkil2uvw_plane(u,v,t,w,invtransmat);
			Vector dir = hkil2uvw_dir(h,k,i,l,transmat);

			manager.add_slipsystem(new SlipSystem{plane,dir,tau,m,name,1,1,1,1,{1,1,1,1,1,1,1,1,1,1,1,1}});

		}
		
	
	}else {
		input_file.close();
		throw runtime_error("load_slip_systems(): slip systems familly not known");

	}

	input_file.close();

}
