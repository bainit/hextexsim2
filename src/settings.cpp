#include <stdexcept>
#include <cstdio>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "../includes/vpsc_classes.h"
#include "../includes/vpsc_core.h"

using namespace std;

Settings* Settings::instance = NULL;

Settings::Settings(double sr, double ds, size_t md, double t, size_t mi, double evm, Matrix& d) : STRAIN_RATE{sr}, DEF_STEP{ds}, MAX_DEF{md}, TOLERANCE{t}, MAXITER{mi}, DEFORMATION{d}, VONMISES_STRAIN{evm}, IDENTITY_MATRIX{1,0,0,0,1,0,0,0,1} {}

void Settings::create_instance(const string& pathname) {

	if(instance == NULL) {

		printf("Loading settings from: %s\n", pathname.c_str());

		ifstream input_file;

		input_file.open(pathname);

		if(!input_file.is_open()) {
			throw runtime_error("Settings::create_instance(): cant open settings file");
		}

		double strain_rate, def_step, tolerance;
		size_t max_def, max_iter;

		string line;
		if(getline(input_file,line)) {
			istringstream input_data(line);
			input_data>>strain_rate;
			input_data>>def_step;
			input_data>>max_def;
			input_data>>tolerance;
			input_data>>max_iter;		

		}else {
			input_file.close();
			throw runtime_error("Settings::create_instance(): cant read settings parameters");
		}

		double u[9];
	
		for(size_t counter = 0;counter<3;counter++) {
			if(getline(input_file,line)) {
				istringstream input_data(line);
				input_data>>u[counter*3+0];
				input_data>>u[counter*3+1];
				input_data>>u[counter*3+2];
			}else {
				input_file.close();
				throw runtime_error("Settings::create_instance(): cant read deformation matrix");
			}

		}

		input_file.close();

		Matrix def{u[0],u[1],u[2], u[3],u[4],u[5], u[6],u[7],u[8]};

	   	double evm = evmises_strain(def);
	       	instance = new Settings{strain_rate, def_step, max_def, tolerance, max_iter, evm, def};
      
	}

}

const Settings& Settings::get_instance(){

	if(instance==NULL) {
		throw runtime_error("const Settings::get_instance(), Settings instance is not created");
	}

	return *instance;

}

Settings::~Settings() {
	
	if(instance!=NULL) {
		delete instance;
	}

}
