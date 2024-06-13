#include <stdexcept>

#include "../includes/vpsc_classes.h"

using namespace std;

Results* Results::instance = NULL;

Results::Results() {}

void Results::create_instance() {

	if(instance == NULL) {
		instance = new Results();
	}

}

Results& Results::get_instance() {
	
	if(instance==NULL) {
		throw runtime_error("Results::get_instance(), Results instance is not created");
	}

	return *instance;
}

void Results::add_to_taylor_factor(double val) {

	taylor_factor.push_back(val);

}

vector<double> const &  Results::get_taylor_factor() const{

	return taylor_factor;

}
