#include <stdexcept>

#include "../includes/vpsc_classes.h"
#include "../includes/vpsc_util.h"
#include "../includes/vpsc_core.h"

using namespace std;

SlipSystemsManager::SlipSystemsManager() {

}

void SlipSystemsManager::add_slipsystem(SlipSystem* system) {
	
	slipsystems.push_back(system);

}

SlipSystem& SlipSystemsManager::get_slipsystem(size_t p) const {

	if(p>=slipsystems.size()) {
	
		throw runtime_error("SlipSystemManager::get_slip_system(), position out of size");

	}

	return *(slipsystems.at(p));

}

size_t SlipSystemsManager::get_size() const {

	return slipsystems.size();

}


SlipSystemsManager::~SlipSystemsManager() {

	for(size_t i=0;i<slipsystems.size();i++) {
		SlipSystem* system = slipsystems.at(i);
		if(system!=NULL) {
			delete(slipsystems.at(i));
		}
	}

}
