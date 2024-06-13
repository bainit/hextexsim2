#include <iostream>
#include <cstdlib>

#ifdef MPI_ENABLED
#include <mpi.h>
#endif

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "./includes/vpsc_classes.h"
#include "./includes/vpsc_core.h"
#include "./includes/vpsc_util.h"

using namespace std;

int main(int argc, char** argv) {

	cout<<"Test for simulations"<<endl;

	#ifdef MPI_ENABLED
    	MPI_Init(NULL, NULL);

	int world_size, world_rank, name_len;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

   	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    	MPI_Get_processor_name(processor_name, &name_len);
	#endif

	time_t t;
	srand((unsigned) time(&t));

	Settings::create_instance("settings_rolling");
	const Settings& settings = Settings::get_instance();	

	SlipSystemsManager manager;
	load_slip_systems("hcp_mg", manager);

	vector<Grain*> grains = load_grains("input_texture_1248g.csv", manager.get_size());

	#ifdef MPI_ENABLED
	double p = grains.size()/world_size;
	size_t parts = (size_t)(p);

	size_t start_grains = world_rank*parts;
	size_t end_grains = (world_rank+1)*parts - 1;
	size_t local_grains_size = end_grains-start_grains;
	printf("Processor %s, rank %d, I will calculate grains from %d to %d\n, total: %d", processor_name, world_rank, start_grains,end_grains, local_grains_size);
	
	for(size_t ng=0;ng<grains.size();ng++) {
		Grain& grain = *(grains.at(ng));
		for(size_t i=0;i<manager.get_size();i++) {
		
			const SlipSystem& system = manager.get_slipsystem(i);
			grain.set_crss_at(i,system.get_crss());	
				
		}
	}
	#endif

	
	cout<<"Number of loaded grains"<<endl;
	cout<<grains.size()<<endl;
	cout<<endl;

	Results::create_instance();
	Results& results = Results::get_instance();
	vector<double> total_taylor_factor_to_remove(5);

	for(size_t ng=0;ng<grains.size();ng++) {
		Grain& grain = *(grains.at(ng));
		for(size_t i=0;i<manager.get_size();i++) {
		
			const SlipSystem& system = manager.get_slipsystem(i);
			grain.set_crss_at(i,system.get_crss());	
				
		}
	}

	Matrix def = settings.DEF_STEP * settings.DEFORMATION;
	double evm = settings.VONMISES_STRAIN * settings.DEF_STEP;

	int status;

	for(size_t step=0;step<settings.MAX_DEF;step++) {

		double local_taylor_factor = 0;

		#pragma omp parallel for
		#ifdef MPI_ENABLED
		for(size_t ng=start_grains;ng<end_grains;ng++) {
		#else
		for(size_t ng=0;ng<grains.size();ng++) {
						
			cout<<"def: "<<step<<", grain: "<<ng<<endl;
		#endif
				
			Grain& grain = *(grains.at(ng));
						
			for(size_t i=0;i<9;i++) {
				gsl_vector_set(grain.get_stress(), i, rand()%100-50);
			}
						
			status = simulations(grain, def, manager, settings.STRAIN_RATE, settings.MAXITER, settings.TOLERANCE);

			vector<double> shear_strains = calculate_shear_strains_helper(grain.get_stress(),grain,manager, settings.STRAIN_RATE);

			Matrix rm = calc_rm_after_def(shear_strains, manager);

			grain.update_rotation(rm);

			local_taylor_factor += (sum_all_shear_strains(shear_strains)/evm);			
		  	
		}

		#ifdef MPI_ENABLED
		results.add_to_taylor_factor(local_taylor_factor/local_grains_size);	
		#else
		results.add_to_taylor_factor(local_taylor_factor/grains.size());
		#endif
	
		#ifdef MPI_ENABLED
		MPI_Barrier(MPI_COMM_WORLD);

		if(world_rank==0) {
			cout<<"Finished def: "<<step<<endl;
			if(step%10==0) {
				string filename = "outputtexture"+to_string(step)+".csv";
				export_grains(grains,filename);
			}
		}
		#else
		cout<<"Finished def: "<<step<<endl;
		if(step%10==0) {
			string filename = "outputtexture"+to_string(step)+".csv";
			export_grains(grains,filename);
		}
		#endif
					
	}

	#ifdef MPI_ENABLED
	MPI_Reduce(results.get_taylor_factor().data(),total_taylor_factor_to_remove.data(), 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	#endif

	#ifdef MPI_ENABLED
	if(world_rank==0) {	
		export_grains(grains,"outtexture.csv");	
		for(size_t i=0;i<total_taylor_factor_to_remove.size();i++) {
			cout<<total_taylor_factor_to_remove.at(i)/world_size<<endl;
		}	
	}

	MPI_Finalize();
	#else
	export_grains(grains,"outtexture.csv");
	vector<double> tf = results.get_taylor_factor();
		for(size_t i=0;i<tf.size();i++) {
			cout<<tf.at(i)<<endl;
		}
	#endif

	return 0;

}
