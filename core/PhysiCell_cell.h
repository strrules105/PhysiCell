/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#ifndef __PhysiCell_cell_h__
#define __PhysiCell_cell_h__

#include "./PhysiCell_custom.h" 

#include "../BioFVM/BioFVM.h"
#include "./PhysiCell_phenotype.h"
#include "./PhysiCell_cell_container.h"
#include "./PhysiCell_constants.h"

#include "../modules/PhysiCell_settings.h" 

#include "./PhysiCell_standard_models.h" 

#include <cstdint>

using namespace BioFVM; 

namespace PhysiCell{
class Cell_Container;

class Cell_Parameters
{
 private:
 public:
	// oxygen values (in mmHg) for critical phenotype changes
	double o2_hypoxic_threshold; // value at which hypoxic signaling starts
	double o2_hypoxic_response; // value at which omics changes are observed 
	double o2_hypoxic_saturation; // value at which hypoxic signalign saturates 
	// o2_hypoxic_saturation < o2_hypoxic_threshold
	
	double o2_proliferation_saturation; // value at which extra o2 does not increase proliferation
	double o2_proliferation_threshold; // value at which o2 is sufficient for proliferation

	double o2_reference; // physioxic reference value, in the linked reference Phenotype
	// o2_proliferation_threshold < o2_reference < o2_proliferation_saturation; 
	
	double o2_necrosis_threshold; // value at which cells start experiencing necrotic death 
	double o2_necrosis_max; // value at which necrosis reaches its maximum rate 
	// o2_necrosis_max < o2_necrosis_threshold
	
	Phenotype* pReference_live_phenotype; // reference live phenotype (typically physioxic) 
	// Phenotype* pReference_necrotic_phenotype; // reference live phenotype (typically physioxic) 

	// necrosis parameters (may evenually be moved into a reference necrotic phenotype 
	double max_necrosis_rate; // deprecate
	int necrosis_type; // deprecate 
	
	Cell_Parameters(); 
}; 

class Cell_Definition
{
 private:
 public: 
	int type; 
	std::string name; 
 
	Microenvironment* pMicroenvironment; 
	
	Cell_Parameters parameters; 
	Custom_Cell_Data custom_data; 
	Cell_Functions functions; 
	Phenotype phenotype; 

	Cell_Definition();  // done 
	Cell_Definition( Cell_Definition& cd ); // copy constructor 
	Cell_Definition& operator=( const Cell_Definition& cd ); // copy assignment 
};

extern Cell_Definition cell_defaults; 

class Cell_State
{
 public:
	std::vector<Cell*> neighbors; // not currently tracked! 
	std::vector<double> orientation;
	double simple_pressure; 
	
	Cell_State(); 
};

class Cell_GPU_UpdateAll_Secretion_Advance{
	private:

	public:
		//Variables used to mimic 'void Secretion::sync_to_microenvironment(Microenvironment* pNew_Microenvironment)' 
		/*Basic_Agent -> Cell -> 'microenvironment' fields*/
		//int cell_phenotype_secretion_pMicroenv_num_of_densities_GPU; //Cell->phenotype->secretion->pMicroenvironment->(*p_density_vectors)[0].size()
		//BioFVM::Microenvironment * cell_microenv_address; //should point to Cell->microenvironment
		//BioFVM::Microenvironment * cell_phenotype_secretion_pMicroenv_address; //should point to Cell->phenotype->secretion->pMicroenvironment
		
		//unsigned long int *default_microenv_address; //should point to BioFVM::default_microenvironment

		//Basic_Agent -> Cell -> 'phenotype' -> 'secretion' -> rates fields
		// uintptr_t cell_phenotype_secretion_secretion_rates_GPU;
		// uintptr_t cell_phenotype_secretion_update_rates_GPU;
		// uintptr_t cell_phenotype_secretion_saturation_densities_GPU;
		// uintptr_t cell_phenotype_secretion_net_export_rates_GPU;

		// //Basic_Agent -> Cell -> rates fields
		// uintptr_t cell_secretion_rates_GPU;
		// uintptr_t cell_update_rates_GPU;
		// uintptr_t cell_saturation_densities_GPU;
		// uintptr_t cell_net_export_rates_GPU;

		//Basic_Agent -> Cell -> 'phenotype' -> 'secretion' -> rates fields
		// double * cell_phenotype_secretion_secretion_rates_GPU;
		// double * cell_phenotype_secretion_update_rates_GPU;
		// double * cell_phenotype_secretion_saturation_densities_GPU;
		// double * cell_phenotype_secretion_net_export_rates_GPU;

		//Basic_Agent -> Cell -> rates fields
		double * cell_secretion_rates_GPU;
		int cell_secretion_rates_GPU_size;
		double * cell_update_rates_GPU;
		int cell_update_rates_GPU_size;
		double * cell_saturation_densities_GPU;
		int cell_saturation_densities_GPU_size;
		double * cell_net_export_rates_GPU;
		int cell_net_export_rates_GPU_size;


		/*---------------*/

		/*Variables used in pCell->set_total_volume( phenotype.volume.total ); */
		double * cell_phenotype_volume_total; //Cell->phenotype.volume.total
		double * volume_GPU;


		/*-----------------------*/

		

		/*GPU Version of void Basic_Agent::simulate_secretion_and_uptake( Microenvironment* pS, double dt )*/
		bool * is_active_GPU;
		bool * volume_is_changed_GPU;

		double * total_extracellular_substrate_change_GPU; //Pointer to array of doubles (total_extracellular_substrate_change)
		int total_extracellular_substrate_change_GPU_size;

		double * pS_current_voxel_index_arr_GPU; //Pointer to array of doubles ((*pS)(current_voxel_index))
		int pS_current_voxel_index_arr_GPU_size;

		double * internalized_substrates_GPU; //Pointer to array of doubles (internalized_substrates)
		int internalized_substrates_GPU_size;

		double * cell_source_sink_solver_temp1_GPU;
		int cell_source_sink_solver_temp1_GPU_size;

		double * cell_source_sink_solver_temp2_GPU;
		int cell_source_sink_solver_temp2_GPU_size;

		double * cell_source_sink_solver_temp_export1_GPU; 
		int cell_source_sink_solver_temp_export1_GPU_size;

		double * cell_source_sink_solver_temp_export2_GPU; 	
		int cell_source_sink_solver_temp_export2_GPU_size;

		//Below is part of Basic_Agent::set_internal_update_constants(dt)
		/*Cell->phenotype->secretion->pMicroenvironment->mesh->voxels   (std::vector<Voxel> voxels) -> volume of each Voxel object*/
		//double ** cell_phenotype_secretion_pMicroenv_mesh_voxels_volumes_GPU;
		//int cell_phenotype_secretion_pMicroenv_mesh_voxels_GPU_size; //the size of std::vector<Voxel> voxels
		//Better: we can probably just store the Voxel[voxel_index].volume of each Cell before going to GPU
		double * cell_phenotype_secretion_pMicroenv_mesh_current_voxel_volume;
		

		Cell_GPU_UpdateAll_Secretion_Advance();
		Cell_GPU_UpdateAll_Secretion_Advance(Cell *cell); 
		~Cell_GPU_UpdateAll_Secretion_Advance();
		void copy_Cell_GPU_to_device();


		void update_device();
		void update_host();

		#pragma acc routine
		void axpy_GPU( double* y, int y_size, double a , double* x, int x_size );

		#pragma acc routine
		void secretion_advance_GPU(double dt, bool default_microenvironment_options_track_internalized_substrates_in_each_agent_GPU);

		#pragma acc routine
		void simulate_secretion_and_uptake_GPU(double dt, bool default_microenvironment_options_track_internalized_substrates_in_each_agent_GPU);

		#pragma acc routine
		void set_internal_uptake_constants_GPU( double dt );

		#pragma acc routine
		void copyValsArr(double *arr1, int arr1_size, double *arr2, int arr2_size);

		#pragma acc routine
		void setArrConst(double *arr1, int arr1_size, double c);

		#pragma acc routine
		void addArrs(double *arr1, int arr1_size, double *arr2, int arr2_size);

		#pragma acc routine
		void multiArrs(double *arr1, int arr1_size, double *arr2, int arr2_size);

		#pragma acc routine
		void subArrs(double *arr1, int arr1_size, double *arr2, int arr2_size);

		#pragma acc routine
		void divArrs(double *arr1, int arr1_size, double *arr2, int arr2_size);

		#pragma acc routine
		void arrAddConst(double *arr, int arr_size, double c);

		#pragma acc routine
		void arrMultiConst(double *arr, int arr_size, double c);

		#pragma acc routine
		void arrDivConst(double *arr, int arr_size, double c);

			

	/*Takes in a vector of Cell pointers and outputs a pointer to an array of 'Cell_GPU_UpdateAll_Secretion_Advance' objects*/
	static Cell_GPU_UpdateAll_Secretion_Advance* create_GPU_Cells_Arr(std::vector<Cell*> *all_cells_ , int all_cells_size_, double dt_);


	
};

class Cell : public Basic_Agent 
{
 private: 
	Cell_Container * container;
	int current_mechanics_voxel_index;
	int updated_current_mechanics_voxel_index; // keeps the updated voxel index for later adjusting of current voxel index
		
 public:

	std::string type_name; 
	bool is_out_of_domain;
	bool is_movable;

	Cell_State state;  //contains 1 pointer member
	Custom_Cell_Data custom_data; //no pointer members
	Cell_Parameters parameters; //1 pointer member
	Cell_Functions functions;  
	Phenotype phenotype; 


	//Methods:
	Cell();
	~Cell();

	/*GPU Methods-------*/
	void testCell();
	void update_device();
	void update_host();



	/*------------------*/
	
	void update_motility_vector( double dt_ );
	void advance_bundled_phenotype_functions( double dt_ ); 
	
	void add_potentials(Cell*);       // Add repulsive and adhesive forces.
	void set_previous_velocity(double xV, double yV, double zV);
	int get_current_mechanics_voxel_index();
	void turn_off_reactions(double); 		  // Turn off all the reactions of the cell
	
	void flag_for_division( void ); // done 
	void flag_for_removal( void ); // done 
	
	void start_death( int death_model_index ); 
	void lyse_cell( void ); 

	Cell* divide( void );
	void die( void );
	void step(double dt);
	bool assign_position(std::vector<double> new_position);
	bool assign_position(double, double, double);
	void set_total_volume(double);
	
	double& get_total_volume(void); // NEW
	
	void set_target_volume(double); 
	void set_target_radius(double); 
	void set_radius(double); 


	
	
	
	
	// mechanics 
	void update_position( double dt ); //
	std::vector<double> displacement; // this should be moved to state, or made private  

	
	void assign_orientation();  // if set_orientaion is defined, uses it to assign the orientation
								// otherwise, it assigns a random orientation to the cell.
	
	void copy_function_pointers(Cell*);
	
	void update_voxel_in_container(void);
	void copy_data(Cell *);
	
	void ingest_cell( Cell* pCell_to_eat ); // for use in predation, e.g., immune cells 

	// I want to eventually deprecate this, by ensuring that 
	// critical BioFVM and PhysiCell data elements are synced when they are needed 
	
	void set_phenotype( Phenotype& phenotype ); // no longer needed?
	void update_radius();
	Cell_Container * get_container();
	
	std::vector<Cell*>& cells_in_my_container( void ); 
	
	void convert_to_cell_definition( Cell_Definition& cd ); 
};

Cell* create_cell( void );  
Cell* create_cell( Cell_Definition& cd );  


void delete_cell( int ); 
void delete_cell( Cell* ); 
void save_all_cells_to_matlab( std::string filename ); 

//function to check if a neighbor voxel contains any cell that can interact with me
bool is_neighbor_voxel(Cell* pCell, std::vector<double> myVoxelCenter, std::vector<double> otherVoxelCenter, int otherVoxelIndex);  


extern std::unordered_map<std::string,Cell_Definition*> cell_definitions_by_name; 
extern std::unordered_map<int,Cell_Definition*> cell_definitions_by_type; 
extern std::vector<Cell_Definition*> cell_definitions_by_index; // works 

void display_cell_definitions( std::ostream& os ); // done 
void build_cell_definitions_maps( void ); // done 

Cell_Definition* find_cell_definition( std::string search_string ); // done 
Cell_Definition* find_cell_definition( int search_type );  

Cell_Definition& get_cell_definition( std::string search_string ); // done 
Cell_Definition& get_cell_definition( int search_type );  

Cell_Definition* initialize_cell_definition_from_pugixml( pugi::xml_node cd_node ); 
void initialize_cell_definitions_from_pugixml( pugi::xml_node root ); 
void initialize_cell_definitions_from_pugixml( void );

extern std::vector<double> (*cell_division_orientation)(void);


};

#endif