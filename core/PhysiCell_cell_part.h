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

#ifndef __PhysiCell_cell_part_h__
#define __PhysiCell_cell_part_h__

#include "./PhysiCell_custom.h" 

#include "../BioFVM/BioFVM.h"
#include "./PhysiCell_phenotype.h"
#include "./PhysiCell_cell_container.h"
#include "./PhysiCell_constants.h"
#include "./PhysiCell_cell.h"

#include "../modules/PhysiCell_settings.h" 

#include "./PhysiCell_standard_models.h" 

#include <cstdint>

class Cell_Part_Definition
{
 private:
 public: 
    //Not sure what the fields should be, Paul should update with these
    
	/*int type; 
	std::string name; 
 
	Microenvironment* pMicroenvironment; 
	
	Cell_Parameters parameters; 
	Custom_Cell_Data custom_data; 
	Cell_Functions functions; 
	Phenotype phenotype; 

	Cell_Part_Definition();  // done 
	Cell_Part_Definition( Cell_Part_Definition& cd ); // copy constructor 
	Cell_Part_Definition& operator=( const Cell_Part_Definition& cd ); // copy assignment */
};

//Cell_Part is a subcellular agent class
class Cell_Part //: public Basic_Agent 
{
    static int Cell_Part_ID_Counter = 0; //keeps track of the cell id (increments each time a new cell is created)

    private:
        

    public:
        //1.unique ID
        int unique_id;

        //2. voxel ID
        int voxel_id; //(where do you live?)

        //3. cell ID
        int cell_id; //(who do you belong to?)

        //4. Cell*
        Cell_* containing_cell; //pointer to containing cell

        //5. Type (integer)
        int type;

        //6. Pointer to an (undefined) cell part definition?
        Cell_Part_Definition* cell_part_definition;

        //7. Position
        double position[3]; // [x,y,z] (double array for easier transfer to/from gpu)

        //8. Velocity, double
        double velocity[3]; //[x,y,z]

        //9. Previous Velcocity (for Adams-Bashforth)
        double previous_velocity[3]; // [x,y,z]

        //10. Volume: need to ask if this is the Volume class or double?
        double total_volume; //based off of basic_agent, total_volume == double volume

        //11.Secretion Rates (vector<double>)
        std::vector<double> secretion_rates;

        //12.Secretion Targets type? need to ask what these refer to
        double secretion_targets;

        //13. Uptake rates
        std::vector<double> uptake_rates;

        //14.Net export rates
        std::vector<double> net_export_rates;

        //15. Internalized Substrates
        std::vector<double> internalized_substrates;

        //16. Equivalent Radius
        double equivalent_radius; //this is 'double radius' within the current Cell from Phenotype.Geometry?

        //17.Repulsion Strength
        double repulsion_strength;

        //18.Cell-cell Adhesion Strength, vector?
        double adhesion_strength;

        //19.Max adhesion distance
        double max_adhesion_distance;

         //20.Vector<Cell_Part*> (vector of mechanically interacting cell parts)
         std::vector<Cell_Part*>;

         //21.Vector<Cell*> (vector of mechanically interacting cells)
         std::vector<Cell*>;

         //22...


        //Methods:

        //Constructs need to place agent in data structures
        Cell_Part();
        Cell_Part(Cell_* cell,int voxel_id);

        update_velocity(); //need mechanics DS, use this to tell cell parts what to interact with

        update_position(); //via velocities (Adams-Bashforth)

        secretion_uptake_export();

        assign_position();

        ~Cell_Part(); //need removal/deconstructor to delist from data structures and connected cells


        //Above are the overleaf instructions-------------------------

		void assign_position(double x, double y, double z);

        void set_velocity(double x, double y, double z);

        void set_previous_velocity(double x, double y, double z);

}

#endif