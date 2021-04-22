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

#include "../BioFVM/BioFVM_agent_container.h"
#include "PhysiCell_constants.h"
#include "../BioFVM/BioFVM_vector.h"
#include "PhysiCell_cell.h"
#include "PhysiCell_cell_part.h"

using namespace BioFVM;

namespace PhysiCell{



Cell_Part::Cell_Part(){

}

Cell_Part::Cell_Part(Cell_* cell,int voxel_id){
	//1. Setting Cell_Part unique ID
	this.unique_id = Cell_Part_ID_Counter;
	Cell_Part_ID_Counter++;

	//2. Voxel ID
	this.voxel_id = voxel_id;

	//3. Cell ID (indicates the cell that this Cell_Part belongs to)
	this.cell_id = (*cell).unique_id;

	//4. Cell* pointer to containing cell
	this.containing_cell = cell;

	//5. Type (integer)
	this.type = cell_defaults.type; //not sure what this should be, just placing 'cell_defaults.type' for now

	//6. Pointer to Cell_Part_Definition
	this.cell_part_definition = NULL;

	//7. Position double[x,y,z]
	this.assign_position(0,0,0); //need to figure out the initial positions

	//8. Velocity double[x,y,z]
	this.set_velocity(0,0,0); //need to figure out the initial velocity

	//9. Previous Velocity double[x,y,z] (for Adams-Bashforth)
	this.set_previous_velocity(0,0,0);

	//10. Volume total
	this.total_volume = 1.0; //I believe this == basic_agent volume

	//11. Secretion rates vector<double>
	this.secretion_rates = new std::vector<double>(0);

	//12. Secretion targets, not sure what this refers to? set_target_radius, target_fluid, target_solid_nuclear

	//13. Uptake rates
	this.uptake_rates = new std::vector<double>(0);

	//14. net export rates
	this.net_export_rates = new std::vector<double>(0);

	//15. internalized substrates
	this.internalized_substrates = new std::vector<double>(0);

	//16. equivalent radius
	this.equivalent_radius = 8.412710547954228; //this is the Geometry.radius (not sure if this is the correct one?)

	//17. Repulsion strength (should this be a vector? does this include BM_adhesion_strength)
	this.repulsion_strength = 10.0; //located in Mechanics::Mechanics()

	//18. Adhesion strength (may be generalized to a vector? should this include BM_repulsion_strength?)
	this.adhesion_strength = 0.4; //cell_cell_adhesion_strength

	//19. Max adhesion distance (I believe this is relative_maximum_adhesion_distance)
	this.max_adhesion_distance = 1.25;

	//20. vector<Cell_Part*> (vector of mechanically interacting cell parts)
	std::vector<Cell_Part*> interacting_cell_parts;

	//21. vector<Cell*> (vector of mechanically interacting cells)
	std::vector<Cell*> interacting_cells;

	//22...
}

Cell_Part::assign_position(double x, double y, double z){
	this.position[0] = x;
	this.position[1] = y;
	this.position[2] = z;
}

Cell_Part::set_velocity(double x, double y, double z){
	this.velocity[0] = x;
	this.velocity[1] = y;
	this.velocity[2] = z;
}

Cell_Part::set_previous_velocity(double x, double y, double z){
	this.previous_velocity[0] = x;
	this.previous_velocity[1] = y;
	this.previous_velocity[2] = z;
}


};
