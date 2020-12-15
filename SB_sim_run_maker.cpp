#include <iostream>
#include <fstream>
#include <array>
#include <cmath> //for fmod function
using namespace std; 

int main() {


	std::ofstream ofs ("SB_sim_min_runner.sh", std::ofstream::out);

	// What do we need to do with the following commands:
	// Remove old test.txt file
	// Make new one with correct energy	"echo <value> > test.txt"
	// Run geant thinger
	// Rename root file
	int const nEmitted =5;
	array<double,nEmitted> emittedEnergies ={40.0,60.0,80.0,100.0,200.0}; 	// Emitted energies we want to run

	int N_iter =20; 
	double E_step = 0.5; 													// Steps for the magnetic seleced energy scan
	for (int k(0);k<nEmitted;k++){
		double E_start = emittedEnergies[k]-N_iter*E_step/2.0; 
		int EnergyEmittedName=emittedEnergies[k]; 							// Change to int to write
		//Change the emitted enrgy:
		ofs << "rm particle_energy.txt" << endl; 
		ofs << "echo " << emittedEnergies[k] << " > particle_energy.txt" << endl; 

		double E_value = 0; 

		for(int i=0; i<N_iter; i++)
		{

			E_value = double(E_start + i*E_step); 
			
			if(fmod(E_value,2.0)!=0.0){										// /!\ Was added to simulate only the 0.5 steps with 2.0 steps already computed. Remove this if statement to do a full scan

				ofs << "rm test.txt" << endl; 
				ofs << "echo " << E_value << " > test.txt" << endl; 
				ofs << "./exampleN02 run2.mac" << endl; 
				ofs << "mv ExN02.root ../ResEtot/SimulationFiles/Epeak"<<EnergyEmittedName<<"_" << ("%g",E_value) << "keV.root" << endl; 
				ofs << endl; 
			}
		}
	}
    
    return 0;
}