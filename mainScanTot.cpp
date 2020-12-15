// To compile: g++ -o analyzer mainScanTot.cpp `root-config --cflags --glibs`
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <cmath>
							// Needed root libraries
#include "TH1.h"			// Histogram
#include "TH2.h"
#include "TFile.h"			// Acces the tree structure: tree,branch,leaves...
#include "TBranch.h"
#include "TTree.h"
#include "TCanvas.h"		// To draw several graphs before printing it to the output root file
#include "TGraph.h"
#include "TGraphErrors.h"	// Graph with error bars
#include "TF1.h"
#include "TAxis.h"			// Change axis

using namespace std; 

int main(int argc, char *argv[])
{
 	

 	//Information contained in tree at generation (primaries) in TreeP
 	//Note that this information is recorded at generation, so it contains all generated particles	
 	//[0] 	Event_num 	Event number
 	//[1] 	Track_num 	Track number in event
 	//[2] 	PID 		Particle ID (-1 = electron, 0 = photon, +1 = positron) (DOES NOT FOLLOW STANDARD GEANT4 CONVENTIONS!)
 	//[4] 	E_kin 		Kinetic energy at generation

	//Information contained in tree at detection in TreeM
	//Note that this information is recorded at detection, so it will not contain all generated particles
	//[0] 	X 		Detected position
	//[1] 	Y		
	//[2] 	Z		
	//[3] 	PX		Detected momentum vector
	//[4] 	PY		
	//[5] 	PZ		
	//[6] 	Xi		Position at generation
	//[7] 	Yi		
	//[8] 	Zi		
	//[9] 	PXi		Momentum vector at generation
	//[10] 	PYi		
	//[11] 	PZi		
	//[12] 	R		sqrt(X^2 + Z^2) (note that this is a bad definition!) 
	//[13] 	Ener	Kinetic energy at detection
	//[14] 	Evt		Detected event number
	//[15] 	EnerP	Kinetic energy at generation
	//[16] 	Bin		Energy set in test.txt for this file


//------------------------------- 1: Initialization-----------------------------------------------------

 	TTree *TreeM; 
 	TTree *TreeP;

	//For the primary tree, make four separate branches, and four variables that will be linked to them. 
 	TBranch *BranchP_Event_num;
 	TBranch *BranchP_Track_num;
 	TBranch *BranchP_PID;
 	TBranch *BranchP_E_kin;

 	double PID = 0.; 
	double E_kin = 0.; 
	double Event_num = 0.; 
	double Track_num = 0.; 

	//For the detected event tree, BranchList is used to fill the parameters in a loop and not 17 individual calls 

 	TBranch *BranchM[17];
 	double data_var_det[17] = {0.}; 		// Create a c array we will fill with the treeM parameters
 	
 	const int N_file=154;					// Number of input files
 	double nTot[N_file]={0};				// Number of hits at the detector, will be filled later
 	const double slitSize=0.8; 				// "radius" [cm] from 0 to 8mm, x2 to have total width

 	
 	int energyScan[N_file]={1,10,20,30,40,50,60,70,80,90,100,    			 // List of available energies
 							110,120,130,140,150,160,170,180,190,200,
 							210,220,230,240,250,260,270,280,290,300,
 							310,320,330,340,350,360,370,380,390,400,
 							410,420,430,440,450,460,470,480,490,500,
 							510,520,530,540,550,560,570,580,590,594,
 							596,598,600,602,604,606,608,610,612,614,
 							616,618,620,622,624,626,628,630,632,634,
 							636,638,640,642,644,646,648,650,652,654,
 							656,658,660,662,664,666,668,670,672,674,
 							676,680,690,700,
 							710,720,730,740,750,760,770,780,790,800,
 							810,820,830,840,850,860,870,880,890,900,
 							910,920,930,940,950,960,970,980,990,1000,
 							1010,1020,1030,1040,1050,1060,1070,1080,1090,1100,
 							1110,1120,1130,1140,1150,1160,1170,1180,1190};

 	double energyScanDouble[N_file]={1,10,20,30,40,50,60,70,80,90,100,     		// List of the same available energies in double (root graphs doesn't like integers)
 							110,120,130,140,150,160,170,180,190,200,
 							210,220,230,240,250,260,270,280,290,300,
 							310,320,330,340,350,360,370,380,390,400,
 							410,420,430,440,450,460,470,480,490,500,
 							510,520,530,540,550,560,570,580,590,594,
 							596,598,600,602,604,606,608,610,612,614,
 							616,618,620,622,624,626,628,630,632,634,
 							636,638,640,642,644,646,648,650,652,654,
 							656,658,660,662,664,666,668,670,672,674,
 							676,680,690,700,
 							710,720,730,740,750,760,770,780,790,800,
 							810,820,830,840,850,860,870,880,890,900,
 							910,920,930,940,950,960,970,980,990,1000,
 							1010,1020,1030,1040,1050,1060,1070,1080,1090,1100,
 							1110,1120,1130,1140,1150,1160,1170,1180,1190}; 		// Tgrah need doubles for x coordinate

/* 	array<string,>energyScan={"1","10","20","30","40","50","60","70","80","90","100",
 							"110","120","130","140","150","160","170","180","190","200",
 							"210","220","230","240","250","260","270","280","290","300",
 							"310","320","330","340","350","360","370","380","390","400",
 							"410","420","430","440","450","460","470","480","490","500",
 							"510","520","530","540","550","560","570","580","590","594",
 							"596","598","600","602","604","606","608","610","612","614",
 							"616","618","620","622","624","626","628","630","632","634",
 							"636","638","640","642","644","646","648","650","652","654",
 							"656","658","660","662","664","666","668","670","672","674",
 							"676","680","690","700",
 							"710","720","730","740","750","760","770","780","790","800",
 							"810","820","830","840","850","860","870","880","890","900",
 							"910","920","930","940","950","960","970","980","990","1000",
 							"1010","1020","1030","1040","1050","1060","1070","1080","1090","1100",
 							"1110","1120","1130","1140","1150","1160","1170","1180","1190"}*/ // Copy in string version if needed




//------------------------------- 2: fill Ekin for each hit of each energy--------------------------------

	vector<vector<double>> E_kinHit={};					// Will be filled for each energy with the Ekin at generation of each particle that hit the detector
 	for(int i(0);i<N_file;i++){ 						// Loop on the selected energies with the magnetic field
 		string fileName("./Measurement/ExN02_"+to_string(energyScan[i])+"keV.root");
 	//cout<<fileName<<endl;
 		TFile *rootfile_in = TFile::Open(fileName.c_str(), "READ");
 		int nTotThisEnergy (0);
 		//Use list of branches for a tree to get the information out. 
 		TreeM = (TTree*)rootfile_in->Get("TreeM");
 		TreeP = (TTree*)rootfile_in->Get("TreeP");

 		for(int j=0; j<17; j++)				//Loop over the 17 variables of treeM to fill the parameters at detection
 		{
 			BranchM[j] = (TBranch*)TreeM->GetListOfBranches()->At(j); 
 			BranchM[j]->SetAddress(&data_var_det[j]); 
 		}

 		
		// Get the branches from the treeP: 
		BranchP_Event_num = (TBranch*)TreeP->GetBranch("Event_num");
		BranchP_Track_num = (TBranch*)TreeP->GetBranch("Track_num");
		BranchP_PID = (TBranch*)TreeP->GetBranch("PID");
		BranchP_E_kin = (TBranch*)TreeP->GetBranch("E_kin");

		// Link the branches to he variables that we made above. 
		BranchP_Event_num->SetAddress(&Event_num);
		BranchP_Track_num->SetAddress(&Track_num);
		BranchP_E_kin->SetAddress(&E_kin); 
		BranchP_PID->SetAddress(&PID); 



		//Get all the energies at generation: Is done with one file combining all the others for simplicity
		

		//Get the particles that reach the detector
		vector<double> E_kinHitThisEnergy={};

	 	for(int j=0; j<TreeM->GetEntries(); j++){	// For every energy:
			for(int q=0; q<17; q++){ 				// Extract needed parameters
			
				BranchM[q]->GetEntry(j); 
			}
			//if particle detected, put Ekin in tab:
			if(sqrt(pow((data_var_det[0]-18.),2)+pow((data_var_det[2]),2)) < 0.8 && abs(data_var_det[0]-18) < slitSize){

				E_kinHitThisEnergy.push_back(data_var_det[15]); 			// Fill the tab
				nTotThisEnergy+=1;											// counter of hits
			}
		}

		nTot[i]=nTotThisEnergy;
		


	//---------- 3: Create root file for each energy to later display the fit-------------------
 	string fileName_energy="./Output/output_plots_energy_"+to_string(energyScan[i])+".root";
	TFile *rootfile_out = new TFile(fileName_energy.c_str(), "RECREATE"); 
 	//Histogram to verify the fits
 	//Syntax for initialization: TH1D("Object name", "title of plot", number of bins, minimum, maximum)
	TH1D *plot_E_kin = new TH1D("plot_E_kin", "plot_E_kin", 81, energyScan[i]-40, energyScan[i]+40); 	// +-20 around the selected energy
	
	for (size_t k (0);k<E_kinHitThisEnergy.size();k++){
		plot_E_kin -> Fill(E_kinHitThisEnergy[k]*1000); 	// Conversion from MeV to keV
	}
	plot_E_kin -> Write();
	rootfile_out->Close();



 	rootfile_in->Close(); 									// Close the input file before going to the next
 	}	// End of energy loop


 	//Write the full energy result
 	string fileNameInAll("./Measurement/sum_file.root");	// Merge of all the others
 	TFile *rootfile_inAll = TFile::Open(fileNameInAll.c_str(), "READ");
 	TreeP = (TTree*)rootfile_inAll->Get("TreeP");
 	BranchP_E_kin = (TBranch*)TreeP->GetBranch("E_kin");
	BranchP_E_kin->SetAddress(&E_kin);

	string fileEnergyAll("./Output/output_EnergyAll.root");
 	TFile *rootfile_outAll = TFile::Open(fileEnergyAll.c_str(), "RECREATE"); 

 	//Syntax for initialization: TH1D("Object name", "title of plot", number of bins, minimum, maximum)
 	int nDiv(1200);
	TH1D *plot_E_kinAll = new TH1D("E_kinAll", "E_kinAll", nDiv, 0.0, 1.2);		// >0 to take only the non empty ones
	cout<<"Nb enteries total: "<<TreeP -> GetEntries()<<endl;
	
	int counterNorm (0);
	int counterTallPeak(0);
	for(size_t j(0); j<TreeP -> GetEntries();j++){							 	// Fill the TH1D with energy at generation
		BranchP_E_kin-> GetEntry(j); 											// Link E_kin(j) to the double E_kin

		
		if(	   abs(E_kin-0.624216) >pow(10,-5) 									// Cut the peaks
			&& abs(E_kin-0.661659) >pow(10,-5) 
			&& abs(E_kin-0.00425)  >pow(10,-5) 
			//&& abs(E_kin-0.0374027)>pow(10,-5)
			&& abs(E_kin-0.031)  >pow(10,-5)
			&& abs(E_kin-0.026)  >pow(10,-5)
			&& abs(E_kin-0.037)  >pow(10,-5)
			&& abs(E_kin-0.03206)  >pow(10,-5)
			&& abs(E_kin-0.656)    >pow(10,-5)
			&& abs(E_kin-0.0045)   >pow(10,-5)){

			plot_E_kinAll -> Fill(E_kin);
			counterNorm++;
			//cout<<E_kin<<endl;
		}
		if(abs(E_kin-0.624216) <pow(10,-5) ){counterTallPeak++;}
		if(j%1000000==0){cout<<j<<endl;}
	}
	cout<<"Normalze counter: "<<counterNorm<<endl;								// To compare with experimental data
	cout<<"Number of generation at the highest peak: "<<counterTallPeak<<endl;

	cout<<"Normalized position:"<<endl;
	for(int k(1);k<=nDiv;k++){
		cout<<plot_E_kinAll->GetBinCenter(k)<< " " <<plot_E_kinAll->GetBinContent(k)/counterNorm<<endl;

	}


	plot_E_kinAll ->Write(); 


	rootfile_outAll->Close();
	rootfile_inAll ->Close();

 	TFile *rootfile_out = new TFile("plot_energy_nTot_scan.root", "RECREATE");	// Plot of the counts on the full spectrum
 	TGraph *plot_energy_nTot_scan = new TGraph(N_file, energyScanDouble, nTot);

	 	plot_energy_nTot_scan->Draw("AL*");
	 	plot_energy_nTot_scan->Write();
	 	rootfile_out->Close();

	
}
