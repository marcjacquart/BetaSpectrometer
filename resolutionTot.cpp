// To compile: g++ -o analyzer resolutionTot.cpp `root-config --cflags --glibs`


#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <cmath>
								// Needed root libraries
#include "TH1.h"				// Histogram
#include "TH2.h"
#include "TFile.h"				// Acces the tree structure: tree,branch,leaves...
#include "TBranch.h"
#include "TTree.h"
#include "TCanvas.h"			// To draw several graphs before printing it to the output root file
#include "TGraph.h"
#include "TGraphErrors.h"		// Graph with error bars
#include "TF1.h"
#include "TAxis.h"				// Change axis
#include "TLegend.h"			// Insert a legend

using namespace std; 

int main(int argc, char *argv[]){

	//Variables:
	//int nSlits=9;																// length of the slitSize tab, easier to create a new variable
	//double slitSize[nSlits]={0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6};		// Slit sizes chosen /!\ It is the radius, x2 for the total slit size
	int nSlits=4;
	double slitSize[nSlits]={0.2,0.3,0.4,0.6};
	string slitSizeName[nSlits]={""};											// Formated values of the slit in cm for the legend
	for(int i(0);i<nSlits;i++){ 												// Not elegant but easy way to display only one decimal
		double slitHere=20*slitSize[i];
		int rest = fmod(slitHere,10.0);
		int decimal = 0.1*slitHere;
		slitSizeName[i]=to_string(decimal)+"."+to_string(rest);
	}

	//int nFiles(20);															// Number of simulation files at each energy: defined later beacause varies (mor at lower energies)
	int nScanEmitted(15);
	
	double emittedEnergy[nScanEmitted]={40.0,60.0,80.0,100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,1100.0,1200.0}; // Vlaues of energy simulated
	double sigma[nSlits][nScanEmitted]={0.0};									// Result tab: sigma of gaussian fit
	double errSigma[nSlits][nScanEmitted]={0.0};								// Error on sigma
	double mean[nSlits][nScanEmitted]={0.0};									// Result tab: difference to the mean of gaussian fit (to point out the magnetic field miscalibration)
	double errMean[nSlits][nScanEmitted]={0.0};									// Error on mean difference
	double nullArray[nScanEmitted]={0.0};										// array of 0.0 to display zero error on x axis on the TGraphErrors


	double sigmaThisEmittedEnergy[nScanEmitted]={0.0};							// Temporary storage array, will be overwritten at each slit size loop



	TTree *TreeP;																// Primary tree: contains parameters at generation
	TTree *TreeM;																// Detected particles tree



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

	bool isSpecial(false);										// special: low energy, has more files than 20, need 35, for more precision arounf the peak.
	double nMaxESlitVal [nSlits][nScanEmitted]={0.0};			// Max counts for one slit, needs this for magnetic efficiency because the max counts are NOT with selected energy=generated energy
	double errNMaxESlitVal [nSlits][nScanEmitted]={0.0};		// Stastistical error
	double nMaxESlitGaus [nSlits][nScanEmitted]={0.0};			// The same with maximum value of the gaussian fit
	for(int slitIndex(0);slitIndex<nSlits;slitIndex++){
		cout<<"Slit size: "<<2*slitSize[slitIndex]<<" cm"<<endl;


		for(int emitIndex(0);emitIndex<nScanEmitted;emitIndex++){
			double energyHere=emittedEnergy[emitIndex];

			int nFiles=0;
			if(energyHere==40.0||energyHere==60.0||energyHere==80.0||energyHere==100.0||energyHere==200.0){nFiles=35;isSpecial=true;}
			else{
				nFiles=20;
				isSpecial=false;
			}
			double nTotThisEmittedEnergy [nFiles]={0}; 															// Is normally integer but TGraph accepts only doubles

			double energiesSelectedThisEmittedEnergy[nFiles]={0.0};												// Will be filled with the 20 or 35 selected energies around the emitted value

			for (int k(0);k<nFiles;k++){																		// 8 files with step 2keV, then 20 with step 0.5 then again 7 files with 2keV steps
				if(isSpecial){ 																					// Tie togeter simulation index
					if(k<8){energiesSelectedThisEmittedEnergy[k]=emittedEnergy[emitIndex]-20+2*k;} 
					else if(k>27){energiesSelectedThisEmittedEnergy[k]=emittedEnergy[emitIndex]-20+2*(k-15);} 
					else {energiesSelectedThisEmittedEnergy[k]=emittedEnergy[emitIndex]-9+0.5*k;} 				// To align indicises
				}


				else{energiesSelectedThisEmittedEnergy[k]=emittedEnergy[emitIndex]-20+2*k;}
				
			}
			for(int i(0);i<nFiles;i++){																			// For all the files for a given emitted energy:
				int emittedEnergyName=emittedEnergy[emitIndex];

				int energiesSelectedThisEmittedEnergyName=energiesSelectedThisEmittedEnergy[i];
				string fileName;																				// Root file with results of the simulation to open
				if(fmod(energiesSelectedThisEmittedEnergy[i],1.0)==0.0){fileName="./SimulationFiles/Epeak"+to_string(emittedEnergyName)+"_"+to_string(energiesSelectedThisEmittedEnergyName)+"keV.root";} // if integer no need to add .5 to the name
			 	else/*case .5*/{fileName="./SimulationFiles/Epeak"+to_string(emittedEnergyName)+"_"+to_string(energiesSelectedThisEmittedEnergyName)+".5keV.root";}
			 	TFile *rootfile_in = TFile::Open(fileName.c_str(), "READ"); 									// Open the root input file for a given selected energy
			 	

			 	//------------Extract parameters for TreeP:
				TreeP = (TTree*)rootfile_in->Get("TreeP");

				//For the primary tree, make four separate branches, and four variables that will be linked to them. 
			 	TBranch *BranchP_Event_num;
			 	TBranch *BranchP_Track_num;
			 	TBranch *BranchP_PID;
			 	TBranch *BranchP_E_kin;

			 	double PID = 0.; 
				double E_kin = 0.; 
				double Event_num = 0.; 
				double Track_num = 0.; 

				//Get the branches from the treeP: 
				BranchP_Event_num = (TBranch*)TreeP->GetBranch("Event_num");
				BranchP_Track_num = (TBranch*)TreeP->GetBranch("Track_num");
				BranchP_PID = (TBranch*)TreeP->GetBranch("PID");
				BranchP_E_kin = (TBranch*)TreeP->GetBranch("E_kin");


				//Link the branches to he variables that we made above. 
				BranchP_Event_num->SetAddress(&Event_num);
				BranchP_Track_num->SetAddress(&Track_num);
				BranchP_E_kin->SetAddress(&E_kin); 
				BranchP_PID->SetAddress(&PID); 

				

				// Same for TreeM:
				TBranch *BranchM[17];
			 	double data_var_det[17] = {0.}; 					// Create a C array, will be filled with the 17 differents parameters available at detection

			 	TreeM = (TTree*)rootfile_in->Get("TreeM");

			 	for(int j=0; j<17; j++){							// Fill the tab from the data of treeM
					BranchM[j] = (TBranch*)TreeM->GetListOfBranches()->At(j); 
					BranchM[j]->SetAddress(&data_var_det[j]); 
				}
				int nTotThisB(0);
				for(int j=0; j<TreeM->GetEntries(); j++){
					for(int q=0; q<3; q++){ 						// Extract needed parameters
					
						BranchM[q]->GetEntry(j); 
					}
					
					//if particle detected, put Ekin in tab:
					if(sqrt(pow((data_var_det[0]-18.),2)+pow((data_var_det[2]),2)) < 0.8 && abs(data_var_det[0]-18) < slitSize[slitIndex]){ // If hits the detector, we give the boudaries corresponding to the chosen slit size
						nTotThisB+=1;	
					}
				}
				nTotThisEmittedEnergy[i]=nTotThisB; 				// Fill the result in the tab before moving to the next file
				rootfile_in->Close();								// Close the input file
			}
			double resultMax=0.0;
			for(int maxIndex(0);maxIndex<nFiles;maxIndex++){		// Find the maximal recoreded value to fit the tab (for the magnetic efficiency)
				if (nTotThisEmittedEnergy[maxIndex]>resultMax){resultMax=nTotThisEmittedEnergy[maxIndex];}
			}
			nMaxESlitVal[slitIndex][emitIndex]=resultMax*0.009923/50000; // 50k emitted particle, 1% from emission angle for performances
			errNMaxESlitVal[slitIndex][emitIndex]=sqrt(resultMax)*0.009923/50000;
			
			//plots:
			string fileName_emittedEnergy="./Output/output_plots_emittedEnergy_"+to_string(emittedEnergy[emitIndex])+"_slit"+to_string(slitSize[slitIndex])+".root";
			TFile *rootfile_out_emitted = new TFile(fileName_emittedEnergy.c_str(), "RECREATE"); // Create new root file to fit and display

			TGraph *plot_energy_emitted = new TGraph(nFiles, energiesSelectedThisEmittedEnergy, nTotThisEmittedEnergy);
			TCanvas *q1 = new TCanvas("q1", "q1", 100, 100, 800, 800);
			q1->cd(); 																			// To be sure we will draw on q1 and not on a previous one
		 	plot_energy_emitted->Draw("AL*");													// Draw the result
		 	auto legend1 = new TLegend();														// Add legend
		 	plot_energy_emitted -> SetTitle("");												// Format title and axis
			plot_energy_emitted -> GetXaxis() -> SetTitle("Energy [keV]");
			plot_energy_emitted -> GetYaxis() -> SetTitle("Events [ ]");

			//----------------------Fit gaussien--------------------------

			TF1 *fit =new TF1("fit","gaus",emittedEnergy[emitIndex]-50,emittedEnergy[emitIndex]+50); 	// Fit around the expected mean value
		 	plot_energy_emitted->Fit("fit","Q"); 
		 	nMaxESlitGaus[slitIndex][emitIndex]=fit ->Eval(fit->GetParameter(1))*0.009923/50000;		// Normalize by number of emission
		 
		 	fit->Draw(); 																				// Draw the fit on q1

		 	sigma[slitIndex][emitIndex]=fit->GetParameter(2); 											// Fill the array for general plot
		 	errSigma[slitIndex][emitIndex]= fit ->GetParError(2);
		 	mean[slitIndex][emitIndex]=fit->GetParameter(1)-emittedEnergy[emitIndex];
		 	errMean[slitIndex][emitIndex]=fit->GetParError(1);

		 	legend1 -> AddEntry(plot_energy_emitted,"Measurements","pe");								// Add legend entery
		 	legend1 -> AddEntry(fit,"Gaussian fit","pe");
		 	legend1 -> Draw();
		 	q1->Write();
			plot_energy_emitted->Write();		
			q1 -> Close();														// Write canevas on output file
		 	rootfile_out_emitted ->Close();
		}	// End loop over emitted energies

	}	// End of loop over slit sizes






	//General plot with sigma at all energies and all slit sizes
	
	TFile *rootfile_all = new TFile("./Output/output_plots_sigmaAllEnergy.root","RECREATE");
	TCanvas *qall = new TCanvas("q_all", "q_all", 100, 100, 800, 800);
	qall->cd();


	//First of the array, the following in the loop

	TGraphErrors *graphAll = new TGraphErrors(nScanEmitted,emittedEnergy,sigma[0],nullArray,errSigma[0]); // Plot the first one with special parameters and axis definition
	graphAll -> SetMarkerColor(4);
	graphAll -> SetLineColor(4);
	graphAll -> SetMarkerStyle(nSlits+1);
	graphAll -> Draw("ALP");

	graphAll -> SetTitle("Energy resolution");
	graphAll -> GetXaxis() -> SetTitle("Energy [keV]");
	graphAll -> GetYaxis() -> SetTitle("sigma [keV]");
	graphAll -> GetYaxis() -> SetRange(0.0,20.0);
	


	auto legend = new TLegend();
	string legendName("Slit size: "+slitSizeName[0]+" cm");
	legend -> AddEntry(graphAll,legendName.c_str(),"pe");


	

	for(int slitIndex(1);slitIndex<nSlits;slitIndex++){ 												// Plot the other on top
		qall->cd();
		TGraphErrors *graphAll = new TGraphErrors(nScanEmitted,emittedEnergy,sigma[slitIndex],nullArray,errSigma[slitIndex]);
		
	 	graphAll -> SetMarkerColor(slitIndex);
	 	graphAll -> SetLineColor(slitIndex);
	 	graphAll -> SetMarkerStyle(slitIndex);


		graphAll ->Draw("SAMELP");


		string legendName("Slit size: "+slitSizeName[slitIndex]+" cm");
		legend -> AddEntry(graphAll,legendName.c_str(),"pe");



	}

	qall->cd();
	legend -> Draw();
	qall ->Write();																						// Write on output file
	qall -> Close();
	rootfile_all -> Close();





//Same for the mean
	
	TFile *rootfile_mean = new TFile("./Output/output_plots_meanAllEnergy.root","RECREATE");
	TCanvas *qmean = new TCanvas("q_all", "q_all", 100, 100, 800, 800);
	qmean->cd();


	//First of the array, the following in the loop

	TGraphErrors *graphMean = new TGraphErrors(nScanEmitted,emittedEnergy,mean[0],nullArray,errMean[0]);
	graphMean -> SetMarkerColor(kBlue);
	graphMean -> SetLineColor(kBlue);
	graphMean -> SetMarkerStyle(0);
	graphMean -> Draw("AL*");

	graphMean -> SetTitle("Mean of detected energy");
	graphMean -> GetXaxis() -> SetTitle("Energy [keV]");
	graphMean -> GetYaxis() -> SetTitle("#Delta#mu [keV]");
	graphMean -> GetYaxis() -> SetRange(0.0,20.0);
	


	auto legendMean = new TLegend();
	legendName="Slit size: "+slitSizeName[0]+" cm";
	legendMean -> AddEntry(graphMean,legendName.c_str(),"pe");


	

	for(int slitIndex(1);slitIndex<nSlits;slitIndex++){ 
		qmean->cd();
		TGraphErrors *graphMean = new TGraphErrors(nScanEmitted,emittedEnergy,mean[slitIndex],nullArray,errMean[slitIndex]);
		
	 	graphMean -> SetMarkerColor(slitIndex);
	 	graphMean -> SetLineColor(slitIndex);
	 	graphMean -> SetMarkerStyle(slitIndex);


		graphMean ->Draw("SAMELP");


		string legendName("Slit size: "+slitSizeName[slitIndex]+" cm");
		legendMean -> AddEntry(graphMean,legendName.c_str(),"pe");

	}

	qmean->cd();
	legendMean -> Draw();
	qmean ->Write();
	qmean ->Close();
	rootfile_mean -> Close();





	// Same for the magnetic efficiency:
	TFile *rootfile_efficiency=new TFile("./Output/output_plots_magneticEfficiency.root","RECREATE");
	TCanvas *q_eff = new TCanvas("q_eff", "q_eff", 100, 100, 800, 800);
	q_eff->cd();
	auto legend_eff = new TLegend();

	// First of the array, the following in the loop

	TGraphErrors *graphEff = new TGraphErrors(nScanEmitted,emittedEnergy,nMaxESlitVal[0],nullArray,errNMaxESlitVal[0]);
	graphEff -> SetMarkerColor(kBlue);
	graphEff -> SetLineColor(kBlue);
	graphEff -> SetMarkerStyle(5);
	graphEff -> Draw("ALP");

	graphEff -> SetTitle("magnetic efficiency");
	graphEff -> GetXaxis() -> SetTitle("Energy [keV]");
	graphEff -> GetYaxis() -> SetTitle("magnetic efficiency [-]");
	graphEff -> GetYaxis() -> SetRange(0.0,20.0);
	graphEff -> GetYaxis() -> SetMaxDigits(3);


	
	legendName="Slit size: "+slitSizeName[0]+" cm";
	legend_eff -> AddEntry(graphMean,legendName.c_str(),"LP");


	

	for(int slitIndex(1);slitIndex<nSlits;slitIndex++){ //put 1 for the GraphAll
		q_eff->cd();
		TGraphErrors *graphEff = new TGraphErrors(nScanEmitted,emittedEnergy,nMaxESlitVal[slitIndex],nullArray,errNMaxESlitVal[slitIndex]);
		
	 	graphEff -> SetMarkerColor(slitIndex);
	 	graphEff -> SetLineColor(slitIndex);
	 	graphEff -> SetMarkerStyle(slitIndex);


		graphEff ->Draw("SAMELP");


		string legendName("Slit size: "+slitSizeName[slitIndex]+" cm");
		legend_eff -> AddEntry(graphEff,legendName.c_str(),"LP");



	}

	q_eff->cd();
	legend_eff -> Draw();
	q_eff ->Write();
	q_eff -> Close();
	rootfile_efficiency -> Close();


} // ENd of the main