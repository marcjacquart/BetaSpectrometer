// To compile: g++ -o analyzer mainScan.cpp `root-config --cflags --glibs`


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
#include "TLegend.h"		// Insert a legend

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

 	TTree *TreeM; 

	//For the detected event tree, BranckList is used to fill the parameters in a loop
 	TBranch *BranchM[17];
 	double data_var_det[17] = {0.}; 			// Create c array to fill the parameters
 	
 	const int N_file=22;						// Number of input files
 	double data_scan[N_file]={0}; 				// 22 files
 	double energy_scan[N_file]={0};

 	const int NSlitScan(81); 					// Number of slit scans
 	double slit_size_tab [NSlitScan]={0}; 		// Tab with result of slit size
 	double sigma[NSlitScan]={0}; 				// Final result of sigma gaussian fit
 	double sigmaError[NSlitScan]={0}; 			// Associated error
 	double mean[NSlitScan]={0}; 				// Mean of the gaussian fit
 	double meanError[NSlitScan]={0};			// Associated error
 	double nullArray[NSlitScan]={0};			// For null error on x axis

 	double sigmaOverMean[NSlitScan]={0};		// Result sigma over mean
 	double sigmaOverMeanError[NSlitScan]={0};
 	double sigmaQuad[NSlitScan]={0}; 			// Result with FWHM for larger slit sizes when gaussian isn't a good fit
 	double sigmaQuadError[NSlitScan]={0};
 	vector<vector<int>> Ntot (NSlitScan,vector<int>(N_file)); // Number of hits for each input file and slit size
 	int counter(0); 							// File index to open

 	//Measurements values:
 	double energyMeasurement[4]={1.6,1.2,0.8,0.4};
 	double resolutionMeasurement[4]={14.75,11.99,11.35,9.17};
 	double errorMeasurement[4]={0.77,0.46,0.40,0.54};




 	for(int i(594);i<637;i+=2){ 				// Loop over the files, 594keV then 2keV steps
 		string fileName("./Measurement/ExN02_"+to_string(i)+"keV.root");
 		TFile *rootfile_in = TFile::Open(fileName.c_str(), "READ");
 		//Use list of branches for a tree to get the information out. 
 		TreeM = (TTree*)rootfile_in->Get("TreeM");
 		for(int j=0; j<17; j++)					// Loop over the 17 parameters to fill the array
 		{
 			BranchM[j] = (TBranch*)TreeM->GetListOfBranches()->At(j); 
 			BranchM[j]->SetAddress(&data_var_det[j]); 
 		}

 		

	 	for(int j=0; j<TreeM->GetEntries(); j++)
		{
			for(int q=0; q<3; q++) 				// We don't need the other parameters for the following analysis, saves time
			{
				BranchM[q]->GetEntry(j); 
			}

			int counter_slit(0);

			for(double slitSize(0.0);slitSize<=0.8001;slitSize+=(0.8/(NSlitScan-1))){
				if (j==1){slit_size_tab[counter_slit]=2*slitSize;} //2 fois pour avoir la grandeur totale
			
				if(sqrt(pow((data_var_det[0]-18.),2)+pow((data_var_det[2]),2)) < 0.8 && abs(data_var_det[0]-18) < slitSize)
					{
					Ntot[counter_slit][counter]+=1;
					}
				counter_slit++;
			}

		}

		energy_scan[counter]=i; 	// End of loop, we fill the result for the corresponding energy
 		counter++;					// Update energy couner to go to the next field
 		rootfile_in->Close(); 		// Close file before going to the next
 	}
 
 	// One plot each slit size:
 	for (int k(0);k<NSlitScan;k++){
 		int slitSize_forName =(0.8/(NSlitScan-1))*k*1000; 	// *1000 to have different names and avoid overw
 		string fileName_slit="./Output/output_plots_slit_"+to_string(slitSize_forName)+".root";
	 	TFile *rootfile_out = new TFile(fileName_slit.c_str(), "RECREATE"); 

	 	TCanvas *q1 =new TCanvas("q1","q1",100,100,800,800);
	 	q1->cd();


		double temp_array[N_file]={0};
		for(int q(0);q<N_file;q++){
			temp_array[q]=Ntot.at(k).at(q);

		}


	 	TGraph *plot_energy_scan = new TGraph(N_file, energy_scan, temp_array);

	 	auto legend1 = new TLegend();
	 	legend1 -> AddEntry(plot_energy_scan,"Measurements","pe");
	 	plot_energy_scan->Draw("AL*");
	 	


		//---------------------- Gaussian fit --------------------------

		TF1 *fit =new TF1("fit","gaus",600,700);
		legend1 -> AddEntry(fit,"Gaussian fit","pe");
	 	plot_energy_scan->Fit("fit","Q"); 	// No need to restrict fit here


	 	double intervalMin=fit->GetParameter(1)-1.5*fit->GetParameter(2);
	 	double intervalMax=fit->GetParameter(1)+1.5*fit->GetParameter(2);






	 	//---------------------- Quadratic fit --------------------------
	 	
	 	cout<<"***fit"<<endl;
	 	TF1 *fitQuad =new TF1("fitQuad","[0]+[1]*x+[2]*x^2",intervalMin,intervalMax);//[0]-[2]*(x-[1])^2
		fitQuad->SetLineColor(kBlue);

	 	plot_energy_scan->Fit("fitQuad","RQ"); 	// R: limited to the range, Q:"quiet" ->don't show the root infos of the fit on terminal
		legend1 -> AddEntry(fitQuad,"Quadratic fit","pe");
		double a = fitQuad -> GetParameter(2);
		double b = fitQuad -> GetParameter(1);
		double c = fitQuad -> GetParameter(0);
		cout<<"a: "<<a<<" b: "<<b<<" c: "<<c<<endl;	 	

	 	double errB=fitQuad -> GetParError(1);	//ax^2+bx+c
		double errA=fitQuad -> GetParError(2);
		double errC=fitQuad -> GetParError(0);
		cout<<"Err a: "<<errA<<" Err b: "<<errB<<" Err c: "<<errC<<endl;


		//error:
		double sqt= sqrt(b*b-4*a*c);	
		cout<<"sqt: "<<sqt<<endl;


		double dFWHMda= abs(((2*a*c/sqt)+sqt)/(sqrt(2)*a*a));			// Check mathematical definition of FWHM error in the report
		cout<<"dFWHMda: "<< dFWHMda <<endl;
		double dFWHMdb=	abs(b/(sqrt(2)*sqt*a));					
		cout<<"dFWHMdb: "<< dFWHMdb <<endl;
		double dFWHMdc=	abs(sqrt(2)/sqt);							
		cout<<"dFWHMdc: "<< dFWHMdc <<endl;

		
		sigmaQuad[k]=abs(sqt/(sqrt(2)*a*2.355));
		cout<<"SigmaQuad: "<<sigmaQuad[k]<<endl;
		cout<<"***fill sigmaQuadError"<<endl;
		cout<<"3 factors: "<<dFWHMda*errA<<" ; "<<dFWHMdb*errB<<" ; "<<dFWHMdc*errC<<endl;
		sigmaQuadError[k]=(dFWHMda*errA+dFWHMdb*errB+dFWHMdc*errC)/2.355;
		cout<<"------->sigmaQuadError "<<k<<": "<<sigmaQuadError[k]<<endl;





	 	TCanvas *q_new = new TCanvas("q_new", "q_new", 100, 100, 800, 800);		// Plot of hits over full spectrum
	 	plot_energy_scan -> SetTitle(("Counts at  "+to_string(slitSize_forName/50)+" mm slit size.").c_str());
		plot_energy_scan -> GetXaxis() -> SetTitle("Selected energy [keV]");
		plot_energy_scan -> GetYaxis() -> SetTitle("Hits [-]");	 	
	 	plot_energy_scan->Draw("AL*");
	 	q_new->cd();
	 	legend1 -> Draw();

	 	fit->Draw("SAMEL");
	 	fitQuad->Draw("SAMEL");

	 	sigma[k]=fit->GetParameter(2); 											// 0: constant, 1: mean , 2: sigma
	 	sigmaError[k]=fit -> GetParError(2);
	 	mean[k]=fit->GetParameter(1); 
	 	meanError[k]=fit -> GetParError(1);
	 	
	 	fit->Write();
	 	q1->Write();
	 	plot_energy_scan->Write();
	 	q_new->Write();
	 	rootfile_out->Close();
	 	cout<<"Fit quadratique "<<k<<": "<< fitQuad ->GetParameter(2)<<" x^2+ "<< fitQuad ->GetParameter(1)<<" x+ "<< fitQuad ->GetParameter(0)<<endl;

	 }


	 for (int m(0);m<NSlitScan-1;m++){ 						// Ugly but easy and fast way to delete first element of array
	 	sigma[m]=sigma[m+1];
	 	sigmaError[m]=sigmaError[m+1];
	 	mean[m]=mean[m+1];
	 	sigmaQuad[m]=sigmaQuad[m+1];
	 	sigmaQuadError[m]=sigmaQuadError[m+1];
	 	meanError[m]=meanError[m+1];
	 	sigmaOverMean[m]=sigma[m]/mean[m]; 					// Fill sigmaOverMean: ugly but easy
	 	sigmaOverMeanError[m]=sigmaOverMean[m]*sqrt(pow((sigmaError[m]/sigma[m]),2)+pow((meanError[m]/mean[m]),2));

	 	slit_size_tab[m]=slit_size_tab[m+1];

	 }
	 sigmaOverMean[NSlitScan]=sigmaOverMean[NSlitScan-1]; 	// For the last element not to be empty
	 sigmaOverMeanError[NSlitScan]=sigmaOverMeanError[NSlitScan-1];
	 
	 //Graph sigma:
	 TFile *rootfile_out_sigma = new TFile("./Output/output_plots_sigma.root", "RECREATE"); 									// Create final graph
	 TGraph *plot_slit_sigma = new TGraphErrors(NSlitScan-1,  slit_size_tab, sigma,nullArray,sigmaError);
	 TCanvas *q2 =new TCanvas("q2","q2",100,100,800,800);
	 TGraphErrors *plot_slit_sigma_quad = new TGraphErrors(NSlitScan-1,  slit_size_tab, sigmaQuad,nullArray,sigmaQuadError);	// N,x,y,err_x,err_y

	 plot_slit_sigma_quad -> SetLineColor(kBlue);																				// Color and shape settings
	 plot_slit_sigma_quad -> SetMarkerColor(kBlue);
	 plot_slit_sigma_quad -> SetMarkerStyle(4);
	 double zeroFour[4]={0.0};

	 TGraphErrors *plot_measurement = new TGraphErrors(4,energyMeasurement,resolutionMeasurement,zeroFour,errorMeasurement);

	 //Axes:
	 plot_slit_sigma -> SetTitle("Energy resolution");
	 plot_slit_sigma -> GetXaxis() -> SetTitle("Slit size [cm]");
	 plot_slit_sigma -> GetYaxis() -> SetTitle("sigma [keV]");
	 plot_slit_sigma -> SetMaximum(16);
	 plot_slit_sigma -> SetMinimum(0);
	 q2->cd();
	 plot_slit_sigma->Draw("AL*");											 	// A: draw axis, L:connect dots with a line, *:markers 
	 plot_slit_sigma_quad ->Draw("SAMELP"); 									// P:Point --> 4 (circle)
	 plot_slit_sigma_quad -> SetLineColor(kBlue);
	 plot_slit_sigma_quad -> SetMarkerColor(kBlue);
	 plot_measurement-> SetMarkerColor(2);
	 plot_measurement -> SetLineColor(2);
	 plot_measurement ->SetMarkerStyle(5);
	 plot_measurement -> Draw("SAMEP");

	 auto legend = new TLegend();
	 legend->AddEntry(plot_measurement,"Measurements","pe");
	 legend->AddEntry(plot_slit_sigma,"Simulations, gaussian fit","pe");
	 legend->AddEntry(plot_slit_sigma_quad,"Simulations, quadratic fit","pe");
	 legend -> Draw();

	 q2->Write();




	 plot_slit_sigma->Write();
	 rootfile_out_sigma->Close();


	 //Graph Mean
	 TFile *rootfile_out_mean = new TFile("./Output/output_plots_mean.root", "RECREATE"); 
	 TGraph *plot_slit_mean = new TGraph(NSlitScan, slit_size_tab, mean);
	 TCanvas *q3 =new TCanvas("q3","q3",100,100,800,800);
	 q3->cd();

	 //Axes:
	 plot_slit_mean -> SetTitle("Gaussian fit pic Energy");
	 plot_slit_mean -> GetXaxis() -> SetTitle("Slit size [cm]");
	 plot_slit_mean -> GetYaxis() -> SetTitle("Mean [keV]");

	 plot_slit_mean->Draw("AL*"); 								// A: draw axis, L:connect dots with a line, *:markers
	 q3->Write();
	 
	 plot_slit_mean->Write();
	 rootfile_out_mean->Close();


	 //Graph Sigma over Mean
	 TFile *rootfile_out_sigmaOverMean = new TFile("./Output/output_plots_sigmaOverMean.root", "RECREATE"); 
	 TGraphErrors *plot_slit_sigmaOverMean = new TGraphErrors(NSlitScan-1, slit_size_tab, sigmaOverMean,nullArray,sigmaOverMeanError); // N,x,y,err_x,err_y
	 TCanvas *q4 =new TCanvas("q4","q4",100,100,800,800); 
	 q4->cd();

	 //Axes:
	 plot_slit_sigmaOverMean -> SetTitle("Gaussian: Sigma over mean");
	 plot_slit_sigmaOverMean -> GetXaxis() -> SetTitle("Slit size [cm]");
	 plot_slit_sigmaOverMean -> GetYaxis() -> SetTitle("Sigma/Mean [ ]");

	 plot_slit_sigmaOverMean->Draw("AL*"); //A: draw axis, L:connect dots with a line, *:markers
	 q4->Write();
	 
	 plot_slit_sigmaOverMean->Write();
	 rootfile_out_sigmaOverMean->Close();
}
