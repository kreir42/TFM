/*
**************************************************************************************************************
*
* February 2023. 
* Takes the TOF*.root files from analysisDAQfiles.C and gets the TOF and Energy spectra and 
* important parameter for the analysis.
*
**************************************************************************************************************
*/

//C,C++ Libraries
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <string.h>

// ROOT libraries
#include "TROOT.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TChain.h"
#include "TObject.h"
#include "TStyle.h"
#include "TLine.h"
#include "TBrowser.h"
#include "TApplication.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TCutG.h"
#include "TString.h"
#include "TObjArray.h"
#include "TPaveLabel.h"
#include "TGraph.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "Riostream.h"
#include <TSystemDirectory.h>

using namespace std;

void process_ToF() {

	gStyle->SetOptStat(0);


	/////////////////////////////Definimos variables///////////////////////////
	char hname[100];

	const Double_t flightPath = 0.5;        // metros
	const Double_t tofFactor = 72.2977*1e3; // ev^1/2*ns/m
	Float_t neutronEnergy;
	Float_t neutronBckg;

	Double_t TotalCurrent;

	const int nDetectors = 3;
	const Float_t detWidth[] = {12.7,12.7,3}; // mm
	Float_t L_eff[nDetectors];
	for (int i=0; i<nDetectors; i++) {
		L_eff[i] = flightPath + detWidth[i]*1.E-3/2;
	}

	Float_t Li2O_dens = 2.6;        // g/cm3
	Float_t Li2O_masaMol = 28.029;  // g/mol
	Float_t nAvogadro = 6.02214e23; // mol
	Float_t thickness[nDetectors];
	for (int i=0; i<nDetectors; i++) {
		thickness[i] = 2*0.18*0.96*(Li2O_dens/Li2O_masaMol)*detWidth[i]*0.1*(nAvogadro/1.e24);  // at/b
	}


	const int ampCutLow[3] = {1040,1420,1130};
	const int ampCutHigh[3] = {1630,2640,1900}; //1912keV_pos3*/
	/*const int ampCutLow[3] = {,,};
	const int ampCutHigh[3] = {,,}; //2000keV_pos3*/



	const float minNeutronTof[3] = {-205.,-195.,-205.}; //1912keV_pos3
	//const float minNeutronTof[3] = {,,}; //2000keV_pos3


	///////////////////////////////Definimos histogramas///////////////////////////////
	//Definimos el axis de nuestros histogramas (BPD)
    int BPD = 500;
    double lowBin = 3.E2;
    double highBin = 1*1E6;
    int nOfBins = BPD*(log10(highBin)-log10(lowBin));
    double nBPDaxis[nOfBins+1]; 

    for (double i=0.; i<=BPD*(log10(highBin) - log10(lowBin)); i++) {
    	nBPDaxis[(int)i] = pow(10.,i/BPD + log10(lowBin));
    }

    const int nBinsEn = 1500;
    const float maxEn = 150.;  //1912keV*/
    /*const int nBinsEn = 2600;
    const float maxEn = 260.;  //2000keV*/



	TH1F* h_amplitude[nDetectors];
	TH1F* h_tof[nDetectors];
	TH1F* h_neutronEnergy[nDetectors];
	TH1F* h_neutronEnergyWeighted[nDetectors];

	TH1F* h_neutronEnergyNoBckg[nDetectors];
	TH1F* h_neutronEnergyWeightedNoBckg[nDetectors];

	TH1F* h_neutronBckg[nDetectors];

	TH1F* h_amplitude_timeCuts[nDetectors];
	TH1F* h_tof_ampCuts[nDetectors];

	for (int i=0; i<nDetectors; i++) {
		sprintf(hname,"h_amplitude_Li-Glass%d",i+1);
		h_amplitude[i] = new TH1F(hname,"",4096,0,4096);

		sprintf(hname,"h_tof_Li-Glass%d",i+1);
		h_tof[i] = new TH1F(hname,"",2000,-1000,1000);

		sprintf(hname,"h_amplitude_timeCuts_Li-Glass%d",i+1);
		h_amplitude_timeCuts[i] = new TH1F(hname,"",4096,0,4096);

		sprintf(hname,"h_tof_ampCuts_Li-Glass%d",i+1);
		h_tof_ampCuts[i] = new TH1F(hname,"",2000,-1000,1000);


		sprintf(hname,"h_neutronEnergy_Li-Glass%d",i+1);
		h_neutronEnergy[i] = new TH1F(hname,"",nBinsEn,0,maxEn);

		sprintf(hname,"h_neutronEnergyWeighted_Li-Glass%d",i+1);
		h_neutronEnergyWeighted[i] = new TH1F(hname,"",nBinsEn,0,maxEn);

		sprintf(hname,"h_neutronEnergyNoBckg_Li-Glass%d",i+1);
		h_neutronEnergyNoBckg[i] = new TH1F(hname,"",nBinsEn,0,maxEn);

		sprintf(hname,"h_neutronEnergyWeightedNoBckg_Li-Glass%d",i+1);
		h_neutronEnergyWeightedNoBckg[i] = new TH1F(hname,"",nBinsEn,0,maxEn);

		sprintf(hname,"h_neutronBckg_Li-Glass%d",i+1);
		h_neutronBckg[i] = new TH1F(hname,"",nBinsEn,0,maxEn);

	}

	TH1F* h_current = new TH1F("h_current","",5400,0,5400);
	TH1F* h_ampCurrent = new TH1F("h_ampCurrent","",5400,0,5400);

	//Li-Glass1 -> Channel 2
	//Li-Glass2 -> Channel 3
	//Li-Glass3 -> Channel 4

	//////////////////////////Importamos XS de 6Li(n,a)3H/////////////////////////////////
	TGraph* graph6Li_na = new TGraph();
	double x,y;
	int xsIndex=0;
	std::ifstream inputfile_6Li("6Li_n_alfa_xs.dat",std::ios::in);
	while(inputfile_6Li >> x >> y) {
		graph6Li_na->SetPoint(xsIndex,x*1e3,y*thickness[1]);
		xsIndex++;
	}
	inputfile_6Li.close();


	//////////////////////////Importamos archivos del DAQ/////////////////////////////////
	TString sPathName = "ToF_PSD/DAQ/BeamON_2000keV_pos3/";
	TString sFileName = "SData_BeamON_2000keV_pos3.root";
	TString sName = sPathName+sFileName;
	TString sTreeName = "Data";
	TChain* t = new TChain(sTreeName);

	cout << "\n";
	cout << "Leyendo archivo -> " << sFileName << " ... " << endl;
	cout << "\n";

	UShort_t nchannel;
	ULong64_t time;
	UShort_t Elong;
	UShort_t Eshort; 
	Double_t tof;

	t->SetBranchAddress("Channel", &nchannel);  
	t->SetBranchAddress("Timestamp", &time);     
	t->SetBranchAddress("Energy", &Elong);       
	t->SetBranchAddress("EnergyShort", &Eshort);       
	t->SetBranchAddress("tof", &tof);
	t->Add(sName.Data());

	for (long long iEntry=0; iEntry < t->GetEntries(); iEntry++) {
		t->GetEntry(iEntry);

		if (nchannel==1) {
			h_current->Fill(time*1.E-12);
			h_ampCurrent->Fill(Elong);
		}

		if (nchannel>1 && nchannel<5) {
			h_amplitude[nchannel-2]->Fill(Elong);
			h_tof[nchannel-2]->Fill(tof);

			if (Elong > ampCutLow[nchannel-2] && Elong < ampCutHigh[nchannel-2]) {
				h_tof_ampCuts[nchannel-2]->Fill(tof);
			}
		}
	}

	///////////////////////////Fiteamos el g-flash/////////////////////////////////////////
	cout << "\n";
	cout << "Fiteando gamma-flash ... " << endl;
	cout << "\n";

	TF1* fit_gFlash[nDetectors];
	float tflash[3];


	int fit_gFlash_low[] = {-230,-220,-225};
	int fit_gFlash_high[] = {-205,-200,-205};  //1912keV_pos3*/

	/*int fit_gFlash_low[] = {,,};
	int fit_gFlash_high[] = {,,};  //2000keV_pos3*/



	for (int i=0; i<nDetectors; i++) {
		sprintf(hname,"Fit%d",i+1);
		fit_gFlash[i] = new TF1(hname,"[0]*exp(-pow(x-[1],2)/(2*pow([2],2)))",fit_gFlash_low[i],fit_gFlash_high[i]);
		//fit_gFlash[i]->SetParLimits(0,100,1000);  //1912keV_pos2&rest
		fit_gFlash[i]->SetParLimits(1,fit_gFlash_low[i],fit_gFlash_high[i]);
		fit_gFlash[i]->SetParLimits(2,1,4);
		h_tof[i]->Fit(fit_gFlash[i],"R");
		tflash[i] = fit_gFlash[i]->GetParameter(1);
	}



	float binEdgeTof[1001];
	float binEdgeEn[nBinsEn+1];
	float corrFactor[nBinsEn];
	float energyAux;

	///////////////////////////Fiteamos el fondo de neutrones de la sala///////////////////////
	cout << "\n";
	cout << "Fiteando fondo de neutrones ... " << endl;
	cout << "\n";

	TF1* fit_bckg[nDetectors];
	for (int i=0; i<nDetectors; i++) {
		sprintf(hname,"FitBckg%d",i+1);
		fit_bckg[i] = new TF1(hname,"pol0(0)",-200,-140); //1912keV
		//fit_bckg[i] = new TF1(hname,"pol0(0)",,); //2000keV

		h_tof_ampCuts[i]->Fit(fit_bckg[i],"R");
	}


	////////////////////////////Calculamos energía de los neutrones///////////////////////////
	cout << "\n";
	cout << "Calculando energía de los neutrones ... " << endl;
	cout << "\n";

	for (long long iEntry=0; iEntry < t->GetEntries(); iEntry++) {
		t->GetEntry(iEntry);

		if (nchannel>1 && nchannel<5) {

			if (Elong > ampCutLow[nchannel-2] && Elong < ampCutHigh[nchannel-2] && tof > minNeutronTof[nchannel-2] && tof < 1000.) {
				neutronEnergy = tofFactor*tofFactor*L_eff[nchannel-2]*L_eff[nchannel-2]/(tof-tflash[nchannel-2])/(tof-tflash[nchannel-2]);
				h_neutronEnergy[nchannel-2]->Fill(neutronEnergy/1000.);
				h_neutronEnergyWeighted[nchannel-2]->Fill(neutronEnergy/1000.,TMath::Sqrt(neutronEnergy/1000.));
				//h_neutronEnergyWeighted[nchannel-2]->Fill( neutronEnergy/1000.,1./graph6Li_na->Eval(neutronEnergy/1000.) );

			}
			else if (Elong > ampCutLow[nchannel-2] && Elong < ampCutHigh[nchannel-2] && tof >-1000 && tof < tflash[nchannel-2]){
				neutronEnergy = tofFactor*tofFactor*L_eff[nchannel-2]*L_eff[nchannel-2]/((tof+2000.)-tflash[nchannel-2])/((tof+2000.)-tflash[nchannel-2]);
				h_neutronEnergy[nchannel-2]->Fill(neutronEnergy/1000.);
				h_neutronEnergyWeighted[nchannel-2]->Fill(neutronEnergy/1000.,TMath::Sqrt(neutronEnergy/1000.));
				//h_neutronEnergyWeighted[nchannel-2]->Fill( neutronEnergy/1000.,1./graph6Li_na->Eval(neutronEnergy/1000.) );

			}
		}
	}



	////////////////////////////Eliminamos el fondo de neutrones///////////////////////////
	for (int i=0; i<nDetectors; i++) {
		for (int j=1; j<nBinsEn; j++) {
			energyAux = tofFactor*L_eff[i]/TMath::Sqrt(h_neutronEnergy[i]->GetBinLowEdge(j)*1000. ) + tflash[i] - (tofFactor*L_eff[i]/TMath::Sqrt(h_neutronEnergy[i]->GetBinLowEdge(j+1)*1000. ) + tflash[i]);
			//cout << energyAux << endl;
			if (energyAux > 1e7) continue;
			h_neutronEnergyNoBckg[i]->SetBinContent(j, h_neutronEnergy[i]->GetBinContent(j) - fit_bckg[i]->GetParameter(0)*energyAux );
			h_neutronEnergyWeightedNoBckg[i]->SetBinContent(j, h_neutronEnergyNoBckg[i]->GetBinContent(j)*TMath::Sqrt(h_neutronEnergyNoBckg[i]->GetBinCenter(j)) );
			//h_neutronEnergyWeightedNoBckg[i]->SetBinContent(j, h_neutronEnergyNoBckg[i]->GetBinContent(j)/graph6Li_na->Eval(h_neutronEnergyNoBckg[i]->GetBinCenter(j)) );
			h_neutronBckg[i]->SetBinContent(j, fit_bckg[i]->GetParameter(0)*energyAux);
		}
	}


	///////////////////////////Escalamos a la corriente///////////////////////////////////////
	cout << "\n";
	cout << "Escalamos histogramas ... " << endl;
	cout << "\n";

	h_current->Scale(0.0001);  //1e-10 C/pulso --> 1e-4 uC/pulso
	h_ampCurrent->Scale(0.0001);
	TotalCurrent = h_ampCurrent->Integral();
	cout << "Total current -> " << TotalCurrent << " uC" << endl;


	/*for (int i=0; i<nDetectors; i++) {
		h_amplitude[i]->Scale(1./TotalCurrent);
		h_tof[i]->Scale(1./TotalCurrent);
		h_tof_ampCuts[i]->Scale(1./TotalCurrent);
		h_neutronEnergy[i]->Scale(1./TotalCurrent);
		h_neutronEnergyWeighted[i]->Scale(1./TotalCurrent);
		h_neutronEnergyNoBckg[i]->Scale(1./TotalCurrent);
		h_neutronEnergyWeightedNoBckg[i]->Scale(1./TotalCurrent);
		if (i==2) {
			h_neutronEnergy[i]->Scale(12.7/3.);
			h_neutronEnergyWeighted[i]->Scale(12.7/3.);
			h_neutronEnergyNoBckg[i]->Scale(12.7/3.);
			h_neutronEnergyWeightedNoBckg[i]->Scale(12.7/3.);
		}
	}*/




	///////////////////////////Pintamos los histogramas/////////////////////////////////////////
	TCanvas* cXS = new TCanvas("cXS","cXS");
	cXS->SetGridx(); cXS->SetGridy();
	cXS->SetLogx(); cXS->SetLogy();
	graph6Li_na->GetXaxis()->SetTitle("Neutron energy (keV)");
	graph6Li_na->GetYaxis()->SetTitle("Cross section (barns)");
	graph6Li_na->Draw("AL");

	
	TCanvas* c_amp = new TCanvas("c_amp");
	c_amp->Divide(1,3);
	for (int i=0; i<nDetectors; i++) {
		c_amp->cd(i+1); gPad->SetGridx(); gPad->SetGridy(); //gPad->SetLogy();
		sprintf(hname,"Li-Glass%d",i+1);
		h_amplitude[i]->SetTitle(hname);
		h_amplitude[i]->GetXaxis()->SetTitle("Amplitude");
		h_amplitude[i]->GetYaxis()->SetTitle("Counts/uC");
		h_amplitude[i]->Draw("HIST");
	}

	TCanvas* c_tof = new TCanvas("c_tof","c_tof");
	c_tof->Divide(1,3);
	for (int i=0; i<nDetectors; i++) {
		c_tof->cd(i+1); gPad->SetGridx(); gPad->SetGridy(); //gPad->SetLogy();
		sprintf(hname,"Li-Glass%d",i+1);
		h_tof[i]->SetTitle(hname);
		h_tof[i]->GetXaxis()->SetTitle("Time (ns)");
		h_tof[i]->GetYaxis()->SetTitle("Counts/uC");
		h_tof[i]->Draw("HIST");
		fit_gFlash[i]->Draw("SAME");
	}

	TCanvas* c_tof_cuts = new TCanvas("c_tof_cuts","c_tof_cuts");
	c_tof_cuts->Divide(1,3);
	for (int i=0; i<nDetectors; i++) {
		c_tof_cuts->cd(i+1); gPad->SetGridx(); gPad->SetGridy(); //gPad->SetLogy();
		sprintf(hname,"Li-Glass%d, amp -> (%d,%d)",i+1,ampCutLow[i],ampCutHigh[i]);
		h_tof_ampCuts[i]->SetTitle(hname);
		h_tof_ampCuts[i]->GetXaxis()->SetTitle("Time (ns)");
		h_tof_ampCuts[i]->GetYaxis()->SetTitle("Counts/uC");
		h_tof_ampCuts[i]->Draw("HIST");
		fit_bckg[i]->Draw("SAME");
	}

	TCanvas* c_nE = new TCanvas("c_nE","c_nE");
	c_nE->Divide(1,3);
	for (int i=0; i<nDetectors; i++) {
		c_nE->cd(i+1); gPad->SetGridx(); gPad->SetGridy(); //gPad->SetLogy();
		sprintf(hname,"Li-Glass%d",i+1);
		h_neutronEnergy[i]->SetTitle(hname);
		h_neutronEnergy[i]->GetXaxis()->SetTitle("Neutron energy (keV)");
		h_neutronEnergy[i]->GetYaxis()->SetTitle("Counts/uC");
		h_neutronEnergy[i]->SetMarkerStyle(20);
		h_neutronEnergy[i]->Draw("histo");
	}

	TCanvas* c_nEW = new TCanvas("c_nEW","c_nEW");
	c_nEW->Divide(1,3);
	for (int i=0; i<nDetectors; i++) {
		c_nEW->cd(i+1); gPad->SetGridx(); gPad->SetGridy(); //gPad->SetLogy();
		sprintf(hname,"Li-Glass%d",i+1);
		h_neutronEnergyWeighted[i]->SetTitle(hname);
		h_neutronEnergyWeighted[i]->GetXaxis()->SetTitle("Neutron energy (keV)");
		h_neutronEnergyWeighted[i]->GetYaxis()->SetTitle("Counts/sqrt(E)/uC (scaled)");
		h_neutronEnergyWeighted[i]->SetMarkerStyle(20);
		h_neutronEnergyWeighted[i]->Draw("histo");
	}

	TCanvas* c_nEB = new TCanvas("c_nEB","c_nEB");
	c_nEB->Divide(1,3);
	for (int i=0; i<nDetectors; i++) {
		c_nEB->cd(i+1); gPad->SetGridx(); gPad->SetGridy(); //gPad->SetLogy();
		sprintf(hname,"Li-Glass%d",i+1);
		h_neutronEnergyNoBckg[i]->SetTitle(hname);
		h_neutronEnergyNoBckg[i]->GetXaxis()->SetTitle("Neutron energy (keV)");
		h_neutronEnergyNoBckg[i]->GetYaxis()->SetTitle("Counts/uC (noBckg)");
		h_neutronEnergyNoBckg[i]->SetMarkerStyle(20);
		h_neutronEnergyNoBckg[i]->SetMinimum(0);
		h_neutronEnergyNoBckg[i]->Draw("histo");
	}

	TCanvas* c_nEWB = new TCanvas("c_nEWB","c_nEWB");
	c_nEWB->Divide(1,3);
	for (int i=0; i<nDetectors; i++) {
		c_nEWB->cd(i+1); gPad->SetGridx(); gPad->SetGridy(); //gPad->SetLogy();
		sprintf(hname,"Li-Glass%d",i+1);
		h_neutronEnergyWeightedNoBckg[i]->SetTitle(hname);
		h_neutronEnergyWeightedNoBckg[i]->GetXaxis()->SetTitle("Neutron energy (keV)");
		h_neutronEnergyWeightedNoBckg[i]->GetYaxis()->SetTitle("Counts/sqrt(E)/uC (scaled&noBckg)");
		h_neutronEnergyWeightedNoBckg[i]->SetMarkerStyle(20);
		h_neutronEnergyWeightedNoBckg[i]->SetMinimum(0);
		h_neutronEnergyWeightedNoBckg[i]->Draw("histo");
	}


	TCanvas* c_B = new TCanvas("c_B","c_B");
	c_B->Divide(1,3);
	for (int i=0; i<nDetectors; i++) {
		c_B->cd(i+1); gPad->SetGridx(); gPad->SetGridy(); gPad->SetLogy();
		sprintf(hname,"Li-Glass%d",i+1);
		h_neutronBckg[i]->SetTitle(hname);
		h_neutronBckg[i]->GetXaxis()->SetTitle("Neutron energy (keV)");
		h_neutronBckg[i]->GetYaxis()->SetTitle("Background");
		h_neutronBckg[i]->SetMarkerStyle(20);
		h_neutronBckg[i]->Draw("histo");
	}



	TCanvas* c_CI = new TCanvas("c_CI","c_CI");
	c_CI->cd();
	h_current->GetXaxis()->SetTitle("Time (s)");
	h_current->GetYaxis()->SetTitle("Current (uA)");
	h_current->Draw("hist");

	TCanvas* c_ampCI = new TCanvas("c_ampCI","c_ampCI");
	c_ampCI->cd();
	h_ampCurrent->GetXaxis()->SetTitle("Amplitude");
	h_ampCurrent->GetYaxis()->SetTitle("Counts");
	h_ampCurrent->Draw("hist");


	//////////////////////////////////Guardamos los resultados////////////////////////////////////////
	cout << "\n";
	cout << "Guardamos los resultados ... " << endl;
	cout << "\n";

	system("mkdir Processed");

	TString fOut_path = "Processed/OUT_scaled_BckgMoved_";
	TFile* fOut = new TFile(fOut_path+sFileName, "RECREATE");
	for (int i=0; i<nDetectors; i++) {
		sprintf(hname,"amp_scaled_LiGlass%d",i+1);
		h_amplitude[i]->Write(hname);

		sprintf(hname,"tof_scaled_LiGlass%d",i+1);
		h_tof[i]->Write(hname);

		sprintf(hname,"tof_scaled_cuts_LiGlass%d",i+1);
		h_tof_ampCuts[i]->Write(hname);

		sprintf(hname,"neutronEnergy_scaled_LiGlass%d",i+1);
		h_neutronEnergy[i]->Write(hname);

		sprintf(hname,"neutronEnergyWeighted_scaled_LiGlass%d",i+1);
		h_neutronEnergyWeighted[i]->Write(hname);

		sprintf(hname,"neutronEnergyNoBckg_scaled_LiGlass%d",i+1);
		h_neutronEnergyNoBckg[i]->Write(hname);

		sprintf(hname,"neutronEnergyWeightedNoBckg_scaled_LiGlass%d",i+1);
		h_neutronEnergyWeightedNoBckg[i]->Write(hname);

		sprintf(hname,"neutronBckg_LiGlass%d",i+1);
		h_neutronBckg[i]->Write(hname);

	}

	h_current->Write("h_current");
	h_ampCurrent->Write("h_ampCurrent");

	fOut->Close();
}
