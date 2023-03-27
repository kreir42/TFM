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
#include <string.h>

// ROOT libraries
#include "TROOT.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// GLOBAL PARAMETERS AND ROUTINES ////////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#define channelsused 8
#define maxTOF 1100
#define minTOF -100
#define maxIntegral 4096
#define binsIntegral 512
#define maxNEne 250
#define binsNEne 250
#define binsPSD 100
#define CoincWindow 1000 // ns

int number_files = 0;

void process_filelist();
void createToF(Char_t filename[500]); //The original ROOT file from CAEN digitizer and crate a new one names "*_ToF"

/************************************************************************************************************
************************************************************************************************************* 
 Runs the createToF routine over a entire project folder                                                     
*************************************************************************************************************
************************************************************************************************************/

void process_filelist(){
	Char_t filepath[500];

	ifstream read("filelist");
	
	while (!read.eof()){
		read >>filepath;
		createToF(filepath);
	}

	cout << "That's all folks!" <<endl;
	return;
}

/**************************************************************************************************************
*************************************************************************************************************** 
 Takes the original ROOT file from CAEN digitizer and crate a new one with ToF and PSD branches named "ToF_*"
***************************************************************************************************************
**************************************************************************************************************/

void createToF(Char_t filename[500]){
	gStyle->SetOptStat(0);
	int i;
	Char_t commandline[1000];
	Char_t newfilename[500];

	sprintf(commandline,"install -Dv %s with_ToF_PSD/%s", filename, filename);
	system(commandline);
	sprintf(newfilename,"with_ToF_PSD/%s", filename);
	
	TFile* f = TFile::Open(newfilename,"Update");


	TTree *Data = (TTree*)f->Get("Data");
	UShort_t nchannel;                                         
	Data->SetBranchAddress("Channel", &nchannel);  
	ULong64_t time;                                          
	Data->SetBranchAddress("Timestamp", &time);     
	UShort_t Elong;                                          
	Data->SetBranchAddress("Energy", &Elong);       
	UShort_t Eshort;                                          
	Data->SetBranchAddress("EnergyShort", &Eshort);       
    Double_t tof;
    TBranch *BTOF = Data->Branch("tof", &tof, "tof/D");
	Double_t psd;
    TBranch *BPSD = Data->Branch("psd", &psd, "psd/D");

//	TH1F *htof[channelsused];
//	TH2F *hamplitudetof[channelsused];
//	Char_t hname[20];
//	for(i=0; i<channelsused; i++)
//	{
//		sprintf(hname,"htof%d",i);
//		htof[i] = new TH1F(hname,";time-of-flight (ns); counts",6000,0,6000.);
//	 	sprintf(hname,"hamplitudetof%d",i);
//	 	hamplitudetof[i]= new TH2F(hname,";time-of-flight (ns); Signal integral (channels);  counts",6000,0,6000.,4192,0,4192);		
// 	}

	int nentries = Data->GetEntries();
	cout<< "Numero señales: "<<nentries<< endl;
	double tref, signaltime, newtof;

	for (i=0; i<nentries; i++) // FOR EACH SIGNAL IN THE RUN	
	{ 	
		Data->GetEntry(i);
		
		if(Elong>0)
			psd=(Double_t) (Elong-Eshort)/Elong;
		else
			psd=0;

		if(nchannel==0)
		{
			//ht0->Fill(time/1.e9);
			tref=time;
			tof=0;
		}
		else
		{
			tof=(Double_t) (time-tref)/1000.;
			if(tof>CoincWindow) // try with next event
			{
				signaltime=time;
				i++;
				Data->GetEntry(i);
				i--;
				if(nchannel==0)
					tof=(Double_t)(signaltime-time)/1000.;
			}
			//cout<<"Filling tof="<<tof<<"---------------"<<endl<<endl;
			//if(abs(tof)<CoincWindow){getchar();}
//			htof[nchannel]->Fill(tof);
//			hamplitudetof[nchannel]->Fill(tof, Elong);
		}
		BTOF->Fill();	
		BPSD->Fill();
	}
	
//	TCanvas *ctof = new TCanvas("ctof","", 1600,900);
//	ctof->Divide(4,2);
//	TCanvas *camplitudetof = new TCanvas("camplitudetof","", 1600,900);
//	camplitudetof->Divide(4,2);
	
//	for(i=0; i<channelsused; i++)
//	{	
//		ctof->cd(i+1);
//		gPad->SetTicks();
//		htof[i]->DrawCopy();
//		camplitudetof->cd(i+1);
//		gPad->SetTicks();
//		gPad->SetLogy(1);
//		hamplitudetof[i]->DrawCopy("zcol");
//	}


	Data->Write("", TObject::kOverwrite);
	
	
	f->Close();
	cout <<"New file with TOF: "<< endl<<newfilename << endl;
	return; 
}
