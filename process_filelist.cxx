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

//Necesarias
#include "ROOT/RDataFrame.hxx"

using namespace std;
using namespace ROOT;

#include "create_output_file.cxx"
#include "activation.cxx"


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
void add_histograms(Char_t filepath[500]);

//process all files inside "input" folder with relative filepaths in "filelist" file
void process_filelist(){
	Char_t filepath[500];
	Char_t new_filepath[500];

	ifstream read("filelist");
	
	signed char tof_flag = 1;
	signed char histo_flag = 1;
	char answer;
	cout << "Add ToF and PSD? (y/n): ";
	cin >> answer;	//TBD:unsafe?
	if(answer=='n'||answer=='N'){
		tof_flag = -1;
	}else if(answer!='y'&&answer!='Y'){
		cout << "Invalid answer. Returning." << endl;
		return;
	}
	cout << "Add histograms? (y/n): ";
	cin >> answer;	//TBD:unsafe?
	if(answer=='n'||answer=='N'){
		histo_flag = -1;
	}else if(answer!='y'&&answer!='Y'){
		cout << "Invalid answer. Returning." << endl;
		return;
	}


	while (!read.eof()){
		read >>filepath;
		if(tof_flag==1){
			createToF(filepath);
		}
		sprintf(new_filepath,"output/%s", filepath);
		if(histo_flag==1){
			add_histograms(new_filepath);
		}
		cout << endl;
	}

	create_output_file();
	activation();

	cout << "That's all folks!" <<endl;
	return;
}

//add histograms to root file. Assumes tree with ToF and PSD
void add_histograms(Char_t filepath[500]){
	EnableImplicitMT();	//multithreading
	RDataFrame d("Data", filepath);

	//Histograms
	auto tadeo_1_time = d.Filter("Channel==2 && Energy>0").Histo1D({"tadeo_1_time", "; Timestamp; Counts", 500, 0, 2500E12}, "Timestamp");
	auto tadeo_1_spectrum = d.Filter("Channel==2 && Energy>0").Histo1D({"tadeo_1_spectrum", "; Energy; Counts", 4096, 0, 4096}, "Energy");
	auto tadeo_2_time = d.Filter("Channel==3 && Energy>0").Histo1D({"tadeo_2_time", "; Timestamp; Counts", 500, 0, 2500E12}, "Timestamp");
	auto tadeo_2_spectrum = d.Filter("Channel==3 && Energy>0").Histo1D({"tadeo_2_spectrum", "; Energy; Counts", 4096, 0, 4096}, "Energy");
	auto monster_2_time = d.Filter("Channel==4 && Energy>0").Histo1D({"monster_time", "; Timestamp; Counts", 500, 0, 2500E12}, "Timestamp");
	auto monster_2_spectrum = d.Filter("Channel==4 && Energy>0").Histo1D({"monster_spectrum", "; Energy; Counts", 4096, 0, 4096}, "Energy");
	auto stylbeno_2_time = d.Filter("Channel==5 && Energy>0").Histo1D({"stylbeno_time", "; Timestamp; Counts", 500, 0, 2500E12}, "Timestamp");
	auto stylbeno_2_spectrum = d.Filter("Channel==5 && Energy>0").Histo1D({"stylbeno_spectrum", "; Energy; Counts", 4096, 0, 4096}, "Energy");
	auto labr_1_time = d.Filter("Channel==6 && Energy>0").Histo1D({"labr_1_time", "; Timestamp; Counts", 500, 0, 2500E12}, "Timestamp");
	auto labr_1_spectrum = d.Filter("Channel==6 && Energy>0").Histo1D({"labr_1_spectrum", "; Energy; Counts", 4096, 0, 4096}, "Energy");
	auto labr_2_time = d.Filter("Channel==7 && Energy>0").Histo1D({"labr_2_time", "; Timestamp; Counts", 500, 0, 2500E12}, "Timestamp");
	auto labr_2_spectrum = d.Filter("Channel==7 && Energy>0").Histo1D({"labr_2_spectrum", "; Energy; Counts", 4096, 0, 4096}, "Energy");

	TFile f(filepath, "UPDATE");
	tadeo_1_time->Write();
	tadeo_1_spectrum->Write();
	tadeo_2_time->Write();
	tadeo_2_spectrum->Write();
	monster_2_time->Write();
	monster_2_spectrum->Write();
	stylbeno_2_time->Write();
	stylbeno_2_spectrum->Write();
	labr_1_time->Write();
	labr_1_spectrum->Write();
	labr_2_time->Write();
	labr_2_spectrum->Write();
	f.Close();

	DisableImplicitMT();	//multithreading
	cout << "Created histograms." <<endl;
}

//Takes the original ROOT file from CAEN digitizer and crate a new one with ToF and PSD branches named "ToF_*"
void createToF(Char_t filename[500]){
	gStyle->SetOptStat(0);
	int i;
	Char_t commandline[1000];
	Char_t newfilename[500];

	sprintf(commandline,"install -Dv input/%s output/%s", filename, filename);
	system(commandline);
	sprintf(newfilename,"output/%s", filename);
	
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
