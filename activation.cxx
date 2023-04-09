#define ACTIVATION_NBINS 1000

static void per_file(Char_t filepath[500], Double_t results[2][4]);

void activation(){
	cout << "Activation" << endl;
	Char_t filepath_1[100] = "output/SData_aAl_J78kV_GVM1808kV_positions2_activacion.root";
	Char_t filepath_2[100] = "output/SData_aAl_J78kV_GVM2312kV_positions2_activacion.root";
	Char_t filepath_3[100] = "output/SData_aAl_J78kV_GVM2810kV_positions2_activacion.root";
	Char_t filepath_4[100] = "output/SData_aAl_J78kV_GVM2810kV_positions2_activacion_20230223.root";

	TFile f("output/output.root", "UPDATE");
	gDirectory->cd("Activation");
	Double_t results[4][2][4];
	Double_t activation_energies[] = {5500, 7000, 8500, 8500};	//keV

	gDirectory->cd("activation_1");
	cout << "activation_1" << endl;
	per_file(filepath_1, results[0]);
	gDirectory->cd("..");

	gDirectory->cd("activation_2");
	cout << "activation_2" << endl;
	per_file(filepath_2, results[1]);
	gDirectory->cd("..");

	gDirectory->cd("activation_3");
	cout << "activation_3" << endl;
	per_file(filepath_3, results[2]);
	gDirectory->cd("..");

	gDirectory->cd("activation_4");
	cout << "activation_4" << endl;
	per_file(filepath_4, results[3]);
	gDirectory->cd("..");

	//graficas
	Double_t x[8];
	Double_t y[8];
	Double_t yerr[8];
	for(short i=0; i<4; i++){
		x[2*i] = activation_energies[i];
		x[2*i+1] = activation_energies[i];
		y[2*i] = results[i][0][2];
		y[2*i+1] = results[i][1][2];
		yerr[2*i] = results[i][0][3];
		yerr[2*i+1] = results[i][1][3];
	}
	TCanvas* myCanvas = new TCanvas("reactions_v_energy");
	TGraph* rectionsvenergy = new TGraphErrors(8, x, y, NULL, yerr);
	rectionsvenergy->SetTitle("(a,n) reactions v a energy;Energy of a (keV);Inferred (a,n)/Number of a");
	rectionsvenergy->SetMarkerStyle(20);
	rectionsvenergy->Draw("ap");
	myCanvas->Write();

	myCanvas->Close();
	f.Close();
}

class unified_fit{
	public:
		TH1F* current_histogram;
		Double_t stepsize;
		unified_fit(TH1F* histo):current_histogram(histo){
			stepsize = current_histogram->GetBinWidth(1);
		}
		Double_t operator()(Double_t* x, Double_t* p){
			Double_t n = 0;	//number of nuclei
			Double_t decay = exp(-p[2]*stepsize);
			Double_t helper_ratio = p[1]/(p[2]*stepsize);
			Double_t I = 0;	//created nuclei per time over lambda
			unsigned long steps = x[0]/stepsize;
			for(unsigned long i=0; i<=steps; i++){
				//n = n*decay + (current_histogram->GetBinContent(i))*p[1];
				I = (current_histogram->GetBinContent(i))*helper_ratio;
				n = I + (n-I)*decay;
			}
			return p[0] + n*(1-decay);
		}
};

static void per_file(Char_t filepath[500], Double_t results[2][4]){
	EnableImplicitMT();	//multithreading
	RDataFrame d("Data", filepath);

	auto integrator_signals = d.Filter("Channel==1");
	Double_t activation_start = integrator_signals.Min("Timestamp").GetValue();	//TBD!:muy ineficiente!!
	Double_t activation_end = integrator_signals.Max("Timestamp").GetValue();
	Double_t measurement_end = d.Filter("(Channel==6||Channel==7) && Energy>0").Max("Timestamp").GetValue();
	ULong64_t number_of_alphas = integrator_signals.Count().GetValue()/(2*1.60217646E-10);

	//histogramas
	auto rise_filter = [&](ULong64_t Timestamp){return Timestamp>=activation_start && Timestamp<=activation_end;};
	auto decay_filter = [&](ULong64_t Timestamp){return Timestamp>activation_end;};
	auto labr_1_filter = d.Filter("Channel==6 && Energy>525 && Energy<650");
	auto labr_2_filter = d.Filter("Channel==7 && Energy>525 && Energy<650");

	auto current_integrator = integrator_signals.Histo1D({"current_integrator", "; Timestamp; Counts", ACTIVATION_NBINS, 0, measurement_end}, "Timestamp");
	auto labr_1 = labr_1_filter.Histo1D({"labr_1", "; Timestamp; Counts", ACTIVATION_NBINS, 0, measurement_end}, "Timestamp");
	auto labr_2 = labr_2_filter.Histo1D({"labr_2", "; Timestamp; Counts", ACTIVATION_NBINS, 0, measurement_end}, "Timestamp");

	auto labr_1_rise = labr_1_filter.Filter(rise_filter, {"Timestamp"}).Histo1D({"labr_1_rise", "; Timestamp; Counts", ACTIVATION_NBINS, activation_start, activation_end}, "Timestamp");
	auto labr_1_decay = labr_1_filter.Filter(decay_filter, {"Timestamp"}).Histo1D({"labr_1_decay", "; Timestamp; Counts", ACTIVATION_NBINS, activation_end, measurement_end}, "Timestamp");
	auto labr_2_rise = labr_2_filter.Filter(rise_filter, {"Timestamp"}).Histo1D({"labr_2_rise", "; Timestamp; Counts", ACTIVATION_NBINS, activation_start, activation_end}, "Timestamp");
	auto labr_2_decay = labr_2_filter.Filter(decay_filter, {"Timestamp"}).Histo1D({"labr_2_decay", "; Timestamp; Counts", ACTIVATION_NBINS, activation_end, measurement_end}, "Timestamp");

	cout << "Rise/decay histograms" << endl;
	current_integrator->Write();
	labr_1->Write();
	labr_2->Write();
	labr_1_rise->Write();
	labr_1_decay->Write();
	labr_2_rise->Write();
	labr_2_decay->Write();

	//fittings
	TCanvas* myCanvas = new TCanvas("myCanvas");
	TFitResultPtr fitresult;

	//unified
	cout << "Unified rise/decay fittings" << endl;
	TF1* unified = new TF1("unified_fit", unified_fit((TH1F*)gDirectory->Get("current_integrator")), 0, measurement_end, 3);
	unified->SetNpx(ACTIVATION_NBINS);
	unified->SetNumberFitPoints(ACTIVATION_NBINS);
	unified->SetParLimits(0, 0, 1E5);
	unified->SetParLimits(1, 0, 1E10);
	unified->SetParLimits(2, 4E-15, 5E-15);
	unified->SetParameters(15, 1E3, 4.62406E-15);
	unified->FixParameter(2, 4.62406E-15);
	unified->SetParNames("Background activity", "current to (a,n)", "Decay constant");

	fitresult = labr_1->Fit("unified_fit", "SLN");
	labr_1->Draw();
	unified->Draw("CSAME");
	myCanvas->SetName("labr_1_unified_fit");
	myCanvas->Write();

	fitresult = labr_2->Fit("unified_fit", "SLN");
	labr_2->Draw();
	unified->Draw("CSAME");
	myCanvas->SetName("labr_2_unified_fit");
	myCanvas->Write();

	//rise
	cout << "Rise fittings" << endl;

	fitresult = labr_1_rise->Fit("unified_fit", "SL");
	results[0][2] = fitresult->Parameter(1)*(activation_end-activation_start)/number_of_alphas;
	results[0][3] = fitresult->ParError(1)*(activation_end-activation_start)/number_of_alphas;
	myCanvas->SetName("labr_1_rise");
	myCanvas->Write();

	fitresult = labr_2_rise->Fit("unified_fit", "SL");
	results[1][2] = fitresult->Parameter(1)*(activation_end-activation_start)/number_of_alphas;
	results[1][3] = fitresult->ParError(1)*(activation_end-activation_start)/number_of_alphas;
	myCanvas->SetName("labr_2_rise");
	myCanvas->Write();

	//decay
	cout << "Decay fittings" << endl;
	TF1* decay = new TF1("decay","[0]+[1]*exp(-[2]*(x[0]-[3]))");
	decay->SetNpx(ACTIVATION_NBINS);
	decay->SetNumberFitPoints(ACTIVATION_NBINS);
	decay->SetParLimits(0, 0, 100);
	decay->SetParLimits(2, 4E-15, 5E-15);
	decay->SetParameters(15, 1000, 4.62406E-15, activation_end);
	decay->FixParameter(3, activation_end);
	decay->FixParameter(2, 4.62406E-15);
	decay->SetParNames("Background activity", "Initial activiy", "Decay constant", "activation_end");

	fitresult = labr_1_decay->Fit("decay", "SL");
	results[0][0] = fitresult->Parameter(1);
	results[0][1] = fitresult->ParError(1);		//TBD:error no completo
	myCanvas->SetName("labr_1_decay");
	myCanvas->Write();

	fitresult = labr_2_decay->Fit("decay", "SL");
	results[1][0] = fitresult->Parameter(1);
	results[1][1] = fitresult->ParError(1);		//TBD:error no completo
	myCanvas->SetName("labr_2_decay");
	myCanvas->Write();

	myCanvas->Close();
	DisableImplicitMT();	//multithreading
}
