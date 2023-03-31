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
		x[2*i+1] = x[i];
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

	auto labr_1_rise = d.Filter("Channel==6 && Energy>525 && Energy<650").Filter(rise_filter, {"Timestamp"}).Histo1D({"labr_1_rise", "; Timestamp; Counts", 500, activation_start, activation_end}, "Timestamp");
	auto labr_1_decay = d.Filter("Channel==6 && Energy>525 && Energy<650").Filter(decay_filter, {"Timestamp"}).Histo1D({"labr_1_decay", "; Timestamp; Counts", 500, activation_end, measurement_end}, "Timestamp");
	auto labr_2_rise = d.Filter("Channel==7 && Energy>525 && Energy<650").Filter(rise_filter, {"Timestamp"}).Histo1D({"labr_2_rise", "; Timestamp; Counts", 500, activation_start, activation_end}, "Timestamp");
	auto labr_2_decay = d.Filter("Channel==7 && Energy>525 && Energy<650").Filter(decay_filter, {"Timestamp"}).Histo1D({"labr_2_decay", "; Timestamp; Counts", 500, activation_end, measurement_end}, "Timestamp");

	cout << "Rise/decay histograms" << endl;
	labr_1_rise->Write();
	labr_1_decay->Write();
	labr_2_rise->Write();
	labr_2_decay->Write();

	//fittings
	TCanvas* myCanvas = new TCanvas("myCanvas");
	TFitResultPtr fitresult;

	//decay
	cout << "Decay fittings" << endl;
	TF1* decay = new TF1("decay","[0]+[1]*exp(-[2]*(x[0]-[3]))");
	decay->SetParLimits(0, 0, 100);
	decay->SetParLimits(2, 4E-15, 5E-15);
	decay->SetParameters(15, 1000, 4.6E-15, activation_end);
	decay->FixParameter(3, activation_end);
	decay->FixParameter(2, 4.6E-15);
	decay->SetParNames("Background activity", "Initial activiy", "Decay constant", "activation_end");

	fitresult = labr_1_decay->Fit("decay", "S");
	results[0][0] = exp(fitresult->Parameter(1)+fitresult->Parameter(2)*activation_end);
	results[0][1] = exp(fitresult->ParError(1));	//TBD:error no completo
	myCanvas->SetName("labr_1_decay");
	myCanvas->Write();

	fitresult = labr_2_decay->Fit("decay", "S");
	results[1][0] = fitresult->Parameter(1);
	results[1][1] = fitresult->ParError(1);	//TBD:error no completo
	myCanvas->SetName("labr_2_decay");
	myCanvas->Write();

	//rise
	cout << "Rise fittings" << endl;
	TF1* rise = new TF1("rise","[0]+[1]*(1-exp(-[2]*(x[0]-[3])))");
	rise->SetParLimits(0, 0, 100);
	rise->SetParLimits(1, 0, 1E5);
	rise->SetParLimits(2, 4E-15, 5E-15);
	rise->SetParameters(15, 5E3, 4.6E-15, activation_start);
	rise->FixParameter(3, activation_start);
	rise->FixParameter(2, 4.6E-15);
	rise->SetParNames("Background activity", "Constant creation", "Decay constant", "activation_start");

	fitresult = labr_1_rise->Fit("rise", "S");
	results[0][2] = fitresult->Parameter(1)*(activation_end-activation_start)/number_of_alphas;
	results[0][3] = fitresult->ParError(1)*(activation_end-activation_start)/number_of_alphas;
	myCanvas->SetName("labr_1_rise");
	myCanvas->Write();

	fitresult = labr_2_rise->Fit("rise", "S");
	results[1][2] = fitresult->Parameter(1)*(activation_end-activation_start)/number_of_alphas;
	results[1][3] = fitresult->ParError(1)*(activation_end-activation_start)/number_of_alphas;
	myCanvas->SetName("labr_2_rise");
	myCanvas->Write();

	myCanvas->Close();
	DisableImplicitMT();	//multithreading
}
