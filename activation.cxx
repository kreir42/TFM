static void per_file(Char_t filepath[500], Double_t* results);

void activation(){
	cout << "Activation" << endl;
	Char_t filepath_1[100] = "output/SData_aAl_J78kV_GVM1808kV_positions2_activacion.root";
	Char_t filepath_2[100] = "output/SData_aAl_J78kV_GVM2312kV_positions2_activacion.root";
	Char_t filepath_3[100] = "output/SData_aAl_J78kV_GVM2810kV_positions2_activacion.root";
	Char_t filepath_4[100] = "output/SData_aAl_J78kV_GVM2810kV_positions2_activacion_20230223.root";

	TFile f("output/output.root", "UPDATE");
	gDirectory->cd("Activation");
	Double_t results[8];
	Double_t activation_energies[8];

	gDirectory->cd("activation_1");
	cout << "activation_1" << endl;
	per_file(filepath_1, results);
	gDirectory->cd("..");

	gDirectory->cd("activation_2");
	cout << "activation_2" << endl;
	per_file(filepath_2, results+1*2);
	gDirectory->cd("..");

	gDirectory->cd("activation_3");
	cout << "activation_3" << endl;
	per_file(filepath_3, results+2*2);
	gDirectory->cd("..");

	gDirectory->cd("activation_4");
	cout << "activation_4" << endl;
	per_file(filepath_4, results+3*2);
	gDirectory->cd("..");

	f.Close();
}

static void per_file(Char_t filepath[500], Double_t* results){
	EnableImplicitMT();	//multithreading
	RDataFrame d("Data", filepath);

	Double_t activation_start = d.Filter("Channel==1").Min("Timestamp").GetValue();	//TBD!:muy ineficiente!!
	Double_t activation_end = d.Filter("Channel==1").Max("Timestamp").GetValue();
	Double_t measurement_end = d.Filter("(Channel==6||Channel==7) && Energy>0").Max("Timestamp").GetValue();

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
	//decay
	TF1* decay = new TF1("decay","pol0(0)+expo(1)");
	decay->SetParameters(15, 14, 4.6E-15);
	decay->SetParNames("Background activity", "Exponential constant", "Decay constant");

	cout << "Decay fittings" << endl;
	labr_1_decay->Fit("decay");
	labr_2_decay->Fit("decay");

	//rise
	TF1* rise = new TF1("rise","[0]+[1]*(1-exp(-[2]*x[0]))");
	rise->SetParameters(15, 1E3, 4.6E-15);
	rise->SetParNames("Background activity", "Constant creation", "Decay constant");

	cout << "Rise fittings" << endl;
	labr_1_rise->Fit("rise");
	labr_2_rise->Fit("rise");

	DisableImplicitMT();	//multithreading
}
