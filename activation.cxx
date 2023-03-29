void activation(){
	Double_t activation1_decay_start = 1350E12;
	Double_t activation1_decay_end = 2240E12;
	Double_t activation2_decay_start = 1440E12;
	Double_t activation2_decay_end = 2305E12;

	TF1* decay = new TF1("decay","pol0(0)+expo(1)");
	decay->SetParameters(250, 14, 4E-15);
	decay->SetParNames("Background activity", "Exponential constant", "Decay constant");

	TFile f("output/output.root", "UPDATE");
	gDirectory->cd("Activation");

	gDirectory->cd("activation_1");
	TH1F* labr_1_time_1 = (TH1F*)gDirectory->Get("labr_1_time");
	labr_1_time_1->Fit("decay", "", "", activation1_decay_start, activation1_decay_end);
	TH1F* labr_2_time_1 = (TH1F*)gDirectory->Get("labr_2_time");
	labr_2_time_1->Fit("decay", "", "", activation1_decay_start, activation1_decay_end);

	f.Close();
}
