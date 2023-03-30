void activation(){
	Double_t activation1_decay_start = 1350E12;
	Double_t activation1_decay_end = 2240E12;
	Double_t activation2_decay_start = 1440E12;
	Double_t activation2_decay_end = 2305E12;
	Double_t activation3_decay_start = 1180E12;
	Double_t activation3_decay_end = 2060E12;
	Double_t activation4_decay_start = 1290E12;
	Double_t activation4_decay_end = 2480E12;

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

	gDirectory->cd("activation_2");
	TH1F* labr_1_time_2 = (TH1F*)gDirectory->Get("labr_1_time");
	labr_1_time_1->Fit("decay", "", "", activation2_decay_start, activation2_decay_end);
	TH1F* labr_2_time_2 = (TH1F*)gDirectory->Get("labr_2_time");
	labr_2_time_1->Fit("decay", "", "", activation2_decay_start, activation2_decay_end);

	gDirectory->cd("activation_3");
	TH1F* labr_1_time_3 = (TH1F*)gDirectory->Get("labr_1_time");
	labr_1_time_1->Fit("decay", "", "", activation3_decay_start, activation3_decay_end);
	TH1F* labr_2_time_3 = (TH1F*)gDirectory->Get("labr_2_time");
	labr_2_time_1->Fit("decay", "", "", activation3_decay_start, activation3_decay_end);

	gDirectory->cd("activation_4");
	TH1F* labr_1_time_4 = (TH1F*)gDirectory->Get("labr_1_time");
	labr_1_time_1->Fit("decay", "same", "", activation4_decay_start, activation4_decay_end);
	TH1F* labr_2_time_4 = (TH1F*)gDirectory->Get("labr_2_time");
	labr_2_time_1->Fit("decay", "", "", activation4_decay_start, activation4_decay_end);

	f.Close();
}
