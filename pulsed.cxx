void pulsed_per_file(char filepath[500]){
	EnableImplicitMT();
	RDataFrame d("Data", filepath);

	auto monster = d.Filter("Channel==4");
	Double_t min_tof = 250;
	Double_t max_tof = 650;
	auto tof_plot = monster.Histo1D({"tof_plot", ";ToF;Counts", 1000, min_tof, max_tof}, "tof");
	auto psd_plot = monster.Histo1D({"psd_plot", ";psd;Counts", 1000, 0, 1}, "psd");
	auto tof_id_plot = monster.Histo2D({"tof_id_plot", ";ToF;PSD;Counts", 100, min_tof, max_tof, 100, 0, 1}, "tof", "psd");
	auto energy_id_plot = monster.Histo2D({"energy_id_plot", ";Energy;PSD;Counts", 4096/4, 0, 4096, 100, 0, 1}, "Energy", "psd");

	auto monster_neutrons = monster.Filter("psd>0.3");
	auto neutron_tof_plot = monster_neutrons.Histo1D({"neutron_tof_plot", ";ToF;Counts", 1000, min_tof, max_tof}, "tof");

	TCanvas* myCanvas = new TCanvas("");

	tof_plot->Write("", TObject::kOverwrite);
	psd_plot->Write("", TObject::kOverwrite);
	tof_id_plot->Draw("COLZ");
	myCanvas->Write("tof_id_plot", TObject::kOverwrite);
	energy_id_plot->Draw("COLZ");
	myCanvas->Write("energy_id_plot", TObject::kOverwrite);

	neutron_tof_plot->Write("", TObject::kOverwrite);

	myCanvas->Close();
	DisableImplicitMT();
}

void pulsed(){
	char pulsed_1[500] = PULSED_1_PATH;
	char pulsed_2[500] = PULSED_2_PATH;
	char pulsed_3[500] = PULSED_3_PATH;
	char pulsed_4[500] = PULSED_4_PATH;
	char pulsed_5[500] = PULSED_5_PATH;

	TFile f("output.root", "UPDATE");
	gDirectory->cd("Pulsed");

	gDirectory->cd("pulsed_1");
	pulsed_per_file(pulsed_1);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_2");
	pulsed_per_file(pulsed_2);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_3");
	pulsed_per_file(pulsed_3);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_4");
	pulsed_per_file(pulsed_4);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_5");
	pulsed_per_file(pulsed_5);
	gDirectory->cd("..");

	gDirectory->cd("..");
	f.Close();
}
