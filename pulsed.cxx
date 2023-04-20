void pulsed_per_file(char filepath[500]){
	EnableImplicitMT();
	RDataFrame d("Data", filepath);

	auto monster = d.Filter("Channel==4");
	Double_t max_tof = monster.Max("tof").GetValue();
	auto tof_plot = monster.Histo1D({"tof_plot", ";ToF;Counts", 1000, 0, max_tof}, "tof");
	auto tof_id_plot = monster.Histo2D({"tof_id_plot", ";ToF;PSD;Counts", 1000, 0, max_tof, 25, 0, 1}, "tof", "psd");
	auto energy_id_plot = monster.Histo2D({"energy_id_plot", ";Energy;PSD;Counts", 4096, 0, 4096, 25, 0, 1}, "Energy", "psd");

	tof_plot->Write("", TObject::kOverwrite);
	tof_id_plot->Write("", TObject::kOverwrite);
	energy_id_plot->Write("", TObject::kOverwrite);
	DisableImplicitMT();
}

void pulsed(){
	char pulsed_1[500] = PULSED_1_PATH;

	TFile f("output.root", "UPDATE");
	gDirectory->cd("Pulsed");

	gDirectory->cd("pulsed_1");
	pulsed_per_file(pulsed_1);
	gDirectory->cd("..");

	gDirectory->cd("..");
	f.Close();
}
