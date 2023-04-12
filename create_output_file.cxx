void create_output_file(){
	TFile activation_1("output/SData_aAl_J78kV_GVM1808kV_positions2_activacion.root", "READ");
	TFile activation_2("output/SData_aAl_J78kV_GVM2312kV_positions2_activacion.root", "READ");
	TFile activation_3("output/SData_aAl_J78kV_GVM2810kV_positions2_activacion.root", "READ");
	TFile activation_4("output/SData_aAl_J78kV_GVM2810kV_positions2_activacion_20230223.root", "READ");

	TFile f("output.root", "RECREATE");

	//Activation measurements
	f.mkdir("Activation");
	gDirectory->cd("Activation");

	gDirectory->mkdir("activation_1");
	gDirectory->cd("activation_1");
	gDirectory->WriteObject(activation_1.Get("labr_1_spectrum"), "labr_1_spectrum");
	gDirectory->WriteObject(activation_1.Get("labr_1_time"), "labr_1_time");
	gDirectory->WriteObject(activation_1.Get("labr_2_spectrum"), "labr_2_spectrum");
	gDirectory->WriteObject(activation_1.Get("labr_2_time"), "labr_2_time");
	gDirectory->cd("..");

	gDirectory->mkdir("activation_2");
	gDirectory->cd("activation_2");
	gDirectory->WriteObject(activation_2.Get("labr_1_spectrum"), "labr_1_spectrum");
	gDirectory->WriteObject(activation_2.Get("labr_1_time"), "labr_1_time");
	gDirectory->WriteObject(activation_2.Get("labr_2_spectrum"), "labr_2_spectrum");
	gDirectory->WriteObject(activation_2.Get("labr_2_time"), "labr_2_time");
	gDirectory->cd("..");

	gDirectory->mkdir("activation_3");
	gDirectory->cd("activation_3");
	gDirectory->WriteObject(activation_3.Get("labr_1_spectrum"), "labr_1_spectrum");
	gDirectory->WriteObject(activation_3.Get("labr_1_time"), "labr_1_time");
	gDirectory->WriteObject(activation_3.Get("labr_2_spectrum"), "labr_2_spectrum");
	gDirectory->WriteObject(activation_3.Get("labr_2_time"), "labr_2_time");
	gDirectory->cd("..");

	gDirectory->mkdir("activation_4");
	gDirectory->cd("activation_4");
	gDirectory->WriteObject(activation_4.Get("labr_1_spectrum"), "labr_1_spectrum");
	gDirectory->WriteObject(activation_4.Get("labr_1_time"), "labr_1_time");
	gDirectory->WriteObject(activation_4.Get("labr_2_spectrum"), "labr_2_spectrum");
	gDirectory->WriteObject(activation_4.Get("labr_2_time"), "labr_2_time");
	gDirectory->cd("..");

	gDirectory->cd("..");

	f.Close();
	activation_1.Close();
	activation_2.Close();
	activation_3.Close();
	activation_4.Close();
	cout << "Archivo output.root creado" << endl;
}
