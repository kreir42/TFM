#define ACTIVATION_NBINS 1000
#define ACT1_AENERGY 5500
#define ACT2_AENERGY 7000
#define ACT3_AENERGY 8500
#define ACT4_AENERGY 8500
#define ACT5_AENERGY 5500
#define ACT6_AENERGY 5500
#define ACT7_AENERGY 8250
#define ACT8_AENERGY 7000
#define ACT9_AENERGY 5500
#define ACT10_AENERGY 7500

#include "peak_activity.cxx"

static void per_file(Char_t filepath[500], Double_t results[2][6]);

void activation_results(){
	//escalado con na22
	//en logbook
	Double_t labr1_sodio_1 = peak_activity("output/SData_LaBr_Na22atTarget_calib_20230223.root", 6, 550, 620);
	cout << "labr1_sodio_1: " << labr1_sodio_1 << endl;
	Double_t labr2_sodio_1 = peak_activity("output/SData_LaBr_Na22atTarget_calib_20230223.root", 7, 560, 640);
	cout << "labr2_sodio_1: " << labr2_sodio_1 << endl;
	//en logbook
	Double_t labr1_sodio_2 = peak_activity("output/SData_LaBr_Na22atTarget_calib_20230418.root", 6, 1180, 1310);
	cout << "labr1_sodio_2: " << labr1_sodio_2 << endl;
	Double_t labr2_sodio_2 = peak_activity("output/SData_LaBr_Na22atTarget_calib_20230418.root", 7, 1280, 1400);
	cout << "labr2_sodio_2: " << labr2_sodio_2 << endl;
	//no en logbook, pone febrero pero probablemente abril
	Double_t labr1_sodio_3 = peak_activity("output/SData_LaBr1y2_Na22atTarget_calib_20230223.root", 6, 1190, 1300);
	cout << "labr1_sodio_3: " << labr1_sodio_3 << endl;
	Double_t labr2_sodio_3 = peak_activity("output/SData_LaBr1y2_Na22atTarget_calib_20230223.root", 7, 1280, 1400);
	cout << "labr2_sodio_3: " << labr2_sodio_3 << endl;

	cout << "Activation results" << endl;
	TFile f("output.root", "UPDATE");
	gDirectory->cd("Activation");

	Double_t exfor_energies[] = {3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500, 5600, 5700, 5800, 5900, 6000, 6100, 6200, 6300, 6400, 6500, 6600, 6700, 6800, 6900, 7000, 7100, 7200, 7300, 7400, 7500, 7600, 7700, 7800, 7900, 8000, 8100, 8200, 8300, 8400, 8500, 8600, 8700, 8800, 8900, 9000, 9100, 9200, 9300, 9400, 9500, 9600, 9700, 9800, 9900};	//TBD:hardcoded, read .txt
	Double_t exfor_data_1[] = {3.147E-09, 5.904E-09, 1.034E-08, 1.655E-08, 2.462E-08, 3.464E-08, 4.686E-08, 6.203E-08, 8.124E-08, 1.073E-07, 1.426E-07, 1.847E-07, 2.306E-07, 2.812E-07, 3.403E-07, 4.150E-07, 5.119E-07, 6.278E-07, 7.555E-07, 8.856E-07, 1.011E-06, 1.150E-06, 1.330E-06, 1.549E-06, 1.797E-06, 2.062E-06, 2.339E-06, 2.651E-06, 3.015E-06, 3.401E-06, 3.774E-06, 4.147E-06, 4.552E-06, 4.999E-06, 5.489E-06, 6.013E-06, 6.562E-06, 7.131E-06, 7.716E-06, 8.319E-06, 8.943E-06, 9.593E-06, 1.027E-05, 1.099E-05, 1.173E-05, 1.252E-05, 1.333E-05, 1.416E-05, 1.502E-05, 1.589E-05, 1.679E-05, 1.771E-05, 1.865E-05, 1.962E-05, 2.062E-05, 2.165E-05, 2.272E-05, 2.383E-05, 2.497E-05, 2.616E-05, 2.740E-05, 2.869E-05, 3.003E-05};	//TBD:hardcoded, read .txt
	Double_t exfor_data_2[] = {3.147E-09, 5.904E-09, 1.034E-08, 1.655E-08, 2.462E-08, 3.464E-08, 4.686E-08, 6.203E-08, 8.124E-08, 1.073E-07, 1.426E-07, 1.847E-07, 2.306E-07, 2.812E-07, 3.403E-07, 4.150E-07, 5.119E-07, 6.278E-07, 7.555E-07, 8.856E-07, 1.011E-06, 1.150E-06, 1.330E-06, 1.549E-06, 1.797E-06, 2.062E-06, 2.339E-06, 2.651E-06, 3.015E-06, 3.401E-06, 3.774E-06, 4.147E-06, 4.552E-06, 4.999E-06, 5.489E-06, 6.013E-06, 6.562E-06, 7.131E-06, 7.716E-06, 8.319E-06, 8.943E-06, 9.593E-06, 1.027E-05, 1.099E-05, 1.173E-05, 1.252E-05, 1.333E-05, 1.416E-05, 1.502E-05, 1.589E-05, 1.679E-05, 1.771E-05, 1.865E-05, 1.962E-05, 2.062E-05, 2.165E-05, 2.272E-05, 2.383E-05, 2.497E-05, 2.616E-05, 2.740E-05, 2.869E-05, 3.003E-05};	//TBD:hardcoded, read .txt
	Double_t exfor_errors_1[63];
	Double_t exfor_errors_2[63];
	for(short i=0; i<63; i++){	//TBD:escalado temporal, números hardcoded
		exfor_data_1[i]*=2.220341468E-5;	//eficiencia feb
		exfor_errors_1[i] = exfor_data_1[i] * 0.05;
		exfor_data_2[i]*=2.642141092E-5;	//eficiencia apr
		exfor_errors_2[i] = exfor_data_2[i] * 0.08;
	}

	//EXFOR data
	TGraph* rectionsvenergy_exfor_feb = new TGraphErrors(63, exfor_energies, exfor_data_1, NULL, exfor_errors_1);	//TBD:hardcoded number
	rectionsvenergy_exfor_feb->SetTitle("(a,n) reactions v a energy;Energy of a (keV);Inferred (a,n)/Number of a");
	rectionsvenergy_exfor_feb->SetFillColorAlpha(kGray, 0.4);
	rectionsvenergy_exfor_feb->SetMarkerStyle(20);
	rectionsvenergy_exfor_feb->SetLineColor(kBlack);
	TGraph* rectionsvenergy_exfor_apr = new TGraphErrors(63, exfor_energies, exfor_data_2, NULL, exfor_errors_2);	//TBD:hardcoded number
	rectionsvenergy_exfor_apr->SetTitle("(a,n) reactions v a energy;Energy of a (keV);Inferred (a,n)/Number of a");
	rectionsvenergy_exfor_apr->SetFillColorAlpha(kGray+1, 0.4);
	rectionsvenergy_exfor_apr->SetMarkerStyle(20);
	rectionsvenergy_exfor_apr->SetLineColor(kBlack);

	//meter resultados en array
	Double_t results[ACTIVATION_N][2][6];
	TTree* tree = (TTree*)gDirectory->Get("activation_results_tree");
	tree->SetBranchAddress("results", results);
	tree->GetEntry(0);

	//graficas
	Double_t activation_energies[] = {ACT1_AENERGY, ACT2_AENERGY, ACT3_AENERGY, ACT4_AENERGY, ACT5_AENERGY, ACT6_AENERGY, ACT7_AENERGY, ACT8_AENERGY, ACT9_AENERGY, ACT10_AENERGY};
	Double_t x1[ACTIVATION_N];
	Double_t x2[ACTIVATION_N];
	Double_t y1[ACTIVATION_N];
	Double_t y2[ACTIVATION_N];
	Double_t yerr1[ACTIVATION_N];
	Double_t yerr2[ACTIVATION_N];
	TCanvas* myCanvas = new TCanvas("30P per alpha sin escalar");

	//unified fit sin escalar
	for(short i=0; i<ACTIVATION_N; i++){
		x1[i] = activation_energies[i];
		x2[i] = activation_energies[i];
		y1[i] = results[i][0][0];
		y2[i] = results[i][1][0];
		yerr1[i] = results[i][0][1];
		yerr2[i] = results[i][1][1];
	}
	TGraph* rectionsvenergy_unified_feb_sin_escalar_1 = new TGraphErrors(4, x1, y1, NULL, yerr1);
	rectionsvenergy_unified_feb_sin_escalar_1->SetMarkerStyle(22);
	TGraph* rectionsvenergy_unified_apr_sin_escalar_1 = new TGraphErrors(3, &x1[4], &y1[4], NULL, yerr1);
	rectionsvenergy_unified_apr_sin_escalar_1->SetMarkerStyle(23);
	TGraph* rectionsvenergy_unified_feb_sin_escalar_2 = new TGraphErrors(4, x2, y2, NULL, yerr2);
	rectionsvenergy_unified_feb_sin_escalar_2->SetMarkerStyle(34);
	TGraph* rectionsvenergy_unified_apr_sin_escalar_2 = new TGraphErrors(3, &x2[4], &y2[4], NULL, yerr2);
	rectionsvenergy_unified_apr_sin_escalar_2->SetMarkerStyle(47);

	//rise fit sin escalar
	for(short i=0; i<ACTIVATION_N; i++){
		x1[i] = activation_energies[i];
		x2[i] = activation_energies[i];
		y1[i] = results[i][0][2];
		y2[i] = results[i][1][2];
		yerr1[i] = results[i][0][3];
		yerr2[i] = results[i][1][3];
	}
	TGraph* rectionsvenergy_rise_feb_sin_escalar_1 = new TGraphErrors(4, x1, y1, NULL, yerr1);
	rectionsvenergy_rise_feb_sin_escalar_1->SetMarkerStyle(22);
	TGraph* rectionsvenergy_rise_apr_sin_escalar_1 = new TGraphErrors(2, &x1[4], &y1[4], NULL, yerr1);
	rectionsvenergy_rise_apr_sin_escalar_1->SetMarkerStyle(23);
	TGraph* rectionsvenergy_rise_feb_sin_escalar_2 = new TGraphErrors(4, x2, y2, NULL, yerr2);
	rectionsvenergy_rise_feb_sin_escalar_2->SetMarkerStyle(34);
	TGraph* rectionsvenergy_rise_apr_sin_escalar_2 = new TGraphErrors(2, &x2[4], &y2[4], NULL, yerr2);
	rectionsvenergy_rise_apr_sin_escalar_2->SetMarkerStyle(47);

	//decay fit sin escalar
	for(short i=0; i<ACTIVATION_N; i++){
		x1[i] = activation_energies[i];
		x2[i] = activation_energies[i];
		y1[i] = results[i][0][4];
		y2[i] = results[i][1][4];
		yerr1[i] = results[i][0][5];
		yerr2[i] = results[i][1][5];
	}
	TGraph* rectionsvenergy_decay_feb_sin_escalar_1 = new TGraphErrors(4, x1, y1, NULL, yerr1);
	rectionsvenergy_decay_feb_sin_escalar_1->SetMarkerStyle(22);
	TGraph* rectionsvenergy_decay_apr_sin_escalar_1 = new TGraphErrors(6, &x1[4], &y1[4], NULL, yerr1);
	rectionsvenergy_decay_apr_sin_escalar_1->SetMarkerStyle(23);
	TGraph* rectionsvenergy_decay_feb_sin_escalar_2 = new TGraphErrors(4, x2, y2, NULL, yerr2);
	rectionsvenergy_decay_feb_sin_escalar_2->SetMarkerStyle(34);
	TGraph* rectionsvenergy_decay_apr_sin_escalar_2 = new TGraphErrors(6, &x2[4], &y2[4], NULL, yerr2);
	rectionsvenergy_decay_apr_sin_escalar_2->SetMarkerStyle(47);

	//todos los resultados antes de los escalados
	TMultiGraph* multigraph_sin_escalar = new TMultiGraph();
	rectionsvenergy_unified_feb_sin_escalar_1->SetMarkerColor(kRed);
	rectionsvenergy_unified_feb_sin_escalar_1->SetTitle("Unified fit, February, LaBr1");
	multigraph_sin_escalar->Add(rectionsvenergy_unified_feb_sin_escalar_1);
	rectionsvenergy_unified_feb_sin_escalar_2->SetMarkerColor(kRed);
	rectionsvenergy_unified_feb_sin_escalar_2->SetTitle("Unified fit, February, LaBr2");
	multigraph_sin_escalar->Add(rectionsvenergy_unified_feb_sin_escalar_2);
	rectionsvenergy_unified_apr_sin_escalar_1->SetMarkerColor(kRed);
	rectionsvenergy_unified_apr_sin_escalar_1->SetTitle("Unified fit, April, LaBr1");
	multigraph_sin_escalar->Add(rectionsvenergy_unified_apr_sin_escalar_1);
	rectionsvenergy_unified_apr_sin_escalar_2->SetMarkerColor(kRed);
	rectionsvenergy_unified_apr_sin_escalar_2->SetTitle("Unified fit, April, LaBr2");
	multigraph_sin_escalar->Add(rectionsvenergy_unified_apr_sin_escalar_2);
	rectionsvenergy_rise_feb_sin_escalar_1->SetMarkerColor(kGreen);
	rectionsvenergy_rise_feb_sin_escalar_1->SetTitle("Rise fit, February, LaBr1");
	multigraph_sin_escalar->Add(rectionsvenergy_rise_feb_sin_escalar_1);
	rectionsvenergy_rise_feb_sin_escalar_2->SetMarkerColor(kGreen);
	rectionsvenergy_rise_feb_sin_escalar_2->SetTitle("Rise fit, February, LaBr2");
	multigraph_sin_escalar->Add(rectionsvenergy_rise_feb_sin_escalar_2);
	rectionsvenergy_rise_apr_sin_escalar_1->SetMarkerColor(kGreen);
	rectionsvenergy_rise_apr_sin_escalar_1->SetTitle("Rise fit, April, LaBr1");
	multigraph_sin_escalar->Add(rectionsvenergy_rise_apr_sin_escalar_1);
	rectionsvenergy_rise_apr_sin_escalar_2->SetMarkerColor(kGreen);
	rectionsvenergy_rise_apr_sin_escalar_2->SetTitle("Rise fit, April, LaBr2");
	multigraph_sin_escalar->Add(rectionsvenergy_rise_apr_sin_escalar_2);
	rectionsvenergy_decay_feb_sin_escalar_1->SetMarkerColor(kBlue);
	rectionsvenergy_decay_feb_sin_escalar_1->SetTitle("Decay fit, February, LaBr1");
	multigraph_sin_escalar->Add(rectionsvenergy_decay_feb_sin_escalar_1);
	rectionsvenergy_decay_feb_sin_escalar_2->SetMarkerColor(kBlue);
	rectionsvenergy_decay_feb_sin_escalar_2->SetTitle("Decay fit, February, LaBr2");
	multigraph_sin_escalar->Add(rectionsvenergy_decay_feb_sin_escalar_2);
	rectionsvenergy_decay_apr_sin_escalar_1->SetMarkerColor(kBlue);
	rectionsvenergy_decay_apr_sin_escalar_1->SetTitle("Decay fit, April, LaBr1");
	multigraph_sin_escalar->Add(rectionsvenergy_decay_apr_sin_escalar_1);
	rectionsvenergy_decay_apr_sin_escalar_2->SetMarkerColor(kBlue);
	rectionsvenergy_decay_apr_sin_escalar_2->SetTitle("Decay fit, April, LaBr2");
	multigraph_sin_escalar->Add(rectionsvenergy_decay_apr_sin_escalar_2);
	multigraph_sin_escalar->SetTitle("30P v a energy;Energy of a (keV);Inferred (a,n)/Number of a");
	multigraph_sin_escalar->Draw("AP");
	myCanvas->BuildLegend();
	myCanvas->Write("", TObject::kOverwrite);

	//escalados
	for(unsigned short j=0; j<6; j++){
		//escalado debido a la mala medida de la carga
		results[3][0][j]*=172/239.1;
		results[3][1][j]*=172/239.1;

		//escalado na22
		results[0][0][j]/=labr1_sodio_1;
		results[0][1][j]/=labr2_sodio_1;
		results[1][0][j]/=labr1_sodio_1;
		results[1][1][j]/=labr2_sodio_1;
		results[2][0][j]/=labr1_sodio_1;
		results[2][1][j]/=labr2_sodio_1;
		results[3][0][j]/=labr1_sodio_1;
		results[3][1][j]/=labr2_sodio_1;
		results[4][0][j]/=labr1_sodio_3;
		results[4][1][j]/=labr2_sodio_3;
		results[5][0][j]/=labr1_sodio_2;
		results[5][1][j]/=labr2_sodio_2;
		results[6][0][j]/=labr1_sodio_3;
		results[6][1][j]/=labr2_sodio_3;
		results[7][0][j]/=labr1_sodio_3;
		results[7][1][j]/=labr2_sodio_3;
		results[8][0][j]/=labr1_sodio_3;
		results[8][1][j]/=labr2_sodio_3;
		results[9][0][j]/=labr1_sodio_3;
		results[9][1][j]/=labr2_sodio_3;
	}

	//unified_fit
	for(short i=0; i<ACTIVATION_N; i++){
		x1[i] = activation_energies[i];
		x2[i] = activation_energies[i];
		y1[i] = results[i][0][0];
		y2[i] = results[i][1][0];
		yerr1[i] = results[i][0][1];
		yerr2[i] = results[i][1][1];
	}
	TGraph* rectionsvenergy_unified_feb_1 = new TGraphErrors(4, x1, y1, NULL, yerr1);
	rectionsvenergy_unified_feb_1->SetTitle("(a,n) reactions v a energy (unified fit, February);Energy of a (keV);Inferred (a,n)/Number of a");
	myCanvas->SetName("reactions_v_energy_unified_feb_1");
	rectionsvenergy_unified_feb_1->SetMarkerStyle(22);
	rectionsvenergy_exfor_feb->Draw("E3");
	rectionsvenergy_exfor_feb->Draw("same LX");
	rectionsvenergy_unified_feb_1->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TGraph* rectionsvenergy_unified_feb_2 = new TGraphErrors(4, x2, y2, NULL, yerr2);
	rectionsvenergy_unified_feb_2->SetTitle("(a,n) reactions v a energy (unified fit, February);Energy of a (keV);Inferred (a,n)/Number of a");
	myCanvas->SetName("reactions_v_energy_unified_feb_2");
	rectionsvenergy_unified_feb_2->SetMarkerStyle(34);
	rectionsvenergy_exfor_feb->Draw("E3");
	rectionsvenergy_exfor_feb->Draw("same LX");
	rectionsvenergy_unified_feb_2->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TGraph* rectionsvenergy_unified_apr_1 = new TGraphErrors(3, &x1[4], &y1[4], NULL, yerr1);
	rectionsvenergy_unified_apr_1->SetTitle("(a,n) reactions v a energy (unified fit, April);Energy of a (keV);Inferred (a,n)/Number of a");
	myCanvas->SetName("reactions_v_energy_unified_apr_1");
	rectionsvenergy_unified_apr_1->SetMarkerStyle(23);
	rectionsvenergy_exfor_apr->Draw("E3");
	rectionsvenergy_exfor_apr->Draw("same LX");
	rectionsvenergy_unified_apr_1->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TGraph* rectionsvenergy_unified_apr_2 = new TGraphErrors(3, &x2[4], &y2[4], NULL, yerr2);
	rectionsvenergy_unified_apr_2->SetTitle("(a,n) reactions v a energy (unified fit, April);Energy of a (keV);Inferred (a,n)/Number of a");
	myCanvas->SetName("reactions_v_energy_unified_apr_2");
	rectionsvenergy_unified_apr_2->SetMarkerStyle(47);
	rectionsvenergy_exfor_apr->Draw("E3");
	rectionsvenergy_exfor_apr->Draw("same LX");
	rectionsvenergy_unified_apr_2->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TMultiGraph* reactions_v_energy_unified = new TMultiGraph();
	reactions_v_energy_unified->Add(rectionsvenergy_exfor_feb, "E3");
	reactions_v_energy_unified->Add(rectionsvenergy_exfor_feb, "LX");
	rectionsvenergy_exfor_feb->SetTitle("EXFOR data, feb");
	reactions_v_energy_unified->Add(rectionsvenergy_exfor_apr, "E3");
	reactions_v_energy_unified->Add(rectionsvenergy_exfor_apr, "LX");
	rectionsvenergy_exfor_apr->SetTitle("EXFOR data, apr");
	rectionsvenergy_unified_feb_1->SetTitle("Unified fit, February, LaBr1");
	reactions_v_energy_unified->Add(rectionsvenergy_unified_feb_1);
	rectionsvenergy_unified_apr_1->SetTitle("Unified fit, April, LaBr1");
	reactions_v_energy_unified->Add(rectionsvenergy_unified_apr_1);
	rectionsvenergy_unified_feb_2->SetTitle("Unified fit, February, LaBr2");
	reactions_v_energy_unified->Add(rectionsvenergy_unified_feb_2);
	rectionsvenergy_unified_apr_2->SetTitle("Unified fit, April, LaBr2");
	reactions_v_energy_unified->Add(rectionsvenergy_unified_apr_2);
	reactions_v_energy_unified->SetTitle("(a,n) reactions v a energy (unified fit);Energy of a (keV);Inferred (a,n)/Number of a");
	reactions_v_energy_unified->Draw("AP");
	myCanvas->BuildLegend();
	myCanvas->SetName("reactions_v_energy_unified");
	myCanvas->Write("", TObject::kOverwrite);

	//rise fit
	for(short i=0; i<ACTIVATION_N; i++){
		x1[i] = activation_energies[i];
		x2[i] = activation_energies[i];
		y1[i] = results[i][0][2];
		y2[i] = results[i][1][2];
		yerr1[i] = results[i][0][3];
		yerr2[i] = results[i][1][3];
	}
	TGraph* rectionsvenergy_rise_feb_1 = new TGraphErrors(4, x1, y1, NULL, yerr1);
	rectionsvenergy_rise_feb_1->SetTitle("(a,n) reactions v a energy (rise fit, February);Energy of a (keV);Inferred (a,n)/Number of a");
	rectionsvenergy_rise_feb_1->SetMarkerStyle(22);
	myCanvas->SetName("reactions_v_energy_rise_feb_1");
	rectionsvenergy_exfor_feb->Draw("E3");
	rectionsvenergy_exfor_feb->Draw("same LX");
	rectionsvenergy_rise_feb_1->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TGraph* rectionsvenergy_rise_feb_2 = new TGraphErrors(4, x2, y2, NULL, yerr2);
	rectionsvenergy_rise_feb_2->SetTitle("(a,n) reactions v a energy (rise fit, February);Energy of a (keV);Inferred (a,n)/Number of a");
	rectionsvenergy_rise_feb_2->SetMarkerStyle(34);
	myCanvas->SetName("reactions_v_energy_rise_feb_2");
	rectionsvenergy_exfor_feb->Draw("E3");
	rectionsvenergy_exfor_feb->Draw("same LX");
	rectionsvenergy_rise_feb_2->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TGraph* rectionsvenergy_rise_apr_1 = new TGraphErrors(2, &x1[4], &y1[4], NULL, yerr1);
	rectionsvenergy_rise_apr_1->SetTitle("(a,n) reactions v a energy (rise fit, April);Energy of a (keV);Inferred (a,n)/Number of a");
	rectionsvenergy_rise_apr_1->SetMarkerStyle(23);
	myCanvas->SetName("reactions_v_energy_rise_apr_1");
	rectionsvenergy_exfor_apr->Draw("E3");
	rectionsvenergy_exfor_apr->Draw("same LX");
	rectionsvenergy_rise_apr_1->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TGraph* rectionsvenergy_rise_apr_2 = new TGraphErrors(2, &x2[4], &y2[4], NULL, yerr2);
	rectionsvenergy_rise_apr_2->SetTitle("(a,n) reactions v a energy (rise fit, April);Energy of a (keV);Inferred (a,n)/Number of a");
	rectionsvenergy_rise_apr_2->SetMarkerStyle(47);
	myCanvas->SetName("reactions_v_energy_rise_apr_2");
	rectionsvenergy_exfor_apr->Draw("E3");
	rectionsvenergy_exfor_apr->Draw("same LX");
	rectionsvenergy_rise_apr_2->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TMultiGraph* reactions_v_energy_rise = new TMultiGraph();
	reactions_v_energy_rise->Add(rectionsvenergy_exfor_feb, "E3");
	reactions_v_energy_rise->Add(rectionsvenergy_exfor_feb, "LX");
	reactions_v_energy_rise->Add(rectionsvenergy_exfor_apr, "E3");
	reactions_v_energy_rise->Add(rectionsvenergy_exfor_apr, "LX");
	rectionsvenergy_rise_feb_1->SetTitle("Rise fit, February, LaBr1");
	reactions_v_energy_rise->Add(rectionsvenergy_rise_feb_1);
	rectionsvenergy_rise_apr_1->SetTitle("Rise fit, April, LaBr1");
	reactions_v_energy_rise->Add(rectionsvenergy_rise_apr_1);
	rectionsvenergy_rise_feb_2->SetTitle("Rise fit, February, LaBr2");
	reactions_v_energy_rise->Add(rectionsvenergy_rise_feb_2);
	rectionsvenergy_rise_apr_2->SetTitle("Rise fit, April, LaBr2");
	reactions_v_energy_rise->Add(rectionsvenergy_rise_apr_2);
	reactions_v_energy_rise->SetTitle("(a,n) reactions v a energy (rise fit);Energy of a (keV);Inferred (a,n)/Number of a");
	reactions_v_energy_rise->Draw("AP");
	myCanvas->BuildLegend();
	myCanvas->SetName("reactions_v_energy_rise");
	myCanvas->Write("", TObject::kOverwrite);

	//decay fit
	for(short i=0; i<ACTIVATION_N; i++){
		x1[i] = activation_energies[i];
		x2[i] = activation_energies[i];
		y1[i] = results[i][0][4];
		y2[i] = results[i][1][4];
		yerr1[i] = results[i][0][5];
		yerr2[i] = results[i][1][5];
	}
	TGraph* rectionsvenergy_decay_feb_1 = new TGraphErrors(4, x1, y1, NULL, yerr1);
	rectionsvenergy_decay_feb_1->SetTitle("(a,n) reactions v a energy (decay fit, February);Energy of a (keV);Inferred (a,n)/Number of a");
	rectionsvenergy_decay_feb_1->SetMarkerStyle(22);
	myCanvas->SetName("reactions_v_energy_decay_feb_1");
	rectionsvenergy_exfor_feb->Draw("E3");
	rectionsvenergy_exfor_feb->Draw("same LX");
	rectionsvenergy_decay_feb_1->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TGraph* rectionsvenergy_decay_feb_2 = new TGraphErrors(4, x2, y2, NULL, yerr2);
	rectionsvenergy_decay_feb_2->SetTitle("(a,n) reactions v a energy (decay fit, February);Energy of a (keV);Inferred (a,n)/Number of a");
	rectionsvenergy_decay_feb_2->SetMarkerStyle(34);
	myCanvas->SetName("reactions_v_energy_decay_feb_2");
	rectionsvenergy_exfor_feb->Draw("E3");
	rectionsvenergy_exfor_feb->Draw("same LX");
	rectionsvenergy_decay_feb_2->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TGraph* rectionsvenergy_decay_apr_1 = new TGraphErrors(6, &x1[4], &y1[4], NULL, yerr1);
	rectionsvenergy_decay_apr_1->SetTitle("(a,n) reactions v a energy (decay it, April);Energy of a (keV);Inferred (a,n)/Number of a");
	rectionsvenergy_decay_apr_1->SetMarkerStyle(23);
	myCanvas->SetName("reactions_v_energy_decay_apr_1");
	rectionsvenergy_exfor_apr->Draw("E3");
	rectionsvenergy_exfor_apr->Draw("same LX");
	rectionsvenergy_decay_apr_1->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TGraph* rectionsvenergy_decay_apr_2 = new TGraphErrors(6, &x2[4], &y2[4], NULL, yerr2);
	rectionsvenergy_decay_apr_2->SetTitle("(a,n) reactions v a energy (decay it, April);Energy of a (keV);Inferred (a,n)/Number of a");
	rectionsvenergy_decay_apr_2->SetMarkerStyle(47);
	myCanvas->SetName("reactions_v_energy_decay_apr_2");
	rectionsvenergy_exfor_apr->Draw("E3");
	rectionsvenergy_exfor_apr->Draw("same LX");
	rectionsvenergy_decay_apr_2->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TMultiGraph* reactions_v_energy_decay = new TMultiGraph();
	reactions_v_energy_decay->Add(rectionsvenergy_exfor_feb, "E3");
	reactions_v_energy_decay->Add(rectionsvenergy_exfor_feb, "LX");
	reactions_v_energy_decay->Add(rectionsvenergy_exfor_apr, "E3");
	reactions_v_energy_decay->Add(rectionsvenergy_exfor_apr, "LX");
	rectionsvenergy_decay_feb_1->SetTitle("Decay fit, February, LaBr1");
	reactions_v_energy_decay->Add(rectionsvenergy_decay_feb_1);
	rectionsvenergy_decay_apr_1->SetTitle("Decay fit, April, LaBr1");
	reactions_v_energy_decay->Add(rectionsvenergy_decay_apr_1);
	rectionsvenergy_decay_feb_2->SetTitle("Decay fit, February, LaBr2");
	reactions_v_energy_decay->Add(rectionsvenergy_decay_feb_2);
	rectionsvenergy_decay_apr_2->SetTitle("Decay fit, April, LaBr2");
	reactions_v_energy_decay->Add(rectionsvenergy_decay_apr_2);
	reactions_v_energy_decay->SetTitle("(a,n) reactions v a energy (decay fit);Energy of a (keV);Inferred (a,n)/Number of a");
	reactions_v_energy_decay->Draw("AP");
	myCanvas->BuildLegend();
	myCanvas->SetName("reactions_v_energy_decay");
	myCanvas->Write("", TObject::kOverwrite);

	//all results
	TMultiGraph* multigraph = new TMultiGraph();
	multigraph->Add(rectionsvenergy_exfor_feb, "E3");
	multigraph->Add(rectionsvenergy_exfor_feb, "LX");
	multigraph->Add(rectionsvenergy_exfor_apr, "E3");
	multigraph->Add(rectionsvenergy_exfor_apr, "LX");
	rectionsvenergy_unified_feb_1->SetMarkerColor(kRed);
	rectionsvenergy_unified_feb_1->SetTitle("Unified fit, February, LaBr1");
	multigraph->Add(rectionsvenergy_unified_feb_1);
	rectionsvenergy_unified_apr_1->SetMarkerColor(kRed);
	rectionsvenergy_unified_apr_1->SetTitle("Unified fit, April, LaBr1");
	multigraph->Add(rectionsvenergy_unified_apr_1);
	rectionsvenergy_rise_feb_1->SetMarkerColor(kGreen);
	rectionsvenergy_rise_feb_1->SetTitle("Rise fit, February, LaBr1");
	multigraph->Add(rectionsvenergy_rise_feb_1);
	rectionsvenergy_rise_apr_1->SetMarkerColor(kGreen);
	rectionsvenergy_rise_apr_1->SetTitle("Rise fit, April, LaBr1");
	multigraph->Add(rectionsvenergy_rise_apr_1);
	rectionsvenergy_decay_feb_1->SetMarkerColor(kBlue);
	rectionsvenergy_decay_feb_1->SetTitle("Decay fit, February, LaBr1");
	multigraph->Add(rectionsvenergy_decay_feb_1);
	rectionsvenergy_decay_apr_1->SetMarkerColor(kBlue);
	rectionsvenergy_decay_apr_1->SetTitle("Decay fit, April, LaBr1");
	multigraph->Add(rectionsvenergy_decay_apr_1);
	rectionsvenergy_unified_feb_2->SetMarkerColor(kRed);
	rectionsvenergy_unified_feb_2->SetTitle("Unified fit, February, LaBr2");
	multigraph->Add(rectionsvenergy_unified_feb_2);
	rectionsvenergy_unified_apr_2->SetMarkerColor(kRed);
	rectionsvenergy_unified_apr_2->SetTitle("Unified fit, April, LaBr2");
	multigraph->Add(rectionsvenergy_unified_apr_2);
	rectionsvenergy_rise_feb_2->SetMarkerColor(kGreen);
	rectionsvenergy_rise_feb_2->SetTitle("Rise fit, February, LaBr2");
	multigraph->Add(rectionsvenergy_rise_feb_2);
	rectionsvenergy_rise_apr_2->SetMarkerColor(kGreen);
	rectionsvenergy_rise_apr_2->SetTitle("Rise fit, April, LaBr2");
	multigraph->Add(rectionsvenergy_rise_apr_2);
	rectionsvenergy_decay_feb_2->SetMarkerColor(kBlue);
	rectionsvenergy_decay_feb_2->SetTitle("Decay fit, February, LaBr2");
	multigraph->Add(rectionsvenergy_decay_feb_2);
	rectionsvenergy_decay_apr_2->SetMarkerColor(kBlue);
	rectionsvenergy_decay_apr_2->SetTitle("Decay fit, April, LaBr2");
	multigraph->Add(rectionsvenergy_decay_apr_2);
	multigraph->SetTitle("(a,n) reactions v a energy;Energy of a (keV);Inferred (a,n)/Number of a");
	multigraph->Draw("AP");
	myCanvas->BuildLegend();
	myCanvas->SetName("reactions_v_energy");
	myCanvas->Write("", TObject::kOverwrite);

	myCanvas->Close();
	gDirectory->cd("..");
	f.Close();
}

Double_t activation_window_low;
Double_t activation_window_high;

void activation(){
	cout << "Activation" << endl;
	cout << endl;
	//Feb
	Char_t filepath_1[100] = "output/SData_aAl_J78kV_GVM1808kV_positions2_activacion.root";
	Char_t filepath_2[100] = "output/SData_aAl_J78kV_GVM2312kV_positions2_activacion.root";
	Char_t filepath_3[100] = "output/SData_aAl_J78kV_GVM2810kV_positions2_activacion.root";

	Char_t filepath_4[100] = "output/SData_aAl_J78kV_GVM2810kV_positions2_activacion_20230223.root";

	//Apr
	Char_t filepath_5[100] =  "output/SData_aAl_J78keV_GVM1808keV_LaBr1_20cm-135deg_LaBr2_20cm135deg_activacion.root";
	Char_t filepath_6[100] =  "output/SData_aAl_J78keV_GVM1808keV_LaBr1_5cmdelante_LaBr2_20cm_activacion.root";

	Char_t filepath_7[100] =  "output/SData_aAl_J78keV_GVM2731kV_LaBr1_20cm-135deg_LaBr2_20cm135deg_activacion_20230418.root";
	Char_t filepath_8[100] =  "output/SData_aAl_J78keV_GVM2310kV_LaBr1_20cm-135deg_LaBr2_20cm135deg_activacion_20230418.root";
	Char_t filepath_9[100] =  "output/SData_aAl_J78keV_GVM1808kV_LaBr1_20cm-135deg_LaBr2_20cm135deg_activacion_20230418.root";
	Char_t filepath_10[100] = "output/SData_aAl_J78keV_GVM2478kV_LaBr1_20cm-135deg_LaBr2_20cm135deg_activacion_20230418.root";

	TFile f("output.root", "UPDATE");
	gDirectory->cd("Activation");
	Double_t results[ACTIVATION_N][2][6];

	activation_window_low=550;
	activation_window_high=650;

	gDirectory->cd("activation_1");
	cout << "activation_1" << endl;
	per_file(filepath_1, results[0]);
	gDirectory->cd("..");

	gDirectory->cd("activation_2");
	cout << "activation_2" << endl;
	per_file(filepath_2, results[1]);
	gDirectory->cd("..");

	activation_window_low=550;
	activation_window_high=750;

	gDirectory->cd("activation_3");
	cout << "activation_3" << endl;
	per_file(filepath_3, results[2]);
	gDirectory->cd("..");

	activation_window_low=550;
	activation_window_high=825;

	gDirectory->cd("activation_4");
	cout << "activation_4" << endl;
	per_file(filepath_4, results[3]);
	gDirectory->cd("..");

	activation_window_low=1200;
	activation_window_high=1450;

	gDirectory->cd("activation_5");
	cout << "activation_5" << endl;
	per_file(filepath_5, results[4]);
	gDirectory->cd("..");

	activation_window_low=1200;
	activation_window_high=1600;

	gDirectory->cd("activation_6");
	cout << "activation_6" << endl;
	per_file(filepath_6, results[5]);
	gDirectory->cd("..");

	activation_window_low=1200;
	activation_window_high=1800;

	gDirectory->cd("activation_7");
	cout << "activation_7" << endl;
	per_file(filepath_7, results[6]);
	gDirectory->cd("..");

	activation_window_low=1200;
	activation_window_high=1700;

	gDirectory->cd("activation_8");
	cout << "activation_8" << endl;
	per_file(filepath_8, results[7]);
	gDirectory->cd("..");

	activation_window_low=1150;
	activation_window_high=1400;

	gDirectory->cd("activation_9");
	cout << "activation_9" << endl;
	per_file(filepath_9, results[8]);
	gDirectory->cd("..");

	activation_window_low=1200;
	activation_window_high=1700;

	gDirectory->cd("activation_10");
	cout << "activation_10" << endl;
	per_file(filepath_10, results[9]);
	gDirectory->cd("..");

	TTree* tree = new TTree("activation_results_tree", "Tree with activation results");
	tree->Branch("results", results, "results[10][2][6]/D");
	tree->SetBranchAddress("results", results);
	tree->Fill();
	tree->Write("", TObject::kOverwrite);

	gDirectory->cd("..");
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
			Double_t n = p[4];	//initial number of nuclei
			Double_t decay = exp(-p[2]*stepsize);
			Double_t helper_ratio = p[1]/p[2];
			Double_t I = 0;	//created nuclei per time over lambda
			unsigned long steps = x[0]/stepsize;
			for(unsigned long i=0; i<=steps; i++){
				I = (current_histogram->GetBinContent(i))*helper_ratio;
				n = I + (n-I)*decay;
			}
			return p[0] + n*(1-decay) + p[3]*current_histogram->GetBinContent(steps);
		}
};

static void per_file(Char_t filepath[500], Double_t results[2][6]){
	cout << "***************" << endl;
	EnableImplicitMT();	//multithreading
	RDataFrame d("Data", filepath);

	Double_t measurement_end = d.Filter("(Channel==6||Channel==7) && Energy>0").Max("Timestamp").GetValue();
	auto integrator_signals = d.Filter("Channel==1");
	auto current_integrator = integrator_signals.Define("t", "Timestamp/1E12").Histo1D({"current_integrator", ";Timestamp (s);Counts", ACTIVATION_NBINS, 0, measurement_end/1E12}, "t");
	Double_t activation_start = integrator_signals.Min("Timestamp").GetValue();	//TBD!:muy ineficiente!!
	Double_t activation_end = integrator_signals.Max("Timestamp").GetValue();
	Double_t current2alpha = 1/(2*1.60217646E-9);
	Double_t number_of_alphas = current_integrator->Integral()*current2alpha;
	if(strcmp(filepath, "output/SData_aAl_J78keV_GVM2310kV_LaBr1_20cm-135deg_LaBr2_20cm135deg_activacion_20230418.root")==0){	//8
		number_of_alphas = current2alpha * 101.1E4;
		activation_start= 364E12;
		activation_end= 830E12;
	}else if(strcmp(filepath, "output/SData_aAl_J78keV_GVM1808kV_LaBr1_20cm-135deg_LaBr2_20cm135deg_activacion_20230418.root")==0){	//9
		number_of_alphas = current2alpha * 67.96E4;
		activation_start= 302E12;
		activation_end= 753E12;
	}else if(strcmp(filepath, "output/SData_aAl_J78keV_GVM2478kV_LaBr1_20cm-135deg_LaBr2_20cm135deg_activacion_20230418.root")==0){	//10
		number_of_alphas = current2alpha * 81.82E4;
		activation_start= 304E12;
		activation_end= 752E12;
	}
	cout << "Number of alphas: " << number_of_alphas << endl;
	cout << "Alphas per second: " << number_of_alphas/(activation_end-activation_start)*1E12 << endl;
	cout << "Coulomb per nanosecond: " << number_of_alphas/(activation_end-activation_start)/current2alpha << endl;

	//histogramas
	auto rise_filter = [&](ULong64_t t){return t>=activation_start && t<=activation_end;};
	auto decay_filter = [&](ULong64_t t){return t>activation_end;};
	auto energy_window = d.Filter("Energy>activation_window_low && Energy<activation_window_high");
	auto labr_1_filter = energy_window.Filter("Channel==6").Define("t", "Timestamp/1E12");
	auto labr_2_filter = energy_window.Filter("Channel==7").Define("t", "Timestamp/1E12");

	auto labr_1 = labr_1_filter.Histo1D({"labr_1", ";Time (s);Counts", ACTIVATION_NBINS, 0, measurement_end/1E12}, "t");
	auto labr_2 = labr_2_filter.Histo1D({"labr_2", ";Time (s);Counts", ACTIVATION_NBINS, 0, measurement_end/1E12}, "t");

	int rise_nbins = (activation_end-activation_start)/(labr_1->GetBinWidth(1)*1E12);
	int decay_nbins = (measurement_end-activation_end)/(labr_1->GetBinWidth(1)*1E12);
	auto labr_1_rise = labr_1_filter.Filter(rise_filter, {"Timestamp"}).Histo1D({"labr_1_rise", ";Time (s);Counts", rise_nbins, activation_start/1E12, activation_end/1E12}, "t");
	auto labr_1_decay = labr_1_filter.Filter(decay_filter, {"Timestamp"}).Histo1D({"labr_1_decay", ";Time (s); Counts", decay_nbins, activation_end/1E12, measurement_end/1E12}, "t");
	auto labr_2_rise = labr_2_filter.Filter(rise_filter, {"Timestamp"}).Histo1D({"labr_2_rise", ";Time (s);Counts", rise_nbins, activation_start/1E12, activation_end/1E12}, "t");
	auto labr_2_decay = labr_2_filter.Filter(decay_filter, {"Timestamp"}).Histo1D({"labr_2_decay", ";Time (s);Counts", decay_nbins, activation_end/1E12, measurement_end/1E12}, "t");

	auto labr1_E_t_plot = d.Filter("Channel==6").Define("t", "Timestamp/1E12").Histo2D({"labr1_E_t_plot", ";Time (s);Energy (channel);Counts", 100, 0, measurement_end/1E12, 100, 0, 4096}, "t", "Energy");
	auto labr2_E_t_plot = d.Filter("Channel==7").Define("t", "Timestamp/1E12").Histo2D({"labr2_E_t_plot", ";Time (s);Energy (channel);Counts", 100, 0, measurement_end/1E12, 100, 0, 4096}, "t", "Energy");

	current_integrator->Write("", TObject::kOverwrite);
	labr_1->Write("", TObject::kOverwrite);
	labr_2->Write("", TObject::kOverwrite);
	labr_1_rise->Write("", TObject::kOverwrite);
	labr_1_decay->Write("", TObject::kOverwrite);
	labr_2_rise->Write("", TObject::kOverwrite);
	labr_2_decay->Write("", TObject::kOverwrite);
	TCanvas* myCanvas = new TCanvas("myCanvas");
	labr1_E_t_plot->Draw("COLZ");
	myCanvas->Write("labr1_E_t_plot", TObject::kOverwrite);
	labr2_E_t_plot->Draw("COLZ");
	myCanvas->Write("labr2_E_t_plot", TObject::kOverwrite);
	labr2_E_t_plot->Draw("COLZ");

	activation_start/=1E12;
	activation_end/=1E12;
	measurement_end/=1E12;

	//fittings
	TFitResultPtr fitresult;

	//unified
	TF1* unified = new TF1("unified_fit", unified_fit((TH1F*)gDirectory->Get("current_integrator")), 0, measurement_end, 5);
	unified->SetNpx(ACTIVATION_NBINS);
	unified->SetNumberFitPoints(ACTIVATION_NBINS);
	unified->SetParLimits(0, 0, 1E5);
	unified->SetParLimits(1, 0, 1E2);
	unified->SetParLimits(2, 4E-3, 5E-3);
	unified->SetParLimits(3, 0, 1E2);
	unified->SetParLimits(4, 0, 1E8);
	unified->SetParameters(15, 1, 4.62406E-3, 0, 0);
	unified->FixParameter(2, 4.62406E-3);
	unified->SetParNames("Background activity", "current to (a,n)", "Decay constant", "extra bg", "initial number of 30P");

	fitresult = labr_1->Fit("unified_fit", "SLEQ");
	results[0][0] = fitresult->Parameter(1)/current2alpha*labr_1->GetBinWidth(1);
	results[0][1] = fitresult->ParError(1)/current2alpha*labr_1->GetBinWidth(1);
	myCanvas->SetName("labr_1_unified_fit");
	myCanvas->Write("", TObject::kOverwrite);

	fitresult = labr_2->Fit("unified_fit", "SLEQ");
	results[1][0] = fitresult->Parameter(1)/current2alpha*labr_1->GetBinWidth(1);
	results[1][1] = fitresult->ParError(1)/current2alpha*labr_1->GetBinWidth(1);
	myCanvas->SetName("labr_2_unified_fit");
	myCanvas->Write("", TObject::kOverwrite);

	cout << "labr1 unified 30P per alpha: " << results[0][0] << endl;
	cout << "labr2 unified 30P per alpha: " << results[1][0] << endl;
	cout << "labr1 unified number of 30P: " << results[0][0]*number_of_alphas << endl;
	cout << "labr2 unified number of 30P: " << results[1][0]*number_of_alphas << endl;

	//rise
	unified->SetNpx(rise_nbins);

	fitresult = labr_1_rise->Fit("unified_fit", "SLEQ");
	results[0][2] = fitresult->Parameter(1)/current2alpha*labr_1->GetBinWidth(1);
	results[0][3] = fitresult->ParError(1)/current2alpha*labr_1->GetBinWidth(1);
	myCanvas->SetName("labr_1_rise_fit");
	myCanvas->Write("", TObject::kOverwrite);

	fitresult = labr_2_rise->Fit("unified_fit", "SLEQ");
	results[1][2] = fitresult->Parameter(1)/current2alpha*labr_1->GetBinWidth(1);
	results[1][3] = fitresult->ParError(1)/current2alpha*labr_1->GetBinWidth(1);
	myCanvas->SetName("labr_2_rise_fit");
	myCanvas->Write("", TObject::kOverwrite);

	cout << "labr1 rise 30P per alpha: " << results[0][2] << endl;
	cout << "labr2 rise 30P per alpha: " << results[1][2] << endl;
	cout << "labr1 rise number of 30P: " << results[0][2]*number_of_alphas << endl;
	cout << "labr2 rise number of 30P: " << results[1][2]*number_of_alphas << endl;

	//decay
	TF1* decay = new TF1("decay","[0]+[1]*exp(-[2]*(x[0]-[3]))");
	decay->SetNpx(decay_nbins);
	decay->SetNumberFitPoints(ACTIVATION_NBINS);
	decay->SetParLimits(0, 0, 100);
	decay->SetParLimits(2, 4E-3, 5E-3);
	decay->SetParameters(15, 1000, 4.62406E-3, activation_end);
	decay->FixParameter(3, activation_end);
	decay->FixParameter(2, 4.62406E-3);
	decay->SetParNames("Background activity", "Initial activiy", "Decay constant", "activation_end");

	fitresult = labr_1_decay->Fit("decay", "SLEQ");
	results[0][4] = fitresult->Parameter(1)*(activation_end-activation_start)/((1-exp(-fitresult->Parameter(2)*(activation_end-activation_start)))*number_of_alphas*labr_1->GetBinWidth(1));
	results[0][5] = fitresult->ParError(1)*(activation_end-activation_start)/((1-exp(-fitresult->Parameter(2)*(activation_end-activation_start)))*number_of_alphas*labr_1->GetBinWidth(1));
	myCanvas->SetName("labr_1_decay_fit");
	myCanvas->Write("", TObject::kOverwrite);

	fitresult = labr_2_decay->Fit("decay", "SLEQ");
	results[1][4] = fitresult->Parameter(1)*(activation_end-activation_start)/((1-exp(-fitresult->Parameter(2)*(activation_end-activation_start)))*number_of_alphas*labr_2->GetBinWidth(1));
	results[1][5] = fitresult->ParError(1)*(activation_end-activation_start)/((1-exp(-fitresult->Parameter(2)*(activation_end-activation_start)))*number_of_alphas*labr_2->GetBinWidth(1));
	myCanvas->SetName("labr_2_decay_fit");
	myCanvas->Write("", TObject::kOverwrite);

	cout << "labr1 decay 30P per alpha: " << results[0][4] << endl;
	cout << "labr2 decay 30P per alpha: " << results[1][4] << endl;
	cout << "labr1 decay number of 30P: " << results[0][4]*number_of_alphas << endl;
	cout << "labr2 decay number of 30P: " << results[1][4]*number_of_alphas << endl;

	myCanvas->Close();
	DisableImplicitMT();	//multithreading
	cout << "---------------" << endl;
}
