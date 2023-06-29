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

#define NA22_511_INTENSITY 1.798
#define NA22_CALIBRATION_ACTIVITY_NOMINAL 53590
#define NA22_CALIBRATION_ACTIVITY_PABLO 85973
#define NA22_ORIGINAL_ACTIVITY NA22_CALIBRATION_ACTIVITY_PABLO
#define CS137_INTENSITY 0.851
#define CS137_ORIGINAL_ACTIVITY 4360

#define PABLO_EFF_10CM 0.0027

#define MARKER_SIZE 2.5
#define ERROR_I 0	//0 unified, 2 rise, 4 decay

#include "peak_activity.cxx"

static void per_file(Char_t filepath[500], Double_t results[2][6]);

void activation_results(){
	//redirige cout a archivo de texto
	std::ofstream outputFile("resultados_activacion.txt");
	std::streambuf* originalBuffer = std::cout.rdbuf();
	std::cout.rdbuf(outputFile.rdbuf());

	Double_t na22_lambda = log(2)/(2.6018*365.2425*24*3600);
	Double_t na22_calibration_activity_feb = NA22_ORIGINAL_ACTIVITY * NA22_511_INTENSITY * exp(-na22_lambda*24*3600*108);
	Double_t na22_calibration_activity_apr = NA22_ORIGINAL_ACTIVITY * NA22_511_INTENSITY * exp(-na22_lambda*24*3600*162);
	Double_t cs137_lambda = log(2)/(30.07*365.2425*24*3600);
	Double_t cs137_calibration_activity = CS137_ORIGINAL_ACTIVITY * CS137_INTENSITY * exp(-cs137_lambda*24*3600*302);	//for apr

	cout << "Calcular emisiones reales de muestras de calibración como calibration activity * intensity * e^(-lambda*t):" << endl;
	cout << "Na22 intensity: " << NA22_511_INTENSITY << endl;
	cout << "Na22 lambda (s^-1): " << na22_lambda << endl;
	cout << "Na22 original activity (decays/s): " << NA22_ORIGINAL_ACTIVITY << endl;
	cout << "Na22 activity Feb (decays/s): " << na22_calibration_activity_feb << endl;
	cout << "Na22 activity Apr (decays/s): " << na22_calibration_activity_apr << endl;
	cout << "--------------" << endl;
	cout << "Cs137 intensity: " << CS137_INTENSITY << endl;
	cout << "Cs137 lambda (s^-1): " << cs137_lambda << endl;
	cout << "Cs137 original activity (decays/s): " << CS137_ORIGINAL_ACTIVITY << endl;
	cout << "Cs137 activity (decays/s): " << cs137_calibration_activity << endl;
	cout << endl;

	//calibración con na22
	cout << "Calcular cuentas/s medidas en runs de calibración:" << endl;
	//en logbook
	Double_t labr1_sodio_1 = peak_activity("output/SData_LaBr_Na22atTarget_calib_20230223.root", 6, 500, 700, "labr1_sodio_1");
	Double_t labr2_sodio_1 = peak_activity("output/SData_LaBr_Na22atTarget_calib_20230223.root", 7, 500, 700, "labr2_sodio_1");
	//en logbook
	Double_t labr1_sodio_2 = peak_activity("output/SData_LaBr_Na22atTarget_calib_20230418.root", 6, 1100, 1400, "labr1_sodio_2");
	Double_t labr2_sodio_2 = peak_activity("output/SData_LaBr_Na22atTarget_calib_20230418.root", 7, 1200, 1450, "labr2_sodio_2");
	//no en logbook, pone febrero pero probablemente abril
	Double_t labr1_sodio_3 = peak_activity("output/SData_LaBr1y2_Na22atTarget_calib_20230223.root", 6, 1100, 1400, "labr1_sodio_3");
	Double_t labr2_sodio_3 = peak_activity("output/SData_LaBr1y2_Na22atTarget_calib_20230223.root", 7, 1150, 1500, "labr2_sodio_3");
	cout << "labr1_sodio_1: " << labr1_sodio_1 << " counts/s" << endl;
	cout << "labr2_sodio_1: " << labr2_sodio_1 << " counts/s" << endl;
	cout << "labr1_sodio_2: " << labr1_sodio_2 << " counts/s" << endl;
	cout << "labr2_sodio_2: " << labr2_sodio_2 << " counts/s" << endl;
	cout << "labr1_sodio_3: " << labr1_sodio_3 << " counts/s" << endl;
	cout << "labr2_sodio_3: " << labr2_sodio_3 << " counts/s" << endl;
	cout << "--------------" << endl;
	Double_t labr1_cesio_1 = peak_activity("output/SData_LaBr_Cs137atTarget_calib_20230223.root", 6, 1600, 1800, "labr1_cesio_1");
	Double_t labr2_cesio_1 = peak_activity("output/SData_LaBr_Cs137atTarget_calib_20230223.root", 7, 1600, 1850, "labr2_cesio_1");
	Double_t labr1_cesio_2 = peak_activity("output/SData_LaBr_Cs137atTarget_calib_20230418.root", 6, 1450, 1750, "labr1_cesio_2");
	Double_t labr2_cesio_2 = peak_activity("output/SData_LaBr_Cs137atTarget_calib_20230418.root", 7, 1650, 1850, "labr2_cesio_2");
	cout << "labr1_cesio_1: " << labr1_cesio_1 << " counts/s" << endl;
	cout << "labr2_cesio_1: " << labr2_cesio_1 << " counts/s" << endl;
	cout << "labr1_cesio_2: " << labr1_cesio_2 << " counts/s" << endl;
	cout << "labr2_cesio_2: " << labr2_cesio_2 << " counts/s" << endl;
	cout << endl;

	TFile f("output.root", "UPDATE");
	gDirectory->cd("Activation");

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
	TGraph* reactionsvenergy_unified_feb_sin_escalar_1 = new TGraphErrors(3, x1, y1, NULL, yerr1);
	reactionsvenergy_unified_feb_sin_escalar_1->SetMarkerStyle(22);
	reactionsvenergy_unified_feb_sin_escalar_1->SetMarkerSize(MARKER_SIZE);
	TGraph* reactionsvenergy_unified_apr_sin_escalar_1 = new TGraphErrors(3, &x1[4], &y1[4], NULL, yerr1);
	reactionsvenergy_unified_apr_sin_escalar_1->SetMarkerStyle(23);
	reactionsvenergy_unified_apr_sin_escalar_1->SetMarkerSize(MARKER_SIZE);
	TGraph* reactionsvenergy_unified_feb_sin_escalar_2 = new TGraphErrors(3, x2, y2, NULL, yerr2);
	reactionsvenergy_unified_feb_sin_escalar_2->SetMarkerStyle(34);
	reactionsvenergy_unified_feb_sin_escalar_2->SetMarkerSize(MARKER_SIZE);
	TGraph* reactionsvenergy_unified_apr_sin_escalar_2 = new TGraphErrors(3, &x2[4], &y2[4], NULL, yerr2);
	reactionsvenergy_unified_apr_sin_escalar_2->SetMarkerStyle(47);
	reactionsvenergy_unified_apr_sin_escalar_2->SetMarkerSize(MARKER_SIZE);

	//rise fit sin escalar
	for(short i=0; i<ACTIVATION_N; i++){
		x1[i] = activation_energies[i];
		x2[i] = activation_energies[i];
		y1[i] = results[i][0][2];
		y2[i] = results[i][1][2];
		yerr1[i] = results[i][0][3];
		yerr2[i] = results[i][1][3];
	}
	TGraph* reactionsvenergy_rise_feb_sin_escalar_1 = new TGraphErrors(3, x1, y1, NULL, yerr1);
	reactionsvenergy_rise_feb_sin_escalar_1->SetMarkerStyle(22);
	reactionsvenergy_rise_feb_sin_escalar_1->SetMarkerSize(MARKER_SIZE);
	TGraph* reactionsvenergy_rise_apr_sin_escalar_1 = new TGraphErrors(3, &x1[4], &y1[4], NULL, yerr1);
	reactionsvenergy_rise_apr_sin_escalar_1->SetMarkerStyle(23);
	reactionsvenergy_rise_apr_sin_escalar_1->SetMarkerSize(MARKER_SIZE);
	TGraph* reactionsvenergy_rise_feb_sin_escalar_2 = new TGraphErrors(3, x2, y2, NULL, yerr2);
	reactionsvenergy_rise_feb_sin_escalar_2->SetMarkerStyle(34);
	reactionsvenergy_rise_feb_sin_escalar_2->SetMarkerSize(MARKER_SIZE);
	TGraph* reactionsvenergy_rise_apr_sin_escalar_2 = new TGraphErrors(3, &x2[4], &y2[4], NULL, yerr2);
	reactionsvenergy_rise_apr_sin_escalar_2->SetMarkerStyle(47);
	reactionsvenergy_rise_apr_sin_escalar_2->SetMarkerSize(MARKER_SIZE);

	//decay fit sin escalar
	for(short i=0; i<ACTIVATION_N; i++){
		x1[i] = activation_energies[i];
		x2[i] = activation_energies[i];
		y1[i] = results[i][0][4];
		y2[i] = results[i][1][4];
		yerr1[i] = results[i][0][5];
		yerr2[i] = results[i][1][5];
	}
	TGraph* reactionsvenergy_decay_feb_sin_escalar_1 = new TGraphErrors(4, x1, y1, NULL, yerr1);
	reactionsvenergy_decay_feb_sin_escalar_1->SetMarkerStyle(22);
	reactionsvenergy_decay_feb_sin_escalar_1->SetMarkerSize(MARKER_SIZE);
	TGraph* reactionsvenergy_decay_apr_sin_escalar_1 = new TGraphErrors(6, &x1[4], &y1[4], NULL, yerr1);
	reactionsvenergy_decay_apr_sin_escalar_1->SetMarkerStyle(23);
	reactionsvenergy_decay_apr_sin_escalar_1->SetMarkerSize(MARKER_SIZE);
	TGraph* reactionsvenergy_decay_feb_sin_escalar_2 = new TGraphErrors(4, x2, y2, NULL, yerr2);
	reactionsvenergy_decay_feb_sin_escalar_2->SetMarkerStyle(34);
	reactionsvenergy_decay_feb_sin_escalar_2->SetMarkerSize(MARKER_SIZE);
	TGraph* reactionsvenergy_decay_apr_sin_escalar_2 = new TGraphErrors(6, &x2[4], &y2[4], NULL, yerr2);
	reactionsvenergy_decay_apr_sin_escalar_2->SetMarkerStyle(47);
	reactionsvenergy_decay_apr_sin_escalar_2->SetMarkerSize(MARKER_SIZE);

	//todos los resultados antes de los escalados
	TMultiGraph* multigraph_sin_escalar = new TMultiGraph();
	reactionsvenergy_unified_feb_sin_escalar_1->SetMarkerColor(kRed);
	reactionsvenergy_unified_feb_sin_escalar_1->SetTitle("Unified fit, February, LaBr1");
	multigraph_sin_escalar->Add(reactionsvenergy_unified_feb_sin_escalar_1);
	reactionsvenergy_unified_feb_sin_escalar_2->SetMarkerColor(kRed);
	reactionsvenergy_unified_feb_sin_escalar_2->SetTitle("Unified fit, February, LaBr2");
	multigraph_sin_escalar->Add(reactionsvenergy_unified_feb_sin_escalar_2);
	reactionsvenergy_unified_apr_sin_escalar_1->SetMarkerColor(kRed);
	reactionsvenergy_unified_apr_sin_escalar_1->SetTitle("Unified fit, April, LaBr1");
	multigraph_sin_escalar->Add(reactionsvenergy_unified_apr_sin_escalar_1);
	reactionsvenergy_unified_apr_sin_escalar_2->SetMarkerColor(kRed);
	reactionsvenergy_unified_apr_sin_escalar_2->SetTitle("Unified fit, April, LaBr2");
	multigraph_sin_escalar->Add(reactionsvenergy_unified_apr_sin_escalar_2);
	reactionsvenergy_rise_feb_sin_escalar_1->SetMarkerColor(kGreen);
	reactionsvenergy_rise_feb_sin_escalar_1->SetTitle("Rise fit, February, LaBr1");
	multigraph_sin_escalar->Add(reactionsvenergy_rise_feb_sin_escalar_1);
	reactionsvenergy_rise_feb_sin_escalar_2->SetMarkerColor(kGreen);
	reactionsvenergy_rise_feb_sin_escalar_2->SetTitle("Rise fit, February, LaBr2");
	multigraph_sin_escalar->Add(reactionsvenergy_rise_feb_sin_escalar_2);
	reactionsvenergy_rise_apr_sin_escalar_1->SetMarkerColor(kGreen);
	reactionsvenergy_rise_apr_sin_escalar_1->SetTitle("Rise fit, April, LaBr1");
	multigraph_sin_escalar->Add(reactionsvenergy_rise_apr_sin_escalar_1);
	reactionsvenergy_rise_apr_sin_escalar_2->SetMarkerColor(kGreen);
	reactionsvenergy_rise_apr_sin_escalar_2->SetTitle("Rise fit, April, LaBr2");
	multigraph_sin_escalar->Add(reactionsvenergy_rise_apr_sin_escalar_2);
	reactionsvenergy_decay_feb_sin_escalar_1->SetMarkerColor(kBlue);
	reactionsvenergy_decay_feb_sin_escalar_1->SetTitle("Decay fit, February, LaBr1");
	multigraph_sin_escalar->Add(reactionsvenergy_decay_feb_sin_escalar_1);
	reactionsvenergy_decay_feb_sin_escalar_2->SetMarkerColor(kBlue);
	reactionsvenergy_decay_feb_sin_escalar_2->SetTitle("Decay fit, February, LaBr2");
	multigraph_sin_escalar->Add(reactionsvenergy_decay_feb_sin_escalar_2);
	reactionsvenergy_decay_apr_sin_escalar_1->SetMarkerColor(kBlue);
	reactionsvenergy_decay_apr_sin_escalar_1->SetTitle("Decay fit, April, LaBr1");
	multigraph_sin_escalar->Add(reactionsvenergy_decay_apr_sin_escalar_1);
	reactionsvenergy_decay_apr_sin_escalar_2->SetMarkerColor(kBlue);
	reactionsvenergy_decay_apr_sin_escalar_2->SetTitle("Decay fit, April, LaBr2");
	multigraph_sin_escalar->Add(reactionsvenergy_decay_apr_sin_escalar_2);
	multigraph_sin_escalar->SetTitle("30P v a energy;Energy of a (keV);Thick target (a,n) yield");
	multigraph_sin_escalar->Draw("AP");
	myCanvas->BuildLegend();
	myCanvas->Write("", TObject::kOverwrite);

	cout << "Calcular eficiencias como cuentas/s en fotopico durante run de calibración dividido por emisiones de la muestra:" << endl;

	Double_t eficiencia_feb_labr1 = labr1_sodio_1 / na22_calibration_activity_feb;
	Double_t eficiencia_feb_labr2 = labr2_sodio_1 / na22_calibration_activity_feb;
	Double_t eficiencia_apr2_labr1 = labr1_sodio_2 / na22_calibration_activity_apr;
	Double_t eficiencia_apr2_labr2 = labr2_sodio_2 / na22_calibration_activity_apr;
	Double_t eficiencia_apr3_labr1 = labr1_sodio_3 / na22_calibration_activity_apr;
	Double_t eficiencia_apr3_labr2 = labr2_sodio_3 / na22_calibration_activity_apr;
	Double_t eficiencia_cesio_feb_labr1 = labr1_cesio_1 / cs137_calibration_activity;
	Double_t eficiencia_cesio_feb_labr2 = labr2_cesio_1 / cs137_calibration_activity;
	Double_t eficiencia_cesio_apr_labr1 = labr1_cesio_2 / cs137_calibration_activity;
	Double_t eficiencia_cesio_apr_labr2 = labr2_cesio_2 / cs137_calibration_activity;

	cout << "Eficiencia sodio labr1 febrero: " << eficiencia_feb_labr1*100 << "%" << endl;
	cout << "Distancia equivalente a Pablo: " << sqrt(100*PABLO_EFF_10CM/eficiencia_feb_labr1) << "cm" << endl;
	cout << "Eficiencia sodio labr2 febrero: " << eficiencia_feb_labr2*100 << "%" << endl;
	cout << "Distancia equivalente a Pablo: " << sqrt(100*PABLO_EFF_10CM/eficiencia_feb_labr2) << "cm" << endl;
	cout << "Eficiencia sodio labr1 abril no logbook: " << eficiencia_apr2_labr1*100 << "%" << endl;
	cout << "Distancia equivalente a Pablo: " << sqrt(100*PABLO_EFF_10CM/eficiencia_apr2_labr1) << "cm" << endl;
	cout << "Eficiencia sodio labr2 abril no logbook: " << eficiencia_apr2_labr2*100 << "%" << endl;
	cout << "Distancia equivalente a Pablo: " << sqrt(100*PABLO_EFF_10CM/eficiencia_apr2_labr2) << "cm" << endl;
	cout << "Eficiencia sodio labr1 abril sí logbook: " << eficiencia_apr3_labr1*100 << "%" << endl;
	cout << "Distancia equivalente a Pablo: " << sqrt(100*PABLO_EFF_10CM/eficiencia_apr3_labr1) << "cm" << endl;
	cout << "Eficiencia sodio labr2 abril sí logbook: " << eficiencia_apr3_labr2*100 << "%" << endl;
	cout << "Distancia equivalente a Pablo: " << sqrt(100*PABLO_EFF_10CM/eficiencia_apr3_labr2) << "cm" << endl;
	cout << "--------------" << endl;
	cout << "Eficiencia cesio labr1 febrero: " << eficiencia_cesio_feb_labr1*100 << "%" << endl;
	cout << "Eficiencia cesio labr2 febrero: " << eficiencia_cesio_feb_labr2*100 << "%" << endl;
	cout << "Eficiencia cesio labr1 abril: " << eficiencia_cesio_apr_labr1*100 << "%" << endl;
	cout << "Eficiencia cesio labr2 abril: " << eficiencia_cesio_apr_labr2*100 << "%" << endl;
	cout << endl;

	//escalados
	for(unsigned short j=0; j<6; j++){
		//escalado debido a la mala medida de la carga
		results[3][0][j]*=172/239.1;
		results[3][1][j]*=172/239.1;

		//escalado na22
		results[0][0][j] /= eficiencia_feb_labr1;
		results[0][1][j] /= eficiencia_feb_labr2;
		results[1][0][j] /= eficiencia_feb_labr1;
		results[1][1][j] /= eficiencia_feb_labr2;
		results[2][0][j] /= eficiencia_feb_labr1;
		results[2][1][j] /= eficiencia_feb_labr2;
		results[3][0][j] /= eficiencia_feb_labr1;
		results[3][1][j] /= eficiencia_feb_labr2;
		results[4][0][j] /= eficiencia_apr3_labr1;
		results[4][1][j] /= eficiencia_apr3_labr2;
		results[5][0][j] /= eficiencia_apr2_labr1;
		results[5][1][j] /= eficiencia_apr2_labr2;
		results[6][0][j] /= eficiencia_apr3_labr1;
		results[6][1][j] /= eficiencia_apr3_labr2;
		results[7][0][j] /= eficiencia_apr3_labr1;
		results[7][1][j] /= eficiencia_apr3_labr2;
		results[8][0][j] /= eficiencia_apr3_labr1;
		results[8][1][j] /= eficiencia_apr3_labr2;
		results[9][0][j] /= eficiencia_apr3_labr1;
		results[9][1][j] /= eficiencia_apr3_labr2;
	}

	Double_t exfor_energies[] = {3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500, 5600, 5700, 5800, 5900, 6000, 6100, 6200, 6300, 6400, 6500, 6600, 6700, 6800, 6900, 7000, 7100, 7200, 7300, 7400, 7500, 7600, 7700, 7800, 7900, 8000, 8100, 8200, 8300, 8400, 8500, 8600, 8700, 8800, 8900, 9000, 9100, 9200, 9300, 9400, 9500, 9600, 9700, 9800, 9900};	//TBD:hardcoded, read .txt
	Double_t exfor_data[] = {3.147E-09, 5.904E-09, 1.034E-08, 1.655E-08, 2.462E-08, 3.464E-08, 4.686E-08, 6.203E-08, 8.124E-08, 1.073E-07, 1.426E-07, 1.847E-07, 2.306E-07, 2.812E-07, 3.403E-07, 4.150E-07, 5.119E-07, 6.278E-07, 7.555E-07, 8.856E-07, 1.011E-06, 1.150E-06, 1.330E-06, 1.549E-06, 1.797E-06, 2.062E-06, 2.339E-06, 2.651E-06, 3.015E-06, 3.401E-06, 3.774E-06, 4.147E-06, 4.552E-06, 4.999E-06, 5.489E-06, 6.013E-06, 6.562E-06, 7.131E-06, 7.716E-06, 8.319E-06, 8.943E-06, 9.593E-06, 1.027E-05, 1.099E-05, 1.173E-05, 1.252E-05, 1.333E-05, 1.416E-05, 1.502E-05, 1.589E-05, 1.679E-05, 1.771E-05, 1.865E-05, 1.962E-05, 2.062E-05, 2.165E-05, 2.272E-05, 2.383E-05, 2.497E-05, 2.616E-05, 2.740E-05, 2.869E-05, 3.003E-05};	//TBD:hardcoded, read .txt
	Double_t exfor_data_1[] = {3.147E-09, 5.904E-09, 1.034E-08, 1.655E-08, 2.462E-08, 3.464E-08, 4.686E-08, 6.203E-08, 8.124E-08, 1.073E-07, 1.426E-07, 1.847E-07, 2.306E-07, 2.812E-07, 3.403E-07, 4.150E-07, 5.119E-07, 6.278E-07, 7.555E-07, 8.856E-07, 1.011E-06, 1.150E-06, 1.330E-06, 1.549E-06, 1.797E-06, 2.062E-06, 2.339E-06, 2.651E-06, 3.015E-06, 3.401E-06, 3.774E-06, 4.147E-06, 4.552E-06, 4.999E-06, 5.489E-06, 6.013E-06, 6.562E-06, 7.131E-06, 7.716E-06, 8.319E-06, 8.943E-06, 9.593E-06, 1.027E-05, 1.099E-05, 1.173E-05, 1.252E-05, 1.333E-05, 1.416E-05, 1.502E-05, 1.589E-05, 1.679E-05, 1.771E-05, 1.865E-05, 1.962E-05, 2.062E-05, 2.165E-05, 2.272E-05, 2.383E-05, 2.497E-05, 2.616E-05, 2.740E-05, 2.869E-05, 3.003E-05};	//TBD:hardcoded, read .txt
	Double_t exfor_data_2[] = {3.147E-09, 5.904E-09, 1.034E-08, 1.655E-08, 2.462E-08, 3.464E-08, 4.686E-08, 6.203E-08, 8.124E-08, 1.073E-07, 1.426E-07, 1.847E-07, 2.306E-07, 2.812E-07, 3.403E-07, 4.150E-07, 5.119E-07, 6.278E-07, 7.555E-07, 8.856E-07, 1.011E-06, 1.150E-06, 1.330E-06, 1.549E-06, 1.797E-06, 2.062E-06, 2.339E-06, 2.651E-06, 3.015E-06, 3.401E-06, 3.774E-06, 4.147E-06, 4.552E-06, 4.999E-06, 5.489E-06, 6.013E-06, 6.562E-06, 7.131E-06, 7.716E-06, 8.319E-06, 8.943E-06, 9.593E-06, 1.027E-05, 1.099E-05, 1.173E-05, 1.252E-05, 1.333E-05, 1.416E-05, 1.502E-05, 1.589E-05, 1.679E-05, 1.771E-05, 1.865E-05, 1.962E-05, 2.062E-05, 2.165E-05, 2.272E-05, 2.383E-05, 2.497E-05, 2.616E-05, 2.740E-05, 2.869E-05, 3.003E-05};	//TBD:hardcoded, read .txt
	Double_t exfor_errors[63];
	Double_t exfor_errors_1[63];
	Double_t exfor_errors_2[63];
	Double_t factor_feb = (results[0][1][4]+results[0][0][4])/(2*exfor_data[18]);
	Double_t error_feb = abs(results[0][1][4]-results[0][0][4])/(exfor_data[18]*factor_feb);
	Double_t factor_apr = (results[8][1][4]+results[8][0][4])/(2*exfor_data[18]);
	Double_t error_apr = abs(results[8][1][4]-results[8][0][4])/(exfor_data[18]*factor_apr);
	cout << "Factor febrero: " << factor_feb << endl;
	cout << "Error febrero: " << error_feb*100 << "%" << endl;
	cout << "Factor abril: " << factor_apr << endl;
	cout << "Error abril: " << error_apr*100 << "%" << endl;
	cout << "Eficiencia Febrero LaBr1 para acuerdo: " << eficiencia_feb_labr1*factor_feb*100 << "%" << endl;
	cout << "Distancia equivalente a Pablo: " << sqrt(100*PABLO_EFF_10CM/(eficiencia_feb_labr1*factor_feb)) << "cm" << endl;
	cout << "Eficiencia Abril3 LaBr2 para acuerdo: " << eficiencia_apr3_labr2*factor_apr*100 << "%" << endl;
	cout << "Distancia equivalente a Pablo: " << sqrt(100*PABLO_EFF_10CM/(eficiencia_apr2_labr2*factor_apr)) << "cm" << endl;
	for(short i=0; i<63; i++){	//TBD:escalado temporal, números hardcoded
		exfor_errors[i] = exfor_data[i] * 0.05;
//		exfor_data_1[i]*=2.220341468E-5*(19.5/16.8)*(19.358/19.5);	//eficiencia feb
		exfor_data_1[i]*=factor_feb;
		exfor_errors_1[i] = exfor_data_1[i] * error_feb;
//		exfor_data_2[i]*=2.642141092E-5*(22.5/20)*(21.847/22.3)*(21.934/21.5958)*(21.3124/22.1447);	//eficiencia apr
		exfor_data_2[i]*=factor_apr;
		exfor_errors_2[i] = exfor_data_2[i] * error_apr;
	}

	//EXFOR data
	TGraph* reactionsvenergy_exfor = new TGraphErrors(63, exfor_energies, exfor_data, NULL, exfor_errors);	//TBD:hardcoded number
	reactionsvenergy_exfor->SetTitle("EXFOR yield data;Energy of a (keV);Thick target (a,n) yield");
	reactionsvenergy_exfor->SetFillColorAlpha(kGray+2, 0.2);
	reactionsvenergy_exfor->SetMarkerStyle(20);
	reactionsvenergy_exfor->SetLineColor(kBlack);
	TGraph* reactionsvenergy_exfor_feb = new TGraphErrors(63, exfor_energies, exfor_data_1, NULL, exfor_errors_1);	//TBD:hardcoded number
	reactionsvenergy_exfor_feb->SetTitle("EXFOR yield data;Energy of a (keV);Thick target (a,n)");
	reactionsvenergy_exfor_feb->SetFillColorAlpha(kGray, 0.4);
	reactionsvenergy_exfor_feb->SetMarkerStyle(20);
	reactionsvenergy_exfor_feb->SetLineColor(kBlack);
	TGraph* reactionsvenergy_exfor_apr = new TGraphErrors(63, exfor_energies, exfor_data_2, NULL, exfor_errors_2);	//TBD:hardcoded number
	reactionsvenergy_exfor_apr->SetTitle("EXFOR yield data;Energy of a (keV);Thick target (a,n)");
	reactionsvenergy_exfor_apr->SetFillColorAlpha(kGray+1, 0.4);
	reactionsvenergy_exfor_apr->SetMarkerStyle(20);
	reactionsvenergy_exfor_apr->SetLineColor(kBlack);

	//medias
	Double_t decay_average[ACTIVATION_N];
	Double_t unified_average[ACTIVATION_N];
	Double_t rise_average[ACTIVATION_N];
	Double_t decay_average_diff[ACTIVATION_N];
	Double_t unified_average_diff[ACTIVATION_N];
	Double_t rise_average_diff[ACTIVATION_N];
	Double_t decay_average_diff_per[ACTIVATION_N];
	Double_t unified_average_diff_per[ACTIVATION_N];
	Double_t rise_average_diff_per[ACTIVATION_N];

	//unified_fit
	cout << "-----UNIFIED-----" << endl;
	for(short i=0; i<ACTIVATION_N; i++){
		x1[i] = activation_energies[i];
		x2[i] = activation_energies[i];
		y1[i] = results[i][0][0];
		y2[i] = results[i][1][0];
		yerr1[i] = results[i][0][1];
		yerr2[i] = results[i][1][1];
		cout << "Activation N " << i+1 << endl;
		cout << "E: " << x1[i] << endl;
		cout << "detector 1: " << y1[i] << endl;
		cout << "detector 2: " << y2[i] << endl;
		unified_average[i] = (y1[i]+y2[i])/2;
		unified_average_diff[i] = abs(y1[i]-y2[i])/2;
		unified_average_diff_per[i] = unified_average_diff[i]/unified_average[i]*100;
	}
	cout << "----- ----- -----" << endl << endl;
	TGraph* reactionsvenergy_unified_feb_1 = new TGraphErrors(3, x1, y1, NULL, yerr1);
	reactionsvenergy_unified_feb_1->SetTitle("(a,n) reactions v a energy (unified fit, February);Energy of a (keV);Thick target (a,n) yield");
	myCanvas->SetName("reactions_v_energy_unified_feb_1");
	reactionsvenergy_unified_feb_1->SetMarkerStyle(22);
	reactionsvenergy_unified_feb_1->SetMarkerSize(MARKER_SIZE);
	reactionsvenergy_unified_feb_1->SetMarkerColor(kRed);
//	reactionsvenergy_exfor_feb->Draw("E3");
//	reactionsvenergy_exfor_feb->Draw("same LX");
	reactionsvenergy_exfor->Draw("same E3");
//	reactionsvenergy_exfor->Draw("same LX");
	reactionsvenergy_unified_feb_1->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TGraph* reactionsvenergy_unified_feb_2 = new TGraphErrors(3, x2, y2, NULL, yerr2);
	reactionsvenergy_unified_feb_2->SetTitle("(a,n) reactions v a energy (unified fit, February);Energy of a (keV);Thick target (a,n) yield");
	myCanvas->SetName("reactions_v_energy_unified_feb_2");
	reactionsvenergy_unified_feb_2->SetMarkerStyle(34);
	reactionsvenergy_unified_feb_2->SetMarkerSize(MARKER_SIZE);
	reactionsvenergy_unified_feb_2->SetMarkerColor(kRed);
//	reactionsvenergy_exfor_feb->Draw("E3");
//	reactionsvenergy_exfor_feb->Draw("same LX");
	reactionsvenergy_exfor->Draw("same E3");
//	reactionsvenergy_exfor->Draw("same LX");
	reactionsvenergy_unified_feb_2->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TGraph* reactionsvenergy_unified_apr_1 = new TGraphErrors(3, &x1[4], &y1[4], NULL, yerr1);
	reactionsvenergy_unified_apr_1->SetTitle("(a,n) reactions v a energy (unified fit, April);Energy of a (keV);Thick target (a,n) yield");
	myCanvas->SetName("reactions_v_energy_unified_apr_1");
	reactionsvenergy_unified_apr_1->SetMarkerStyle(23);
	reactionsvenergy_unified_apr_1->SetMarkerSize(MARKER_SIZE);
	reactionsvenergy_unified_apr_1->SetMarkerColor(kBlue);
//	reactionsvenergy_exfor_apr->Draw("E3");
//	reactionsvenergy_exfor_apr->Draw("same LX");
	reactionsvenergy_exfor->Draw("same E3");
//	reactionsvenergy_exfor->Draw("same LX");
	reactionsvenergy_unified_apr_1->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TGraph* reactionsvenergy_unified_apr_2 = new TGraphErrors(3, &x2[4], &y2[4], NULL, yerr2);
	reactionsvenergy_unified_apr_2->SetTitle("(a,n) reactions v a energy (unified fit, April);Energy of a (keV);Thick target (a,n) yield");
	myCanvas->SetName("reactions_v_energy_unified_apr_2");
	reactionsvenergy_unified_apr_2->SetMarkerStyle(47);
	reactionsvenergy_unified_apr_2->SetMarkerSize(MARKER_SIZE);
	reactionsvenergy_unified_apr_2->SetMarkerColor(kBlue);
//	reactionsvenergy_exfor_apr->Draw("E3");
//	reactionsvenergy_exfor_apr->Draw("same LX");
	reactionsvenergy_exfor->Draw("same E3");
//	reactionsvenergy_exfor->Draw("same LX");
	reactionsvenergy_unified_apr_2->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TMultiGraph* reactions_v_energy_unified = new TMultiGraph();
//	reactions_v_energy_unified->Add(reactionsvenergy_exfor_feb, "E3");
//	reactions_v_energy_unified->Add(reactionsvenergy_exfor_feb, "LX");
	reactionsvenergy_exfor_feb->SetTitle("EXFOR data, feb");
//	reactions_v_energy_unified->Add(reactionsvenergy_exfor_apr, "E3");
//	reactions_v_energy_unified->Add(reactionsvenergy_exfor_apr, "LX");
	reactionsvenergy_exfor_apr->SetTitle("EXFOR data, apr");
	reactions_v_energy_unified->Add(reactionsvenergy_exfor, "E3");
//	reactions_v_energy_unified->Add(reactionsvenergy_exfor, "LX");
	reactionsvenergy_exfor_apr->SetTitle("EXFOR data");
	reactionsvenergy_unified_feb_1->SetTitle("Unified fit, February, LaBr1");
	reactions_v_energy_unified->Add(reactionsvenergy_unified_feb_1);
	reactionsvenergy_unified_apr_1->SetTitle("Unified fit, April, LaBr1");
	reactions_v_energy_unified->Add(reactionsvenergy_unified_apr_1);
	reactionsvenergy_unified_feb_2->SetTitle("Unified fit, February, LaBr2");
	reactions_v_energy_unified->Add(reactionsvenergy_unified_feb_2);
	reactionsvenergy_unified_apr_2->SetTitle("Unified fit, April, LaBr2");
	reactions_v_energy_unified->Add(reactionsvenergy_unified_apr_2);
	reactions_v_energy_unified->SetTitle("(a,n) reactions v a energy (unified fit);Energy of a (keV);Thick target (a,n) yield");
	reactions_v_energy_unified->Draw("AP");
	myCanvas->BuildLegend();
	myCanvas->SetName("reactions_v_energy_unified");
	myCanvas->Write("", TObject::kOverwrite);

	//rise fit
	cout << "-----RISE-----" << endl;
	for(short i=0; i<ACTIVATION_N; i++){
		x1[i] = activation_energies[i];
		x2[i] = activation_energies[i];
		y1[i] = results[i][0][2];
		y2[i] = results[i][1][2];
		yerr1[i] = results[i][0][3];
		yerr2[i] = results[i][1][3];
		cout << "Activation N " << i+1 << endl;
		cout << "E: " << x1[i] << endl;
		cout << "detector 1: " << y1[i] << endl;
		cout << "detector 2: " << y2[i] << endl;
		rise_average[i] = (y1[i]+y2[i])/2;
		rise_average_diff[i] = abs(y1[i]-y2[i])/2;
		rise_average_diff_per[i] = rise_average_diff[i]/rise_average[i]*100;
	}
	cout << "----- ----- -----" << endl << endl;
	TGraph* reactionsvenergy_rise_feb_1 = new TGraphErrors(3, x1, y1, NULL, yerr1);
	reactionsvenergy_rise_feb_1->SetTitle("(a,n) reactions v a energy (rise fit, February);Energy of a (keV);Thick target (a,n) yield");
	reactionsvenergy_rise_feb_1->SetMarkerStyle(22);
	reactionsvenergy_rise_feb_1->SetMarkerSize(MARKER_SIZE);
	reactionsvenergy_rise_feb_1->SetMarkerColor(kRed);
	myCanvas->SetName("reactions_v_energy_rise_feb_1");
//	reactionsvenergy_exfor_feb->Draw("E3");
//	reactionsvenergy_exfor_feb->Draw("same LX");
	reactionsvenergy_exfor->Draw("same E3");
//	reactionsvenergy_exfor->Draw("same LX");
	reactionsvenergy_rise_feb_1->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TGraph* reactionsvenergy_rise_feb_2 = new TGraphErrors(3, x2, y2, NULL, yerr2);
	reactionsvenergy_rise_feb_2->SetTitle("(a,n) reactions v a energy (rise fit, February);Energy of a (keV);Thick target (a,n) yield");
	reactionsvenergy_rise_feb_2->SetMarkerStyle(34);
	reactionsvenergy_rise_feb_2->SetMarkerSize(MARKER_SIZE);
	reactionsvenergy_rise_feb_2->SetMarkerColor(kRed);
	myCanvas->SetName("reactions_v_energy_rise_feb_2");
//	reactionsvenergy_exfor_feb->Draw("E3");
//	reactionsvenergy_exfor_feb->Draw("same LX");
	reactionsvenergy_exfor->Draw("same E3");
//	reactionsvenergy_exfor->Draw("same LX");
	reactionsvenergy_rise_feb_2->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TGraph* reactionsvenergy_rise_apr_1 = new TGraphErrors(2, &x1[4], &y1[4], NULL, yerr1);
	reactionsvenergy_rise_apr_1->SetTitle("(a,n) reactions v a energy (rise fit, April);Energy of a (keV);Thick target (a,n) yield");
	reactionsvenergy_rise_apr_1->SetMarkerStyle(23);
	reactionsvenergy_rise_apr_1->SetMarkerSize(MARKER_SIZE);
	reactionsvenergy_rise_apr_1->SetMarkerColor(kBlue);
	myCanvas->SetName("reactions_v_energy_rise_apr_1");
//	reactionsvenergy_exfor_apr->Draw("E3");
//	reactionsvenergy_exfor_apr->Draw("same LX");
	reactionsvenergy_exfor->Draw("same E3");
//	reactionsvenergy_exfor->Draw("same LX");
	reactionsvenergy_rise_apr_1->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TGraph* reactionsvenergy_rise_apr_2 = new TGraphErrors(2, &x2[4], &y2[4], NULL, yerr2);
	reactionsvenergy_rise_apr_2->SetTitle("(a,n) reactions v a energy (rise fit, April);Energy of a (keV);Thick target (a,n) yield");
	reactionsvenergy_rise_apr_2->SetMarkerStyle(47);
	reactionsvenergy_rise_apr_2->SetMarkerSize(MARKER_SIZE);
	reactionsvenergy_rise_apr_2->SetMarkerColor(kBlue);
	myCanvas->SetName("reactions_v_energy_rise_apr_2");
//	reactionsvenergy_exfor_apr->Draw("E3");
//	reactionsvenergy_exfor_apr->Draw("same LX");
	reactionsvenergy_exfor->Draw("same E3");
//	reactionsvenergy_exfor->Draw("same LX");
	reactionsvenergy_rise_apr_2->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TMultiGraph* reactions_v_energy_rise = new TMultiGraph();
//	reactions_v_energy_rise->Add(reactionsvenergy_exfor_feb, "E3");
//	reactions_v_energy_rise->Add(reactionsvenergy_exfor_feb, "LX");
//	reactions_v_energy_rise->Add(reactionsvenergy_exfor_apr, "E3");
//	reactions_v_energy_rise->Add(reactionsvenergy_exfor_apr, "LX");
	reactions_v_energy_rise->Add(reactionsvenergy_exfor, "E3");
//	reactions_v_energy_rise->Add(reactionsvenergy_exfor, "LX");
	reactionsvenergy_rise_feb_1->SetTitle("Rise fit, February, LaBr1");
	reactions_v_energy_rise->Add(reactionsvenergy_rise_feb_1);
	reactionsvenergy_rise_apr_1->SetTitle("Rise fit, April, LaBr1");
	reactions_v_energy_rise->Add(reactionsvenergy_rise_apr_1);
	reactionsvenergy_rise_feb_2->SetTitle("Rise fit, February, LaBr2");
	reactions_v_energy_rise->Add(reactionsvenergy_rise_feb_2);
	reactionsvenergy_rise_apr_2->SetTitle("Rise fit, April, LaBr2");
	reactions_v_energy_rise->Add(reactionsvenergy_rise_apr_2);
	reactions_v_energy_rise->SetTitle("(a,n) reactions v a energy (rise fit);Energy of a (keV);Thick target (a,n) yield");
	reactions_v_energy_rise->Draw("AP");
	myCanvas->BuildLegend();
	myCanvas->SetName("reactions_v_energy_rise");
	myCanvas->Write("", TObject::kOverwrite);

	//decay fit
	cout << "-----DECAY-----" << endl;
	for(short i=0; i<ACTIVATION_N; i++){
		x1[i] = activation_energies[i];
		x2[i] = activation_energies[i];
		y1[i] = results[i][0][4];
		y2[i] = results[i][1][4];
		yerr1[i] = results[i][0][5];
		yerr2[i] = results[i][1][5];
		cout << "Activation N " << i+1 << endl;
		cout << "E: " << x1[i] << endl;
		cout << "detector 1: " << y1[i] << endl;
		cout << "detector 2: " << y2[i] << endl;
		decay_average[i] = (y1[i]+y2[i])/2;
		decay_average_diff[i] = abs(y1[i]-y2[i])/2;
		decay_average_diff_per[i] = decay_average_diff[i]/decay_average[i]*100;
	}
	cout << "----- ----- -----" << endl << endl;
	TGraph* reactionsvenergy_decay_feb_1 = new TGraphErrors(4, x1, y1, NULL, yerr1);
	reactionsvenergy_decay_feb_1->SetTitle("(a,n) reactions v a energy (decay fit, February);Energy of a (keV);Thick target (a,n) yield");
	reactionsvenergy_decay_feb_1->SetMarkerStyle(22);
	reactionsvenergy_decay_feb_1->SetMarkerSize(MARKER_SIZE);
	reactionsvenergy_decay_feb_1->SetMarkerColor(kRed);
	myCanvas->SetName("reactions_v_energy_decay_feb_1");
//	reactionsvenergy_exfor_feb->Draw("E3");
//	reactionsvenergy_exfor_feb->Draw("same LX");
	reactionsvenergy_exfor->Draw("same E3");
//	reactionsvenergy_exfor->Draw("same LX");
	reactionsvenergy_decay_feb_1->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TGraph* reactionsvenergy_decay_feb_2 = new TGraphErrors(4, x2, y2, NULL, yerr2);
	reactionsvenergy_decay_feb_2->SetTitle("(a,n) reactions v a energy (decay fit, February);Energy of a (keV);Thick target (a,n) yield");
	reactionsvenergy_decay_feb_2->SetMarkerStyle(34);
	reactionsvenergy_decay_feb_2->SetMarkerSize(MARKER_SIZE);
	reactionsvenergy_decay_feb_2->SetMarkerColor(kRed);
	myCanvas->SetName("reactions_v_energy_decay_feb_2");
//	reactionsvenergy_exfor_feb->Draw("E3");
//	reactionsvenergy_exfor_feb->Draw("same LX");
	reactionsvenergy_exfor->Draw("same E3");
//	reactionsvenergy_exfor->Draw("same LX");
	reactionsvenergy_decay_feb_2->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TGraph* reactionsvenergy_decay_apr_1 = new TGraphErrors(6, &x1[4], &y1[4], NULL, yerr1);
	reactionsvenergy_decay_apr_1->SetTitle("(a,n) reactions v a energy (decay it, April);Energy of a (keV);Thick target (a,n) yield");
	reactionsvenergy_decay_apr_1->SetMarkerStyle(23);
	reactionsvenergy_decay_apr_1->SetMarkerSize(MARKER_SIZE);
	reactionsvenergy_decay_apr_1->SetMarkerColor(kBlue);
	myCanvas->SetName("reactions_v_energy_decay_apr_1");
//	reactionsvenergy_exfor_apr->Draw("E3");
//	reactionsvenergy_exfor_apr->Draw("same LX");
	reactionsvenergy_exfor->Draw("same E3");
//	reactionsvenergy_exfor->Draw("same LX");
	reactionsvenergy_decay_apr_1->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TGraph* reactionsvenergy_decay_apr_2 = new TGraphErrors(6, &x2[4], &y2[4], NULL, yerr2);
	reactionsvenergy_decay_apr_2->SetTitle("(a,n) reactions v a energy (decay it, April);Energy of a (keV);Thick target (a,n) yield");
	reactionsvenergy_decay_apr_2->SetMarkerStyle(47);
	reactionsvenergy_decay_apr_2->SetMarkerSize(MARKER_SIZE);
	reactionsvenergy_decay_apr_2->SetMarkerColor(kBlue);
	myCanvas->SetName("reactions_v_energy_decay_apr_2");
//	reactionsvenergy_exfor_apr->Draw("E3");
//	reactionsvenergy_exfor_apr->Draw("same LX");
	reactionsvenergy_exfor->Draw("same E3");
//	reactionsvenergy_exfor->Draw("same LX");
	reactionsvenergy_decay_apr_2->Draw("same p");
	myCanvas->Write("", TObject::kOverwrite);

	TMultiGraph* reactions_v_energy_decay = new TMultiGraph();
//	reactions_v_energy_decay->Add(reactionsvenergy_exfor_feb, "E3");
//	reactions_v_energy_decay->Add(reactionsvenergy_exfor_feb, "LX");
//	reactions_v_energy_decay->Add(reactionsvenergy_exfor_apr, "E3");
//	reactions_v_energy_decay->Add(reactionsvenergy_exfor_apr, "LX");
	reactions_v_energy_decay->Add(reactionsvenergy_exfor, "E3");
//	reactions_v_energy_decay->Add(reactionsvenergy_exfor, "LX");
	reactionsvenergy_decay_feb_1->SetTitle("Decay fit, February, LaBr1");
	reactions_v_energy_decay->Add(reactionsvenergy_decay_feb_1);
	reactionsvenergy_decay_apr_1->SetTitle("Decay fit, April, LaBr1");
	reactions_v_energy_decay->Add(reactionsvenergy_decay_apr_1);
	reactionsvenergy_decay_feb_2->SetTitle("Decay fit, February, LaBr2");
	reactions_v_energy_decay->Add(reactionsvenergy_decay_feb_2);
	reactionsvenergy_decay_apr_2->SetTitle("Decay fit, April, LaBr2");
	reactions_v_energy_decay->Add(reactionsvenergy_decay_apr_2);
	reactions_v_energy_decay->SetTitle("(a,n) reactions v a energy (decay fit);Energy of a (keV);Thick target (a,n) yield");
	reactions_v_energy_decay->Draw("AP");
	myCanvas->BuildLegend();
	myCanvas->SetName("reactions_v_energy_decay");
	myCanvas->Write("", TObject::kOverwrite);

	//all results
	TMultiGraph* multigraph = new TMultiGraph();
//	multigraph->Add(reactionsvenergy_exfor_feb, "E3");
//	multigraph->Add(reactionsvenergy_exfor_feb, "LX");
//	multigraph->Add(reactionsvenergy_exfor_apr, "E3");
//	multigraph->Add(reactionsvenergy_exfor_apr, "LX");
	multigraph->Add(reactionsvenergy_exfor, "E3");
//	multigraph->Add(reactionsvenergy_exfor, "LX");
	reactionsvenergy_unified_feb_1->SetMarkerColor(kRed);
	reactionsvenergy_unified_feb_1->SetTitle("Unified fit, February, LaBr1");
	multigraph->Add(reactionsvenergy_unified_feb_1);
	reactionsvenergy_unified_apr_1->SetMarkerColor(kRed);
	reactionsvenergy_unified_apr_1->SetTitle("Unified fit, April, LaBr1");
	multigraph->Add(reactionsvenergy_unified_apr_1);
	reactionsvenergy_rise_feb_1->SetMarkerColor(kGreen);
	reactionsvenergy_rise_feb_1->SetTitle("Rise fit, February, LaBr1");
	multigraph->Add(reactionsvenergy_rise_feb_1);
	reactionsvenergy_rise_apr_1->SetMarkerColor(kGreen);
	reactionsvenergy_rise_apr_1->SetTitle("Rise fit, April, LaBr1");
	multigraph->Add(reactionsvenergy_rise_apr_1);
	reactionsvenergy_decay_feb_1->SetMarkerColor(kBlue);
	reactionsvenergy_decay_feb_1->SetTitle("Decay fit, February, LaBr1");
	multigraph->Add(reactionsvenergy_decay_feb_1);
	reactionsvenergy_decay_apr_1->SetMarkerColor(kBlue);
	reactionsvenergy_decay_apr_1->SetTitle("Decay fit, April, LaBr1");
	multigraph->Add(reactionsvenergy_decay_apr_1);
	reactionsvenergy_unified_feb_2->SetMarkerColor(kRed);
	reactionsvenergy_unified_feb_2->SetTitle("Unified fit, February, LaBr2");
	multigraph->Add(reactionsvenergy_unified_feb_2);
	reactionsvenergy_unified_apr_2->SetMarkerColor(kRed);
	reactionsvenergy_unified_apr_2->SetTitle("Unified fit, April, LaBr2");
	multigraph->Add(reactionsvenergy_unified_apr_2);
	reactionsvenergy_rise_feb_2->SetMarkerColor(kGreen);
	reactionsvenergy_rise_feb_2->SetTitle("Rise fit, February, LaBr2");
	multigraph->Add(reactionsvenergy_rise_feb_2);
	reactionsvenergy_rise_apr_2->SetMarkerColor(kGreen);
	reactionsvenergy_rise_apr_2->SetTitle("Rise fit, April, LaBr2");
	multigraph->Add(reactionsvenergy_rise_apr_2);
	reactionsvenergy_decay_feb_2->SetMarkerColor(kBlue);
	reactionsvenergy_decay_feb_2->SetTitle("Decay fit, February, LaBr2");
	multigraph->Add(reactionsvenergy_decay_feb_2);
	reactionsvenergy_decay_apr_2->SetMarkerColor(kBlue);
	reactionsvenergy_decay_apr_2->SetTitle("Decay fit, April, LaBr2");
	multigraph->Add(reactionsvenergy_decay_apr_2);
	multigraph->SetTitle("(a,n) reactions v a energy;Energy of a (keV);Thick target (a,n) yield");
	multigraph->Draw("AP");
	myCanvas->BuildLegend();
	myCanvas->SetName("reactions_v_energy");
	myCanvas->Write("", TObject::kOverwrite);

	//errores
	Double_t abs_errors_feb_labr1[4];	//errores absolutos respecto a exfor
	Double_t abs_errors_feb_labr2[4];
	Double_t abs_errors_apr_labr1[6];
	Double_t abs_errors_apr_labr2[6];
	Double_t per_errors_feb_labr1[4];	//errores porcentuales respecto a exfor
	Double_t per_errors_feb_labr2[4];
	Double_t per_errors_apr_labr1[6];
	Double_t per_errors_apr_labr2[6];
	Double_t rel_abs_errors_feb[4];		//diferencias absolutas entre detectores
	Double_t rel_abs_errors_apr[6];
	Double_t rel_per_errors_feb[4];		//diferencias porcentuales entre detectores
	Double_t rel_per_errors_apr[6];
	for(short i=0; i<4; i++){
		abs_errors_feb_labr1[i] = results[i][0][ERROR_I] - reactionsvenergy_exfor->Eval(activation_energies[i]);
		abs_errors_feb_labr2[i] = results[i][1][ERROR_I] - reactionsvenergy_exfor->Eval(activation_energies[i]);
		per_errors_feb_labr1[i] = abs_errors_feb_labr1[i]/reactionsvenergy_exfor->Eval(activation_energies[i]) * 100;
		per_errors_feb_labr2[i] = abs_errors_feb_labr2[i]/reactionsvenergy_exfor->Eval(activation_energies[i]) * 100;
		rel_abs_errors_feb[i] = results[i][0][ERROR_I] - results[i][1][ERROR_I];
		rel_per_errors_feb[i] = rel_abs_errors_feb[i] / ((results[i][0][ERROR_I]+results[i][1][ERROR_I])/2) * 100;
	}
	for(short i=0; i<6; i++){
		abs_errors_apr_labr1[i] = results[4+i][0][ERROR_I] - reactionsvenergy_exfor->Eval(activation_energies[4+i]);
		abs_errors_apr_labr2[i] = results[4+i][1][ERROR_I] - reactionsvenergy_exfor->Eval(activation_energies[4+i]);
		per_errors_apr_labr1[i] = abs_errors_apr_labr1[i]/reactionsvenergy_exfor->Eval(activation_energies[4+i]) * 100;
		per_errors_apr_labr2[i] = abs_errors_apr_labr2[i]/reactionsvenergy_exfor->Eval(activation_energies[4+i]) * 100;
		rel_abs_errors_apr[i] = results[4+i][0][ERROR_I] - results[4+i][1][ERROR_I];
		rel_per_errors_apr[i] = rel_abs_errors_apr[i] / ((results[4+i][0][ERROR_I]+results[4+i][1][ERROR_I])/2) * 100;
	}
	TGraph* abs_errors_feb_labr1_graph = new TGraph(4, activation_energies, abs_errors_feb_labr1);
	abs_errors_feb_labr1_graph->SetTitle("Absolute error, February, LaBr1");
	abs_errors_feb_labr1_graph->SetMarkerColor(kRed);
	abs_errors_feb_labr1_graph->SetMarkerStyle(22);
	abs_errors_feb_labr1_graph->SetMarkerSize(MARKER_SIZE);
	TGraph* abs_errors_feb_labr2_graph = new TGraph(4, activation_energies, abs_errors_feb_labr2);
	abs_errors_feb_labr2_graph->SetTitle("Absolute error, February, LaBr2");
	abs_errors_feb_labr2_graph->SetMarkerColor(kRed);
	abs_errors_feb_labr2_graph->SetMarkerStyle(34);
	abs_errors_feb_labr2_graph->SetMarkerSize(MARKER_SIZE);
	TGraph* abs_errors_apr_labr1_graph = new TGraph(6, &activation_energies[4], abs_errors_apr_labr1);
	abs_errors_apr_labr1_graph->SetTitle("Absolute error, April, LaBr1");
	abs_errors_apr_labr1_graph->SetMarkerColor(kBlue);
	abs_errors_apr_labr1_graph->SetMarkerStyle(23);
	abs_errors_apr_labr1_graph->SetMarkerSize(MARKER_SIZE);
	TGraph* abs_errors_apr_labr2_graph = new TGraph(6, &activation_energies[4], abs_errors_apr_labr2);
	abs_errors_apr_labr2_graph->SetTitle("Absolute error, April, LaBr2");
	abs_errors_apr_labr2_graph->SetMarkerColor(kBlue);
	abs_errors_apr_labr2_graph->SetMarkerStyle(47);
	abs_errors_apr_labr2_graph->SetMarkerSize(MARKER_SIZE);
	TGraph* per_errors_feb_labr1_graph = new TGraph(4, activation_energies, per_errors_feb_labr1);
	per_errors_feb_labr1_graph->SetTitle("Relative error, February, LaBr1");
	per_errors_feb_labr1_graph->SetMarkerColor(kRed);
	per_errors_feb_labr1_graph->SetMarkerStyle(22);
	per_errors_feb_labr1_graph->SetMarkerSize(MARKER_SIZE);
	TGraph* per_errors_feb_labr2_graph = new TGraph(4, activation_energies, per_errors_feb_labr2);
	per_errors_feb_labr2_graph->SetTitle("Relative error, February, LaBr2");
	per_errors_feb_labr2_graph->SetMarkerColor(kRed);
	per_errors_feb_labr2_graph->SetMarkerStyle(34);
	per_errors_feb_labr2_graph->SetMarkerSize(MARKER_SIZE);
	TGraph* per_errors_apr_labr1_graph = new TGraph(6, &activation_energies[4], per_errors_apr_labr1);
	per_errors_apr_labr1_graph->SetTitle("Relative error, April, LaBr1");
	per_errors_apr_labr1_graph->SetMarkerColor(kBlue);
	per_errors_apr_labr1_graph->SetMarkerStyle(23);
	per_errors_apr_labr1_graph->SetMarkerSize(MARKER_SIZE);
	TGraph* per_errors_apr_labr2_graph = new TGraph(6, &activation_energies[4], per_errors_apr_labr2);
	per_errors_apr_labr2_graph->SetTitle("Relative error, April, LaBr2");
	per_errors_apr_labr2_graph->SetMarkerColor(kBlue);
	per_errors_apr_labr2_graph->SetMarkerStyle(47);
	per_errors_apr_labr2_graph->SetMarkerSize(MARKER_SIZE);
	TGraph* rel_abs_errors_feb_graph = new TGraph(4, activation_energies, rel_abs_errors_feb);
	rel_abs_errors_feb_graph->SetTitle("Absolute difference, February");
	rel_abs_errors_feb_graph->SetMarkerColor(kRed);
	rel_abs_errors_feb_graph->SetMarkerStyle(20);
	rel_abs_errors_feb_graph->SetMarkerSize(MARKER_SIZE);
	TGraph* rel_abs_errors_apr_graph = new TGraph(6, &activation_energies[4], rel_abs_errors_apr);
	rel_abs_errors_apr_graph->SetTitle("Absolute difference, April");
	rel_abs_errors_apr_graph->SetMarkerColor(kBlue);
	rel_abs_errors_apr_graph->SetMarkerStyle(21);
	rel_abs_errors_apr_graph->SetMarkerSize(MARKER_SIZE);
	TGraph* rel_per_errors_feb_graph = new TGraph(4, activation_energies, rel_per_errors_feb);
	rel_per_errors_feb_graph->SetTitle("Percentage difference, February");
	rel_per_errors_feb_graph->SetMarkerColor(kRed);
	rel_per_errors_feb_graph->SetMarkerStyle(20);
	rel_per_errors_feb_graph->SetMarkerSize(MARKER_SIZE);
	TGraph* rel_per_errors_apr_graph = new TGraph(6, &activation_energies[4], rel_per_errors_apr);
	rel_per_errors_apr_graph->SetTitle("Percentage difference, April");
	rel_per_errors_apr_graph->SetMarkerColor(kBlue);
	rel_per_errors_apr_graph->SetMarkerStyle(21);
	rel_per_errors_apr_graph->SetMarkerSize(MARKER_SIZE);

	TMultiGraph* abs_errors = new TMultiGraph();
	abs_errors->Add(abs_errors_feb_labr1_graph);
	abs_errors->Add(abs_errors_feb_labr2_graph);
	abs_errors->Add(abs_errors_apr_labr1_graph);
	abs_errors->Add(abs_errors_apr_labr2_graph);
	abs_errors->SetTitle("Absolute errors with respect to EXFOR data");
	abs_errors->Draw("AP");
	myCanvas->BuildLegend();
	myCanvas->Write("abs_errors", TObject::kOverwrite);

	TMultiGraph* per_errors = new TMultiGraph();
	per_errors->Add(per_errors_feb_labr1_graph);
	per_errors->Add(per_errors_feb_labr2_graph);
	per_errors->Add(per_errors_apr_labr1_graph);
	per_errors->Add(per_errors_apr_labr2_graph);
	per_errors->SetTitle("Relative errors with respect to EXFOR data");
	per_errors->Draw("AP");
	myCanvas->BuildLegend();
	myCanvas->Write("per_errors", TObject::kOverwrite);

	TMultiGraph* rel_abs_errors = new TMultiGraph();
	rel_abs_errors->Add(rel_abs_errors_feb_graph);
	rel_abs_errors->Add(rel_abs_errors_apr_graph);
	rel_abs_errors->SetTitle("Absolute errors with respect to average between Labr 1 and 2");
	rel_abs_errors->Draw("AP");
	myCanvas->BuildLegend();
	myCanvas->Write("rel_abs_errors", TObject::kOverwrite);

	TMultiGraph* rel_per_errors = new TMultiGraph();
	rel_per_errors->Add(rel_per_errors_feb_graph);
	rel_per_errors->Add(rel_per_errors_apr_graph);
	rel_per_errors->SetTitle("Relative errors with respect to average between Labr 1 and 2");
	rel_per_errors->Draw("AP");
	myCanvas->BuildLegend();
	myCanvas->Write("rel_per_errors", TObject::kOverwrite);

	unified_average[3] = unified_average_diff[3] = unified_average_diff_per[3] = 0;
	unified_average[5] = unified_average_diff[5] = unified_average_diff_per[5] = 0;
	unified_average[7] = unified_average_diff[7] = unified_average_diff_per[7] = 0;
	unified_average[8] = unified_average_diff[8] = unified_average_diff_per[8] = 0;
	unified_average[9] = unified_average_diff[9] = unified_average_diff_per[9] = 0;
	rise_average[3] = rise_average_diff[3] = rise_average_diff_per[3] = 0;
	rise_average[5] = rise_average_diff[5] = rise_average_diff_per[5] = 0;
	rise_average[7] = rise_average_diff[7] = rise_average_diff_per[7] = 0;
	rise_average[8] = rise_average_diff[8] = rise_average_diff_per[8] = 0;
	rise_average[9] = rise_average_diff[9] = rise_average_diff_per[9] = 0;
	decay_average[5] = decay_average_diff[5] = decay_average_diff_per[5] = 0;
	cout << endl << "----- RESULTS -----" << endl;
	Double_t exfor_result[ACTIVATION_N];
	Double_t decay_exfor_diff[ACTIVATION_N];
	Double_t unified_exfor_diff[ACTIVATION_N];
	Double_t rise_exfor_diff[ACTIVATION_N];
	for(short i=0; i<ACTIVATION_N; i++){
		cout << "Activation N " << i+1 << endl;
		cout << "Activation E: " << activation_energies[i] << endl;
		exfor_result[i] = reactionsvenergy_exfor->Eval(activation_energies[i]);
		cout << "exfor: " << exfor_result[i] << endl;
		decay_exfor_diff[i] = (decay_average[i] - exfor_result[i])/exfor_result[i]*100;
		unified_exfor_diff[i] = (unified_average[i] - exfor_result[i])/exfor_result[i]*100;
		rise_exfor_diff[i] = (rise_average[i] - exfor_result[i])/exfor_result[i]*100;
		cout << "decay average : " << decay_average[i] << " (" << decay_average_diff_per[i] << "%) | " << decay_exfor_diff[i] << "%" << endl;
		cout << "unified average : " << unified_average[i] << " (" << unified_average_diff_per[i] << "%) | " << unified_exfor_diff[i] << "%" << endl;
		cout << "rise average : " << rise_average[i] << " (" << rise_average_diff_per[i] << "%) | " << rise_exfor_diff[i] << "%" << endl;
		cout << endl;
	}

	TMultiGraph* final_exfor_diffs_multigraph = new TMultiGraph();
	TGraph* unified_exfor_diffs_feb = new TGraph(4, activation_energies, unified_exfor_diff);
	unified_exfor_diffs_feb->SetTitle("Unified diffs, feb");
	unified_exfor_diffs_feb->SetMarkerColor(kRed);
	unified_exfor_diffs_feb->SetMarkerStyle(34);
	unified_exfor_diffs_feb->SetMarkerSize(MARKER_SIZE);
	TGraph* unified_exfor_diffs_apr = new TGraph(6, &activation_energies[4], &unified_exfor_diff[4]);
	unified_exfor_diffs_apr->SetTitle("Unified diffs, apr");
	unified_exfor_diffs_apr->SetMarkerColor(kRed);
	unified_exfor_diffs_apr->SetMarkerStyle(47);
	unified_exfor_diffs_apr->SetMarkerSize(MARKER_SIZE);
	TGraph* rise_exfor_diffs_feb = new TGraph(4, activation_energies, rise_exfor_diff);
	rise_exfor_diffs_feb->SetTitle("Rise diffs, feb");
	rise_exfor_diffs_feb->SetMarkerColor(kGreen);
	rise_exfor_diffs_feb->SetMarkerStyle(34);
	rise_exfor_diffs_feb->SetMarkerSize(MARKER_SIZE);
	TGraph* rise_exfor_diffs_apr = new TGraph(6, &activation_energies[4], &rise_exfor_diff[4]);
	rise_exfor_diffs_apr->SetTitle("Rise diffs, apr");
	rise_exfor_diffs_apr->SetMarkerColor(kGreen);
	rise_exfor_diffs_apr->SetMarkerStyle(47);
	rise_exfor_diffs_apr->SetMarkerSize(MARKER_SIZE);
	TGraph* decay_exfor_diffs_feb = new TGraph(4, activation_energies, decay_exfor_diff);
	decay_exfor_diffs_feb->SetTitle("Decay diffs, feb");
	decay_exfor_diffs_feb->SetMarkerColor(kBlue);
	decay_exfor_diffs_feb->SetMarkerStyle(34);
	decay_exfor_diffs_feb->SetMarkerSize(MARKER_SIZE);
	TGraph* decay_exfor_diffs_apr = new TGraph(6, &activation_energies[4], &decay_exfor_diff[4]);
	decay_exfor_diffs_apr->SetTitle("Decay diffs, apr");
	decay_exfor_diffs_apr->SetMarkerColor(kBlue);
	decay_exfor_diffs_apr->SetMarkerStyle(47);
	decay_exfor_diffs_apr->SetMarkerSize(MARKER_SIZE);
	final_exfor_diffs_multigraph->Add(unified_exfor_diffs_feb);
	final_exfor_diffs_multigraph->Add(rise_exfor_diffs_feb);
	final_exfor_diffs_multigraph->Add(decay_exfor_diffs_feb);
	final_exfor_diffs_multigraph->Add(unified_exfor_diffs_apr);
	final_exfor_diffs_multigraph->Add(rise_exfor_diffs_apr);
	final_exfor_diffs_multigraph->Add(decay_exfor_diffs_apr);
	final_exfor_diffs_multigraph->SetTitle("Relative error with respect to EXFOR data");
	final_exfor_diffs_multigraph->Draw("AP");
	myCanvas->BuildLegend();
	myCanvas->Write("exfor_diffs", TObject::kOverwrite);

	myCanvas->Close();
	gDirectory->cd("..");

	std::cout.rdbuf(originalBuffer);	//restaurar cout
	outputFile.close();
	f.Close();
}

Double_t activation_window_low;
Double_t activation_window_high;

void activation(){
	//redirige cout a archivo de texto
	std::ofstream outputFile("fitting_activacion.txt");
	std::streambuf* originalBuffer = std::cout.rdbuf();
	std::cout.rdbuf(outputFile.rdbuf());

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
	std::cout.rdbuf(originalBuffer);	//restaurar cout
	outputFile.close();
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
	if(strcmp(filepath, "output/SData_aAl_J78keV_GVM2731kV_LaBr1_20cm-135deg_LaBr2_20cm135deg_activacion_20230418.root")==0){	//exception for activation 7
//		activation_start = 300.899E12;
//		activation_end = 753.207E12;
		integrator_signals = d.Filter("Channel==1 && Timestamp>299E12 && Timestamp<755E12");
	}

	auto current_integrator = integrator_signals.Define("t", "Timestamp/1E12").Histo1D({"current_integrator", ";Timestamp (s);Counts", ACTIVATION_NBINS, 0, measurement_end/1E12}, "t");
	Double_t activation_start = integrator_signals.Min("Timestamp").GetValue();	//TBD!:muy ineficiente!!
	Double_t activation_end = integrator_signals.Max("Timestamp").GetValue();
	Double_t current2alpha = 1/(2*1.60217646E-9); //10^-10C to number of alphas
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
	Double_t activation_time = activation_end-activation_start;
	cout << "*Number of alphas: " << number_of_alphas << endl;
	cout << "Alphas per second: " << number_of_alphas/(activation_time/1E12) << endl;
//	cout << "Current integrator integral (nC): " << current_integrator->Integral()/10 << endl;
//	cout << "current2alpha: " << current2alpha << endl;
//	cout << "current (nA): " << current_integrator->Integral()/(activation_time/1E12)/10 << endl;
	cout << "current (nA): " << (number_of_alphas/current2alpha)/(activation_time/1E12)/10 << endl;

	//histogramas
	auto rise_filter = [&](ULong64_t t){return t>=activation_start && t<=activation_end;};
	auto decay_filter = [&](ULong64_t t){return t>activation_end;};
	auto energy_window = d.Filter("Energy>activation_window_low && Energy<activation_window_high");
	auto labr_1_filter = energy_window.Filter("Channel==6").Define("t", "Timestamp/1E12");
	auto labr_2_filter = energy_window.Filter("Channel==7").Define("t", "Timestamp/1E12");

	auto labr_1 = labr_1_filter.Histo1D({"labr_1", ";Time (s);Counts", ACTIVATION_NBINS, 0, measurement_end/1E12}, "t");
	auto labr_2 = labr_2_filter.Histo1D({"labr_2", ";Time (s);Counts", ACTIVATION_NBINS, 0, measurement_end/1E12}, "t");

	int rise_nbins = activation_time/(labr_1->GetBinWidth(1)*1E12);
	int decay_nbins = (measurement_end-activation_end)/(labr_1->GetBinWidth(1)*1E12);
//	cout << "rise_nbins: " << rise_nbins << endl;
//	cout << "decay_nbins: " << decay_nbins << endl;
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
	activation_time/=1E12;
	measurement_end/=1E12;
//	cout << "Activation start: " << activation_start << "s" << endl;
//	cout << "Activation end: " << activation_end << "s" << endl;
	cout << "Activation time: " << activation_time << "s" << endl;
//	cout << "Measurement end: " << measurement_end << "s" << endl;
//	cout << "activation_window_low: " << activation_window_low << endl;
//	cout << "activation_window_high: " << activation_window_high << endl;

	Double_t labr1_binwidth = labr_1->GetBinWidth(1);
	Double_t labr2_binwidth = labr_2->GetBinWidth(1);
	Double_t labr1_decay_binwidth = labr_1_decay->GetBinWidth(1);
	Double_t labr2_decay_binwidth = labr_2_decay->GetBinWidth(1);
	cout << "labr1 binwidth: " << labr1_binwidth << endl;
//	cout << "labr2 binwidth: " << labr2_binwidth << endl;
	cout << "labr1 decay binwidth: " << labr1_decay_binwidth << endl;
//	cout << "labr2 decay binwidth: " << labr2_decay_binwidth << endl;

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
	results[0][0] = fitresult->Parameter(1)/current2alpha/1.99;
	results[0][1] = fitresult->ParError(1)/current2alpha/1.99;
	myCanvas->SetName("labr_1_unified_fit");
	myCanvas->Write("", TObject::kOverwrite);

	fitresult = labr_2->Fit("unified_fit", "SLEQ");
	results[1][0] = fitresult->Parameter(1)/current2alpha/1.99;
	results[1][1] = fitresult->ParError(1)/current2alpha/1.99;
	myCanvas->SetName("labr_2_unified_fit");
	myCanvas->Write("", TObject::kOverwrite);

	cout << "labr1 unified 30P per alpha: " << results[0][0] << endl;
	cout << "labr2 unified 30P per alpha: " << results[1][0] << endl;
//	cout << "labr1 unified number of 30P: " << results[0][0]*number_of_alphas << endl;
//	cout << "labr2 unified number of 30P: " << results[1][0]*number_of_alphas << endl;

	//rise
	unified->SetNpx(rise_nbins);

	fitresult = labr_1_rise->Fit("unified_fit", "SLEQ");
	results[0][2] = fitresult->Parameter(1)/current2alpha/1.99;
	results[0][3] = fitresult->ParError(1)/current2alpha/1.99;
	myCanvas->SetName("labr_1_rise_fit");
	myCanvas->Write("", TObject::kOverwrite);

	fitresult = labr_2_rise->Fit("unified_fit", "SLEQ");
	results[1][2] = fitresult->Parameter(1)/current2alpha/1.99;
	results[1][3] = fitresult->ParError(1)/current2alpha/1.99;
	myCanvas->SetName("labr_2_rise_fit");
	myCanvas->Write("", TObject::kOverwrite);

	cout << "labr1 rise 30P per alpha: " << results[0][2] << endl;
	cout << "labr2 rise 30P per alpha: " << results[1][2] << endl;
//	cout << "labr1 rise number of 30P: " << results[0][2]*number_of_alphas << endl;
//	cout << "labr2 rise number of 30P: " << results[1][2]*number_of_alphas << endl;

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
	Double_t decay_factor = 1/(1-exp(-fitresult->Parameter(2)*activation_time));
	results[0][4] = fitresult->Parameter(1)*decay_factor/(labr1_decay_binwidth*number_of_alphas*fitresult->Parameter(2)*1.99);
	results[0][5] = fitresult->ParError(1)*decay_factor/(labr1_decay_binwidth*number_of_alphas*fitresult->Parameter(2)*1.99);
	myCanvas->SetName("labr_1_decay_fit");
	myCanvas->Write("", TObject::kOverwrite);
	cout << "Decay factor: " << decay_factor << endl << endl;
	cout << "LaBr1 initial activity (/s): " << fitresult->Parameter(1)/labr1_decay_binwidth << endl;

	fitresult = labr_2_decay->Fit("decay", "SLEQ");
	results[1][4] = fitresult->Parameter(1)*decay_factor/(labr2_decay_binwidth*number_of_alphas*fitresult->Parameter(2)*1.99);
	results[1][5] = fitresult->ParError(1)*decay_factor/(labr2_decay_binwidth*number_of_alphas*fitresult->Parameter(2)*1.99);
	myCanvas->SetName("labr_2_decay_fit");
	myCanvas->Write("", TObject::kOverwrite);
	cout << "LaBr2 initial activity (/s): " << fitresult->Parameter(1)/labr2_decay_binwidth << endl;

	cout << "LaBr1 decay 30P per alpha: " << results[0][4] << endl;
	cout << "LaBr1 decay 30P per alpha error: " << results[0][4] << endl;
	cout << "LaBr2 decay 30P per alpha: " << results[1][4] << endl;
	cout << "LaBr2 decay 30P per alpha error: " << results[1][4] << endl;

	myCanvas->Close();
	DisableImplicitMT();	//multithreading
	cout << "---------------" << endl;
}
