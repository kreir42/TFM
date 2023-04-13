#define ACTIVATION_NBINS 1000

static void per_file(Char_t filepath[500], Double_t results[2][6]);

void activation(){
	cout << "Activation" << endl;
	Char_t filepath_1[100] = "output/SData_aAl_J78kV_GVM1808kV_positions2_activacion.root";
	Char_t filepath_2[100] = "output/SData_aAl_J78kV_GVM2312kV_positions2_activacion.root";
	Char_t filepath_3[100] = "output/SData_aAl_J78kV_GVM2810kV_positions2_activacion.root";
	Char_t filepath_4[100] = "output/SData_aAl_J78kV_GVM2810kV_positions2_activacion_20230223.root";

	TFile f("output.root", "UPDATE");
	gDirectory->cd("Activation");
	Double_t results[4][2][6];
	Double_t activation_energies[] = {5500, 7000, 8500, 8500};	//keV

	Double_t exfor_energies[] = {3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500, 5600, 5700, 5800, 5900, 6000, 6100, 6200, 6300, 6400, 6500, 6600, 6700, 6800, 6900, 7000, 7100, 7200, 7300, 7400, 7500, 7600, 7700, 7800, 7900, 8000, 8100, 8200, 8300, 8400, 8500, 8600, 8700, 8800, 8900, 9000, 9100, 9200, 9300, 9400, 9500, 9600, 9700, 9800, 9900};	//TBD:hardcoded, read .txt

	Double_t exfor_data[] = {3.147E-09, 5.904E-09, 1.034E-08, 1.655E-08, 2.462E-08, 3.464E-08, 4.686E-08, 6.203E-08, 8.124E-08, 1.073E-07, 1.426E-07, 1.847E-07, 2.306E-07, 2.812E-07, 3.403E-07, 4.150E-07, 5.119E-07, 6.278E-07, 7.555E-07, 8.856E-07, 1.011E-06, 1.150E-06, 1.330E-06, 1.549E-06, 1.797E-06, 2.062E-06, 2.339E-06, 2.651E-06, 3.015E-06, 3.401E-06, 3.774E-06, 4.147E-06, 4.552E-06, 4.999E-06, 5.489E-06, 6.013E-06, 6.562E-06, 7.131E-06, 7.716E-06, 8.319E-06, 8.943E-06, 9.593E-06, 1.027E-05, 1.099E-05, 1.173E-05, 1.252E-05, 1.333E-05, 1.416E-05, 1.502E-05, 1.589E-05, 1.679E-05, 1.771E-05, 1.865E-05, 1.962E-05, 2.062E-05, 2.165E-05, 2.272E-05, 2.383E-05, 2.497E-05, 2.616E-05, 2.740E-05, 2.869E-05, 3.003E-05};	//TBD:hardcoded, read .txt

	for(short i=0; i<63; i++){	//TBD:escalado temporal, números hardcoded
		exfor_data[i]*=8.71921676E-5;
	}

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
	TCanvas* myCanvas = new TCanvas("reactions_v_energy_unified");

	//EXFOR data
	x[0] = activation_energies[0];
	x[1] = activation_energies[1];
	x[2] = activation_energies[2];
	TGraph* rectionsvenergy_exfor = new TGraphErrors(63, exfor_energies, exfor_data, NULL, NULL);	//TBD:hardcoded number
	rectionsvenergy_exfor->SetTitle("(a,n) reactions v a energy;Energy of a (keV);Inferred (a,n)/Number of a");
	rectionsvenergy_exfor->SetMarkerStyle(20);

	//unified_fit
	for(short i=0; i<4; i++){
		x[2*i] = activation_energies[i];
		x[2*i+1] = activation_energies[i];
		y[2*i] = results[i][0][0];
		y[2*i+1] = results[i][1][0];
		yerr[2*i] = results[i][0][1];
		yerr[2*i+1] = results[i][1][1];
	}
	TGraph* rectionsvenergy_unified = new TGraphErrors(8, x, y, NULL, yerr);
	rectionsvenergy_unified->SetTitle("(a,n) reactions v a energy (unified fit);Energy of a (keV);Inferred (a,n)/Number of a");
	rectionsvenergy_unified->SetMarkerStyle(20);
	rectionsvenergy_unified->Draw("ap");
	rectionsvenergy_exfor->Draw("same");
	myCanvas->Write();

	//rise fit
	for(short i=0; i<4; i++){
		x[2*i] = activation_energies[i];
		x[2*i+1] = activation_energies[i];
		y[2*i] = results[i][0][2];
		y[2*i+1] = results[i][1][2];
		yerr[2*i] = results[i][0][3];
		yerr[2*i+1] = results[i][1][3];
	}
	TGraph* rectionsvenergy_rise = new TGraphErrors(8, x, y, NULL, yerr);
	rectionsvenergy_rise->SetTitle("(a,n) reactions v a energy (rise fit);Energy of a (keV);Inferred (a,n)/Number of a");
	rectionsvenergy_rise->SetMarkerStyle(20);
	myCanvas->SetName("reactions_v_energy_rise");
	rectionsvenergy_rise->Draw("ap");
	rectionsvenergy_exfor->Draw("same");
	myCanvas->Write();

	//decay fit
	for(short i=0; i<4; i++){
		x[2*i] = activation_energies[i];
		x[2*i+1] = activation_energies[i];
		y[2*i] = results[i][0][4];
		y[2*i+1] = results[i][1][4];
		yerr[2*i] = results[i][0][5];
		yerr[2*i+1] = results[i][1][5];
	}
	TGraph* rectionsvenergy_decay = new TGraphErrors(8, x, y, NULL, yerr);
	rectionsvenergy_decay->SetTitle("(a,n) reactions v a energy (decay fit);Energy of a (keV);Inferred (a,n)/Number of a");
	rectionsvenergy_decay->SetMarkerStyle(20);
	myCanvas->SetName("reactions_v_energy_decay");
	rectionsvenergy_decay->Draw("ap");
	rectionsvenergy_exfor->Draw("same");
	myCanvas->Write();

	//all results
	TMultiGraph* multigraph = new TMultiGraph();
	rectionsvenergy_exfor->SetMarkerColor(kBlack);
	rectionsvenergy_exfor->SetTitle("EXFOR data");
	multigraph->Add(rectionsvenergy_exfor,"c");
	rectionsvenergy_unified->SetMarkerColor(kRed);
	rectionsvenergy_unified->SetTitle("Unified fit");
	multigraph->Add(rectionsvenergy_unified);
	rectionsvenergy_rise->SetMarkerColor(kGreen);
	rectionsvenergy_rise->SetTitle("Rise fit");
	multigraph->Add(rectionsvenergy_rise);
	rectionsvenergy_decay->SetMarkerColor(kBlue);
	rectionsvenergy_decay->SetTitle("Decay fit");
	multigraph->Add(rectionsvenergy_decay);
	myCanvas->SetLogy();
	multigraph->SetTitle("(a,n) reactions v a energy;Energy of a (keV);Inferred (a,n)/Number of a");
	multigraph->Draw("AP");
	myCanvas->BuildLegend();
	myCanvas->SetName("reactions_v_energy");
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
	EnableImplicitMT();	//multithreading
	RDataFrame d("Data", filepath);

	auto integrator_signals = d.Filter("Channel==1");
	Double_t activation_start = integrator_signals.Min("Timestamp").GetValue();	//TBD!:muy ineficiente!!
	Double_t activation_end = integrator_signals.Max("Timestamp").GetValue();
	Double_t measurement_end = d.Filter("(Channel==6||Channel==7) && Energy>0").Max("Timestamp").GetValue();
	Double_t current2alpha = 1/(2*1.60217646E-10);
	Double_t number_of_alphas = integrator_signals.Count().GetValue()*current2alpha;

	//histogramas
	auto rise_filter = [&](ULong64_t t){return t>=activation_start && t<=activation_end;};
	auto decay_filter = [&](ULong64_t t){return t>activation_end;};
	auto labr_1_filter = d.Filter("Channel==6 && Energy>525 && Energy<750").Define("t", "Timestamp/1E12");
	auto labr_2_filter = d.Filter("Channel==7 && Energy>525 && Energy<750").Define("t", "Timestamp/1E12");

	auto current_integrator = integrator_signals.Define("t", "Timestamp/1E12").Histo1D({"current_integrator", ";Timestamp (s);Counts", ACTIVATION_NBINS, 0, measurement_end/1E12}, "t");
	auto labr_1 = labr_1_filter.Histo1D({"labr_1", ";Time (s);Counts", ACTIVATION_NBINS, 0, measurement_end/1E12}, "t");
	auto labr_2 = labr_2_filter.Histo1D({"labr_2", ";Time (s);Counts", ACTIVATION_NBINS, 0, measurement_end/1E12}, "t");

	int rise_nbins = (activation_end-activation_start)/(labr_1->GetBinWidth(1)*1E12);
	int decay_nbins = (measurement_end-activation_end)/(labr_1->GetBinWidth(1)*1E12);
	auto labr_1_rise = labr_1_filter.Filter(rise_filter, {"Timestamp"}).Histo1D({"labr_1_rise", ";Time (s);Counts", rise_nbins, activation_start/1E12, activation_end/1E12}, "t");
	auto labr_1_decay = labr_1_filter.Filter(decay_filter, {"Timestamp"}).Histo1D({"labr_1_decay", ";Time (s); Counts", decay_nbins, activation_end/1E12, measurement_end/1E12}, "t");
	auto labr_2_rise = labr_2_filter.Filter(rise_filter, {"Timestamp"}).Histo1D({"labr_2_rise", ";Time (s);Counts", rise_nbins, activation_start/1E12, activation_end/1E12}, "t");
	auto labr_2_decay = labr_2_filter.Filter(decay_filter, {"Timestamp"}).Histo1D({"labr_2_decay", ";Time (s);Counts", decay_nbins, activation_end/1E12, measurement_end/1E12}, "t");

	cout << "Rise/decay histograms" << endl;
	current_integrator->Write();
	labr_1->Write();
	labr_2->Write();
	labr_1_rise->Write();
	labr_1_decay->Write();
	labr_2_rise->Write();
	labr_2_decay->Write();

	activation_start/=1E12;
	activation_end/=1E12;
	measurement_end/=1E12;

	//fittings
	TCanvas* myCanvas = new TCanvas("myCanvas");
	TFitResultPtr fitresult;

	//unified
	cout << "Unified rise/decay fittings" << endl;
	TF1* unified = new TF1("unified_fit", unified_fit((TH1F*)gDirectory->Get("current_integrator")), 0, measurement_end, 4);
	unified->SetNpx(ACTIVATION_NBINS);
	unified->SetNumberFitPoints(ACTIVATION_NBINS);
	unified->SetParLimits(0, 0, 1E5);
	unified->SetParLimits(1, 0, 1E2);
	unified->SetParLimits(2, 4E-3, 5E-3);
	unified->SetParLimits(3, 0, 1E2);
	unified->SetParameters(15, 1, 4.62406E-3, 0);
	unified->FixParameter(2, 4.62406E-3);
	unified->SetParNames("Background activity", "current to (a,n)", "Decay constant", "extra bg");

	fitresult = labr_1->Fit("unified_fit", "SLE");
	results[0][0] = fitresult->Parameter(1)/current2alpha*labr_1->GetBinWidth(1);
	results[0][1] = fitresult->ParError(1)/current2alpha*labr_1->GetBinWidth(1);
	myCanvas->SetName("labr_1_unified_fit");
	myCanvas->Write();

	fitresult = labr_2->Fit("unified_fit", "SLE");
	results[1][0] = fitresult->Parameter(1)/current2alpha*labr_1->GetBinWidth(1);
	results[1][1] = fitresult->ParError(1)/current2alpha*labr_1->GetBinWidth(1);
	myCanvas->SetName("labr_2_unified_fit");
	myCanvas->Write();

	//rise
	cout << "Rise fittings" << endl;
	unified->SetNpx(rise_nbins);

	fitresult = labr_1_rise->Fit("unified_fit", "SLE");
	results[0][2] = fitresult->Parameter(1)/current2alpha*labr_1->GetBinWidth(1);
	results[0][3] = fitresult->ParError(1)/current2alpha*labr_1->GetBinWidth(1);
	myCanvas->SetName("labr_1_rise");
	myCanvas->Write();

	fitresult = labr_2_rise->Fit("unified_fit", "SLE");
	results[1][2] = fitresult->Parameter(1)/current2alpha*labr_1->GetBinWidth(1);
	results[1][3] = fitresult->ParError(1)/current2alpha*labr_1->GetBinWidth(1);
	myCanvas->SetName("labr_2_rise");
	myCanvas->Write();

	//decay
	cout << "Decay fittings" << endl;
	TF1* decay = new TF1("decay","[0]+[1]*exp(-[2]*(x[0]-[3]))");
	decay->SetNpx(decay_nbins);
	decay->SetNumberFitPoints(ACTIVATION_NBINS);
	decay->SetParLimits(0, 0, 100);
	decay->SetParLimits(2, 4E-3, 5E-3);
	decay->SetParameters(15, 1000, 4.62406E-3, activation_end);
	decay->FixParameter(3, activation_end);
	decay->FixParameter(2, 4.62406E-3);
	decay->SetParNames("Background activity", "Initial activiy", "Decay constant", "activation_end");

	fitresult = labr_1_decay->Fit("decay", "SLE");
	results[0][4] = fitresult->Parameter(1)*(activation_end-activation_start)/((1-exp(-fitresult->Parameter(2)*(activation_end-activation_start)))*number_of_alphas*labr_1->GetBinWidth(1));
	results[0][5] = fitresult->ParError(1)*(activation_end-activation_start)/((1-exp(-fitresult->Parameter(2)*(activation_end-activation_start)))*number_of_alphas*labr_1->GetBinWidth(1));
	myCanvas->SetName("labr_1_decay");
	myCanvas->Write();

	fitresult = labr_2_decay->Fit("decay", "SLE");
	results[1][4] = fitresult->Parameter(1)*(activation_end-activation_start)/((1-exp(-fitresult->Parameter(2)*(activation_end-activation_start)))*number_of_alphas*labr_2->GetBinWidth(1));
	results[1][5] = fitresult->ParError(1)*(activation_end-activation_start)/((1-exp(-fitresult->Parameter(2)*(activation_end-activation_start)))*number_of_alphas*labr_2->GetBinWidth(1));
	myCanvas->SetName("labr_2_decay");
	myCanvas->Write();

	myCanvas->Close();
	DisableImplicitMT();	//multithreading
}
