#define C 299792458
#define NEUTRON_MASS 1.67492749804E-27
#define J_TO_KEV (1/1.602177E-16)
#define TOF_TO_S (1E-9)

#define MIN_TOF 250
#define MAX_TOF 650
#define PULSE_ENERGY_FILTER 0
#define GAMMA_FLASH_BINS_N 350
#define NEUTRON_RESPONSE_BINS_N 500
#define PULSE_FIT_PARAMS_N 10	//must be even; if changed, must also change tree size definition
#define MAX_PARAM_TOF 150
#define MIN_PARAM_TOF 0

#define GMIN_TOF -15
#define GMAX_TOF +100
#define NMIN_TOF 20
#define NMAX_TOF 150


#define PULSED1_GMIN	480
#define PULSED1_GCENTER	488.096
#define PULSED1_GMAX	520
#define PULSED1_NMIN	510
#define PULSED1_NMAX	600

#define PULSED2_GMIN	478
#define PULSED2_GCENTER	484.855
#define PULSED2_GMAX	520
#define PULSED2_NMIN	500
#define PULSED2_NMAX	600

#define PULSED3_GMIN	364
#define PULSED3_GCENTER	370.548
#define PULSED3_GMAX	400
#define PULSED3_NMIN	380
#define PULSED3_NMAX	480

#define PULSED4_GMIN	285
#define PULSED4_GCENTER	294.508
#define PULSED4_GMAX	320
#define PULSED4_NMIN	300
#define PULSED4_NMAX	400

#define PULSED5_GMIN	290
#define PULSED5_GCENTER	297.949
#define PULSED5_GMAX	350
#define PULSED5_NMIN	330
#define PULSED5_NMAX	550

class pulse_fit_functor{
	private:
		TH1F* gamma_histo;
		Double_t gamma_stepsize;
	public:
		Double_t gamma_min;
		Double_t neutron_min;
		Double_t neutron_max;
		Double_t min_param_tof;
		Double_t max_param_tof;
		pulse_fit_functor(TH1F* gamma_histo, Double_t gamma_min, Double_t neutron_min, Double_t neutron_max, Double_t min_param_tof, Double_t max_param_tof):gamma_histo(gamma_histo),gamma_min(gamma_min),neutron_min(neutron_min),neutron_max(neutron_max),min_param_tof(min_param_tof),max_param_tof(max_param_tof){
			gamma_stepsize = gamma_histo->GetBinWidth(1);
		}
		Double_t operator()(Double_t* x, Double_t* p){
			Double_t sum = 0;
			Short_t param_index = 0;
			Double_t tof = 0;
			for(Short_t i=0; i<GAMMA_FLASH_BINS_N; i++){
				tof = x[0] - (gamma_min+(0.5+i)*gamma_stepsize);
				if(tof<min_param_tof || tof>max_param_tof){
					continue;
				}
				param_index = 1 + (tof-min_param_tof)/max_param_tof*PULSE_FIT_PARAMS_N;
				if(param_index>0 && param_index<=PULSE_FIT_PARAMS_N){
					sum += gamma_histo->GetBinContent(i)>p[PULSE_FIT_PARAMS_N+1] ? (gamma_histo->GetBinContent(i)-p[PULSE_FIT_PARAMS_N+1]) * p[param_index] : 0;
				}
			}
			return p[0] + sum;
		}
};

void pulsed_per_file(char filepath[500], Double_t gammaflash_min, Double_t gammaflash_center, Double_t gammaflash_max, Double_t neutronresponse_min, Double_t neutronresponse_max, Double_t max_param, Double_t gamma_background, Double_t distance, Double_t max_energy){
	EnableImplicitMT();
	RDataFrame d("Data", filepath);

	Double_t gcenter = gammaflash_center-distance/C/TOF_TO_S;
	auto monster = d.Filter("Channel==4").Define("newtof",[gcenter](const RVec<Double_t>& tof) {return tof-gcenter;},{"tof"}).Define("neutron_energy",[distance](const RVec<Double_t>& newtof) {RVec<Double_t> v=distance/(newtof*TOF_TO_S);return NEUTRON_MASS/2*v*v*J_TO_KEV;},{"newtof"});
	auto tof_plot = monster.Histo1D({"tof_plot", ";ToF;Counts", 1000, MIN_TOF, MAX_TOF}, "tof");
	auto psd_plot = monster.Histo1D({"psd_plot", ";psd;Counts", 1000, 0, 1}, "psd");
	auto tof_id_plot = monster.Histo2D({"tof_id_plot", ";ToF;PSD;Counts", 100, MIN_TOF-gcenter, MAX_TOF-gcenter, 100, 0, 1}, "newtof", "psd");
	auto energy_id_plot = monster.Histo2D({"energy_id_plot", ";Energy;PSD;Counts", 4096/4, 0, 4096, 100, 0, 1}, "Energy", "psd");
	Double_t y1 = 0.3;
	Double_t y2 = 0.42;
	Double_t x1 = 2100;
	Double_t x2 = 2600;
	Double_t m = (y2-y1)/(x2-x1);
	Double_t n = y1-x1*m;
	auto monster_efiltered = monster.Filter("Energy>PULSE_ENERGY_FILTER");
	auto monster_gammas = monster_efiltered.Filter([m, n, y1](Double_t psd, UShort_t Energy){return psd<y1 || psd<m*Energy+n;},{"psd", "Energy"});
	auto gamma_tof_plot = monster_gammas.Histo1D({"gamma_tof_plot", ";ToF;Counts", 1000, MIN_TOF-gcenter, MAX_TOF-gcenter}, "newtof");
	auto monster_neutrons = monster_efiltered.Filter([m, n, y1](Double_t psd, UShort_t Energy){return psd>=y1 && psd>=m*Energy+n;},{"psd", "Energy"});
	auto neutron_tof_plot = monster_neutrons.Histo1D({"neutron_tof_plot", ";ToF;Counts", 1000, MIN_TOF-gcenter, MAX_TOF-gcenter}, "newtof");
	auto neutron_energy_plot = monster_neutrons.Histo1D({"neutron_energy_plot", ";Energy (keV);Counts", 1000, 0, 10000}, "neutron_energy");

	auto gamma_flash = monster_gammas.Histo1D({"gamma_flash", ";ToF;Counts", GAMMA_FLASH_BINS_N, GMIN_TOF, GMAX_TOF}, "newtof");
	auto neutron_response = monster_neutrons.Histo1D({"neutron_response", ";ToF;Counts", NEUTRON_RESPONSE_BINS_N, NMIN_TOF*distance, NMAX_TOF*distance}, "newtof");
	Double_t vmin = distance/(NMAX_TOF*TOF_TO_S);
	Double_t vmax = distance/(NMIN_TOF*TOF_TO_S);
	auto neutron_response_energy = monster_neutrons.Histo1D({"neutron_response_energy", ";Energy (keV);Counts", 1000, NEUTRON_MASS/2*vmin*vmin*J_TO_KEV, NEUTRON_MASS/2*vmax*vmax*J_TO_KEV}, "neutron_energy");

	TCanvas* myCanvas = new TCanvas("");

	tof_plot->Write("", TObject::kOverwrite);
	psd_plot->Write("", TObject::kOverwrite);
	tof_id_plot->Draw("COLZ");
	myCanvas->Write("tof_id_plot", TObject::kOverwrite);
	energy_id_plot->Draw("COLZ");
	TLine* psd_line_1 = new TLine(0, y1, x1, y1);
	TLine* psd_line_2 = new TLine(x1, y1, (1-n)/m, 1);
	psd_line_1->SetLineStyle(3);
	psd_line_2->SetLineStyle(3);
	psd_line_1->Draw("same");
	psd_line_2->Draw("same");
	myCanvas->Write("energy_id_plot", TObject::kOverwrite);

	gamma_tof_plot->Write("", TObject::kOverwrite);
	neutron_tof_plot->Write("", TObject::kOverwrite);
	neutron_energy_plot->Write("", TObject::kOverwrite);

	gamma_flash->Write("", TObject::kOverwrite);
	neutron_response->Write("", TObject::kOverwrite);
	neutron_response_energy->Write("", TObject::kOverwrite);

	//fit
	pulse_fit_functor pulse_fit_obj = pulse_fit_functor((TH1F*)gDirectory->Get("gamma_flash"),GMIN_TOF,NMIN_TOF*distance,NMAX_TOF*distance,MIN_PARAM_TOF*distance,MAX_PARAM_TOF*distance);
	TF1* pulse_fit = new TF1("pulse_fit", pulse_fit_obj, NMIN_TOF*distance, NMAX_TOF*distance, PULSE_FIT_PARAMS_N+2);
	pulse_fit->SetNpx(NEUTRON_RESPONSE_BINS_N);
	pulse_fit->SetNumberFitPoints(NEUTRON_RESPONSE_BINS_N);
	pulse_fit->SetParLimits(0, 0, 500);
	for(UShort_t i=1; i<PULSE_FIT_PARAMS_N+1; i++){
		pulse_fit->SetParLimits(i, 0, max_param);
		pulse_fit->SetParameter(i, 0);
	}
	pulse_fit->SetParLimits(PULSE_FIT_PARAMS_N+1, 0, 500);
	pulse_fit->SetParameter(PULSE_FIT_PARAMS_N+1, gamma_background);
//	pulse_fit->FixParameter(PULSE_FIT_PARAMS_N+1, gamma_background);
	TFitResultPtr fitresult = neutron_response->Fit("pulse_fit", "SLE");
	myCanvas->Write("pulse_fit_plot", TObject::kOverwrite);

	//guardar resultados en tree
	Double_t results[PULSE_FIT_PARAMS_N+2][2];
	for(UShort_t i=0; i<PULSE_FIT_PARAMS_N; i++){
		results[i][0] = fitresult->Parameter(i);
		results[i][1] = fitresult->ParError(i);
	}
	TTree* results_tree = new TTree("results_tree", "Tree with pulsed results");
	results_tree->Branch("results", results, "results[12][2]/D");	//PULSE_FIT_PARAMS_N
	results_tree->SetBranchAddress("results", results);
	results_tree->Fill();
	results_tree->Write("", TObject::kOverwrite);

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
	pulsed_per_file(pulsed_1, PULSED1_GMIN, PULSED1_GCENTER, PULSED1_GMAX, PULSED1_NMIN, PULSED1_NMAX, 1E-4, 20, 1, 5500);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_2");
	pulsed_per_file(pulsed_2, PULSED2_GMIN, PULSED2_GCENTER, PULSED2_GMAX, PULSED2_NMIN, PULSED2_NMAX, 1E-3, 20, 1, 5500);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_3");
	pulsed_per_file(pulsed_3, PULSED3_GMIN, PULSED3_GCENTER, PULSED3_GMAX, PULSED3_NMIN, PULSED3_NMAX, 1E-3, 40, 1, 7000);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_4");
	pulsed_per_file(pulsed_4, PULSED4_GMIN, PULSED4_GCENTER, PULSED4_GMAX, PULSED4_NMIN, PULSED4_NMAX, 1E-2, 10, 1, 8250);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_5");
	pulsed_per_file(pulsed_5, PULSED5_GMIN, PULSED5_GCENTER, PULSED5_GMAX, PULSED5_NMIN, PULSED5_NMAX, 1E-1, 35, 2, 8250);
	gDirectory->cd("..");

	gDirectory->cd("..");
	f.Close();
}

TGraph* pulsed_results_per_file(Double_t g_min, Double_t g_center, Double_t g_max, Double_t n_min, Double_t n_max, Double_t distance, Double_t alpha_energy, Double_t max_energy){
	Double_t results[PULSE_FIT_PARAMS_N+2][2];
	TTree* tree = (TTree*)gDirectory->Get("results_tree");
	tree->SetBranchAddress("results", results);
	tree->GetEntry(0);
	Double_t p[PULSE_FIT_PARAMS_N+2];
	Double_t p_err[PULSE_FIT_PARAMS_N+2];
	for(UShort_t i=0; i<PULSE_FIT_PARAMS_N; i++){
		p[i] = results[i+1][0];
		p_err[i] = results[i+1][1];
	}

	Double_t x[PULSE_FIT_PARAMS_N];
	Double_t y[PULSE_FIT_PARAMS_N];
	Double_t y_err[PULSE_FIT_PARAMS_N];
	Double_t paramwidth = (MAX_PARAM_TOF-MIN_PARAM_TOF)*distance/PULSE_FIT_PARAMS_N * TOF_TO_S;
	Double_t first_param_seconds = MIN_PARAM_TOF*distance*TOF_TO_S + distance/C;
	cout << "first_param_seconds: " << first_param_seconds << endl;
	for(UShort_t i=0; i<PULSE_FIT_PARAMS_N; i++){
		x[i] = first_param_seconds + paramwidth*(i+0.5);
		y[i] = p[i];
		y_err[i] = p_err[i];
	}

	TCanvas* myCanvas = new TCanvas("");
	TGraph* cross_section_result = new TGraphErrors(PULSE_FIT_PARAMS_N, x, y, NULL, y_err);
	cross_section_result->SetTitle("Fit results;ToF (s);Counts");
	cross_section_result->SetMarkerStyle(21);
	cross_section_result->Draw("alp");
	myCanvas->Write("pulse_fit_results", TObject::kOverwrite);

	TH1D* gamma_flash = (TH1D*)gDirectory->Get("gamma_flash");
	Double_t gammas_n = gamma_flash->Integral(0, gamma_flash->GetNbinsX());
	//histograma "delta de dirac" con una sola cuenta
	TH1D* delta_histogram = new TH1D("centered dirac delta", "Dirac delta;ToF;Counts", GAMMA_FLASH_BINS_N, g_min, g_max);
	delta_histogram->Fill(g_center, gammas_n);
	delta_histogram->Write("dirac_delta", TObject::kOverwrite);
	//functor with delta
	pulse_fit_functor pulse_functor = pulse_fit_functor((TH1F*)gDirectory->Get("dirac_delta"),g_min,n_min,n_max,MIN_PARAM_TOF*distance,MAX_PARAM_TOF*distance);
	Double_t x2[200];
	Double_t y2[200];
	for(UShort_t i=0; i<200; i++){
		x2[i] = n_min + (n_max-n_min)/200*(i+0.5);
		y2[i] = pulse_functor(&x2[i], p);
	}
	TGraph* result_for_delta = new TGraphErrors(200, x2, y2, NULL, NULL);
	result_for_delta->SetTitle("delta; ToF (ns);Counts per gamma");
	result_for_delta->SetMarkerStyle(21);
	result_for_delta->Draw("alp");
	((TH1D*)gDirectory->Get("neutron_response"))->Draw("same");
	myCanvas->Write("response_to_delta", TObject::kOverwrite);
	//functor with gamma flash
	pulse_fit_functor pulse_functor_2 = pulse_fit_functor((TH1F*)gDirectory->Get("gamma_flash"),g_min,n_min,n_max,MIN_PARAM_TOF*distance,MAX_PARAM_TOF*distance);
	for(UShort_t i=0; i<200; i++){
		x2[i] = n_min + (n_max-n_min)/200*(i+0.5);
		y2[i] = pulse_functor_2(&x2[i], p);
	}
	TGraph* result_for_gflash = new TGraphErrors(200, x2, y2, NULL, NULL);
	result_for_gflash->SetTitle("delta; ToF (ns);Counts per gamma");
	result_for_gflash->SetMarkerStyle(21);
	result_for_gflash->Draw("alp");
	((TH1D*)gDirectory->Get("neutron_response"))->Draw("same");
	myCanvas->Write("response_to_gflash", TObject::kOverwrite);

	for(UShort_t i=0; i<PULSE_FIT_PARAMS_N; i++){
		x[i] = first_param_seconds + paramwidth*(i+0.5);
		y[i] = p[i]*gammas_n;
		y_err[i] = p_err[i]*gammas_n;
	}
	TGraph* cross_section_result_delta = new TGraphErrors(PULSE_FIT_PARAMS_N, x, y, NULL, y_err);
	cross_section_result_delta->SetTitle("delta;ToF (s);Counts");
	cross_section_result_delta->SetMarkerStyle(21);
	cross_section_result_delta->Draw("alp");
	((TH1D*)gDirectory->Get("neutron_response"))->Draw("same");
	myCanvas->Write("neutron_response+delta", TObject::kOverwrite);

	//energy
	Double_t v;
	for(UShort_t i=0; i<PULSE_FIT_PARAMS_N; i++){
		v = distance/(first_param_seconds + paramwidth*(i+0.5));
		x[i] = (NEUTRON_MASS/2*v*v)*J_TO_KEV;
		y[i] = p[i];
		y_err[i] = p_err[i];
	}
	cout << "min tof energy: " << x[0] << endl;
	cout << "max tof energy: " << x[PULSE_FIT_PARAMS_N-1] << endl;
	TGraph* energy_result= new TGraphErrors(PULSE_FIT_PARAMS_N, x, y, NULL, y_err);
	energy_result->SetTitle("Energy result;Energy (keV);Counts");
	energy_result->SetMarkerStyle(21);
	energy_result->Draw("alp");
	Double_t max_y = *max_element(y,y+PULSE_FIT_PARAMS_N-1);
	TLine* energy_line = new TLine(alpha_energy-2642.41, -0.1*max_y, alpha_energy-2642.41, 1.1*max_y);
	energy_line->SetLineStyle(3);
	energy_line->Draw("same");
	myCanvas->Write("energy_result", TObject::kOverwrite);
	myCanvas->Close();
	cout << "------------" << endl;
	return energy_result;
}

void pulsed_results(){
	TFile f("output.root", "UPDATE");
	gDirectory->cd("Pulsed");

	gDirectory->cd("pulsed_1");
	TGraph* energy_1 = pulsed_results_per_file(PULSED1_GMIN,PULSED1_GCENTER,PULSED1_GMAX,PULSED1_NMIN,PULSED1_NMAX, 1, 5500, 5500);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_2");
	TGraph* energy_2 = pulsed_results_per_file(PULSED2_GMIN,PULSED2_GCENTER,PULSED2_GMAX,PULSED2_NMIN,PULSED2_NMAX, 1, 5500, 5500);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_3");
	TGraph* energy_3 = pulsed_results_per_file(PULSED3_GMIN,PULSED3_GCENTER,PULSED3_GMAX,PULSED3_NMIN,PULSED3_NMAX, 1, 7000, 7000);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_4");
	TGraph* energy_4 = pulsed_results_per_file(PULSED4_GMIN,PULSED4_GCENTER,PULSED4_GMAX,PULSED4_NMIN,PULSED4_NMAX, 1, 8250, 8250);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_5");
	TGraph* energy_5 = pulsed_results_per_file(PULSED5_GMIN,PULSED5_GCENTER,PULSED5_GMAX,PULSED5_NMIN,PULSED5_NMAX, 2, 8250, 8250);
	gDirectory->cd("..");

	TCanvas* myCanvas = new TCanvas("");
	TMultiGraph* energy_results = new TMultiGraph();
	energy_1->SetTitle("5500keV");
	energy_1->SetMarkerColor(kRed);
	energy_1->SetLineColor(kRed);
	energy_results->Add(energy_1);
	energy_2->SetTitle("5500keV");
	energy_2->SetMarkerColor(kBlue);
	energy_2->SetLineColor(kBlue);
	energy_results->Add(energy_2);
	energy_3->SetTitle("7000keV");
	energy_3->SetMarkerColor(kViolet);
	energy_3->SetLineColor(kViolet);
	energy_results->Add(energy_3);
	energy_4->SetTitle("8250keV");
	energy_4->SetMarkerColor(kGreen);
	energy_4->SetLineColor(kGreen);
	energy_results->Add(energy_4);
	energy_5->SetTitle("8250keV, 2m");
	energy_5->SetMarkerColor(kBlack);
	energy_5->SetLineColor(kBlack);
	energy_results->Add(energy_5);

	energy_results->Draw("ALP");
	myCanvas->BuildLegend();
	myCanvas->Write("energy_results", TObject::kOverwrite);
	myCanvas->Close();

	gDirectory->cd("..");
	f.Close();
}
