#define C 299792458
#define NEUTRON_MASS 1.67492749804E-27
#define J_TO_KEV (1/1.602177E-16)
#define TOF_TO_S (1E-9)

#define MIN_TOF 0
#define MAX_TOF 1500
#define PULSE_ENERGY_FILTER 69	//482/700 keV/channel
#define GAMMA_FLASH_BINS_N 400
#define NEUTRON_RESPONSE_BINS_N 600
#define PULSE_NPX 600
#define PULSE_FIT_PARAMS_N 60	//must be even; if changed, must also change tree size definition
#define MAX_PARAM_E 10000
#define MIN_PARAM_E 400
#define MIN_EFF 15

#define GMIN_TOF -15
#define GMAX_TOF +100
#define NMIN_TOF 1
#define NMAX_TOF 301

#define PULSED1_PARMAX 1
#define PULSED2_PARMAX 1
#define PULSED3_PARMAX 1
#define PULSED4_PARMAX 1
#define PULSED5_PARMAX 1

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

Double_t tof_to_energy(Double_t tof, Double_t distance){	//returns keV, input in ns
	return 0.5*NEUTRON_MASS*distance*distance/(tof*tof*TOF_TO_S*TOF_TO_S)*J_TO_KEV;
}
Double_t energy_to_tof(Double_t energy, Double_t distance){	//returns ns, input in keV
	return sqrt(NEUTRON_MASS*0.5/(energy/J_TO_KEV))*distance/TOF_TO_S;
}

class pulse_fit_functor{
	private:
		TH1D* gamma_histo;
		Double_t gamma_stepsize;
	public:
		Double_t gamma_min;
		Double_t distance;
		Double_t min_param_tofdiff;
		Double_t max_param_tofdiff;
		pulse_fit_functor(TH1D* gamma_histo, Double_t gamma_min, Double_t distance, Double_t min_param_tof, Double_t max_param_tof):gamma_histo(gamma_histo),gamma_min(gamma_min),distance(distance){
			gamma_stepsize = gamma_histo->GetBinWidth(1);
			min_param_tofdiff = min_param_tof-distance/C/TOF_TO_S;
			max_param_tofdiff = max_param_tof-distance/C/TOF_TO_S;
		}
		Double_t operator()(Double_t* x, Double_t* p){
			Double_t sum = 0;
			Short_t param_index = 0;
			Double_t tofdiff = 0;
			for(Short_t i=0; i<GAMMA_FLASH_BINS_N; i++){
				tofdiff = x[0] - (gamma_min+(0.5+i)*gamma_stepsize);
				if(tofdiff<min_param_tofdiff || tofdiff>max_param_tofdiff){
					continue;
				}
				param_index = 1 + (tofdiff-min_param_tofdiff)/(max_param_tofdiff-min_param_tofdiff)*PULSE_FIT_PARAMS_N;
				if(param_index>0 && param_index<=PULSE_FIT_PARAMS_N){
					sum += gamma_histo->GetBinContent(i+1)>p[PULSE_FIT_PARAMS_N+1] ? (gamma_histo->GetBinContent(i+1)-p[PULSE_FIT_PARAMS_N+1]) * p[param_index] : 0;
				}
			}
			return p[0] + sum;
		}
};

void pulsed_per_file(char filepath[500], Double_t gammaflash_min, Double_t gammaflash_center, Double_t gammaflash_max, Double_t neutronresponse_min, Double_t neutronresponse_max, Double_t max_param, Double_t gamma_background, Double_t distance, Double_t max_energy){
	EnableImplicitMT();
	RDataFrame d("Data", filepath);

	Double_t min_param_tof = energy_to_tof(MAX_PARAM_E, distance);
	Double_t max_param_tof = energy_to_tof(MIN_PARAM_E, distance);

	Double_t gcenter = gammaflash_center-distance/C/TOF_TO_S;
	auto monster = d.Filter("Channel==4").Define("newtof",[gcenter](const RVec<Double_t>& tof) {return tof-gcenter;},{"tof"}).Define("neutron_energy",[distance](const RVec<Double_t>& newtof) {Double_t v=distance/(newtof[0]*TOF_TO_S);return NEUTRON_MASS/2*v*v*J_TO_KEV;},{"newtof"});
	auto tof_plot = monster.Histo1D({"tof_plot", ";ToF (ns);Counts", 3000, MIN_TOF, MAX_TOF}, "tof");
	auto psd_plot = monster.Histo1D({"psd_plot", ";psd;Counts", 1000, 0, 1}, "psd");
	auto tof_id_plot = monster.Histo2D({"tof_id_plot", ";ToF (ns);PSD;Counts", 100, MIN_TOF-gcenter, MAX_TOF-gcenter, 100, 0, 1}, "newtof", "psd");
	auto energy_id_plot = monster.Histo2D({"energy_id_plot", ";Energy;PSD;Counts", 4096/4, 0, 4096, 100, 0, 1}, "Energy", "psd");
	Double_t y1 = 0.3;
	Double_t y2 = 0.42;
	Double_t x1 = 2100;
	Double_t x2 = 2600;
	Double_t m = (y2-y1)/(x2-x1);
	Double_t n = y1-x1*m;
	auto monster_efiltered = monster.Filter("Energy>PULSE_ENERGY_FILTER");
	auto monster_gammas = monster_efiltered.Filter([m, n, y1](Double_t psd, UShort_t Energy){return psd<y1 || psd<m*Energy+n;},{"psd", "Energy"});
	auto gamma_tof_plot = monster_gammas.Histo1D({"gamma_tof_plot", ";ToF (ns);Counts", 1000, MIN_TOF-gcenter, MAX_TOF-gcenter}, "newtof");
	auto monster_neutrons = monster_efiltered.Filter([m, n, y1](Double_t psd, UShort_t Energy){return psd>=y1 && psd>=m*Energy+n;},{"psd", "Energy"});
	auto monster_neutrons_efiltered = monster_neutrons.Filter([](Double_t neutron_energy){return neutron_energy>MIN_PARAM_E;},{"neutron_energy"});
	auto neutron_tof_plot = monster_neutrons_efiltered.Histo1D({"neutron_tof_plot", ";ToF (ns);Counts", 1000, MIN_TOF-gcenter, MAX_TOF-gcenter}, "newtof");
	auto neutron_tof_plot_raw = monster_neutrons.Histo1D({"neutron_tof_plot_raw", ";ToF (ns);Counts", 1000, MIN_TOF-gcenter, MAX_TOF-gcenter}, "newtof");
	auto neutron_energy_plot = monster_neutrons_efiltered.Histo1D({"neutron_energy_plot", ";Energy (keV);Counts", 1000, 0, 10000}, "neutron_energy");
	auto neutron_energy_plot_raw = monster_neutrons.Histo1D({"neutron_energy_plot_raw", ";Energy (keV);Counts", 1000, 0, 10000}, "neutron_energy");

	auto gamma_flash = monster_gammas.Histo1D({"gamma_flash", ";ToF (ns);Counts", GAMMA_FLASH_BINS_N, GMIN_TOF, GMAX_TOF}, "newtof");
	auto neutron_response = monster_neutrons_efiltered.Histo1D({"neutron_response", ";ToF (ns);Counts", NEUTRON_RESPONSE_BINS_N, NMIN_TOF*distance, max_param_tof+distance/C/TOF_TO_S}, "newtof");
	auto neutron_response_raw = monster_neutrons.Histo1D({"neutron_response_raw", ";ToF (ns);Counts", NEUTRON_RESPONSE_BINS_N, NMIN_TOF*distance, max_param_tof+distance/C/TOF_TO_S}, "newtof");
	auto neutron_response_energy = monster_neutrons_efiltered.Histo1D({"neutron_response_energy", ";Energy (keV);Counts", 1000, 0, 10000}, "neutron_energy");
	auto neutron_response_energy_raw = monster_neutrons.Histo1D({"neutron_response_energy_raw", ";Energy (keV);Counts", 1000, 0, 10000}, "neutron_energy");

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
	neutron_tof_plot_raw->Write("", TObject::kOverwrite);
	neutron_energy_plot->Write("", TObject::kOverwrite);
	neutron_energy_plot_raw->Write("", TObject::kOverwrite);

	gamma_flash->Write("", TObject::kOverwrite);
	neutron_response->Write("", TObject::kOverwrite);
	neutron_response_raw->Write("", TObject::kOverwrite);
	neutron_response_energy->Write("", TObject::kOverwrite);
	neutron_response_energy_raw->Write("", TObject::kOverwrite);

	//fit
	pulse_fit_functor pulse_fit_obj = pulse_fit_functor((TH1D*)gDirectory->Get("gamma_flash"),GMIN_TOF,distance,min_param_tof,max_param_tof);
	TF1* pulse_fit = new TF1("pulse_fit", pulse_fit_obj, NMIN_TOF*distance, max_param_tof+distance/C/TOF_TO_S, PULSE_FIT_PARAMS_N+2);
	pulse_fit->SetNpx(PULSE_NPX);
	pulse_fit->SetNumberFitPoints(PULSE_NPX);
	pulse_fit->SetParLimits(0, 0, 10);
	for(UShort_t i=1; i<PULSE_FIT_PARAMS_N+1; i++){
		pulse_fit->SetParLimits(i, 0, max_param);
		pulse_fit->SetParameter(i, 0);
	}
	pulse_fit->SetParLimits(PULSE_FIT_PARAMS_N+1, 0.65*gamma_background, gamma_background*1.35);
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
	results_tree->Branch("results", results, "results[62][2]/D");	//PULSE_FIT_PARAMS_N
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
//	pulsed_per_file(pulsed_1, PULSED1_GMIN, PULSED1_GCENTER, PULSED1_GMAX, PULSED1_NMIN, PULSED1_NMAX, PULSED1_PARMAX, 115, 1, 5500);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_2");
	pulsed_per_file(pulsed_2, PULSED2_GMIN, PULSED2_GCENTER, PULSED2_GMAX, PULSED2_NMIN, PULSED2_NMAX, PULSED2_PARMAX, 120, 1, 5500);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_3");
//	pulsed_per_file(pulsed_3, PULSED3_GMIN, PULSED3_GCENTER, PULSED3_GMAX, PULSED3_NMIN, PULSED3_NMAX, PULSED3_PARMAX, 155, 1, 7000);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_4");
//	pulsed_per_file(pulsed_4, PULSED4_GMIN, PULSED4_GCENTER, PULSED4_GMAX, PULSED4_NMIN, PULSED4_NMAX, PULSED4_PARMAX, 45, 1, 8250);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_5");
//	pulsed_per_file(pulsed_5, PULSED5_GMIN, PULSED5_GCENTER, PULSED5_GMAX, PULSED5_NMIN, PULSED5_NMAX, PULSED5_PARMAX, 135, 2, 8250);
	gDirectory->cd("..");

	gDirectory->cd("..");
	f.Close();
}

TGraph* pulsed_results_per_file(Double_t g_min, Double_t g_center, Double_t g_max, Double_t n_min, Double_t n_max, Double_t distance, Double_t alpha_energy, Double_t max_energy){
	Double_t min_param_tof = energy_to_tof(MAX_PARAM_E, distance);
	Double_t max_param_tof = energy_to_tof(MIN_PARAM_E, distance);
	Double_t paramwidth_tof = (max_param_tof-min_param_tof)/PULSE_FIT_PARAMS_N;

	Double_t results[PULSE_FIT_PARAMS_N+2][2];
	TTree* tree = (TTree*)gDirectory->Get("results_tree");
	tree->SetBranchAddress("results", results);
	tree->GetEntry(0);
	Double_t p[PULSE_FIT_PARAMS_N];
	Double_t p_err[PULSE_FIT_PARAMS_N];
	for(UShort_t i=0; i<PULSE_FIT_PARAMS_N; i++){
		p[i] = results[i+1][0];
		p_err[i] = results[i+1][1];
	}

	cout << "first param tof: " << energy_to_tof(MAX_PARAM_E,distance)<< endl;
	cout << "last param tof: " << energy_to_tof(MIN_PARAM_E,distance) << endl;

	Double_t x[PULSE_FIT_PARAMS_N];
	Double_t y[PULSE_FIT_PARAMS_N];
	Double_t y_err[PULSE_FIT_PARAMS_N];
	TH1D* gamma_flash = (TH1D*)gDirectory->Get("gamma_flash");
	Double_t gammas_n = gamma_flash->Integral(0, gamma_flash->GetNbinsX());
	TCanvas* myCanvas = new TCanvas("");
	//represent results with neutron response
	for(UShort_t i=0; i<PULSE_FIT_PARAMS_N; i++){
		x[i] = paramwidth_tof*i+min_param_tof;
		y[i] = p[i]*gammas_n;
		y_err[i] = p_err[i]*gammas_n;
	}
	TGraph* cross_section_result_delta = new TGraphErrors(PULSE_FIT_PARAMS_N, x, y, NULL, y_err);
	cross_section_result_delta->SetTitle("delta;ToF (ns);Counts");
	cross_section_result_delta->SetMarkerStyle(21);
	cross_section_result_delta->Draw("alp");
	((TH1D*)gDirectory->Get("neutron_response"))->Draw("same");
	myCanvas->Write("neutron_response+delta", TObject::kOverwrite);

	TGraph* cross_section_result = new TGraphErrors(PULSE_FIT_PARAMS_N, x, y, NULL, y_err);
	cross_section_result->SetTitle("Fit results;ToF (ns);Counts");
	cross_section_result->SetMarkerStyle(21);
	cross_section_result->Draw("alp");
	myCanvas->Write("pulse_fit_results", TObject::kOverwrite);

	//histograma "delta de dirac" con una sola cuenta
	TH1D* delta_histogram = new TH1D("centered dirac delta", "Dirac delta;ToF;Counts", GAMMA_FLASH_BINS_N, GMIN_TOF, GMAX_TOF);
	delta_histogram->Fill(distance/C/TOF_TO_S, gammas_n);
	delta_histogram->Write("dirac_delta", TObject::kOverwrite);
	//functor with delta
	pulse_fit_functor pulse_functor = pulse_fit_functor((TH1D*)gDirectory->Get("dirac_delta"),GMIN_TOF,distance,min_param_tof,max_param_tof);
	Double_t x2[800];
	Double_t y2[800];
	for(UShort_t i=0; i<800; i++){
		x2[i] = NMIN_TOF + (NMAX_TOF-NMIN_TOF)/800*i;
		y2[i] = pulse_functor(&x2[i], p);
	}
	TGraph* result_for_delta = new TGraphErrors(800, x2, y2, NULL, NULL);
	result_for_delta->SetTitle("delta; ToF (ns);Counts per gamma");
	result_for_delta->SetMarkerStyle(21);
	result_for_delta->Draw("alp");
	((TH1D*)gDirectory->Get("neutron_response"))->Draw("same");
	myCanvas->Write("response_to_delta", TObject::kOverwrite);
	//functor with gamma flash
	pulse_fit_functor pulse_functor_2 = pulse_fit_functor((TH1D*)gDirectory->Get("gamma_flash"),GMIN_TOF,distance,min_param_tof,max_param_tof);
	for(UShort_t i=0; i<800; i++){
		y2[i] = pulse_functor_2(&x2[i], p);
	}
	TGraph* result_for_gflash = new TGraphErrors(800, x2, y2, NULL, NULL);
	result_for_gflash->SetTitle("delta; ToF (ns);Counts per gamma");
	result_for_gflash->SetMarkerStyle(21);
	result_for_gflash->Draw("alp");
	((TH1D*)gDirectory->Get("neutron_response"))->Draw("same");
	myCanvas->Write("response_to_gflash", TObject::kOverwrite);


	//scale because of efficiency
	Double_t eff_energy[19] = {50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000};
	Double_t eff_30[19] = {0.0159, 0.0136, 0.0155, 0.0129, 40.8395, 53.697, 57.0457, 57.6266, 57.3337, 57.3179, 57.4631, 52.4621, 48.26, 45.2324, 42.2129, 40.6158, 38.2269, 35.8656, 33.607};
	Double_t eff_100[19] = {0.0142, 0.0128, 0.0146, 0.0121, 0.0098, 0.0086, 0.0087, 1.6957, 11.2538, 21.4897, 30.1175, 41.4082, 40.7344, 39.4135, 37.5031, 36.5453, 34.8368, 33.0715, 31.1554};
	Double_t eff_250[19] = {0.0134, 0.0116, 0.0135, 0.0117, 0.009, 0.008, 0.0082, 0.0077, 0.0081, 0.0064, 0.0058, 12.7651, 27.1103, 30.7212, 31.0754, 30.9491, 30.2988, 29.3758, 28.042};
	Double_t eff = 0;
	Double_t energy = 0;
	for(unsigned short i=0; i<PULSE_FIT_PARAMS_N; i++){
		energy = tof_to_energy(min_param_tof+paramwidth_tof*i,distance);
		if(energy<=eff_energy[0]){
			eff = energy*eff_100[0]/eff_energy[0];
		}else if (energy <= eff_energy[18]){
			for(unsigned short j=1; j<19; j++){
				if(energy<=eff_energy[j] && energy>eff_energy[j-1]){
					eff = eff_100[j-1] + (energy-eff_energy[j-1])*(eff_100[j]-eff_100[j-1])/(eff_energy[j]-eff_energy[j-1]);
				}
			}
		}else{
			eff = eff_100[18];
		}
		if(eff<MIN_EFF){
			p[i] = 0;
			p_err[i] = 0;
		}else{
			p[i] /= eff/100;
			p_err[i] /= eff/100;
		}
	}


	for(UShort_t i=0; i<PULSE_FIT_PARAMS_N; i++){
		y[i] = p[i];
		y_err[i] = p_err[i];
	}

	//energy
	Double_t v;
	for(UShort_t i=0; i<PULSE_FIT_PARAMS_N; i++){
		x[i] = tof_to_energy(min_param_tof+paramwidth_tof*i,distance);
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
	TH1D* energysimple_1 = (TH1D*)gDirectory->Get("neutron_response_energy");
	Double_t ngammas_1 = ((TH1D*)gDirectory->Get("gamma_flash"))->Integral();
	gDirectory->cd("..");

	gDirectory->cd("pulsed_2");
	TGraph* energy_2 = pulsed_results_per_file(PULSED2_GMIN,PULSED2_GCENTER,PULSED2_GMAX,PULSED2_NMIN,PULSED2_NMAX, 1, 5500, 5500);
	TH1D* energysimple_2 = (TH1D*)gDirectory->Get("neutron_response_energy");
	Double_t ngammas_2 = ((TH1D*)gDirectory->Get("gamma_flash"))->Integral();
	gDirectory->cd("..");

	gDirectory->cd("pulsed_3");
	TGraph* energy_3 = pulsed_results_per_file(PULSED3_GMIN,PULSED3_GCENTER,PULSED3_GMAX,PULSED3_NMIN,PULSED3_NMAX, 1, 7000, 7000);
	TH1D* energysimple_3 = (TH1D*)gDirectory->Get("neutron_response_energy");
	Double_t ngammas_3 = ((TH1D*)gDirectory->Get("gamma_flash"))->Integral();
	gDirectory->cd("..");

	gDirectory->cd("pulsed_4");
	TGraph* energy_4 = pulsed_results_per_file(PULSED4_GMIN,PULSED4_GCENTER,PULSED4_GMAX,PULSED4_NMIN,PULSED4_NMAX, 1, 8250, 8250);
	TH1D* energysimple_4 = (TH1D*)gDirectory->Get("neutron_response_energy");
	Double_t ngammas_4 = ((TH1D*)gDirectory->Get("gamma_flash"))->Integral();
	gDirectory->cd("..");

	gDirectory->cd("pulsed_5");
	TGraph* energy_5 = pulsed_results_per_file(PULSED5_GMIN,PULSED5_GCENTER,PULSED5_GMAX,PULSED5_NMIN,PULSED5_NMAX, 2, 8250, 8250);
	TH1D* energysimple_5 = (TH1D*)gDirectory->Get("neutron_response_energy");
	Double_t ngammas_5 = ((TH1D*)gDirectory->Get("gamma_flash"))->Integral();
	gDirectory->cd("..");

	//divide by number of alphas
	Double_t charge_to_alpha = 1/(2*1.60217646E-13);	//microcoulomb to alpha
	energysimple_1->Scale(	1/(1.26*charge_to_alpha));
	energy_1->Scale(	1/(1.26*charge_to_alpha));
	energysimple_2->Scale(	1/(7.08*charge_to_alpha));
	energy_2->Scale(	1/(7.08*charge_to_alpha));
	energysimple_3->Scale(	1/(5.07*charge_to_alpha));
	energy_3->Scale(	1/(5.07*charge_to_alpha));
	energysimple_4->Scale(	1/(1.03*charge_to_alpha));
	energy_4->Scale(	1/(1.03*charge_to_alpha));
	energysimple_5->Scale(	1/(4.20*charge_to_alpha));
	energy_5->Scale(	1/(4.20*charge_to_alpha));

	//multiply by number of gammas
	energy_1->Scale(	ngammas_1);
	energy_2->Scale(	ngammas_2);
	energy_3->Scale(	ngammas_3);
	energy_4->Scale(	ngammas_4);
	energy_5->Scale(	ngammas_5);

	//scale beacause of distance
	Double_t detector_surface = M_PI*0.10*0.10; 	//m2
	energysimple_1->Scale(	4*M_PI*1*1/detector_surface);
	energy_1->Scale(	4*M_PI*1*1/detector_surface);
	energysimple_2->Scale(	4*M_PI*1*1/detector_surface);
	energy_2->Scale(	4*M_PI*1*1/detector_surface);
	energysimple_3->Scale(	4*M_PI*1*1/detector_surface);
	energy_3->Scale(	4*M_PI*1*1/detector_surface);
	energysimple_4->Scale(	4*M_PI*1*1/detector_surface);
	energy_4->Scale(	4*M_PI*1*1/detector_surface);
	energysimple_5->Scale(	4*M_PI*2*2/detector_surface);
	energy_5->Scale(	4*M_PI*2*2/detector_surface);

	//scale because of efficiency
	Double_t eff_energy[19] = {50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000};
	Double_t eff_30[19] = {0.0159, 0.0136, 0.0155, 0.0129, 40.8395, 53.697, 57.0457, 57.6266, 57.3337, 57.3179, 57.4631, 52.4621, 48.26, 45.2324, 42.2129, 40.6158, 38.2269, 35.8656, 33.607};
	Double_t eff_100[19] = {0.0142, 0.0128, 0.0146, 0.0121, 0.0098, 0.0086, 0.0087, 1.6957, 11.2538, 21.4897, 30.1175, 41.4082, 40.7344, 39.4135, 37.5031, 36.5453, 34.8368, 33.0715, 31.1554};
	Double_t eff_250[19] = {0.0134, 0.0116, 0.0135, 0.0117, 0.009, 0.008, 0.0082, 0.0077, 0.0081, 0.0064, 0.0058, 12.7651, 27.1103, 30.7212, 31.0754, 30.9491, 30.2988, 29.3758, 28.042};
	Double_t eff = 0;
	Double_t energy = 0;
	for(short i=0; i<1000; i++){
		energy += 10;	//10000keV max entre 1000 bines
		if(energy<=eff_energy[0]){
			eff = energy*eff_100[0]/eff_energy[0];
		}else if(energy<=eff_energy[18]){
			for(unsigned short j=1; j<19; j++){
				if(energy<=eff_energy[j] && energy>eff_energy[j-1]){
					eff = eff_100[j-1] + (energy-eff_energy[j-1])*(eff_100[j]-eff_100[j-1])/(eff_energy[j]-eff_energy[j-1]);
				}
			}
		}else{
			eff = eff_100[18];
		}
		if(eff<MIN_EFF){
			energysimple_1->SetBinContent(i+1, 0);
			energysimple_2->SetBinContent(i+1, 0);
			energysimple_3->SetBinContent(i+1, 0);
			energysimple_4->SetBinContent(i+1, 0);
			energysimple_5->SetBinContent(i+1, 0);
		}else{
			energysimple_1->SetBinContent(i+1, energysimple_1->GetBinContent(i+1)*100/eff);
			energysimple_2->SetBinContent(i+1, energysimple_2->GetBinContent(i+1)*100/eff);
			energysimple_3->SetBinContent(i+1, energysimple_3->GetBinContent(i+1)*100/eff);
			energysimple_4->SetBinContent(i+1, energysimple_4->GetBinContent(i+1)*100/eff);
			energysimple_5->SetBinContent(i+1, energysimple_5->GetBinContent(i+1)*100/eff);
		}
	}

	TCanvas* myCanvas = new TCanvas("");
	TGraph* efficiency_curve = new TGraph(19,eff_energy,eff_100);
	efficiency_curve->SetTitle("Efficiency curve;Energy (keV);Efficiency (%)");
	efficiency_curve->Draw("ALP");
	myCanvas->Write("efficiency_curve", TObject::kOverwrite);

	//Jacobs data
	Double_t jacobs_energies[119];
	for(short i=0; i<119; i++){
		jacobs_energies[i] =20 + i*20;
	}
	Double_t jacobs_yield_5500[119] = {6.62E-10,6.62E-10,6.62E-10,6.62E-10,6.62E-10,9.06E-10,9.06E-10,9.06E-10,9.06E-10,9.06E-10,4.36E-10,4.36E-10,4.36E-10,4.36E-10,4.36E-10,1.86E-10,3.30E-10,3.77E-10,3.55E-10,4.15E-10,3.98E-10,3.12E-10,1.70E-10,2.93E-10,4.71E-10,5.76E-10,5.35E-10,6.75E-10,6.97E-10,7.54E-10,7.29E-10,6.80E-10,8.96E-10,9.40E-10,6.88E-10,6.25E-10,6.18E-10,6.47E-10,6.02E-10,3.80E-10,2.97E-10,2.95E-10,1.90E-10,2.49E-10,2.56E-10,2.93E-10,3.45E-10,3.71E-10,3.48E-10,3.07E-10,3.86E-10,4.19E-10,5.04E-10,5.20E-10,5.21E-10,5.92E-10,7.03E-10,6.72E-10,5.21E-10,5.28E-10,5.49E-10,4.92E-10,5.20E-10,5.63E-10,5.69E-10,6.23E-10,6.21E-10,5.84E-10,5.53E-10,5.25E-10,6.29E-10,7.23E-10,7.46E-10,8.63E-10,9.34E-10,8.70E-10,7.99E-10,6.51E-10,5.21E-10,5.28E-10,5.87E-10,4.93E-10,5.12E-10,5.13E-10,5.18E-10,4.65E-10,3.88E-10,3.96E-10,4.61E-10,4.98E-10,5.47E-10,6.20E-10,6.62E-10,6.76E-10,6.57E-10,5.49E-10,4.85E-10,4.94E-10,5.12E-10,5.60E-10,5.57E-10,5.66E-10,6.48E-10,7.16E-10,7.07E-10,7.35E-10,7.13E-10,8.31E-10,1.05E-09,1.15E-09,1.02E-09,9.18E-10,8.21E-10,7.75E-10,6.30E-10,5.34E-10,4.71E-10,4.17E-10,3.97E-10};

	TGraph* jacobs_5500 = new TGraphErrors(119, jacobs_energies, jacobs_yield_5500, NULL, NULL);
	jacobs_5500->SetTitle("Jacobs data, 5500keV at 60deg");
	jacobs_5500->SetLineColor(kBlack);
	jacobs_5500->SetMarkerColor(kBlack);



//	---------------------------------------------------

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
	energy_4->SetMarkerColor(kGreen+1);
	energy_4->SetLineColor(kGreen+1);
	energy_results->Add(energy_4);
	energy_5->SetTitle("8250keV, 2m");
	energy_5->SetMarkerColor(kGreen+0.2);
	energy_5->SetLineColor(kGreen+0.2);
	energy_results->Add(energy_5);
	energy_results->Add(jacobs_5500);

	myCanvas->SetTitle(";Energy (keV);Neutrons per alpha");
	energy_results->Draw("ALP");
	myCanvas->BuildLegend();
	myCanvas->Write("energy_results", TObject::kOverwrite);

//	------------------------------------

	energysimple_1->SetTitle("5500keV");
	energysimple_1->SetLineColor(kRed);
	energysimple_1->Draw("HISTO");
	energysimple_2->SetTitle("5500keV");
	energysimple_2->SetLineColor(kBlue);
	energysimple_2->Draw("same HISTO");
	energysimple_3->SetTitle("7000keV");
	energysimple_3->SetLineColor(kViolet);
	energysimple_3->Draw("same HISTO");
	energysimple_4->SetTitle("8250keV");
	energysimple_4->SetLineColor(kGreen+1);
	energysimple_4->Draw("same HISTO");
	energysimple_5->SetTitle("8250keV, 2m");
	energysimple_5->SetLineColor(kGreen+0.2);
	energysimple_5->Draw("same HISTO");
	jacobs_5500->Draw("same");

	myCanvas->BuildLegend();
	myCanvas->SetTitle(";Energy (keV);Neutrons per alpha");
	myCanvas->Write("energysimple", TObject::kOverwrite);

//	------------------------------------

	energy_results->Draw("ALP");
	myCanvas->BuildLegend();
	energysimple_1->Draw("same HISTO");
	energysimple_2->Draw("same HISTO");
	energysimple_3->Draw("same HISTO");
	energysimple_4->Draw("same HISTO");
	energysimple_5->Draw("same HISTO");
	jacobs_5500->Draw("same");
	TLine* energy_line_1 = new TLine(2280, 0, 2280, 1);
	TLine* energy_line_2 = new TLine(3640, 0, 3640, 1);
	TLine* energy_line_3 = new TLine(4760, 0, 4760, 1);
	energy_line_1->SetLineStyle(3);
	energy_line_1->SetLineColor(kRed);
	energy_line_2->SetLineStyle(3);
	energy_line_2->SetLineColor(kViolet);
	energy_line_3->SetLineStyle(3);
	energy_line_3->SetLineColor(kGreen+1);
	energy_line_1->Draw("same");
	energy_line_2->Draw("same");
	energy_line_3->Draw("same");
	myCanvas->SetTitle(";Energy (keV);Neutrons per alpha");
	myCanvas->Write("energysimple_and_deconvolution", TObject::kOverwrite);

	myCanvas->Close();

	gDirectory->cd("..");
	f.Close();
}
