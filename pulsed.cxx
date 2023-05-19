#define MIN_TOF 250
#define MAX_TOF 650
#define GAMMA_FLASH_BINS_N 350
#define NEUTRON_RESPONSE_BINS_N 500
#define PULSE_FIT_PARAMS_N 50

class pulse_fit_functor{
	private:
		TH1F* gamma_histo;
		Double_t gamma_stepsize;
		TH1F* neutron_histo;
		Double_t neutron_stepsize;
		Double_t param_width;
	public:
		Double_t gamma_min;
		Double_t gamma_max;
		Double_t neutron_min;
		Double_t neutron_max;
		pulse_fit_functor(TH1F* gamma_histo, TH1F* neutron_histo):gamma_histo(gamma_histo),neutron_histo(neutron_histo){
			gamma_stepsize = gamma_histo->GetBinWidth(1);
			neutron_stepsize = neutron_histo->GetBinWidth(1);
			param_width = (neutron_max-neutron_min)/PULSE_FIT_PARAMS_N;
		}
		Double_t operator()(Double_t* x, Double_t* p){
			Double_t sum = 0;
			Short_t param_index = 0;
			for(Short_t i=0; i<GAMMA_FLASH_BINS_N; i++){
				param_index = 1 + (x[0] - gamma_stepsize*i - neutron_min)/param_width;
				if(param_index>0 && param_index<PULSE_FIT_PARAMS_N){
					sum += (gamma_histo->GetBinContent(i)-p[PULSE_FIT_PARAMS_N+1]) * p[param_index];
				}
			}
			return p[0] + sum;
		}
};

void pulsed_per_file(char filepath[500], Double_t gammaflash_min, Double_t gammaflash_max, Double_t neutronresponse_min, Double_t neutronresponse_max, Double_t max_param, Double_t gamma_background){
	EnableImplicitMT();
	RDataFrame d("Data", filepath);

	auto monster = d.Filter("Channel==4");
	auto tof_plot = monster.Histo1D({"tof_plot", ";ToF;Counts", 1000, MIN_TOF, MAX_TOF}, "tof");
	auto psd_plot = monster.Histo1D({"psd_plot", ";psd;Counts", 1000, 0, 1}, "psd");
	auto tof_id_plot = monster.Histo2D({"tof_id_plot", ";ToF;PSD;Counts", 100, MIN_TOF, MAX_TOF, 100, 0, 1}, "tof", "psd");
	auto energy_id_plot = monster.Histo2D({"energy_id_plot", ";Energy;PSD;Counts", 4096/4, 0, 4096, 100, 0, 1}, "Energy", "psd");
	Double_t y1 = 0.3;
	Double_t y2 = 0.42;
	Double_t x1 = 2100;
	Double_t x2 = 2600;
	Double_t m = (y2-y1)/(x2-x1);
	Double_t n = y1-x1*m;
	auto monster_gammas = monster.Filter([m, n, y1](Double_t psd, UShort_t Energy){return psd<y1 || psd<m*Energy+n;},{"psd", "Energy"});
	auto gamma_tof_plot = monster_gammas.Histo1D({"gamma_tof_plot", ";ToF;Counts", 1000, MIN_TOF, MAX_TOF}, "tof");
	auto monster_neutrons = monster.Filter([m, n, y1](Double_t psd, UShort_t Energy){return psd>=y1 || psd>=m*Energy+n;},{"psd", "Energy"});
	auto neutron_tof_plot = monster_neutrons.Histo1D({"neutron_tof_plot", ";ToF;Counts", 1000, MIN_TOF, MAX_TOF}, "tof");

	auto gamma_flash = monster_gammas.Histo1D({"gamma_flash", ";ToF;Counts", GAMMA_FLASH_BINS_N, gammaflash_min, gammaflash_max}, "tof");
	auto neutron_response = monster_neutrons.Histo1D({"neutron_response", ";ToF;Counts", NEUTRON_RESPONSE_BINS_N, neutronresponse_min, neutronresponse_max}, "tof");

	TCanvas* myCanvas = new TCanvas("");

	tof_plot->Write("", TObject::kOverwrite);
	psd_plot->Write("", TObject::kOverwrite);
	tof_id_plot->Draw("COLZ");
	TLine* psd_line_1 = new TLine(0, y1, x1, y1);
	TLine* psd_line_2 = new TLine(x1, y1, (1-n)/m, 1);
	psd_line_1->SetLineStyle(3);
	psd_line_2->SetLineStyle(3);
	psd_line_1->Draw("same");
	psd_line_2->Draw("same");
	myCanvas->Write("tof_id_plot", TObject::kOverwrite);
	energy_id_plot->Draw("COLZ");
	myCanvas->Write("energy_id_plot", TObject::kOverwrite);

	gamma_tof_plot->Write("", TObject::kOverwrite);
	neutron_tof_plot->Write("", TObject::kOverwrite);

	gamma_flash->Write("", TObject::kOverwrite);
	neutron_response->Write("", TObject::kOverwrite);

	//fit
	pulse_fit_functor pulse_fit_obj = pulse_fit_functor((TH1F*)gDirectory->Get("gamma_flash"), (TH1F*)gDirectory->Get("neutron_response"));
	pulse_fit_obj.gamma_min = gammaflash_min;
	pulse_fit_obj.gamma_max = gammaflash_max;
	pulse_fit_obj.neutron_min = neutronresponse_min;
	pulse_fit_obj.neutron_max = neutronresponse_max;
	TF1* pulse_fit = new TF1("pulse_fit", pulse_fit_obj, neutronresponse_min, neutronresponse_max, PULSE_FIT_PARAMS_N+2);
	pulse_fit->SetNpx(NEUTRON_RESPONSE_BINS_N);
	pulse_fit->SetNumberFitPoints(NEUTRON_RESPONSE_BINS_N);
	pulse_fit->SetParLimits(0, 0, 50);
	for(UShort_t i=1; i<PULSE_FIT_PARAMS_N+1; i++){
		pulse_fit->SetParLimits(i, 0, max_param);
	}
	pulse_fit->SetParLimits(PULSE_FIT_PARAMS_N+1, 0, 500);
	pulse_fit->SetParameter(PULSE_FIT_PARAMS_N+1, gamma_background);
	pulse_fit->FixParameter(PULSE_FIT_PARAMS_N+1, gamma_background);
	TFitResultPtr fitresult = neutron_response->Fit("pulse_fit", "SLE");
	myCanvas->Write("pulse_fit_plot", TObject::kOverwrite);

	//guardar resultados en tree
	Double_t results[PULSE_FIT_PARAMS_N][2];
	for(UShort_t i=0; i<PULSE_FIT_PARAMS_N; i++){
		results[i][0] = fitresult->Parameter(i+1);
		results[i][1] = fitresult->ParError(i+1);
	}
	TTree* results_tree = new TTree("results_tree", "Tree with pulsed results");
	results_tree->Branch("results", results, "results[50][2]/D");	//PULSE_FIT_PARAMS_N
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
	pulsed_per_file(pulsed_1, 482, 505, 510, 600, 1E-5, 20);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_2");
	pulsed_per_file(pulsed_2, 479, 500, 510, 600, 1E-3, 20);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_3");
	pulsed_per_file(pulsed_3, 365, 390, 390, 500, 1E-2, 40);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_4");
	pulsed_per_file(pulsed_4, 288, 310, 310, 400, 1E-2, 10);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_5");
	pulsed_per_file(pulsed_5, 292, 310, 330, 500, 1E-1, 30);
	gDirectory->cd("..");

	gDirectory->cd("..");
	f.Close();
}

TGraph* pulsed_results_per_file(Double_t g_min, Double_t g_center, Double_t g_max, Double_t n_min, Double_t n_max, Double_t distance, Double_t alpha_energy){
	Double_t results[PULSE_FIT_PARAMS_N][2];
	TTree* tree = (TTree*)gDirectory->Get("results_tree");
	tree->SetBranchAddress("results", results);
	tree->GetEntry(0);

	Double_t tof_to_seconds = 1E-9;
	Double_t c = 299792458;	//speed of light
	Double_t neutron_mass = 1.67492749804E-27;
	Double_t J_to_keV = 1/1.602177E-16;

	Double_t x[PULSE_FIT_PARAMS_N];
	Double_t y[PULSE_FIT_PARAMS_N];
	Double_t y_err[PULSE_FIT_PARAMS_N];
	Double_t base = n_min-g_min + distance/c;
	Double_t paramwidth = (n_max-n_min)/PULSE_FIT_PARAMS_N;
	for(UShort_t i=0; i<PULSE_FIT_PARAMS_N; i++){
		x[i] = (base + paramwidth*i)*tof_to_seconds;
		y[i] = results[i][0];
		y_err[i] = results[i][1];
	}

	TCanvas* myCanvas = new TCanvas("");
	TGraph* cross_section_result = new TGraphErrors(PULSE_FIT_PARAMS_N, x, y, NULL, y_err);
	cross_section_result->SetTitle("Fit results;ToF (s);Counts");
	cross_section_result->SetMarkerStyle(21);
	cross_section_result->Draw("alp");
	myCanvas->Write("pulse_fit_results", TObject::kOverwrite);

	TH1D* gamma_flash = (TH1D*)gDirectory->Get("gamma_flash");
	Double_t gammas_n = gamma_flash->Integral(0, gamma_flash->GetNbinsX());
	for(UShort_t i=0; i<PULSE_FIT_PARAMS_N; i++){
		x[i] = (n_min + g_center-g_min + paramwidth*i);
		y[i] *= gammas_n;
		y_err[i] *= gammas_n;
	}
	TGraph* cross_section_result_delta = new TGraphErrors(PULSE_FIT_PARAMS_N, x, y, NULL, y_err);
	cross_section_result_delta->SetTitle("delta;ToF (ns);Counts");
	cross_section_result_delta->SetMarkerStyle(21);
	cross_section_result_delta->Draw("alp");
	((TH1D*)gDirectory->Get("neutron_response"))->Draw("same");
	myCanvas->Write("neutron_response+delta", TObject::kOverwrite);

	Double_t v;
	for(UShort_t i=0; i<PULSE_FIT_PARAMS_N; i++){
		x[i] = base + paramwidth*i;
		v = distance/((base + paramwidth*i)*tof_to_seconds);
		x[i] = (neutron_mass/2*v*v)*J_to_keV;
		y[i] = results[i][0];
		y_err[i] = results[i][1];
	}
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
	return energy_result;
}

void pulsed_results(){
	TFile f("output.root", "UPDATE");
	gDirectory->cd("Pulsed");

	gDirectory->cd("pulsed_1");
	TGraph* energy_1 = pulsed_results_per_file(482, 488.100, 505, 510, 600, 1, 5500);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_2");
	TGraph* energy_2 = pulsed_results_per_file(479, 485.950, 500, 510, 600, 1, 5500);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_3");
	TGraph* energy_3 = pulsed_results_per_file(365, 360.637, 390, 390, 500, 1, 7000);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_4");
	TGraph* energy_4 = pulsed_results_per_file(288, 294.508, 310, 310, 400, 1, 8250);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_5");
	TGraph* energy_5 = pulsed_results_per_file(292, 297.897, 310, 330, 500, 2, 8250);
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
