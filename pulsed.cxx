#define MIN_TOF 250
#define MAX_TOF 650
#define GAMMA_FLASH_BINS_N 500
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
					sum += gamma_histo->GetBinContent(i) * p[param_index];
				}
			}
			return p[0] + sum;
		}
};

void pulsed_per_file(char filepath[500], Double_t gammaflash_min, Double_t gammaflash_max, Double_t neutronresponse_min, Double_t neutronresponse_max, Double_t max_param){
	EnableImplicitMT();
	RDataFrame d("Data", filepath);

	auto monster = d.Filter("Channel==4");
	auto tof_plot = monster.Histo1D({"tof_plot", ";ToF;Counts", 1000, MIN_TOF, MAX_TOF}, "tof");
	auto psd_plot = monster.Histo1D({"psd_plot", ";psd;Counts", 1000, 0, 1}, "psd");
	auto tof_id_plot = monster.Histo2D({"tof_id_plot", ";ToF;PSD;Counts", 100, MIN_TOF, MAX_TOF, 100, 0, 1}, "tof", "psd");
	auto energy_id_plot = monster.Histo2D({"energy_id_plot", ";Energy;PSD;Counts", 4096/4, 0, 4096, 100, 0, 1}, "Energy", "psd");

	auto monster_gammas = monster.Filter("psd<0.3");
	auto gamma_tof_plot = monster_gammas.Histo1D({"gamma_tof_plot", ";ToF;Counts", 1000, MIN_TOF, MAX_TOF}, "tof");
	auto monster_neutrons = monster.Filter("psd>0.3");
	auto neutron_tof_plot = monster_neutrons.Histo1D({"neutron_tof_plot", ";ToF;Counts", 1000, MIN_TOF, MAX_TOF}, "tof");

	auto gamma_flash = monster_gammas.Histo1D({"gamma_flash", ";ToF;Counts", GAMMA_FLASH_BINS_N, gammaflash_min, gammaflash_max}, "tof");
	auto neutron_response = monster_neutrons.Histo1D({"neutron_response", ";ToF;Counts", NEUTRON_RESPONSE_BINS_N, neutronresponse_min, neutronresponse_max}, "tof");

	TCanvas* myCanvas = new TCanvas("");

	tof_plot->Write("", TObject::kOverwrite);
	psd_plot->Write("", TObject::kOverwrite);
	tof_id_plot->Draw("COLZ");
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
	TF1* pulse_fit = new TF1("pulse_fit", pulse_fit_obj, neutronresponse_min, neutronresponse_max, PULSE_FIT_PARAMS_N+1);
	pulse_fit->SetNpx(200);
	pulse_fit->SetNumberFitPoints(200);
	pulse_fit->SetParLimits(0, 0, 50);
	for(UShort_t i=1; i<PULSE_FIT_PARAMS_N+1; i++){
		pulse_fit->SetParLimits(i, 0, max_param);
	}
	TFitResultPtr fitresult = neutron_response->Fit("pulse_fit", "SLE");
	myCanvas->Write("pulse_fit_plot", TObject::kOverwrite);

	Float_t x[PULSE_FIT_PARAMS_N];
	Float_t y[PULSE_FIT_PARAMS_N];
	Float_t y_err[PULSE_FIT_PARAMS_N];
	for(UShort_t i=0; i<PULSE_FIT_PARAMS_N; i++){
		x[i] = i;	//TBD
		y[i] = fitresult->Parameter(i+1);
		y_err[i] = fitresult->ParError(i+1);
	}
	TGraph* cross_section_result = new TGraphErrors(PULSE_FIT_PARAMS_N, x, y, NULL, y_err);
	cross_section_result->SetTitle("Fit results;ToF (arbitrary);Cross section (arbitrary)");
	cross_section_result->SetMarkerStyle(21);
	cross_section_result->Draw("alp");
	myCanvas->Write("pulse_fit_results", TObject::kOverwrite);

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
	pulsed_per_file(pulsed_1, 482, 505, 510, 600, 1E-4);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_2");
	pulsed_per_file(pulsed_2, 479, 500, 510, 600, 1E-3);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_3");
	pulsed_per_file(pulsed_3, 365, 390, 390, 500, 1E-3);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_4");
	pulsed_per_file(pulsed_4, 288, 310, 310, 400, 1E-3);
	gDirectory->cd("..");

	gDirectory->cd("pulsed_5");
	pulsed_per_file(pulsed_5, 292, 310, 330, 500, 1E-1);
	gDirectory->cd("..");

	gDirectory->cd("..");
	f.Close();
}

void pulsed_results(){
}
