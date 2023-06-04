#include <cmath>

Double_t peak_activity(const TString filepath, unsigned short channel, unsigned int low, unsigned int high, const TString name){
	TFile f(filepath, "READ");

	TH1D* spectrum;
	TH1D* time;
	switch(channel){
		case 6:
			spectrum = (TH1D*)f.Get("labr_1_spectrum");
			time = (TH1D*)f.Get("labr_1_time");
			break;
		case 7:
			spectrum = (TH1D*)f.Get("labr_2_spectrum");
			time = (TH1D*)f.Get("labr_2_time");
			break;
	}

	TF1* fitf = new TF1("gaussplus","[0]+[1]*x[0]+gaus(2)", low, high);
//	fitf->SetNpx(100);
//	fitf->SetNumberFitPoints(100);
//	fitf->SetParLimits(0, 0, 1E6);
//	fitf->SetParLimits(1, -1E6, 1E6);
//	fitf->SetParLimits(2, -1E6, 1E6);
//	fitf->SetParLimits(3, -1E6, 1E6);
//	fitf->SetParLimits(4, 0, 1E6);
	fitf->SetParameters(0, 0, (high-low)/2, (high+low)/2, 1);
	fitf->SetParNames("Constant BG", "BG linear", "gauss a", "gauss b", "gauss c");
	// gauss(x) = a * exp(-0.5*((x-b)/c)**2)
	// integral(gauss(x))(x0,x1) = sqrt(2*pi)*a*|c|

	TCanvas* myCanvas = new TCanvas("peak activity fit to pol1+gaus");
	TFitResultPtr fitfresult = spectrum->Fit("gaussplus", "RSQ");
	std::stringstream resultnamestream;
	resultnamestream << "output/" << name << "_fit.root";
	TString resultname = resultnamestream.str();
	myCanvas->SaveAs(resultname);
	myCanvas->Close();

	Double_t result = sqrt(2*M_PI) * fitfresult->Parameter(2) * abs(fitfresult->Parameter(4));
//	cout << "Result from integral: " << spectrum->Integral(low, high) << endl;
//	cout << "Result from gaussian: " << result << endl;
	result/= time->FindLastBinAbove(0, 1, -1)*time->GetBinWidth(1)/1E12;

	f.Close();
	return result;
}
