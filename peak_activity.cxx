Double_t peak_activity(const Char_t filepath[500], unsigned short channel, unsigned int low, unsigned int high){
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

	Double_t result = spectrum->Integral(low, high);
	result/= time->FindLastBinAbove(0, 1, -1)*time->GetBinWidth(1)/1E12;

	f.Close();
	return result;
}
