Double_t peak_activity(const Char_t filepath[500], unsigned short channel, unsigned int low, unsigned int high){
	TFile f(filepath, "READ");

	TH1D* histo;
	switch(channel){
		case 6:
			histo = (TH1D*)f.Get("labr_1_spectrum");
			break;
		case 7:
			histo = (TH1D*)f.Get("labr_2_spectrum");
			break;
	}

	Double_t result = histo->Integral(low, high);
	result/= histo->FindLastBinAbove(0, 1, -1)*histo->GetBinWidth(1)/1E12;

	f.Close();
	return result;
}
