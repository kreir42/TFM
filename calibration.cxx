void calibrate_energy(Char_t filepath[500], short number_of_peaks, int* left, int* center, int* right, Double_t* energy, int* channel){
}

void energy_calibration(){
	Char_t filepaths[6][500]={"output/SData_LaBr_Cs137atTarget_calib_20230223.root",
		"output/SData_LaBr_Na22atTarget_calib_20230223.root",
		"output/SData_StylbD2_Cf252calibration_20230222_11am.root",
		"output/SData_Tadeo1_Cs137_Tadeo2_Na22_stlbd2_Cf252_calib_20230223.root",
		"output/SData_Tadeo1_Na22_Tadeo2_Cf252_stlbd2_Cs137_calib_20230223.root",
		"output/SData_Tadeo1_nada_Tadeo2_Cs137_stlbd2_nada_calib_20230223.root"};
	Double_t data[]={560, 600, 640, 511, 7,		//TBD:Double_t->int implicit conversion!
	1450, 1500, 1550, 1274.53, 7,
//	1650, 1725, 1800, 0.0, 7,
	550, 585, 630, 511, 6,
	1430, 1490, 1550, 1274.53, 6,
//	1650, 1710, 1760, 0.0, 6,
	0, 0, 0, 0.0, 0};
	short numbers_of_peaks[]={2,2};
	Char_t* filepath;
	unsigned short n = 0;
	//for(short i=0; i<6; i++){
	//	filepath = filepaths[i];
	//	calibrate_energy(filepath, numbers_of_peaks[i], data[n], data[n+1], data[n+2], data[n+3], data[n+4]);
	//	n+=numbers_of_peaks[i];
}
