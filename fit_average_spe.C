Double_t fit_simple_spe(Double_t *x, Double_t *par){

	Double_t fit_val = 0;
	fit_val = par[0]*TMath::Gaus(x[0],par[1],par[2])+par[3]*exp(-(x[0]-par[4])/par[5]) + par[6]*(TMath::Gaus(x[0],par[7],par[8])+TMath::Gaus(x[0],par[9]*par[7],par[10]*par[8]))+par[11]*(TMath::Gaus(x[0],2*par[7],sqrt(2)*par[8])+TMath::Gaus(x[0],(1+par[9])*par[7],sqrt(1+par[10])*par[8]));
	return fit_val;
}


void fit_average_spe(){

	TFile *f = new TFile("LED_SPE_Results.root","READ");

	TProfile *hAvgLUX = (TProfile*) f->Get("hAvgLUX");
	TProfile *hAvgWB = (TProfile*) f->Get("hAvgWB");
	TProfile *hAvgHM = (TProfile*) f->Get("hAvgHM");
	TProfile *hAvgETEL = (TProfile*) f->Get("hAvgETEL");


	TF1 *fit_spe = new TF1("fit_spe",fit_simple_spe,0,6,12);
	Double_t start_pars[13]={1.,0.05,0.04,0.02,-7000000,1000000,0.05,1.0,0.5,0.2,0.3,0.0001};
	std::string par_names[13]={"A_{ped}","#mu_{ped}","#sigma_{ped}","B","a","#tau","A_{1PE}","#mu","#sigma","f_{#mu}","f_{#sigma}","A_{2PE}"};
	fit_spe->SetParName(0,"A_{ped}");
	fit_spe->SetParName(1,"#mu_{ped}");
	fit_spe->SetParName(2,"#sigma_{ped}");
	fit_spe->SetParName(3,"B");
	fit_spe->SetParName(4,"a");
	fit_spe->SetParName(5,"#tau");
	fit_spe->SetParName(6,"A_{1PE}");
	fit_spe->SetParName(7,"#mu");
	fit_spe->SetParName(8,"#sigma");
	fit_spe->SetParName(9,"f_{#mu}");
	fit_spe->SetParName(10,"f_{#sigma}");
	fit_spe->SetParName(11,"A_{2PE}");
	//fit_spe->SetParNames(par_names);
	
	fit_spe->SetParLimits(0,0.5,100);
	fit_spe->SetParLimits(1,0.01,0.2);
	fit_spe->SetParLimits(2,0.001,0.05);
	fit_spe->SetParLimits(6,0.001,0.1);
	fit_spe->SetParLimits(7,0.5,2);
	fit_spe->SetParLimits(8,0.2,1);
	fit_spe->SetParLimits(9,0.01,1.);
	fit_spe->SetParLimits(10,0.01,1.);
	fit_spe->SetParLimits(11,0.0001,0.01);

	//hAvgLUX->Fit("fit_spe","R");
	//hAvgWB->Fit("fit_spe","R");
	//hAvgHM->Fit("fit_spe","R");
	hAvgETEL->Fit("fit_spe","R");


}
