Double_t Scale(Double_t x, Double_t factor)
{
  Double_t v;
  v = factor * x ; // "linear scaling" function example
  return v;
}


void ScaleAxis(TAxis *a, Double_t Factor)
{
  if (!a) return; // just a precaution
  if (a->GetXbins()->GetSize())
    {
      // an axis with variable bins
      // note: bins must remain in increasing order, hence the "Scale"
      // function must be strictly (monotonically) increasing
      TArrayD X(*(a->GetXbins()));
      for(Int_t i = 0; i < X.GetSize(); i++) X[i] = Scale(X[i],Factor);
      a->Set((X.GetSize() - 1), X.GetArray()); // new Xbins
    }
  else
    {
      // an axis with fix bins
      // note: we modify Xmin and Xmax only, hence the "Scale" function
      // must be linear (and Xmax must remain greater than Xmin)
      a->Set( a->GetNbins(),
              Scale(a->GetXmin(),Factor), // new Xmin
              Scale(a->GetXmax(),Factor )); // new Xmax
    }
  return;
}

void ScaleXaxis(TH1 *h, Double_t Factor)
{
  if (!h) return; // just a precaution
  ScaleAxis(h->GetXaxis(), Factor);
  return;
}

void BinLogY(TH2*h) 
{

   TAxis *axis = h->GetYaxis(); 
   int bins = axis->GetNbins();

   Axis_t from = axis->GetXmin();
   Axis_t to = axis->GetXmax();
   Axis_t width = (to - from) / bins;
   Axis_t *new_bins = new Axis_t[bins + 1];

   for (int i = 0; i <= bins; i++) {
     new_bins[i] = TMath::Power(10, from + i * width);
   } 
   axis->Set(bins, new_bins); 
   delete[] new_bins; 
}


TH2F *Average(std::vector<TH1F*> hists,std::string suffix) {
   std::stringstream ss_name, ss_title;
   ss_name << "hAll_"<<suffix;
   ss_title << suffix <<" SPE distributions";
   TH2F *hAvg = new TH2F(ss_name.str().c_str(),ss_title.str().c_str(),200,0,6.,300,-6,0.5);
   BinLogY(hAvg);
   hAvg->GetXaxis()->SetTitle("charge [p.e.]");
   hAvg->GetYaxis()->SetTitle("Probability");
   for (size_t i=0;i < hists.size(); ++i) {
      for (int bin = 0; bin < hists[i]->GetNbinsX(); ++bin) {
         double xbin = hists[i]->GetBinCenter(bin+1);
         int resultBin = hAvg->FindBin(xbin, hists[i]->GetBinContent(bin+1));
         int binContent = hAvg->GetBinContent(resultBin) + 1;
         hAvg->SetBinContent(resultBin,binContent);
      }
   }
   return hAvg;
}


void plot_spe_proposed(){

	TFile *f = new TFile("ProposedVoltagesV3R1301S0_0V_PMTStability_Run0.root","READ");
	//TTree *t = (TTree*) f->Get("fit_parameters");

	int detkey;
	int voltage;
	int found_spe_peak;
	double fit_chi2;
	double ped_mean;
	double ped_sigma;
	double spe_firstgamma_scaling;
	double spe_firstgamma_shape;
	double spe_firstgamma_mean;
	double spe_secondgamma_scaling;
	double spe_secondgamma_mean_scaling;
	double spe_secondgamma_shape_scaling;
	double gain;
	double spe_firstgamma_gain;
/*
	t->SetBranchAddress("DetectorKey",&detkey);
	t->SetBranchAddress("Voltage",&voltage);
	t->SetBranchAddress("found_spe_peak",&found_spe_peak);
	t->SetBranchAddress("fit_chi2",&fit_chi2);
	t->SetBranchAddress("ped_mean",&ped_mean);
	t->SetBranchAddress("ped_sigma",&ped_sigma);
	t->SetBranchAddress("spe_firstgamma_scaling",&spe_firstgamma_scaling);
	t->SetBranchAddress("spe_firstgamma_mean",&spe_firstgamma_mean);
	t->SetBranchAddress("spe_firstgamma_shape",&spe_firstgamma_shape);
	t->SetBranchAddress("spe_secondgamma_scaling",&spe_secondgamma_scaling);
	t->SetBranchAddress("spe_secondgamma_mean_scaling",&spe_secondgamma_mean_scaling);
	t->SetBranchAddress("spe_secondgamma_shape_scaling",&spe_secondgamma_shape_scaling);
	t->SetBranchAddress("gain",&gain);
	t->SetBranchAddress("spe_firstgamma_gain",&spe_firstgamma_gain);
*/
	//int entries = t->GetEntries();

	TFile *fout = new TFile("LED_SPE_Results_Proposed.root","RECREATE");
	TTree *tout = new TTree("led_gains","LED gains");
	TDirectoryFile *dir_gain = new TDirectoryFile("dir_gain","Gain directory");
	TDirectoryFile *dir_spe = new TDirectoryFile("dir_spe","S.P.E. directory");
	TDirectoryFile *dir_shape = new TDirectoryFile("dir_shape","Shape directory");
	TDirectoryFile *dir_summary = new TDirectoryFile("dir_summary","Summary directory");
	TDirectoryFile *dir_calibration = new TDirectoryFile("dir_calibration","Calibration directory");
	TDirectoryFile *dir_calibration_spe = new TDirectoryFile("dir_calibration_spe","SPE Calibration directory");

	ifstream file_hv("ANNIE_HV_Chankey.txt");
	int detKey;
	int hv;
	std::string type;
	std::map<int,int> map_chkey_hv;
	std::map<int,std::string> map_chkey_type;
	std::map<int,TH1F*> h_gains_pmts;
	std::map<int,TH1F*> h_spe_pmts;
	std::map<int,TH1F*> h_shape_pmts;
	while (!file_hv.eof()){
		file_hv >> detKey >> type >> hv;
		std::cout <<"Detkey: "<<detKey<<", Type: "<<type<<", HV: "<<hv<<std::endl;
		map_chkey_hv.emplace(detKey,hv);
		map_chkey_type.emplace(detKey,type);
		std::stringstream ss_hist_gain;
		std::stringstream ss_hist_spe;
		std::stringstream ss_hist_shape;
		ss_hist_gain << "h_gain_"<<detKey;
		ss_hist_spe << "h_spe_"<<detKey;
		ss_hist_shape << "h_shape_"<<detKey;
		TH1F *h_gain = new TH1F(ss_hist_gain.str().c_str(),ss_hist_gain.str().c_str(),100,0,1.0*pow(10,7));
		TH1F *h_spe = new TH1F(ss_hist_spe.str().c_str(),ss_hist_spe.str().c_str(),100,0,10);
		TH1F *h_shape = new TH1F(ss_hist_shape.str().c_str(),ss_hist_shape.str().c_str(),100,0,10);
		h_gains_pmts.emplace(detKey,h_gain);
		h_spe_pmts.emplace(detKey,h_spe);
		h_shape_pmts.emplace(detKey,h_shape);
		if (file_hv.eof()) break;
	}
	file_hv.close();

	ifstream file_gain("ANNIE_Gains.txt");
	int tmp_detkey;
	double tmp_gain;
	std::map<int,double> map_chkey_gain;
	while (!file_gain.eof()){
		file_gain >> tmp_detkey >> tmp_gain;
		std::cout <<"Detkey: "<<tmp_detkey<<", Gain: "<<tmp_gain<<std::endl;
		map_chkey_gain.emplace(tmp_detkey,tmp_gain);
		if (file_gain.eof()) break;
	}
	file_gain.close();

	ifstream file_bad("bad_spe.txt");
	int tmp_det;
	std::vector<int> bad_spe;
	while (!file_bad.eof()){
		file_bad >> tmp_det;
		std::cout <<"Bad detkey: "<<tmp_det<<std::endl;
		bad_spe.push_back(tmp_det);
		if (file_bad.eof()) break;
	}
	file_bad.close();


	std::string pmt_type;
	
	TH1F *h_gains_etel = new TH1F("h_gains_etel","ETEL gain distribution",100,0,1.0*pow(10,7));
	TH1F *h_gains_lux = new TH1F("h_gains_lux","LUX gain distribution",100,0,1.0*pow(10,7));
	TH1F *h_gains_hm = new TH1F("h_gains_hm","Hamamatsu gain distribution",100,0,1.0*pow(10,7));
	TH1F *h_gains_wb = new TH1F("h_gains_wb","Watchboy gain distribution",100,0,1.0*pow(10,7));

	TH1F *h_spe_etel = new TH1F("h_spe_etel","ETEL SPE pC mean distribution",100,0,10);
	TH1F *h_spe_lux = new TH1F("h_spe_lux","ETEL SPE pC mean distribution",100,0,10);
	TH1F *h_spe_hm = new TH1F("h_spe_hm","ETEL SPE pC mean distribution",100,0,10);
	TH1F *h_spe_wb = new TH1F("h_spe_wb","ETEL SPE pC mean distribution",100,0,10);
	
	TH1F *h_shape_etel = new TH1F("h_shape_etel","ETEL SPE shape distribution",100,0,10);
	TH1F *h_shape_lux = new TH1F("h_shape_lux","ETEL SPE shape distribution",100,0,10);
	TH1F *h_shape_hm = new TH1F("h_shape_hm","ETEL SPE shape distribution",100,0,10);
	TH1F *h_shape_wb = new TH1F("h_shape_wb","ETEL SPE shape distribution",100,0,10);
	
	TH1F *h_spe_conv_etel = new TH1F("h_spe_conv_etel","ETEL SPE p.e. mean distribution",100,0,10);
	TH1F *h_spe_conv_lux = new TH1F("h_spe_conv_lux","ETEL SPE p.e. mean distribution",100,0,10);
	TH1F *h_spe_conv_hm = new TH1F("h_spe_conv_hm","ETEL SPE p.e. mean distribution",100,0,10);
	TH1F *h_spe_conv_wb = new TH1F("h_spe_conv_wb","ETEL SPE p.e. mean distribution",100,0,10);

	TH2F *h_spe_gain = new TH2F("h_spe_gain","SPE vs. gain distribution",50,0,10,50,0,1.0*pow(10,7));

	std::vector<TH1F*> good_spe_LUX;
	std::vector<TH1F*> good_spe_ETEL;
	std::vector<TH1F*> good_spe_HM;
	std::vector<TH1F*> good_spe_WB;
/*
	for (int i_entry = 0; i_entry < entries; i_entry++){

		t->GetEntry(i_entry);

		if (map_chkey_hv.count(detkey)>0){
			if (voltage == map_chkey_hv.at(detkey)){
				h_gains_pmts.at(detkey)->Fill(gain);
				h_spe_pmts.at(detkey)->Fill(spe_firstgamma_mean);
				h_spe_gain->Fill(spe_firstgamma_mean,gain);
				h_shape_pmts.at(detkey)->Fill(spe_firstgamma_shape);
			}
		}	
	}
*/
	dir_gain->cd();
	for (std::map<int,TH1F*>::iterator it = h_gains_pmts.begin(); it!= h_gains_pmts.end(); it++){
		int detKey = it->first;
		TH1F *htemp = it->second;
		double mean = htemp->GetMean();
		if (mean>0){
			if (map_chkey_type.at(detKey)=="ETEL") h_gains_etel->Fill(mean);
			else if (map_chkey_type.at(detKey)=="LUX") h_gains_lux->Fill(mean);
			else if (map_chkey_type.at(detKey)=="HM") h_gains_hm->Fill(mean);
			else if (map_chkey_type.at(detKey)=="WB") h_gains_wb->Fill(mean);
		}
		htemp->Write();
	}
        
	dir_shape->cd();
	for (std::map<int,TH1F*>::iterator it = h_shape_pmts.begin(); it!= h_shape_pmts.end(); it++){
		int detKey = it->first;
		TH1F *htemp = it->second;
		double mean = htemp->GetMean();
		if (mean>0){
			if (map_chkey_type.at(detKey)=="ETEL") h_shape_etel->Fill(mean);
			else if (map_chkey_type.at(detKey)=="LUX") h_shape_lux->Fill(mean);
			else if (map_chkey_type.at(detKey)=="HM") h_shape_hm->Fill(mean);
			else if (map_chkey_type.at(detKey)=="WB") h_shape_wb->Fill(mean);
		}
		htemp->Write();
	}

	dir_spe->cd();
	for (std::map<int,TH1F*>::iterator it = h_spe_pmts.begin(); it!= h_spe_pmts.end(); it++){
		int detKey = it->first;
		TH1F *htemp = it->second;
		double mean = htemp->GetMean();
		double gain = h_gains_pmts.at(detKey)->GetMean();
		if (mean>0){
			if (map_chkey_type.at(detKey)=="ETEL") {
				h_spe_etel->Fill(mean);
				h_spe_conv_etel->Fill(mean/(gain*1.6*pow(10,-7)));
			}
			else if (map_chkey_type.at(detKey)=="LUX") {
				h_spe_lux->Fill(mean);
				h_spe_conv_lux->Fill(mean/(gain*1.6*pow(10,-7)));
			}
			else if (map_chkey_type.at(detKey)=="HM"){
				h_spe_hm->Fill(mean);
				h_spe_conv_hm->Fill(mean/(gain*1.6*pow(10,-7)));
			}
			else if (map_chkey_type.at(detKey)=="WB") {
				h_spe_wb->Fill(mean);
				h_spe_conv_wb->Fill(mean/(gain*1.6*pow(10,-7)));
			}
		}
		htemp->Write();
	}

	dir_summary->cd();
	h_gains_etel->Write();
	h_gains_lux->Write();
	h_gains_hm->Write();
	h_gains_wb->Write();
	h_shape_etel->Write();
	h_shape_lux->Write();
	h_shape_hm->Write();
	h_shape_wb->Write();
	h_spe_etel->Write();
	h_spe_conv_etel->Write();
	h_spe_lux->Write();
	h_spe_conv_lux->Write();
	h_spe_hm->Write();
	h_spe_conv_hm->Write();
	h_spe_wb->Write();
	h_spe_conv_wb->Write();

	dir_calibration->cd();
	TCanvas *c = new TCanvas("c","Canvas",900,600);
	TCanvas *c_spe = new TCanvas("c_spe","Canvas SPE",900,600);
	c->Divide(4,2);
	c_spe->Divide(4,2);
	bool first=true;
	bool first_spe=true;
	int counter=0;
	int counter_spe=0;
	for (std::map<int,int>::iterator it = map_chkey_hv.begin(); it!= map_chkey_hv.end(); it++){
		int chkey = it->first;
		int hv = it->second;
		std::stringstream ss_hist;
		//ss_hist << "hist_charge_"<<chkey<<"_"<<hv<<"V";
		ss_hist << "hist_charge_"<<chkey;
		if (f->GetListOfKeys()->Contains(ss_hist.str().c_str())){
			TH1F *htemp = (TH1F*) f->Get(ss_hist.str().c_str());
			std::stringstream ss_title;
			ss_title << "Hit charges for detkey "<<chkey<<" ("<<hv<<" V)";
			htemp->SetTitle(ss_title.str().c_str());
			std::stringstream ss_name;
			ss_name << "fit_ch"<<chkey<<"_hv"<<hv;
			htemp->SetName(ss_name.str().c_str());
			c->cd();
			if (counter == 8) counter=0;
			if (counter==0) {c->Clear();c->Divide(4,2);}
			c->cd(counter+1);
			if (htemp->GetXaxis()->GetXmax()<0.1) ScaleXaxis(htemp,1000);
			double max_value = htemp->GetMaximum();
			if (max_value > 2.) htemp->Scale(1./max_value);
			htemp->GetXaxis()->SetTitle("charge [pC]");
			htemp->Draw();
			gPad->SetLogy();
			double g = map_chkey_gain.at(chkey);
			double conv = 0.001/g;
			double min = htemp->GetMinimum();
			TLine *lgain = new TLine(1./conv,min,1./conv,1);
			lgain->SetLineWidth(1);
			lgain->SetLineColor(kBlack);
			lgain->SetLineStyle(2);
			lgain->Draw("same");
			if (counter==7 || chkey == 463){
			if (first) {c->Print("LED_SPE_Fits_Proposed.pdf(","pdf"); first = false;}
			else if (chkey == 463) c->Print("LED_SPE_Fits_Proposed.pdf)","pdf");
			else c->Print("LED_SPE_Fits_Proposed.pdf","pdf");
			}
			counter++;
			htemp->Write();

			//Create SPE version of histogram
			TH1F *htemp_spe = (TH1F*) htemp->Clone();
			ss_title.str("");
			ss_name.str("");
			ss_title << "P.E. distribution for detkey "<<chkey<<" ("<<hv<<" V)";
			ss_name << "fit_ch"<<chkey<<"_hv"<<hv<<"_pe";
			htemp_spe->SetName(ss_name.str().c_str());
			htemp_spe->SetTitle(ss_title.str().c_str());
			//double g = map_chkey_gain.at(chkey);
			//double conv = 0.001/g;
			c_spe->cd();
 			if (counter_spe == 8) counter_spe=0;
                        if (counter_spe == 0) {c_spe->Clear();c_spe->Divide(4,2);}
                        c_spe->cd(counter_spe+1);
			htemp_spe->GetListOfFunctions()->Remove(htemp_spe->GetFunction("full_fit_func"));
			htemp_spe->GetXaxis()->SetTitle("charge [p.e.]");
                        htemp_spe->Draw();
			double max_value_spe = htemp_spe->GetMaximum();
			if (max_value_spe > 2.) htemp_spe->Scale(1./max_value_spe);
			ScaleXaxis(htemp_spe,conv);
			htemp_spe->ResetStats();
			htemp_spe->Draw();
                        gPad->SetLogy();
                        if (counter_spe==7 || chkey == 463){
                        if (first_spe) {c_spe->Print("LED_SPE_Fits_normed_Proposed.pdf(","pdf"); first_spe = false;}
                        else if (chkey == 463) c_spe->Print("LED_SPE_Fits_normed_Proposed.pdf)","pdf");
                        else c_spe->Print("LED_SPE_Fits_normed_Proposed.pdf","pdf");
			}
			counter_spe++;
			dir_calibration_spe->cd();
			htemp_spe->Write();
			dir_calibration->cd();
			if (map_chkey_type.at(chkey)=="LUX" && std::find(bad_spe.begin(),bad_spe.end(),chkey)==bad_spe.end()){
				good_spe_LUX.push_back(htemp_spe);
			}
			if (map_chkey_type.at(chkey)=="ETEL" && std::find(bad_spe.begin(),bad_spe.end(),chkey)==bad_spe.end()){
				good_spe_ETEL.push_back(htemp_spe);
			}
			if (map_chkey_type.at(chkey)=="HM" && std::find(bad_spe.begin(),bad_spe.end(),chkey)==bad_spe.end()){
				good_spe_HM.push_back(htemp_spe);
			}
			if (map_chkey_type.at(chkey)=="WB" && std::find(bad_spe.begin(),bad_spe.end(),chkey)==bad_spe.end()){
				good_spe_WB.push_back(htemp_spe);
			}
		} else {
			std::stringstream ss_temp;
			ss_temp << "deapfit_inputs/Input_"<<hv<<"V.root";
			TFile *ftemp = new TFile(ss_temp.str().c_str(),"READ");
			std::stringstream ss_hist_temp;
			ss_hist_temp << "hist_charge_"<<chkey;
			TH1F *hcharge = (TH1F*) ftemp->Get(ss_hist_temp.str().c_str());
			std::stringstream ss_new;
			ss_new << "input_ch"<<chkey<<"_hv"<<hv;
			hcharge->SetName(ss_new.str().c_str());
			std::stringstream ss_title;
			ss_title << "Hit charges for detkey "<<chkey<<" ("<<hv<<" V), NOT FIT";
			hcharge->SetTitle(ss_title.str().c_str());
			hcharge->GetXaxis()->SetRangeUser(0,0.006);
			ScaleXaxis(hcharge,1000);
			double max_value_charge = hcharge->GetMaximum();
			if (max_value_charge > 2.) hcharge->Scale(1./max_value_charge);
			hcharge->GetXaxis()->SetTitle("charge [pC]");
			c->cd();
			if (counter==8) counter=0;
			if (counter==0) {c->Clear(); c->Divide(4,2);}
			c->cd(counter+1);
			hcharge->GetListOfFunctions()->Remove(hcharge->GetFunction("total"));
			hcharge->Draw();
			gPad->SetLogy();
			double g = map_chkey_gain.at(chkey);
			double conv = 0.001/g;
			double min = hcharge->GetMinimum();
			TLine *lgain = new TLine(1./conv,min,1./conv,1);
			lgain->SetLineWidth(1);
			lgain->SetLineColor(kBlack);
			lgain->SetLineStyle(2);
			lgain->Draw("same");
			if (counter==7 || chkey == 463){
			if (first) {c->Print("LED_SPE_Fits.pdf(","pdf"); first = false;}
			else if (chkey == 463) c->Print("LED_SPE_Fits.pdf)","pdf");
			else c->Print("LED_SPE_Fits.pdf","pdf");
			}
			counter++;
			fout->cd();
			dir_calibration->cd();
			hcharge->Write();
			
                        //Create SPE version of histogram
                        TH1F *hcharge_spe = (TH1F*) hcharge->Clone();
                        ss_title.str("");
                        ss_new.str("");
                        ss_title << "P.E. distribution for detkey "<<chkey<<" ("<<hv<<" V) NOT FIT";
                        ss_new << "input_ch"<<chkey<<"_hv"<<hv<<"_pe";
                        hcharge_spe->SetName(ss_new.str().c_str());
                        hcharge_spe->SetTitle(ss_title.str().c_str());
                        //double g = map_chkey_gain.at(chkey);
                        //double conv = 0.001/g;
                        c_spe->cd();
                        if (counter_spe == 8) counter_spe=0;
                        if (counter_spe == 0) {c_spe->Clear();c_spe->Divide(4,2);}
                        c_spe->cd(counter_spe+1);
                        hcharge_spe->GetListOfFunctions()->Remove(hcharge_spe->GetFunction("total"));
                        hcharge_spe->GetXaxis()->SetTitle("charge [p.e.]");
                        hcharge_spe->Draw();
                        ScaleXaxis(hcharge_spe,conv);
                        hcharge_spe->ResetStats();
                        hcharge_spe->Draw();
                        double max_value_spe = hcharge_spe->GetMaximum();
                        if (max_value_spe > 2.) hcharge_spe->Scale(1./max_value_spe);
                        gPad->SetLogy();
                        if (counter_spe==7 || chkey == 463){
                        if (first_spe) {c_spe->Print("LED_SPE_Fits_normed.pdf(","pdf"); first_spe = false;}
                        else if (chkey == 463) c_spe->Print("LED_SPE_Fits_normed.pdf)","pdf");
                        else c_spe->Print("LED_SPE_Fits_normed.pdf","pdf");
                        }
                        counter_spe++;
                        dir_calibration_spe->cd();
                        hcharge_spe->Write();
                        dir_calibration->cd();

			if (map_chkey_type.at(chkey)=="LUX" && std::find(bad_spe.begin(),bad_spe.end(),chkey)==bad_spe.end()){
				good_spe_LUX.push_back(hcharge_spe);
			}
			if (map_chkey_type.at(chkey)=="ETEL" && std::find(bad_spe.begin(),bad_spe.end(),chkey)==bad_spe.end()){
				good_spe_ETEL.push_back(hcharge_spe);
			}
			if (map_chkey_type.at(chkey)=="HM" && std::find(bad_spe.begin(),bad_spe.end(),chkey)==bad_spe.end()){
				good_spe_HM.push_back(hcharge_spe);
			}
			if (map_chkey_type.at(chkey)=="WB" && std::find(bad_spe.begin(),bad_spe.end(),chkey)==bad_spe.end()){
				good_spe_WB.push_back(hcharge_spe);
			}
			std::cout <<"Did not find LED charge histogram for chankey "<<chkey<<" at "<<hv << " V! Omit"<<std::endl;
		}
	}

   gStyle->SetPalette(52); //Use Grayscale palette

   TCanvas *c_avg_all = new TCanvas("c_avg_all","Average SPE distributions");
   c_avg_all->Divide(2,2);
   c_avg_all->SetLogy();

   //Show average LUX distribution
   TCanvas *c_avg_lux = new TCanvas("c_avg_lux","Average LUX SPE Distribution");
   //Compute average
   std::string lux = "LUX";
   TH2F *hAvg = Average(good_spe_LUX,lux);
   TProfile *hProf = hAvg->ProfileX("hAvgLUX");

   //Turn off stats
   hAvg->SetStats(kFALSE);

   //Draw the 2D hist with the average prof on top.
   hProf->SetLineColor(kRed);
   hProf->SetMarkerColor(kRed);
   hAvg->Draw("COLZ");
   hProf->Draw("SAME");

   //Use a log scale in Z
   gPad->SetLogz();
   gPad->SetLogy();
   gPad->Update();

   fout->cd();
   c_avg_lux->Write("c_avg_lux");
   c_avg_lux->Print("LUX_Average_SPE.pdf","pdf");
  
   c_avg_all->cd(1);
   hAvg->SetStats(0);
   hAvg->Draw("COLZ");
   hProf->Draw("SAME");
   gPad->SetLogz();
   gPad->SetLogy();

   //Show average ETEL distribution
   TCanvas *c_avg_etel = new TCanvas("c_avg_etel","Average ETEL SPE Distribution");
   //Compute average
   std::string etel = "ETEL";
   TH2F *hAvgETEL = Average(good_spe_ETEL,etel);
   TProfile *hProfETEL = hAvgETEL->ProfileX("hAvgETEL");

   //Turn off stats
   hAvgETEL->SetStats(kFALSE);

   //Draw the 2D hist with the average prof on top.
   hProfETEL->SetLineColor(kRed);
   hProfETEL->SetMarkerColor(kRed);
   hAvgETEL->Draw("COLZ");
   hProfETEL->Draw("SAME");

   //Use a log scale in Z
   gPad->SetLogz();
   gPad->SetLogy();
   gPad->Update();

   fout->cd();
   c_avg_etel->Write("c_avg_etel");
   c_avg_etel->Print("ETEL_Average_SPE.pdf","pdf");
   
   c_avg_all->cd(2);
   hAvgETEL->Draw("COLZ");
   hProfETEL->Draw("SAME");
   gPad->SetLogz();
   gPad->SetLogy();

   //Show average ETEL distribution
   TCanvas *c_avg_hm = new TCanvas("c_avg_hm","Average HM SPE Distribution");
   //Compute average
   std::string hm = "HM";
   TH2F *hAvgHM = Average(good_spe_HM,hm);
   TProfile *hProfHM = hAvgHM->ProfileX("hAvgHM");

   //Turn off stats
   hAvgHM->SetStats(kFALSE);

   //Draw the 2D hist with the average prof on top.
   hProfHM->SetLineColor(kRed);
   hProfHM->SetMarkerColor(kRed);
   hAvgHM->Draw("COLZ");
   hProfHM->Draw("SAME");

   //Use a log scale in Z
   gPad->SetLogz();
   gPad->SetLogy();
   gPad->Update();

   fout->cd();
   c_avg_hm->Write("c_avg_hm");
   c_avg_hm->Print("HM_Average_SPE.pdf","pdf");

   c_avg_all->cd(3);
   hAvgHM->Draw("COLZ");
   hProfHM->Draw("SAME");
   gPad->SetLogz();
   gPad->SetLogy();


   //Show average ETEL distribution
   TCanvas *c_avg_wb = new TCanvas("c_avg_wb","Average WB SPE Distribution");
   //Compute average
   std::string wb = "WB";
   TH2F *hAvgWB = Average(good_spe_WB,wb);
   TProfile *hProfWB = hAvgWB->ProfileX("hAvgWB");

   //Turn off stats
   hAvgWB->SetStats(kFALSE);

   //Draw the 2D hist with the average prof on top.
   hProfWB->SetLineColor(kRed);
   hProfWB->SetMarkerColor(kRed);
   hAvgWB->Draw("COLZ");
   hProfWB->Draw("SAME");

   //Use a log scale in Z
   gPad->SetLogz();
   gPad->SetLogy();
   gPad->Update();

   fout->cd();
   c_avg_wb->Write("c_avg_wb");
   c_avg_wb->Print("WB_Average_SPE.pdf","pdf");

   c_avg_all->cd(4);
   hAvgWB->Draw("COLZ");
   hProfWB->Draw("SAME");
   gPad->SetLogz();
   gPad->SetLogy();

   TCanvas *c_avg_profiles = new TCanvas("c_avg_profiles","Average SPE distributions",900,600);
   c_avg_profiles->cd();
   hProf->SetLineColor(kBlack);
   hProfWB->SetLineColor(kRed);
   hProfHM->SetLineColor(kBlue);
   hProfETEL->SetLineColor(kOrange);
   hProf->SetMarkerColor(kBlack);
   hProfWB->SetMarkerColor(kRed);
   hProfHM->SetMarkerColor(kBlue);
   hProfETEL->SetMarkerColor(kOrange);
   hProf->SetLineWidth(2);
   hProfWB->SetLineWidth(2);
   hProfHM->SetLineWidth(2);
   hProfETEL->SetLineWidth(2);
   hProf->Draw();
   hProfWB->Draw("same");
   hProfHM->Draw("same");
   hProfETEL->Draw("same");
   c_avg_profiles->SetLogy();

   c_avg_all->Print("Average_SPE_Distributions_AllTypes.pdf","pdf");
   c_avg_profiles->Print("Average_SPE_Profiles_AllTypes.pdf","pdf");

   fout->cd();
   c_avg_all->Write("c_avg_all");
   c_avg_profiles->Write("c_avg_profiles");
   hProf->Write();
   hProfWB->Write();
   hProfHM->Write();
   hProfETEL->Write();   


}
