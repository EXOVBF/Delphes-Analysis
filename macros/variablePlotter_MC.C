int variablePlotter_MC()
{
    gROOT->Reset();
    gStyle->SetOptStat("");
    //---
    int I_buffer;
    double F_buffer;
    string S_buffer;    
    //---
    vector<string> cut_number;
    vector<string> cut_value;
    vector<string> samples; 
    vector<int> s_color;
    vector<int> s_line;
    vector<double> norm;
    //---
    string cuts;
    float tot_bkg=0;
    //---
    TString x_title("");
    TString y_title("");
    TH1F* histo_cmssw;
    TGraph* data_point = new TGraph();
    //-----load cuts-----
    ifstream cut_list("config/cut.list");
    while(cut_list >> S_buffer)
    {
	if(strcmp(S_buffer.c_str(), "#") == 0)
	    getline(cut_list, S_buffer);
	if(strcmp(S_buffer.c_str(), "C") == 0)
	{
	    cut_list >> S_buffer;
	    cut_number.push_back(S_buffer);
	    cut_list >> S_buffer;
	    getline(cut_list, S_buffer);
	    getline(cut_list, S_buffer);
	    cut_value.push_back(S_buffer);
	}
    }
    cut_list.close();
    //-----config-----
    ifstream plot_conf("config/macro_plot.conf", ios::in);
    //---get samples
    plot_conf >> S_buffer;
    while(strcmp(S_buffer.c_str(), "!") != 0)
    {
	if(strcmp(S_buffer.c_str(), "#") == 0)
	    getline(plot_conf, S_buffer);
	else
	{
	    samples.push_back(S_buffer);
	    plot_conf >> I_buffer;
	    s_color.push_back(I_buffer);
	    plot_conf >> I_buffer;
	    s_line.push_back(I_buffer);
	    plot_conf >> F_buffer;
	    norm.push_back(F_buffer);
	}
	plot_conf >> S_buffer;
    }
    TFile* sfile[10];
    TTree* stree[10];
    for(int iSample=0; iSample<samples.size(); iSample++)
    {
	sfile[iSample] = TFile::Open("light_ntuples/"+TString(samples.at(iSample))+".root", "r");
	stree[iSample] = (TTree*)sfile[iSample]->Get("light_tree");
    }
    //---get cuts
    plot_conf >> S_buffer;
    while(strcmp(S_buffer.c_str(), "!") != 0)
    {
	if(strcmp(S_buffer.c_str(), "#") == 0)
	    getline(plot_conf, S_buffer);
	else if(strcmp(S_buffer.c_str(), "C") == 0)
	{
	    plot_conf >> S_buffer;
	    cuts = S_buffer;
	}
	plot_conf >> S_buffer;
    }
    TString all_cuts = "evt_weight*(";
    //TString all_cuts =" ";
    for(int iCut=0; iCut<cut_value.size(); iCut++)
    {
	if(TString(cuts).Contains(cut_number.at(iCut)))
	{
	    if(all_cuts == "evt_weight*(")
		//if(all_cuts == " ")
		all_cuts += cut_value.at(iCut);
	    else
		all_cuts += " && " + cut_value.at(iCut);
	}
    }
    all_cuts += ")";
    cout << all_cuts.Data() << endl;
    //-----plot all variables-----
    while(plot_conf >> S_buffer)
    {
	TH1F* histo1D[10];
	TH2F* histo2D[10];
	if(strcmp(S_buffer.c_str(), "#") == 0)
	    getline(plot_conf, S_buffer);
	if(strcmp(S_buffer.c_str(), "1D") == 0)
	{
	    int nbins;
	    float low, high;
	    string unit, step, var_name;
	    plot_conf >> S_buffer;
	    plot_conf >> nbins >> low >> high;
	    plot_conf >> unit >> step;
	    var_name = S_buffer;
	    x_title = TString(S_buffer+" ["+unit+"]");
	    y_title = TString("counts/"+step+unit);
	    THStack* stack = new THStack(S_buffer.c_str(), S_buffer.c_str());
	    histo_errors = new TH1F("errors", "errors", nbins, low, high);
	    histo_errors->Sumw2();
	    for(int iSample=0; iSample<samples.size(); iSample++)
	    {
		sfile[iSample]->cd();
		char name[100];
		sprintf(name, "%s", (S_buffer+"_"+samples.at(iSample)).c_str());
		histo1D[iSample] = new TH1F(name, name, nbins, low, high);
		histo1D[iSample]->SetFillColor(s_color.at(iSample));
		if(s_line.at(iSample) == 0)
		    histo1D[iSample]->SetLineColor(s_color.at(iSample));
		else
		    histo1D[iSample]->SetLineColor(kBlack);
		char var[100];
		sprintf(var, "%s>>%s", S_buffer.c_str(), name);
		stree[iSample]->Draw(var, all_cuts.Data(), "goff");
		histo1D[iSample]->Sumw2();
		histo1D[iSample]->Scale(norm.at(iSample)/histo1D[iSample]->GetEntries());
		if(samples.at(iSample) != "data" && 
		   TString(samples.at(iSample)).Contains("CMSSW") != 1)
		{
		    stack->Add(histo1D[iSample]);
		    histo_errors->Add(histo1D[iSample]);
		    tot_bkg += histo1D[iSample]->GetEntries()*norm.at(iSample);
		}
		else
		{
		    histo_cmssw = histo1D[iSample];
		}
		cout << "events for " << samples.at(iSample) 
		     << " sample: " << histo1D[iSample]->GetEntries()*norm.at(iSample) << endl;
	    }
	    cout << "total background: " << tot_bkg << endl;    
	    //-----Draw-----
	    histo_errors->SetLineColor(kBlack);
	    histo_errors->SetMarkerStyle(20);
	    histo_errors->SetMarkerSize(0.7);
	    histo_cmssw->SetFillColor(kGreen);
	    histo_cmssw->Draw("hist");
	    histo_cmssw->SetTitle(var_name.c_str());
	    histo_cmssw->GetXaxis()->SetTitle(x_title);
	    histo_cmssw->GetYaxis()->SetTitle(y_title);
	    histo_cmssw->SetMinimum(0);
	    if(histo_errors->GetMaximum() > histo_cmssw->GetMaximum())
	    {
		histo_cmssw->SetMaximum((histo_errors->GetBinError(histo_errors->GetMaximumBin())+
					 histo_errors->GetMaximum())*1.1);
	    }
	    else
	    {
		histo_cmssw->SetMaximum((histo_cmssw->GetBinError(histo_cmssw->GetMaximumBin())+
					 histo_cmssw->GetMaximum())*1.1);
	    }
	    gPad->Update();
	    histo_errors->Draw("E1same");
	    //TLegend* lg = new TLegend(0.15, 0.6, 0.30, 0.85);
	    TLegend* lg = new TLegend(0.7, 0.6, 0.85, 0.85);
	    lg->SetBorderSize(0);
	    lg->SetFillStyle(0);
	    lg->AddEntry(histo_cmssw, "CMSSW", "f");
	    lg->AddEntry(histo_errors, "Delphes", "pl");
	    lg->Draw("same");
	    char plot_name[20];
	    sprintf(plot_name, "plot/%s.png", S_buffer.c_str());
	    gPad->Update();
	    c1->Print(plot_name, "png");
	    /*
	    ofstream weights("tmp/delphes_evt_weights_mu.txt", ios::out);
	    for(int iBin=1; iBin<=histo_cmssw->GetNbinsX(); iBin++)
	    {
		if(histo_cmssw->GetBinContent(iBin) != 0)
		{
		    weights << histo_cmssw->GetBinLowEdge(iBin) << "  " 
			    << (histo_errors->GetBinContent(iBin)/
				histo_cmssw->GetBinContent(iBin))
			    << endl;
		}
		else
		{
		    weights << histo_cmssw->GetBinLowEdge(iBin) << "  1" << endl;
		}
	    }
	    weights.close();*/
	    //---reset >> next variable---
	    histo_errors->Delete();
	    stack->Delete(); 
	    if(histo_cmssw)
	      histo_cmssw->Delete();
	}
    }
    //-----Exit-----
    plot_conf.close();
    return 0;
}
