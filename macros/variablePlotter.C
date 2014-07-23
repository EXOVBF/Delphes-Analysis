int variablePlotter()
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
    TH1F* histo_data;
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
    //TString all_cuts = "evt_weight*(";
    TString all_cuts = " ";
    for(int iCut=0; iCut<cut_value.size(); iCut++)
    {
	if(TString(cuts).Contains(cut_number.at(iCut)))
	{
	    //if(all_cuts == "evt_weight*(")
	    if(all_cuts == " ")
		all_cuts += cut_value.at(iCut);
	    else
		all_cuts += " && " + cut_value.at(iCut);
	}
    }
    //all_cuts += ")";
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
	    string unit, step;
	    plot_conf >> S_buffer;
	    plot_conf >> nbins >> low >> high;
	    plot_conf >> unit >> step;
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
		    histo_data = histo1D[iSample];
		cout << "events for " << samples.at(iSample) 
		     << " sample: " << histo1D[iSample]->GetEntries()*norm.at(iSample) << endl;
	    }
	    cout << "total background: " << tot_bkg << endl;    
	    //-----Draw-----
	    if(histo_data)
	    {
		cout << "bkg/data ratio: " << tot_bkg/histo_data->GetEntries() << endl;
		for(int iBin=1; iBin<histo_data->GetNbinsX(); iBin++)
		{
		    data_point->SetPoint(iBin-1, 
					 histo_data->GetBinCenter(iBin),
					 histo_data->GetBinContent(iBin));
		}
		data_point->SetMarkerColor(kBlack);
		data_point->SetMarkerStyle(20);
		data_point->SetMarkerSize(1);
	    }
	    stack->Draw("hist");
	    stack->GetXaxis()->SetTitle(x_title);
	    stack->GetYaxis()->SetTitle(y_title);
	    stack->SetMinimum(0);
	    stack->SetMaximum(histo_errors->GetBinError(histo_errors->GetMaximumBin())+
			      histo_errors->GetMaximum());
	    gPad->Update();
	    histo_errors->SetLineColor(kBlack);
	    histo_errors->Draw("E1X0same");
	    if(histo_data)
		data_point->Draw("sameP");
	    char plot_name[20];
	    sprintf(plot_name, "plot/%s.png", S_buffer.c_str());
	    c1->Print(plot_name, "png");
	    histo_errors->Delete();
	    stack->Delete();
	    if(histo_data)
		histo_data->Delete();
	}
    }
    //-----Exit-----
    c1->Delete();
    plot_conf.close();
    return 0;
}
