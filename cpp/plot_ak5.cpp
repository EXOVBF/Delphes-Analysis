/*****************************************************************************************

    c++ -O2 -lm `root-config --cflags --glibs` -o makePlots makePlots.cpp 
    
    or exexute compile.sh makePlots

/****************************************************************************************/

#include<vector>
#include<iterator>
#include<string>
#include<iostream>
#include<fstream>
#include<cmath>
#include<algorithm>

#include "TChain.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "../include/delphes_tree.h"

using namespace std ;

//****************************************************************************************

float DeltaPhi (float phi1, float phi2)
{
    float delta_phi = TMath::Abs(phi1 - phi2);
    if (delta_phi > 2*TMath::Pi()) 
	delta_phi -= 2*TMath::Pi();
    if (delta_phi > TMath::Pi() && delta_phi < 2*TMath::Pi()) 
	delta_phi = 2*TMath::Pi() - delta_phi;
    return delta_phi;
}

//****************************************************************************************

float DeltaR (float eta1, float eta2, float phi1, float phi2)
{
    float d_phi = DeltaPhi(phi1, phi2);
    float d_eta = TMath::Abs(eta1 - eta2);

    return TMath::Sqrt(pow(d_eta,2)+pow(d_phi,2));
}

//****************************************************************************************

int CreateHistos1D (TString suffix, map<TString, TH1F *> & histos, 
                    vector<string>* var, vector<int>* nbin, vector<float>* min, vector<float>* max) 
{
    //---set sample colors
    int sample_color=0;
    if( suffix.Contains("tt") == 1 ) sample_color=3;
    if( suffix.Contains("W") == 1 ) sample_color=2;
    if( suffix.Contains("Z") == 1) sample_color=5;
    //---create histogram
    for(int iVar=0; iVar<var->size(); iVar++)
    {
        if( suffix.Contains("Signal")==0 && ((TString)(var->at(iVar))).Contains("MWW")==1) continue;
        TH1F * h_tmp = new TH1F (TString (var->at(iVar)) + "__"+suffix, var->at(iVar) + "__"+suffix, 
				 nbin->at(iVar), min->at(iVar), max->at(iVar)); 
        h_tmp->SetFillColor(sample_color);
        histos[(var->at(iVar)).c_str()] = h_tmp;
    }
    return histos.size();
}

//****************************************************************************************

int CreateHistos2D (TString suffix, map<TString, TH2F *> & histos, 
                    vector<string>* var, vector<int>* nbin_x, vector<float>* min_x, vector<float>* max_x,
		    vector<int>* nbin_y, vector<float>* min_y, vector<float>* max_y)
{
    for(int iVar=0; iVar<var->size(); iVar++)
    {
        if( strcmp(suffix.Data(),"W4jets")==0 && ((TString)(var->at(iVar))).Contains("WJ")==1) continue;
        if( strcmp(suffix.Data(),"Z")==0 && ((TString)(var->at(iVar))).Contains("WJ")==1) continue;
        TH2F * h_tmp = new TH2F (TString (var->at(iVar)) + "__"+suffix, var->at(iVar) + "__"+suffix, 
				 nbin_x->at(iVar), min_x->at(iVar), max_x->at(iVar), 
				 nbin_y->at(iVar), min_y->at(iVar), max_y->at(iVar)); 
        histos[(var->at(iVar)).c_str()] = h_tmp;
    }
    return histos.size();
}

//****************************************************************************************

void SaveHistos1D (TFile * outfile, map<TString, TH1F *> histos, float cross_section = 1.)
{
    outfile->cd () ;
    for (map<TString, TH1F *>::iterator iMap = histos.begin () ;
         iMap != histos.end () ; ++iMap)
    {
        //iMap->second->Scale (cross_section / iMap->second->Integral ()) ;
        iMap->second->Write () ;
    }
    return ;
}

//****************************************************************************************

void SaveHistos2D (TFile * outfile, map<TString, TH2F *> histos, float cross_section = 1.)
{
    outfile->cd () ;
    for (map<TString, TH2F *>::iterator iMap = histos.begin () ;
         iMap != histos.end () ; ++iMap)
    {
        iMap->second->Write () ;
    }
    return ;
}

//*****************************************************************************************

int readDataset (TString sampleName, vector<TString> datasetBaseName, map<TString, TH1F *> histos_1D, map<TString, TH2F *> histos_2D) 
{
//-----------------Definitions------------------------------------------------------------

    //---vertex
    float X_ver_tmp=0, Y_ver_tmp=0, Z_ver_tmp=0;
    float d_xy_tmp=0;
    int S_ver_tmp=0;
    //---leptons
    int n_lept_tmp=0, n_lepl_tmp=0, good_lep_tmp=0, nu_pz_comp_type=0;
    //---quarks
    int q1_tmp=0, q2_tmp=0;
    float q1_pt_tmp=-1, q2_pt_tmp=-1;
    vector<float> std_vect_tmp;
    vector<float>::iterator max_tmp;
    //---jets
    float deltaR_jets_tmp=0;
    int n_CA8_tmp=0, n_ak5_tmp=0;
    vector<int> CA8_tmp, ak5_tmp;
    //---Ws
    float deltaR_Wl=0, deltaR_Wlv=0, deltaR_Wj=0, deltaR_Wjgen=0, W_lep_tmp=0, W_had_tmp=0;
    //---counters
    int countPSEvents=0, countBSEvents=0;
    
    TChain * ch = new TChain ("Delphes") ;
    for(int i=0; i<datasetBaseName.size(); i++)
    {
        ch->Add(datasetBaseName.at(i)+"*root") ;
    }
    cout << "read " << ch->GetEntries() << " events in " << ch->GetNtrees() 
         << " files of " << sampleName.Data() << " sample\n" ;
    delphes_tree DT ;
    if(sampleName.Contains("Signal")) DT.Init(ch, true);
    else DT.Init(ch, false);

//-----------------Events loop------------------------------------------------------------
 
    for (int iEvent = 0 ; iEvent < ch->GetEntries() ; iEvent++)
    {
        ch->GetEntry(iEvent);
        if (iEvent % 100000 == 0) cout << "reading event number " << iEvent << " / " << ch->GetEntries() << "\n" ;

        //---Reset---
        X_ver_tmp = 0;
        Y_ver_tmp = 0;
        Z_ver_tmp = 0;
        S_ver_tmp = 0;
        n_lept_tmp = 0;
        n_lepl_tmp = 0;
        n_CA8_tmp = 0;
        n_ak5_tmp = 0;
	q1_tmp = -1;
	q2_tmp = -1;
        q1_pt_tmp = -1;
        q2_pt_tmp = -1;
        std_vect_tmp.clear();
        CA8_tmp.clear();
        ak5_tmp.clear();
        //------------------Preselection--------------------------------------------------
        //---PV selection---
        for(int iVertex=0; iVertex<DT.nPV->size(); iVertex++)
        {
            //---choose highets pt^2 vertex
            if( DT.sum_pt_square->at(iVertex) > S_ver_tmp )
            {
                X_ver_tmp = DT.vertex_X->at(iVertex);
                Y_ver_tmp = DT.vertex_Y->at(iVertex);
                Z_ver_tmp = DT.vertex_Z->at(iVertex);
            }
        }
        //---Selec leptons from PV---        
        for(int iLep=0; iLep<DT.lep_number->at(0); iLep++)
        {        
            d_xy_tmp = TMath::Sqrt(pow(X_ver_tmp - DT.lep_X_vertex->at(iLep),2) + pow(X_ver_tmp - DT.lep_X_vertex->at(iLep),2));
            //---tight electron
            if( abs(DT.lep_flv->at(iLep)) == 11 && abs(DT.lep_Z_vertex->at(iLep)-Z_ver_tmp) < 0.1 && DT.lep_pt->at(iLep) > 35 && 
                (DT.lep_eta->at(iLep) < 1.4442 || DT.lep_eta->at(iLep) > 1.566) && DT.lep_isolation->at(iLep)*DT.lep_pt->at(iLep) < 5 &&
                d_xy_tmp < 0.2 && DT.MET->at(0) > 65) 
            {
                good_lep_tmp = iLep;
                n_lept_tmp++;
            }
            //---loose electron
            else if( abs(DT.lep_flv->at(iLep)) == 11 && abs(DT.lep_Z_vertex->at(iLep)-Z_ver_tmp) < 0.1 && DT.lep_pt->at(iLep) > 20 && 
		     DT.lep_eta->at(iLep) < 2.5 && DT.lep_isolation->at(iLep)*DT.lep_pt->at(iLep) < 5 && d_xy_tmp < 0.2) 
            {
                n_lepl_tmp++;
            }
            //---tight muon
            if( abs(DT.lep_flv->at(iLep)) == 13 && abs(DT.lep_Z_vertex->at(iLep)-Z_ver_tmp) < 5 && DT.lep_pt->at(iLep) > 50 && 
		DT.lep_eta->at(iLep) < 2.1 && DT.lep_isolation->at(iLep) < 0.1 && d_xy_tmp < 2 && DT.MET->at(0) > 50)
            {
                good_lep_tmp = iLep;
                n_lept_tmp++;
            }
            //---loose muon
            else if( abs(DT.lep_flv->at(iLep)) == 13 && abs(DT.lep_Z_vertex->at(iLep)-Z_ver_tmp) < 5 && DT.lep_pt->at(iLep) > 10 && 
                     DT.lep_eta->at(iLep) < 2.5 && DT.lep_isolation->at(iLep) < 0.1 && d_xy_tmp < 2)
            {
                n_lepl_tmp++;
            }
        }
	if(n_lepl_tmp > 0 || n_lept_tmp != 1)
	    continue;
        //---select good ak5---
	float pt_sum=1000;
//	ak5_tmp.push_back(0);
//	ak5_tmp.push_back(0);
	TLorentzVector jj_Wh_4vect_tmp, j1_Wh_4vect_tmp, j2_Wh_4vect_tmp;
        for(int iak5=0; iak5<DT.number_jet_ak5->at(0); iak5++)
        {
	    if(fabs(DT.jet_ak5_eta->at(iak5)) > 4.7)
		continue;
//            j1_Wh_4vect_tmp.SetPtEtaPhiM(DT.jet_ak5_pt->at(iak5), DT.jet_ak5_eta->at(iak5),
//					 DT.jet_ak5_phi->at(iak5), DT.jet_ak5_mass_pruned->at(iak5));
/*
	    for(int i=iak5+1; i<DT.number_jet_ak5->at(0); i++)
            {
		if(fabs(DT.jet_ak5_eta->at(i)) > 2.4)
		    continue;
	        j2_Wh_4vect_tmp.SetPtEtaPhiM(DT.jet_ak5_pt->at(i), DT.jet_ak5_eta->at(i),
					      DT.jet_ak5_phi->at(i), DT.jet_ak5_mass_pruned->at(i));
                jj_Wh_4vect_tmp = j1_Wh_4vect_tmp + j2_Wh_4vect_tmp;
		if(TMath::Abs(jj_Wh_4vect_tmp.M()-80) < pt_sum)
		{
		    ak5_tmp.at(0) = iak5;
		    ak5_tmp.at(1) = i;
		    pt_sum = TMath::Abs(jj_Wh_4vect_tmp.M()-80);
		    n_ak5_tmp = 2;
		    }
	    }*/
	    if(n_ak5_tmp < 2 && DT.jet_ak5_btag->at(iak5) == 0)
	    {
		ak5_tmp.push_back(iak5);
		n_ak5_tmp++;
	    }
	}
	//---order quark wrt pt
	vector<int> ord_quark_tmp;
	ord_quark_tmp.push_back(0);
	for(int iPart=1; iPart<DT.lhe_p_pt->size(); iPart++)
	{
	    bool insert=0;
	    for(int iq=0; iq<ord_quark_tmp.size(); iq++)
	    {
		if(DT.lhe_p_pt->at(iPart) > DT.lhe_p_pt->at(ord_quark_tmp.at(iq)))
		{
		    ord_quark_tmp.insert(ord_quark_tmp.begin()+iq, iPart);
		    insert = 1;
		    break;
		}
	    }
	    if(insert == 0)
		ord_quark_tmp.push_back(iPart);
	}
	for(int iPart=0; iPart<ord_quark_tmp.size(); iPart++)
	{
	    if(q1_tmp == -1 && DT.lhe_p_from_W->at(ord_quark_tmp.at(iPart)) == 1)
		q1_tmp = iPart; 
	    else if(q2_tmp == -1 && DT.lhe_p_from_W->at(ord_quark_tmp.at(iPart)) == 1)
	    {
		q2_tmp = iPart; 
		break;
	    }
	}
	if(n_ak5_tmp < 2 || q1_tmp == -1 || q2_tmp == -1)
	    continue;
	j1_Wh_4vect_tmp.SetPtEtaPhiM(DT.jet_ak5_pt->at(ak5_tmp.at(0)), DT.jet_ak5_eta->at(ak5_tmp.at(0)),
				     DT.jet_ak5_phi->at(ak5_tmp.at(0)), DT.jet_ak5_mass_pruned->at(ak5_tmp.at(0)));
	j2_Wh_4vect_tmp.SetPtEtaPhiM(DT.jet_ak5_pt->at(ak5_tmp.at(1)), DT.jet_ak5_eta->at(ak5_tmp.at(1)),
				     DT.jet_ak5_phi->at(ak5_tmp.at(1)), DT.jet_ak5_mass_pruned->at(ak5_tmp.at(1)));
	jj_Wh_4vect_tmp = j1_Wh_4vect_tmp + j2_Wh_4vect_tmp;
	//---Plots
	//---quark-jet matching
	if(histos_2D["ps_deltaR_j1q1_vs_j2q2"])
	    histos_2D["ps_deltaR_j1q1_vs_j2q2"]->Fill(DeltaR(j2_Wh_4vect_tmp.Eta(), DT.lhe_p_eta->at(q2_tmp),
							     j2_Wh_4vect_tmp.Phi(), DT.lhe_p_phi->at(q2_tmp)),
						      DeltaR(j1_Wh_4vect_tmp.Eta(), DT.lhe_p_eta->at(q1_tmp),
							     j1_Wh_4vect_tmp.Phi(), DT.lhe_p_phi->at(q1_tmp)));
	//---jets
	if(histos_1D["ps_W_j1_pt"])
	    histos_1D["ps_W_j1_pt"]->Fill(j1_Wh_4vect_tmp.Pt());
	if(histos_1D["ps_W_j1_eta"])
	    histos_1D["ps_W_j1_eta"]->Fill(j1_Wh_4vect_tmp.Eta());
	if(histos_1D["ps_W_j1_phi"])
	    histos_1D["ps_W_j1_phi"]->Fill(j1_Wh_4vect_tmp.Phi());
	if(histos_1D["ps_W_j2_pt"])
	    histos_1D["ps_W_j2_pt"]->Fill(j2_Wh_4vect_tmp.Pt());
	if(histos_1D["ps_W_j2_eta"])
	    histos_1D["ps_W_j2_eta"]->Fill(j2_Wh_4vect_tmp.Eta());
	if(histos_1D["ps_W_j2_phi"])
	    histos_1D["ps_W_j2_phi"]->Fill(j2_Wh_4vect_tmp.Phi());
	if(histos_1D["ps_W_jj_deltaR"])
	    histos_1D["ps_W_jj_deltaR"]->Fill(DeltaR(j1_Wh_4vect_tmp.Eta(),
						     j2_Wh_4vect_tmp.Eta(),
						     j1_Wh_4vect_tmp.Phi(),
						     j2_Wh_4vect_tmp.Phi()));
	if(histos_1D["ps_W_jj_M_pruned"])
	   histos_1D["ps_W_jj_M_pruned"]->Fill(jj_Wh_4vect_tmp.M());
	//---quarks
	if(histos_1D["ps_W_q1_pt"])
	    histos_1D["ps_W_q1_pt"]->Fill(DT.lhe_p_pt->at(q1_tmp));
	if(histos_1D["ps_W_q1_eta"])
	    histos_1D["ps_W_q1_eta"]->Fill(DT.lhe_p_eta->at(q1_tmp));
	if(histos_1D["ps_W_q1_phi"])
	    histos_1D["ps_W_q1_phi"]->Fill(DT.lhe_p_phi->at(q1_tmp));
	if(histos_1D["ps_W_q2_pt"])
	    histos_1D["ps_W_q2_pt"]->Fill(DT.lhe_p_pt->at(q2_tmp));
	if(histos_1D["ps_W_q2_eta"])
	    histos_1D["ps_W_q2_eta"]->Fill(DT.lhe_p_eta->at(q2_tmp));
	if(histos_1D["ps_W_q2_phi"])
	    histos_1D["ps_W_q2_phi"]->Fill(DT.lhe_p_phi->at(q2_tmp));
	if(histos_1D["ps_W_qq_deltaR"])
	   histos_1D["ps_W_qq_deltaR"]->Fill(DeltaR(DT.lhe_p_eta->at(q1_tmp),
						    DT.lhe_p_eta->at(q2_tmp),
						    DT.lhe_p_phi->at(q1_tmp),
						    DT.lhe_p_phi->at(q2_tmp))); 

	//---gen lep
	if(histos_1D["gen_lep_pt"])
	    histos_1D["gen_lep_pt"]->Fill(DT.lhe_lep_pt->at(0));
	//---reco lep
	if(histos_1D["ps_l_pt"])
	    histos_1D["ps_l_pt"]->Fill(DT.lep_pt->at(good_lep_tmp));
    }    
    cout << "number of events that passed the preselection:  " << countPSEvents << endl;    
    return 0 ;
}


//****************************************************************************************
//main

int main (int argc, char* argv[]) 
{
//-----------------Config files check-----------------------------------------------------
    TString configFile1;
    TString configFile2;
    TString outFile = "output.root";
    if(argc < 3)
    {
        cout << "ERROR: config files !!! " << endl;
        return 0;
    }
    configFile1 = argv[1];
    configFile2 = argv[2];
    cout << "CONFIG variables file:  " << configFile1.Data() << endl;
    cout << "CONFIG samples file:  " << configFile2.Data() << endl;
    if(argc > 3) outFile = argv[3];
    cout << "OUTPUT file: " << outFile.Data() << endl;
    
//-----------------Definitions------------------------------------------------------------
    //---variables
    string buffer_type;
    int buffer_nbin=0, buffer_nbin_x=0, buffer_nbin_y=0;
    float buffer_min=0, buffer_max=0, buffer_min_x=0, buffer_max_x=0, buffer_min_y=0, buffer_max_y=0;
    vector<int> nbin, nbin_x, nbin_y; 
    vector<float> min, max, min_x, max_x, min_y, max_y;
    string buffer_var_1d, buffer_var_2d;
    vector<string> var_1d, var_2d;
    //---samples
    string sample_path;
    string buffer_sample_dir, buffer_sample_name;
    vector<string> sample_dir, sample_name;
//-----------------Read Config files------------------------------------------------------
    //---variables
    ifstream config_file_1 (configFile1.Data(), ios::in);
    while(config_file_1 >> buffer_type)
    {   
        if(strcmp(buffer_type.c_str(),"#") == 0)
        {
            getline(config_file_1, buffer_type);
        }
        else if(strcmp(buffer_type.c_str(),"1D") == 0)
        {
            config_file_1 >> buffer_var_1d >> buffer_nbin >> buffer_min >> buffer_max;
            var_1d.push_back(buffer_var_1d);
            nbin.push_back(buffer_nbin);
            min.push_back(buffer_min);
            max.push_back(buffer_max);
        }
        else if(strcmp(buffer_type.c_str(),"2D") == 0)
        {
            config_file_1 >> buffer_var_2d >> buffer_nbin_x >> buffer_min_x >> buffer_max_x 
                          >> buffer_nbin_y >> buffer_min_y >> buffer_max_y;
            var_2d.push_back(buffer_var_2d);
            nbin_x.push_back(buffer_nbin_x);
            min_x.push_back(buffer_min_x);
            max_x.push_back(buffer_max_x);
            nbin_y.push_back(buffer_nbin_y);
            min_y.push_back(buffer_min_y);
            max_y.push_back(buffer_max_y);
        }
    }
    config_file_1.close();
    //---samples
    ifstream config_file_2 (configFile2.Data(), ios::in);
    while(config_file_2 >> buffer_type)
    {   
        if(strcmp(buffer_type.c_str(),"#") == 0)
        {
            getline(config_file_2, buffer_type);
        }
        else if(strcmp(buffer_type.c_str(),"P") == 0)
        {
            config_file_2 >> sample_path;
        }
        else if(strcmp(buffer_type.c_str(),"D") == 0)
        {
            config_file_2 >> buffer_sample_name >> buffer_sample_dir;
            sample_dir.push_back(buffer_sample_dir);
            sample_name.push_back(buffer_sample_name);
        }
    }
    config_file_2.close();

//-----------------Make Plots-------------------------------------------------------------
    TFile * outfile = new TFile (outFile.Data(), "recreate") ;
    for(int iSample=0; iSample<sample_dir.size(); iSample++)
    {
        TString sampleName(sample_name.at(iSample));
        TString samplePath(sample_path);
        map<TString, TH1F*> histos_1D_tmp;
        map<TString, TH2F*> histos_2D_tmp;
        vector<TString> sample_set;
        CreateHistos1D(sampleName, histos_1D_tmp, &var_1d, &nbin, &min, &max);
        CreateHistos2D(sampleName, histos_2D_tmp, &var_2d, &nbin_x, &min_x, &max_x, &nbin_y, &min_y, &max_y);
        samplePath += sample_dir.at(iSample);
        sample_set.push_back(samplePath);
        //---call to the analyzer function
        cout << "Analyzing:  " << sample_set.at(0).Data() << endl;
        readDataset(sampleName, sample_set, histos_1D_tmp, histos_2D_tmp);
        //---call to the plotter
        SaveHistos1D (outfile, histos_1D_tmp);
        SaveHistos2D (outfile, histos_2D_tmp);
    }
    outfile->Close () ;

    return 0 ;
}
