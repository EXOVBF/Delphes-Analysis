/***************************************************************************************** 

    c++ -O2 -lm `root-config --cflags --glibs` -o convertCMSSWTree convertCMSSWTree.cpp 
    
    or exexute compile.sh cpp/convertCMSSWTree.cpp

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

#include "../include/data_tree.h"
#include "../include/light_tree.h"

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

//*****************************************************************************************

int main(int argc, char* argv[])
{
//-----------------Definitions------------------------------------------------------------
 
    if(argc < 4)
    {
	cout << "ERROR: missing one or more input files" << endl
	     << "USAGE: convertCMSSWTree input.root output.root weights.txt" << endl;
	return 0;
    }
    TFile* inFile = TFile::Open(argv[1], "read");
    data_tree* otree = new data_tree();
    otree->Init((TTree*)inFile->Get("otree"));
    TFile* outFile = TFile::Open(argv[2], "recreate");
    outFile->cd();
    TTree* LT = new TTree("light_tree", "light_tree");
    InitLightTree(LT);
    //---inport delphes' weights file
    vector<float> wgt_l_edge, wgt_value;
    float wgt_buffer=0;
    ifstream weights(argv[3], ios::in);
    while(weights >> wgt_buffer)
    {
	wgt_l_edge.push_back(wgt_buffer);
	weights >> wgt_buffer;
	wgt_value.push_back(wgt_buffer);
    }
    for(int iEntry=0; iEntry<otree->GetEntries(); iEntry++)
    {
	if(iEntry % 100000 == 0)
	    cout << "read " << iEntry << " / " << otree->GetEntries() << endl;
	otree->GetEntry(iEntry);
	//-----apply the same preselection of Delphes----
	//---skip 0,1 tag jet event
	if(otree->vbf_maxpt_j1_pt == 0 || otree->vbf_maxpt_j2_pt == 0)
	    continue;
	//---CA8_jet eta check
	if(fabs(otree->ungroomed_jet_eta) > 2.4)
	    continue;
	//-----Store reco variables----
	//---vertex
	nPV = otree->nPV;
	//---cmssw weight
	//evt_weight = otree->eff_and_pu_Weight*otree->btag_weight;
	evt_weight = 1;
	//---delphes weight
	int iBin=0;       
	while(iBin<wgt_l_edge.size() && nPV >= wgt_l_edge.at(iBin))
	{
	    iBin++;
	}
	evt_weight = evt_weight*wgt_value.at(iBin-1);
	//---leptons
	lep_pt = otree->l_pt;
	lep_eta = otree->l_eta;
	lep_phi = otree->l_phi;
	//---MET
	MET = otree->pfMET;
	MET_phi = otree->pfMET_Phi;
	//---leptonic W reco
	lv_mass = otree->mass_lv_subj_type0;
	lv_pt = otree->v_pt;
	lv_eta = otree->v_eta;
	lv_phi = otree->v_phi;
	lv_Mt = otree->v_mt;
	//FAKE value
	lv_delta_R = DeltaR(lep_eta, 0,
			    lep_phi, 0);
	//---Wl closestjet
	lv_closestjet_mass = otree->mass_leptonic_closerjet;
	//---CA8 jet (hadronic W reco)
	CA8_jet_pt = otree->ungroomed_jet_pt;
	CA8_jet_eta = otree->ungroomed_jet_eta;
	CA8_jet_phi = otree->ungroomed_jet_phi;
	CA8_jet_mass = otree->jet_mass_pr;
	CA8_jet_t2t1 = otree->jet_tau2tau1;
	CA8_jet_t3t2 = 1;
	//---CA8 closestjet
	CA8_closestjet_mass = otree->mass_ungroomedjet_closerjet; 
	//---object separation
	lv_J_delta_phi = otree->deltaphi_Vca8jet;
	MET_J_delta_phi = otree->deltaphi_METca8jet;
	l_J_delta_R = otree->deltaR_lca8jet;
	//---mlvJ
	lvJ_mass = otree->mass_lvj_type0;
	//---vbf jets
	vbf_jet1_pt = otree->vbf_maxpt_j1_pt;
	vbf_jet1_eta = otree->vbf_maxpt_j1_eta;
	vbf_jet1_phi = otree->vbf_maxpt_j1_phi;
	vbf_jet1_mass = otree->vbf_maxpt_j1_m;
	if(otree->vbf_maxpt_j1_bDiscriminatorCSV < 0.679)
	    vbf_jet1_btag = 0;
	else
	    vbf_jet1_btag = 1;
	vbf_jet2_pt = otree->vbf_maxpt_j2_pt;
	vbf_jet2_eta = otree->vbf_maxpt_j2_eta;
	vbf_jet2_phi = otree->vbf_maxpt_j2_phi;
	vbf_jet2_mass = otree->vbf_maxpt_j2_m;
	if(otree->vbf_maxpt_j2_bDiscriminatorCSV < 0.679)
	    vbf_jet2_btag = 0;
	else
	    vbf_jet2_btag = 1;
	vbf_jj_mass = otree->vbf_maxpt_jj_m;
	vbf_jj_delta_eta = fabs(vbf_jet1_eta - vbf_jet2_eta);
	vbf_jj_delta_phi = DeltaPhi(vbf_jet1_phi, vbf_jet2_phi);
	vbf_jj_delta_R = DeltaR(vbf_jet1_eta, vbf_jet2_eta,
				vbf_jet1_phi, vbf_jet2_phi);

	LT->Fill();
    }
    cout << "Total passed events: " << LT->GetEntriesFast() << endl;
    LT->Write();
    outFile->Close();

    return 0 ;
}


