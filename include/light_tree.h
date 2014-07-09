#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

using namespace std ;

//-----Declaration of leaf types
//---global---
float scale=1;
//---gen---
//---Ws 
float gen_Wh_pt=0;
float gen_Wh_eta=0;
float gen_Wh_phi=0;
float gen_Wh_mass=0;
float gen_Wl_pt=0;
float gen_Wl_eta=0;
float gen_Wl_phi=0;
float gen_Wl_mass=0;
float gen_X_mass=0;
//---lep
float gen_lep_pt=0;
float gen_lep_eta=0;
float gen_lep_phi=0;
float gen_nu_pt=0;
float gen_nu_eta=0;
float gen_nu_phi=0;
//---W's quark
float gen_W_q1_pt=0;
float gen_W_q1_eta=0;
float gen_W_q1_phi=0;
float gen_W_q2_pt=0;
float gen_W_q2_eta=0;
float gen_W_q2_phi=0;
//---vbf quark
float gen_vbf_q1_pt=0;
float gen_vbf_q1_eta=0;
float gen_vbf_q1_phi=0;
float gen_vbf_q1_btag=0;
float gen_vbf_q2_pt=0;
float gen_vbf_q2_eta=0;
float gen_vbf_q2_phi=0;
float gen_vbf_q2_btag=0;
float gen_vbf_qq_mass=0;
//---reco----
//---lep 
float lep_pt=0;
float lep_eta=0;
float lep_phi=0;
float MET=0;
float MET_phi=0;
//---leptonic W
float lv_mass=0;
float lv_pt=0;
float lv_eta=0;
float lv_phi=0;
float lv_delta_R=0;
float lv_closerjet_mass=0;
//---CA8 jet
float CA8_jet_pt=0;
float CA8_jet_eta=0;
float CA8_jet_phi=0;
float CA8_jet_mass=0;
float CA8_jet_t2t1=0;
float CA8_jet_t3t2=0;
float CA8_closerjet_mass=0;
//---vbf jets
float vbf_jet1_pt=0;
float vbf_jet1_eta=0;
float vbf_jet1_phi=0;
float vbf_jet1_mass=0;
float vbf_jet1_btag=0;
float vbf_jet2_pt=0;
float vbf_jet2_eta=0;
float vbf_jet2_phi=0;
float vbf_jet2_mass=0;
float vbf_jet2_btag=0;
float vbf_jj_mass=0;
float vbf_jj_delta_R=0;
//---lvj
float lvJ_mass=0; 

//****************************************************************************************

void InitLightTree (TTree* fTree)
{
    //---global
    fTree->Branhc("scale", &scale, "scale/F");
    //---gen
    fTree->Branch("gen_Wh_pt", &gen_Wh_pt, "gen_Wh_pt/F");
    fTree->Branch("gen_Wh_eta", &gen_Wh_eta, "gen_Wh_eta/F");
    fTree->Branch("gen_Wh_phi", &gen_Wh_phi, "gen_Wh_phi/F");
    fTree->Branch("gen_Wh_mass", &gen_Wh_mass, "gen_Wh_mass/F");
    fTree->Branch("gen_Wl_pt", &gen_Wl_pt, "gen_Wl_pt/F");
    fTree->Branch("gen_Wl_eta", &gen_Wl_eta, "gen_Wl_eta/F");
    fTree->Branch("gen_Wl_phi", &gen_Wl_phi, "gen_Wl_phi/F");
    fTree->Branch("gen_Wl_mass", &gen_Wl_mass, "gen_Wl_mass/F");
    fTree->Branch("gen_X_mass", &gen_X_mass, "gen_X_mass/F");
    fTree->Branch("gen_lep_pt", &gen_lep_pt, "gen_lep_pt/F");
    fTree->Branch("gen_lep_eta", &gen_lep_eta, "gen_lep_eta/F");
    fTree->Branch("gen_lep_phi", &gen_lep_phi, "gen_lep_phi/F"); 
    fTree->Branch("gen_nu_pt", &gen_nu_pt, "gen_nu_pt/F");
    fTree->Branch("gen_nu_eta", &gen_nu_eta, "gen_nu_eta/F");
    fTree->Branch("gen_nu_phi", &gen_nu_phi, "gen_nu_phi/F"); 
    fTree->Branch("gen_W_q1_pt", &gen_W_q1_pt, "gen_W_q1_pt/F");
    fTree->Branch("gen_W_q1_eta", &gen_W_q1_eta, "gen_W_q1_eta/F");
    fTree->Branch("gen_W_q1_phi", &gen_W_q1_phi, "gen_W_q1_phi/F"); 
    fTree->Branch("gen_W_q2_pt", &gen_W_q2_pt, "gen_W_q2_pt/F");
    fTree->Branch("gen_W_q2_eta", &gen_W_q2_eta, "gen_W_q2_eta/F");
    fTree->Branch("gen_W_q2_phi", &gen_W_q2_phi, "gen_W_q2_phi/F"); 
    fTree->Branch("gen_vbf_q1_pt", &gen_vbf_q1_pt, "gen_vbf_q1_pt/F");
    fTree->Branch("gen_vbf_q1_eta", &gen_vbf_q1_eta, "gen_vbf_q1_eta/F");
    fTree->Branch("gen_vbf_q1_phi", &gen_vbf_q1_phi, "gen_vbf_q1_phi/F");
    fTree->Branch("gen_vbf_q1_btag", &gen_vbf_q1_btag, "gen_vbf_q1_btag/F"); 
    fTree->Branch("gen_vbf_q2_pt", &gen_vbf_q2_pt, "gen_vbf_q2_pt/F");
    fTree->Branch("gen_vbf_q2_eta", &gen_vbf_q2_eta, "gen_vbf_q2_eta/F");
    fTree->Branch("gen_vbf_q2_phi", &gen_vbf_q2_phi, "gen_vbf_q2_phi/F"); 
    fTree->Branch("gen_vbf_q2_btag", &gen_vbf_q2_btag, "gen_vbf_q2_btag/F");
    fTree->Branch("gen_vbf_qq_mass", &gen_vbf_qq_mass, "gen_vbf_qq_mass/F"); 
    //---reco
    fTree->Branch("lep_pt", &lep_pt, "lep_pt/F");
    fTree->Branch("lep_eta", &lep_eta, "lep_eta/F");
    fTree->Branch("lep_phi", &lep_phi, "lep_phi/F");
    fTree->Branch("MET", &MET, "MET/F");
    fTree->Branch("MET_phi", &MET_phi, "MET_phi/F");
    fTree->Branch("lv_mass", &lv_mass, "lv_mass/F");
    fTree->Branch("lv_pt", &lv_pt, "lv_pt/F");
    fTree->Branch("lv_eta", &lv_eta, "lv_eta/F");
    fTree->Branch("lv_phi", &lv_phi, "lv_phi/F");
    fTree->Branch("lv_deltaR", &lv_delta_R, "lv_delta_R/F");
    fTree->Branch("lv_closerjet_mass", &lv_closerjet_mass, "lv_closerjet_mass/F");
    fTree->Branch("CA8_jet_pt", &CA8_jet_pt, "CA8_jet_pt/F");
    fTree->Branch("CA8_jet_eta", &CA8_jet_eta, "CA8_jet_eta/F");
    fTree->Branch("CA8_jet_phi", &CA8_jet_phi, "CA8_jet_phi/F");
    fTree->Branch("CA8_jet_mass", &CA8_jet_mass, "CA8_jet_mass/F");
    fTree->Branch("CA8_jet_t2t1", &CA8_jet_t2t1, "CA8_jet_t2t1/F");
    fTree->Branch("CA8_jet_t3t2", &CA8_jet_t3t2, "CA8_jet_t3t2/F");
    fTree->Branch("CA8_closerjet_mass", &CA8_closerjet_mass, "CA8_closerjet_mass/F");
    fTree->Branch("vbf_jet1_pt", &vbf_jet1_pt, "vbf_jet1_pt/F");
    fTree->Branch("vbf_jet1_eta", &vbf_jet1_eta, "vbf_jet1_eta/F");
    fTree->Branch("vbf_jet1_phi", &vbf_jet1_phi, "vbf_jet1_phi/F");
    fTree->Branch("vbf_jet1_mass", &vbf_jet1_mass, "vbf_jet1_mass/F");
    fTree->Branch("vbf_jet1_btag", &vbf_jet1_btag, "vbf_jet1_btag/F");
    fTree->Branch("vbf_jet2_pt", &vbf_jet2_pt, "vbf_jet2_pt/F");
    fTree->Branch("vbf_jet2_eta", &vbf_jet2_eta, "vbf_jet2_eta/F");
    fTree->Branch("vbf_jet2_phi", &vbf_jet2_phi, "vbf_jet2_phi/F");
    fTree->Branch("vbf_jet2_mass", &vbf_jet2_mass, "vbf_jet2_mass/F");
    fTree->Branch("vbf_jet2_btag", &vbf_jet2_btag, "vbf_jet2_btag/F");
    fTree->Branch("vbf_jj_mass", &vbf_jj_mass, "vbf_jj_mass/F");
    fTree->Branch("vbf_jj_delta_R", &vbf_jj_delta_R, "vbf_jj_delta_R/F");
    fTree->Branch("lvJ_mass", &lvJ_mass, "lvJ_mass/F"); 
}

