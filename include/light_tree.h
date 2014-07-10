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
float gen_vbf_qq_delta_phi=0;
float gen_vbf_qq_delta_R=0;
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
float vbf_jj_delta_phi=0;
float vbf_jj_delta_R=0;
//---lvj
float lvJ_mass=0; 


//-----Declaration of branches 
//---global---
TBranch* b_scale;
//---gen---
//---Ws 
TBranch* b_gen_Wh_pt;
TBranch* b_gen_Wh_eta;
TBranch* b_gen_Wh_phi;
TBranch* b_gen_Wh_mass;
TBranch* b_gen_Wl_pt;
TBranch* b_gen_Wl_eta;
TBranch* b_gen_Wl_phi;
TBranch* b_gen_Wl_mass;
TBranch* b_gen_X_mass;
//---lep
TBranch* b_gen_lep_pt;
TBranch* b_gen_lep_eta;
TBranch* b_gen_lep_phi;
TBranch* b_gen_nu_pt;
TBranch* b_gen_nu_eta;
TBranch* b_gen_nu_phi;
//---W's quark
TBranch* b_gen_W_q1_pt;
TBranch* b_gen_W_q1_eta;
TBranch* b_gen_W_q1_phi;
TBranch* b_gen_W_q2_pt;
TBranch* b_gen_W_q2_eta;
TBranch* b_gen_W_q2_phi;
//---vbf quark
TBranch* b_gen_vbf_q1_pt;
TBranch* b_gen_vbf_q1_eta;
TBranch* b_gen_vbf_q1_phi;
TBranch* b_gen_vbf_q1_btag;
TBranch* b_gen_vbf_q2_pt;
TBranch* b_gen_vbf_q2_eta;
TBranch* b_gen_vbf_q2_phi;
TBranch* b_gen_vbf_q2_btag;
TBranch* b_gen_vbf_qq_mass;
TBranch* b_gen_vbf_qq_delta_phi;
TBranch* b_gen_vbf_qq_delta_R;
//---reco----
//---lep 
TBranch* b_lep_pt;
TBranch* b_lep_eta;
TBranch* b_lep_phi;
TBranch* b_MET;
TBranch* b_MET_phi;
//---leptonic W
TBranch* b_lv_mass;
TBranch* b_lv_pt;
TBranch* b_lv_eta;
TBranch* b_lv_phi;
TBranch* b_lv_delta_R;
TBranch* b_lv_closerjet_mass;
//---CA8 jet
TBranch* b_CA8_jet_pt;
TBranch* b_CA8_jet_eta;
TBranch* b_CA8_jet_phi;
TBranch* b_CA8_jet_mass;
TBranch* b_CA8_jet_t2t1;
TBranch* b_CA8_jet_t3t2;
TBranch* b_CA8_closerjet_mass;
//---vbf jets
TBranch* b_vbf_jet1_pt;
TBranch* b_vbf_jet1_eta;
TBranch* b_vbf_jet1_phi;
TBranch* b_vbf_jet1_mass;
TBranch* b_vbf_jet1_btag;
TBranch* b_vbf_jet2_pt;
TBranch* b_vbf_jet2_eta;
TBranch* b_vbf_jet2_phi;
TBranch* b_vbf_jet2_mass;
TBranch* b_vbf_jet2_btag;
TBranch* b_vbf_jj_mass;
TBranch* b_vbf_jj_delta_phi;
TBranch* b_vbf_jj_delta_R;
//---lvj
TBranch* b_lvJ_mass; 

//****************************************************************************************
 
void InitLightTree (TTree* fTree)
{
    //---global
    fTree->Branch("scale", &scale, "scale/F");
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
    fTree->Branch("gen_vbf_qq_delta_phi", &gen_vbf_qq_delta_phi, "gen_vbf_qq_delta_phi/F"); 
    fTree->Branch("gen_vbf_qq_delta_R", &gen_vbf_qq_delta_R, "gen_vbf_qq_delta_R/F"); 
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
    fTree->Branch("vbf_jj_delta_phi", &vbf_jj_delta_phi, "vbf_jj_delta_phi/F");
    fTree->Branch("vbf_jj_delta_R", &vbf_jj_delta_R, "vbf_jj_delta_R/F");
    fTree->Branch("lvJ_mass", &lvJ_mass, "lvJ_mass/F"); 
}

//****************************************************************************************
  
void SetLightTree (TTree* fTree)
{
    //---global
    fTree->SetBranchAddress("scale", &scale, &b_scale);
    //---gen
    fTree->SetBranchAddress("gen_Wh_pt", &gen_Wh_pt, &b_gen_Wh_pt);
    fTree->SetBranchAddress("gen_Wh_eta", &gen_Wh_eta, &b_gen_Wh_eta);
    fTree->SetBranchAddress("gen_Wh_phi", &gen_Wh_phi, &b_gen_Wh_phi);
    fTree->SetBranchAddress("gen_Wh_mass", &gen_Wh_mass, &b_gen_Wh_mass);
    fTree->SetBranchAddress("gen_Wl_pt", &gen_Wl_pt, &b_gen_Wl_pt);
    fTree->SetBranchAddress("gen_Wl_eta", &gen_Wl_eta, &b_gen_Wl_eta);
    fTree->SetBranchAddress("gen_Wl_phi", &gen_Wl_phi, &b_gen_Wl_phi);
    fTree->SetBranchAddress("gen_Wl_mass", &gen_Wl_mass, &b_gen_Wl_mass);
    fTree->SetBranchAddress("gen_X_mass", &gen_X_mass, &b_gen_X_mass);
    fTree->SetBranchAddress("gen_lep_pt", &gen_lep_pt, &b_gen_lep_pt);
    fTree->SetBranchAddress("gen_lep_eta", &gen_lep_eta, &b_gen_lep_eta);
    fTree->SetBranchAddress("gen_lep_phi", &gen_lep_phi, &b_gen_lep_phi);
    fTree->SetBranchAddress("gen_nu_pt", &gen_nu_pt, &b_gen_nu_pt);
    fTree->SetBranchAddress("gen_nu_eta", &gen_nu_eta, &b_gen_nu_eta);
    fTree->SetBranchAddress("gen_nu_phi", &gen_nu_phi, &b_gen_nu_phi);
    fTree->SetBranchAddress("gen_W_q1_pt", &gen_W_q1_pt, &b_gen_W_q1_pt);
    fTree->SetBranchAddress("gen_W_q1_eta", &gen_W_q1_eta, &b_gen_W_q1_eta);
    fTree->SetBranchAddress("gen_W_q1_phi", &gen_W_q1_phi, &b_gen_W_q1_phi);
    fTree->SetBranchAddress("gen_W_q2_pt", &gen_W_q2_pt, &b_gen_W_q2_pt);
    fTree->SetBranchAddress("gen_W_q2_eta", &gen_W_q2_eta, &b_gen_W_q2_eta);
    fTree->SetBranchAddress("gen_W_q2_phi", &gen_W_q2_phi, &b_gen_W_q2_phi);
    fTree->SetBranchAddress("gen_vbf_q1_pt", &gen_vbf_q1_pt, &b_gen_vbf_q1_pt);
    fTree->SetBranchAddress("gen_vbf_q1_eta", &gen_vbf_q1_eta, &b_gen_vbf_q1_eta);
    fTree->SetBranchAddress("gen_vbf_q1_phi", &gen_vbf_q1_phi, &b_gen_vbf_q1_phi);
    fTree->SetBranchAddress("gen_vbf_q1_btag", &gen_vbf_q1_btag, &b_gen_vbf_q1_btag);
    fTree->SetBranchAddress("gen_vbf_q2_pt", &gen_vbf_q2_pt, &b_gen_vbf_q2_pt);
    fTree->SetBranchAddress("gen_vbf_q2_eta", &gen_vbf_q2_eta, &b_gen_vbf_q2_eta);
    fTree->SetBranchAddress("gen_vbf_q2_phi", &gen_vbf_q2_phi, &b_gen_vbf_q2_phi);
    fTree->SetBranchAddress("gen_vbf_q2_btag", &gen_vbf_q2_btag, &b_gen_vbf_q2_btag);
    fTree->SetBranchAddress("gen_vbf_qq_mass", &gen_vbf_qq_mass, &b_gen_vbf_qq_mass);
    fTree->SetBranchAddress("gen_vbf_qq_delta_phi", &gen_vbf_qq_delta_phi, &b_gen_vbf_qq_delta_phi);
    fTree->SetBranchAddress("gen_vbf_qq_delta_R", &gen_vbf_qq_delta_R, &b_gen_vbf_qq_delta_R);
    //---reco    
    fTree->SetBranchAddress("lep_pt", &lep_pt, &b_lep_pt);
    fTree->SetBranchAddress("lep_eta", &lep_eta, &b_lep_eta);
    fTree->SetBranchAddress("lep_phi", &lep_phi, &b_lep_phi);
    fTree->SetBranchAddress("MET", &MET, &b_MET);
    fTree->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
    fTree->SetBranchAddress("lv_mass", &lv_mass, &b_lv_mass);
    fTree->SetBranchAddress("lv_pt", &lv_pt, &b_lv_mass);
    fTree->SetBranchAddress("lv_eta", &lv_eta, &b_lv_eta);
    fTree->SetBranchAddress("lv_phi", &lv_phi, &b_lv_phi);
    fTree->SetBranchAddress("lv_deltaR", &lv_delta_R, &b_lv_delta_R);
    fTree->SetBranchAddress("lv_closerjet_mass", &lv_closerjet_mass, &b_lv_closerjet_mass);
    fTree->SetBranchAddress("CA8_jet_pt", &CA8_jet_pt, &b_CA8_jet_pt);
    fTree->SetBranchAddress("CA8_jet_eta", &CA8_jet_eta, &b_CA8_jet_eta);
    fTree->SetBranchAddress("CA8_jet_phi", &CA8_jet_phi, &b_CA8_jet_phi);
    fTree->SetBranchAddress("CA8_jet_mass", &CA8_jet_mass, &b_CA8_jet_mass);
    fTree->SetBranchAddress("CA8_jet_t2t1", &CA8_jet_t2t1, &b_CA8_jet_t2t1);
    fTree->SetBranchAddress("CA8_jet_t3t2", &CA8_jet_t3t2, &b_CA8_jet_t3t2);
    fTree->SetBranchAddress("CA8_closerjet_mass", &CA8_closerjet_mass, &b_CA8_closerjet_mass);
    fTree->SetBranchAddress("vbf_jet1_pt", &vbf_jet1_pt, &b_vbf_jet1_pt);
    fTree->SetBranchAddress("vbf_jet1_eta", &vbf_jet1_eta, &b_vbf_jet1_eta);
    fTree->SetBranchAddress("vbf_jet1_phi", &vbf_jet1_phi, &b_vbf_jet1_phi);
    fTree->SetBranchAddress("vbf_jet1_mass", &vbf_jet1_mass, &b_vbf_jet1_mass);
    fTree->SetBranchAddress("vbf_jet1_btag", &vbf_jet1_btag, &b_vbf_jet1_btag);
    fTree->SetBranchAddress("vbf_jet2_pt", &vbf_jet2_pt, &b_vbf_jet2_pt);
    fTree->SetBranchAddress("vbf_jet2_eta", &vbf_jet2_eta, &b_vbf_jet2_eta);
    fTree->SetBranchAddress("vbf_jet2_phi", &vbf_jet2_phi, &b_vbf_jet2_phi);
    fTree->SetBranchAddress("vbf_jet2_mass", &vbf_jet2_mass, &b_vbf_jet2_mass);
    fTree->SetBranchAddress("vbf_jet2_btag", &vbf_jet2_btag, &b_vbf_jet2_btag);
    fTree->SetBranchAddress("vbf_jj_mass", &vbf_jj_mass, &b_vbf_jj_mass);
    fTree->SetBranchAddress("vbf_jj_delta_phi", &vbf_jj_delta_phi, &b_vbf_jj_delta_phi);
    fTree->SetBranchAddress("vbf_jj_delta_R", &vbf_jj_delta_R, &b_vbf_jj_delta_R);
    fTree->SetBranchAddress("lvJ_mass", &lvJ_mass, &b_lvJ_mass);
}

