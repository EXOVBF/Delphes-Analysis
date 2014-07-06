#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

using namespace std ;

TTree *fChain ;

//-----Declaration of leaf types
//---gen---
//---Ws 
float gen_pt_W_had;
float gen_eta_W_had;
float gen_phi_W_had;
float gen_mass_W_had;
float gen_pt_W_lep;
float gen_eta_W_lep;
float gen_phi_W_lep;
float gen_mass_W_lep;
float gen_mass_X;
//---lep
float gen_lep_pt;
float gen_lep_eta;
float gen_lep_phi;
float gen_nu_pt;
float gen_nu_eta;
float gen_nu_phi;
//---W's quark
float gen_W_q1_pt;
float gen_W_q1_eta;
float gen_W_q1_phi;
float gen_W_q2_pt;
float gen_W_q2_eta;
float gen_W_q2_phi;
//---vbf quark
float gen_vbf_q1_pt;
float gen_vbf_q1_eta;
float gen_vbf_q1_phi;
float gen_vbf_q2_pt;
float gen_vbf_q2_eta;
float gen_vbf_q2_phi;
float gen_vbf_qq_mass;
//---reco----
//---lep
float lep_pt;
float lep_eta;
float lep_phi;
float MET;
float MET_phi;
//---leptonic W
float lv_mass;
float lv_pt;
float lv_eta;
float lv_phi;
float lv_delta_R;
float lv_closerjet_mass;
//---CA8 jet
float CA8_jet_pt;
float CA8_jet_eta;
float CA8_jet_phi;
float CA8_jet_mass;
float CA8_jet_t2t1;
float CA8_jet_t3t2;
float CA8_closerjet_mass;
//---vbf jets
float vbf_jet1_pt;
float vbf_jet1_eta;
float vbf_jet1_phi;
float vbf_jet1_mass;
float vbf_jet2_pt;
float vbf_jet2_eta;
float vbf_jet2_phi;
float vbf_jet2_mass;
float vbf_jj_mass;
float vbf_jj_delta_R;
//---lvj
float lvJ_mass; 

//****************************************************************************************

void InitEasyTree (TTree * fChain)
{

    //-----Set object pointer
    //---gen
    gen_pt_W_had=0;
    gen_eta_W_had=0;
    gen_phi_W_had=0;
    gen_mass_W_had=0;
    gen_pt_W_lep=0;
    gen_eta_W_lep=0;
    gen_phi_W_lep=0;
    gen_mass_W_lep=0;
    gen_mass_X=0;
    gen_lep_pt=0;
    gen_lep_eta=0;
    gen_lep_phi=0;
    gen_nu_pt=0;
    gen_nu_eta=0;
    gen_nu_phi=0;
    gen_W_q1_pt=0;
    gen_W_q1_eta=0;
    gen_W_q1_phi=0;
    gen_W_q2_pt=0;
    gen_W_q2_eta=0;
    gen_W_q2_phi=0;
    gen_vbf_q1_pt=0;
    gen_vbf_q1_eta=0;
    gen_vbf_q1_phi=0;
    gen_vbf_q2_pt=0;
    gen_vbf_q2_eta=0;
    gen_vbf_q2_phi=0;
    gen_vbf_qq_mass=0;
    //---reco
    lep_pt=0;
    lep_eta=0;
    lep_phi=0;
    MET=0;
    MET_phi=0;
    lv_mass=0;
    lv_pt=0;
    lv_eta=0;
    lv_phi=0;
    lv_delta_R=0;
    lv_closerjet_mass=0;
    CA8_jet_pt=0;
    CA8_jet_eta=0;
    CA8_jet_phi=0;
    CA8_jet_mass=0;
    CA8_jet_t2t1=0;
    CA8_jet_t3t2=0;
    CA8_closerjet_mass=0;
    vbf_jet1_pt=0;
    vbf_jet1_eta=0;
    vbf_jet1_phi=0;
    vbf_jet1_mass=0;
    vbf_jet2_pt=0;
    vbf_jet2_eta=0;
    vbf_jet2_phi=0;
    vbf_jet2_mass=0;
    vbf_jj_mass=0;
    vbf_jj_delta_R=0;
    lvJ_mass=0; 

    // Set branch addresses and branch pointers
    fChain->SetMakeClass (1) ;

    //---gen
    fChain->Branch("gen_pt_W_had", &gen_pt_W_had, "gen_pt_W_had/F");
    fChain->Branch("gen_eta_W_had", &gen_eta_W_had, "gen_eta_W_had/F");
    fChain->Branch("gen_phi_W_had", &gen_phi_W_had, "gen_phi_W_had/F");
    fChain->Branch("gen_mass_W_had", &gen_mass_W_had, "gen_mass_W_had/F");
    fChain->Branch("gen_pt_W_lep", &gen_pt_W_lep, "gen_pt_W_lep/F");
    fChain->Branch("gen_eta_W_lep", &gen_eta_W_lep, "gen_eta_W_lep/F");
    fChain->Branch("gen_phi_W_lep", &gen_phi_W_lep, "gen_phi_W_lep/F");
    fChain->Branch("gen_mass_W_lep", &gen_mass_W_lep, "gen_mass_W_lep/F");
    fChain->Branch("gen_mass_X", &gen_mass_X, "gen_mass_X/F");
    fChain->Branch("gen_lep_pt", &gen_lep_pt, "gen_lep_pt/F");
    fChain->Branch("gen_lep_eta", &gen_lep_eta, "gen_lep_eta/F");
    fChain->Branch("gen_lep_phi", &gen_lep_phi, "gen_lep_phi/F"); 
    fChain->Branch("gen_nu_pt", &gen_nu_pt, "gen_nu_pt/F");
    fChain->Branch("gen_nu_eta", &gen_nu_eta, "gen_nu_eta/F");
    fChain->Branch("gen_nu_phi", &gen_nu_phi, "gen_nu_phi/F"); 
    fChain->Branch("gen_W_q1_pt", &gen_W_q1_pt, "gen_W_q1_pt/F");
    fChain->Branch("gen_W_q1_eta", &gen_W_q1_eta, "gen_W_q1_eta/F");
    fChain->Branch("gen_W_q1_phi", &gen_W_q1_phi, "gen_W_q1_phi/F"); 
    fChain->Branch("gen_W_q2_pt", &gen_W_q2_pt, "gen_W_q2_pt/F");
    fChain->Branch("gen_W_q2_eta", &gen_W_q2_eta, "gen_W_q2_eta/F");
    fChain->Branch("gen_W_q2_phi", &gen_W_q2_phi, "gen_W_q2_phi/F"); 
    fChain->Branch("gen_vbf_q1_pt", &gen_vbf_q1_pt, "gen_vbf_q1_pt/F");
    fChain->Branch("gen_vbf_q1_eta", &gen_vbf_q1_eta, "gen_vbf_q1_eta/F");
    fChain->Branch("gen_vbf_q1_phi", &gen_vbf_q1_phi, "gen_vbf_q1_phi/F"); 
    fChain->Branch("gen_vbf_q2_pt", &gen_vbf_q2_pt, "gen_vbf_q2_pt/F");
    fChain->Branch("gen_vbf_q2_eta", &gen_vbf_q2_eta, "gen_vbf_q2_eta/F");
    fChain->Branch("gen_vbf_q2_phi", &gen_vbf_q2_phi, "gen_vbf_q2_phi/F"); 
    fChain->Branch("gen_vbf_qq_mass", &gen_vbf_qq_mass, "gen_vbf_qq_mass/F"); 
    //---reco
    fChain->Branch("lep_pt", &lep_pt, "lep_pt/F");
    fChain->Branch("lep_eta", &lep_eta, "lep_eta/F");
    fChain->Branch("lep_phi", &lep_phi, "lep_phi/F");
    fChain->Branch("MET", &MET, "MET/F");
    fChain->Branch("MET_phi", &MET_phi, "MET_phi/F");
    fChain->Branch("lv_mass", &lv_mass, "lv_mass/F");
    fChain->Branch("lv_pt", &lv_pt, "lv_pt/F");
    fChain->Branch("lv_eta", &lv_eta, "lv_eta/F");
    fChain->Branch("lv_phi", &lv_phi, "lv_phi/F");
    fChain->Branch("lv_deltaR", &lv_delta_R, "lv_delta_R/F");
    fChain->Branch("lv_closerjet_mass", &lv_closerjet_mass, "lv_closerjet_mass/F");
    fChain->Branch("CA8_jet_pt", &CA8_jet_pt, "CA8_jet_pt/F");
    fChain->Branch("CA8_jet_eta", &CA8_jet_eta, "CA8_jet_eta/F");
    fChain->Branch("CA8_jet_phi", &CA8_jet_phi, "CA8_jet_phi/F");
    fChain->Branch("CA8_jet_mass", &CA8_jet_mass, "CA8_jet_mass/F");
    fChain->Branch("CA8_jet_t2t1", &CA8_jet_t2t1, "CA8_jet_t2t1/F");
    fChain->Branch("CA8_jet_t3t2", &CA8_jet_t3t2, "CA8_jet_t3t2/F");
    fChain->Branch("CA8_closerjet_mass", &CA8_closerjet_mass, "CA8_closerjet_mass/F");
    fChain->Branch("vbf_jet1_pt", &vbf_jet1_pt, "vbf_jet1_pt/F");
    fChain->Branch("vbf_jet1_eta", &vbf_jet1_eta, "vbf_jet1_eta/F");
    fChain->Branch("vbf_jet1_phi", &vbf_jet1_phi, "vbf_jet1_phi/F");
    fChain->Branch("vbf_jet1_mass", &vbf_jet1_mass, "vbf_jet1_mass/F");
    fChain->Branch("vbf_jet2_pt", &vbf_jet2_pt, "vbf_jet2_pt/F");
    fChain->Branch("vbf_jet2_eta", &vbf_jet2_eta, "vbf_jet2_eta/F");
    fChain->Branch("vbf_jet2_phi", &vbf_jet2_phi, "vbf_jet2_phi/F");
    fChain->Branch("vbf_jet2_mass", &vbf_jet2_mass, "vbf_jet2_mass/F");
    fChain->Branch("vbf_jj_mass", &vbf_jj_mass, "vbf_jj_mass/F");
    fChain->Branch("vbf_jj_delta_R", &vbf_jj_delta_R, "vbf_jj_delta_R/F");
    fChain->Branch("lvJ_mass", &lvJ_mass, "lvJ_mass/F"); 
}

