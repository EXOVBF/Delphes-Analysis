#ifndef data_tree_h
#define data_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class data_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Float_t         nPV;
   Int_t           event;
   Int_t           event_lumi;
   Int_t           event_runNo;
   Int_t           issignal;
   Int_t           numberJetBin;
   Int_t           numberJetBin2;
   Int_t           numberJetBin3;
   Int_t           numberJetBin4;
   Int_t           numberJetBinGen;
   Int_t           numberJetBinGen2;
   Int_t           numberJetBinGen3;
   Int_t           numberJetBinGen4;
   Float_t         totalEventWeight;
   Float_t         eff_and_pu_Weight;
   Float_t         wSampleWeight;
   Float_t         event_weight;
   Float_t         btag_weight;
   Float_t         btag_weight_up;
   Float_t         btag_weight_dn;
   Float_t         btag_weight_dn_up;
   Float_t         btag_weight_up_dn;
   Float_t         interference_Weight_H600;
   Float_t         interference_Weight_H700;
   Float_t         interference_Weight_H800;
   Float_t         interference_Weight_H900;
   Float_t         interference_Weight_H1000;
   Float_t         cps_Weight_H600;
   Float_t         cps_Weight_H700;
   Float_t         cps_Weight_H800;
   Float_t         cps_Weight_H900;
   Float_t         cps_Weight_H1000;
   Float_t         bsmReweight_cPrime01_brNew00;
   Float_t         bsmReweight_cPrime01_brNew01;
   Float_t         bsmReweight_cPrime01_brNew02;
   Float_t         bsmReweight_cPrime01_brNew03;
   Float_t         bsmReweight_cPrime01_brNew04;
   Float_t         bsmReweight_cPrime01_brNew05;
   Float_t         bsmReweight_cPrime02_brNew00;
   Float_t         bsmReweight_cPrime02_brNew01;
   Float_t         bsmReweight_cPrime02_brNew02;
   Float_t         bsmReweight_cPrime02_brNew03;
   Float_t         bsmReweight_cPrime02_brNew04;
   Float_t         bsmReweight_cPrime02_brNew05;
   Float_t         bsmReweight_cPrime03_brNew00;
   Float_t         bsmReweight_cPrime03_brNew01;
   Float_t         bsmReweight_cPrime03_brNew02;
   Float_t         bsmReweight_cPrime03_brNew03;
   Float_t         bsmReweight_cPrime03_brNew04;
   Float_t         bsmReweight_cPrime03_brNew05;
   Float_t         bsmReweight_cPrime04_brNew00;
   Float_t         bsmReweight_cPrime04_brNew01;
   Float_t         bsmReweight_cPrime04_brNew02;
   Float_t         bsmReweight_cPrime04_brNew03;
   Float_t         bsmReweight_cPrime04_brNew04;
   Float_t         bsmReweight_cPrime04_brNew05;
   Float_t         bsmReweight_cPrime05_brNew00;
   Float_t         bsmReweight_cPrime05_brNew01;
   Float_t         bsmReweight_cPrime05_brNew02;
   Float_t         bsmReweight_cPrime05_brNew03;
   Float_t         bsmReweight_cPrime05_brNew04;
   Float_t         bsmReweight_cPrime05_brNew05;
   Float_t         bsmReweight_cPrime06_brNew00;
   Float_t         bsmReweight_cPrime06_brNew01;
   Float_t         bsmReweight_cPrime06_brNew02;
   Float_t         bsmReweight_cPrime06_brNew03;
   Float_t         bsmReweight_cPrime06_brNew04;
   Float_t         bsmReweight_cPrime06_brNew05;
   Float_t         bsmReweight_cPrime07_brNew00;
   Float_t         bsmReweight_cPrime07_brNew01;
   Float_t         bsmReweight_cPrime07_brNew02;
   Float_t         bsmReweight_cPrime07_brNew03;
   Float_t         bsmReweight_cPrime07_brNew04;
   Float_t         bsmReweight_cPrime07_brNew05;
   Float_t         bsmReweight_cPrime08_brNew00;
   Float_t         bsmReweight_cPrime08_brNew01;
   Float_t         bsmReweight_cPrime08_brNew02;
   Float_t         bsmReweight_cPrime08_brNew03;
   Float_t         bsmReweight_cPrime08_brNew04;
   Float_t         bsmReweight_cPrime08_brNew05;
   Float_t         bsmReweight_cPrime09_brNew00;
   Float_t         bsmReweight_cPrime09_brNew01;
   Float_t         bsmReweight_cPrime09_brNew02;
   Float_t         bsmReweight_cPrime09_brNew03;
   Float_t         bsmReweight_cPrime09_brNew04;
   Float_t         bsmReweight_cPrime09_brNew05;
   Float_t         bsmReweight_cPrime10_brNew00;
   Float_t         bsmReweight_cPrime10_brNew01;
   Float_t         bsmReweight_cPrime10_brNew02;
   Float_t         bsmReweight_cPrime10_brNew03;
   Float_t         bsmReweight_cPrime10_brNew04;
   Float_t         bsmReweight_cPrime10_brNew05;
   Float_t         mass_lvj_type0_met;
   Float_t         mass_lvj_type2_met;
   Float_t         mass_lvj_type0;
   Float_t         mass_lvj_type2;
   Float_t         mass_lvj_type0_met_jes_up;
   Float_t         mass_lvj_type2_met_jes_up;
   Float_t         mass_lvj_type0_met_jes_dn;
   Float_t         mass_lvj_type2_met_jes_dn;
   Float_t         mass_lvj_type0_met_jer;
   Float_t         mass_lvj_type2_met_jer;
   Float_t         mass_lvj_type0_met_jer_up;
   Float_t         mass_lvj_type2_met_jer_up;
   Float_t         mass_lvj_type0_met_jer_dn;
   Float_t         mass_lvj_type2_met_jer_dn;
   Float_t         mass_lvj_type0_met_lep_scale_up;
   Float_t         mass_lvj_type2_met_lep_scale_up;
   Float_t         mass_lvj_type0_met_lep_scale_dn;
   Float_t         mass_lvj_type2_met_lep_scale_dn;
   Float_t         mass_lvj_type0_met_lep_res;
   Float_t         mass_lvj_type2_met_lep_res;
   Float_t         mass_lv_subj_type0_met;
   Float_t         mass_lv_subj_type2_met;
   Float_t         mass_lv_subj_type0;
   Float_t         mass_lv_subj_type2;
   Float_t         mass_ungroomedjet_vbf_j1;
   Float_t         mass_ungroomedjet_vbf_j2;
   Float_t         mass_ungroomedjet_vbf_j1_pr;
   Float_t         mass_ungroomedjet_vbf_j2_pr;
   Float_t         mass_ungroomedjet_closerjet;
   Float_t         mass_ungroomedjet_closerjet_pr;
   Float_t         mass_leptonic_closerjet;
   Float_t         l_pt;
   Float_t         l_eta;
   Float_t         l_charge;
   Float_t         l_phi;
   Float_t         l_pt_scale_up;
   Float_t         l_eta_scale_up;
   Float_t         l_phi_scale_up;
   Float_t         l_pt_scale_dn;
   Float_t         l_eta_scale_dn;
   Float_t         l_phi_scale_dn;
   Float_t         l_pt_res;
   Float_t         l_eta_res;
   Float_t         l_phi_res;
   Float_t         pfMET;
   Float_t         pfMET_Phi;
   Float_t         pfMET_jes_up;
   Float_t         pfMET_Phi_jes_up;
   Float_t         pfMET_jes_dn;
   Float_t         pfMET_Phi_jes_dn;
   Float_t         pfMET_jer_;
   Float_t         pfMET_Phi_jer_;
   Float_t         pfMET_jer_up;
   Float_t         pfMET_Phi_jer_up;
   Float_t         pfMET_jer_dn;
   Float_t         pfMET_Phi_jer_dn;
   Float_t         pfMET_lep_scale_up;
   Float_t         pfMET_Phi_lep_scale_up;
   Float_t         pfMET_lep_scale_dn;
   Float_t         pfMET_Phi_lep_scale_dn;
   Float_t         pfMET_lep_res;
   Float_t         pfMET_Phi_lep_res;
   Float_t         nu_pz_type0;
   Float_t         nu_pz_type2;
   Float_t         nu_pz_type0_met;
   Float_t         nu_pz_type2_met;
   Float_t         W_pz_type0;
   Float_t         W_pz_type2;
   Float_t         W_pz_type0_met;
   Float_t         W_pz_type2_met;
   Float_t         nu_pz_gen;
   Float_t         W_pz_gen;
   Float_t         W_pt_gen;
   Float_t         v_pt;
   Float_t         v_mt;
   Float_t         v_eta;
   Float_t         v_phi;
   Float_t         ungroomed_jet_eta;
   Float_t         ungroomed_jet_phi;
   Float_t         ungroomed_jet_pt;
   Float_t         ungroomed_jet_e;
   Float_t         ungroomed_gen_jet_eta;
   Float_t         ungroomed_gen_jet_phi;
   Float_t         ungroomed_gen_jet_pt;
   Float_t         ungroomed_gen_jet_e;
   Float_t         jet_mass_pr;
   Float_t         jet_pt_pr;
   Float_t         jet_charge;
   Float_t         jet_charge_k05;
   Float_t         jet_charge_k07;
   Float_t         jet_charge_k10;
   Float_t         gen_jet_mass_pr;
   Float_t         gen_jet_pt_pr;
   Float_t         jet_grsens_ft;
   Float_t         jet_grsens_tr;
   Float_t         jet_massdrop_pr;
   Float_t         jet_qjetvol;
   Float_t         gen_jet_grsens_ft;
   Float_t         gen_jet_grsens_tr;
   Float_t         gen_jet_massdrop_pr;
   Float_t         gen_jet_qjetvol;
   Float_t         jet_tau2tau1;
   Float_t         jet_tau2tau1_exkT;
   Float_t         jet_tau2tau1_pr;
   Float_t         jet_GeneralizedECF;
   Float_t         gen_jet_tau2tau1;
   Float_t         gen_jet_tau2tau1_exkT;
   Float_t         gen_jet_tau2tau1_pr;
   Float_t         gen_jet_GeneralizedECF;
   Float_t         jet_jetconstituents;
   Float_t         gen_jet_jetconstituents;
   Float_t         jet_rcore4;
   Float_t         jet_rcore5;
   Float_t         jet_rcore6;
   Float_t         jet_rcore7;
   Float_t         gen_jet_rcore4;
   Float_t         gen_jet_rcore5;
   Float_t         gen_jet_rcore6;
   Float_t         gen_jet_rcore7;
   Float_t         jet_pt1frac;
   Float_t         jet_pt2frac;
   Float_t         jet_sjdr;
   Float_t         j_jecfactor_up;
   Float_t         j_jecfactor_dn;
   Float_t         jet_mass_pr_jes_up;
   Float_t         jet_mass_pr_jes_dn;
   Float_t         jet_mass_pr_jer;
   Float_t         jet_mass_pr_jer_up;
   Float_t         jet_mass_pr_jer_dn;
   Float_t         ungroomed_jet_pt_jes_dn;
   Float_t         ungroomed_jet_pt_jes_up;
   Float_t         ungroomed_jet_pt_jer;
   Float_t         ungroomed_jet_pt_jer_dn;
   Float_t         ungroomed_jet_pt_jer_up;
   Float_t         jet_planarlow04;
   Float_t         jet_planarlow05;
   Float_t         jet_planarlow06;
   Float_t         jet_planarlow07;
   Float_t         vbf_maxpt_jj_m;
   Float_t         vbf_maxpt_jj_pt;
   Float_t         vbf_maxpt_jj_eta;
   Float_t         vbf_maxpt_jj_phi;
   Float_t         vbf_maxpt_j1_m;
   Float_t         vbf_maxpt_j1_pt;
   Float_t         vbf_maxpt_j1_eta;
   Float_t         vbf_maxpt_j1_phi;
   Float_t         vbf_maxpt_j2_m;
   Float_t         vbf_maxpt_j2_pt;
   Float_t         vbf_maxpt_j2_eta;
   Float_t         vbf_maxpt_j2_phi;
   Float_t         vbf_maxpt_jj_m_gen;
   Float_t         vbf_maxpt_jj_pt_gen;
   Float_t         vbf_maxpt_jj_eta_gen;
   Float_t         vbf_maxpt_jj_phi_gen;
   Float_t         vbf_maxpt_j1_m_gen;
   Float_t         vbf_maxpt_j1_pt_gen;
   Float_t         vbf_maxpt_j1_eta_gen;
   Float_t         vbf_maxpt_j1_phi_gen;
   Float_t         vbf_maxpt_j2_m_gen;
   Float_t         vbf_maxpt_j2_pt_gen;
   Float_t         vbf_maxpt_j2_eta_gen;
   Float_t         vbf_maxpt_j2_phi_gen;
   Float_t         vbf_maxpt_j1_m_jes_up;
   Float_t         vbf_maxpt_j1_pt_jes_up;
   Float_t         vbf_maxpt_j1_eta_jes_up;
   Float_t         vbf_maxpt_j1_phi_jes_up;
   Float_t         vbf_maxpt_j1_m_jes_dn;
   Float_t         vbf_maxpt_j1_pt_jes_dn;
   Float_t         vbf_maxpt_j1_eta_jes_dn;
   Float_t         vbf_maxpt_j1_phi_jes_dn;
   Float_t         vbf_maxpt_j1_m_jer;
   Float_t         vbf_maxpt_j1_pt_jer;
   Float_t         vbf_maxpt_j1_eta_jer;
   Float_t         vbf_maxpt_j1_phi_jer;
   Float_t         vbf_maxpt_j1_m_jer_up;
   Float_t         vbf_maxpt_j1_pt_jer_up;
   Float_t         vbf_maxpt_j1_eta_jer_up;
   Float_t         vbf_maxpt_j1_phi_jer_up;
   Float_t         vbf_maxpt_j1_m_jer_dn;
   Float_t         vbf_maxpt_j1_pt_jer_dn;
   Float_t         vbf_maxpt_j1_eta_jer_dn;
   Float_t         vbf_maxpt_j1_phi_jer_dn;
   Float_t         vbf_maxpt_j2_m_jes_up;
   Float_t         vbf_maxpt_j2_pt_jes_up;
   Float_t         vbf_maxpt_j2_eta_jes_up;
   Float_t         vbf_maxpt_j2_phi_jes_up;
   Float_t         vbf_maxpt_j2_m_jes_dn;
   Float_t         vbf_maxpt_j2_pt_jes_dn;
   Float_t         vbf_maxpt_j2_eta_jes_dn;
   Float_t         vbf_maxpt_j2_phi_jes_dn;
   Float_t         vbf_maxpt_j2_m_jer;
   Float_t         vbf_maxpt_j2_pt_jer;
   Float_t         vbf_maxpt_j2_eta_jer;
   Float_t         vbf_maxpt_j2_phi_jer;
   Float_t         vbf_maxpt_j2_m_jer_up;
   Float_t         vbf_maxpt_j2_pt_jer_up;
   Float_t         vbf_maxpt_j2_eta_jer_up;
   Float_t         vbf_maxpt_j2_phi_jer_up;
   Float_t         vbf_maxpt_j2_m_jer_dn;
   Float_t         vbf_maxpt_j2_pt_jer_dn;
   Float_t         vbf_maxpt_j2_eta_jer_dn;
   Float_t         vbf_maxpt_j2_phi_jer_dn;
   Float_t         vbf_maxpt_j1_QGLikelihood;
   Float_t         vbf_maxpt_j2_QGLikelihood;
   Int_t           vbf_maxpt_j1_isPileUpMedium;
   Int_t           vbf_maxpt_j2_isPileUpMedium;
   Int_t           vbf_maxpt_j1_isPileUpTight;
   Int_t           vbf_maxpt_j2_isPileUpTight;
   Float_t         vbf_maxpt_j1_bDiscriminatorCSV;
   Float_t         vbf_maxpt_j2_bDiscriminatorCSV;
   Float_t         vbf_maxpt_j1_bDiscriminatorCSV_gen;
   Float_t         vbf_maxpt_j2_bDiscriminatorCSV_gen;
   Float_t         nbjets_csvl_veto;
   Float_t         nbjets_csvm_veto;
   Float_t         nbjets_csvt_veto;
   Float_t         nbjets_ssvhem_veto;
   Float_t         nbjets_csvl_veto_cleaned;
   Float_t         nbjets_csvm_veto_cleaned;
   Float_t         nbjets_csvt_veto_cleaned;
   Float_t         nbjets_ssvhem_veto_cleaned;
   Float_t         njets;
   Float_t         deltaR_lca8jet;
   Float_t         deltaphi_METca8jet;
   Float_t         deltaphi_Vca8jet;
   Float_t         deltaphi_METca8jet_met;
   Float_t         deltaphi_Vca8jet_met;
   Float_t         genHMass;
   Float_t         genHphi;
   Float_t         genHeta;
   Float_t         genHpt;
   Float_t         genTagQuark1W;
   Float_t         genTagQuark1phi;
   Float_t         genTagQuark1eta;
   Float_t         genTagQuark1pt;
   Float_t         genTagQuarkE;
   Float_t         genTagQuark2phi;
   Float_t         genTagQuark2eta;
   Float_t         genTagQuark2pt;
   Float_t         ttb_nak5_same;
   Float_t         ttb_nak5_same_csvl;
   Float_t         ttb_nak5_same_csvm;
   Float_t         ttb_nak5_same_csvt;
   Float_t         ttb_nak5_oppo;
   Float_t         ttb_nak5_oppo_csvl;
   Float_t         ttb_nak5_oppo_csvm;
   Float_t         ttb_nak5_oppo_csvt;
   Float_t         ttb_nak5_oppoveto;
   Float_t         ttb_nak5_oppoveto_csvl;
   Float_t         ttb_nak5_oppoveto_csvm;
   Float_t         ttb_nak5_oppoveto_csvt;
   Float_t         ttb_nak5_sameveto;
   Float_t         ttb_nak5_sameveto_csvl;
   Float_t         ttb_nak5_sameveto_csvm;
   Float_t         ttb_nak5_sameveto_csvt;
   Float_t         ttb_ht;
   Float_t         ttb_ca8_mass_pr;
   Float_t         ttb_ca8_charge;
   Float_t         ttb_ca8_charge_k05;
   Float_t         ttb_ca8_charge_k07;
   Float_t         ttb_ca8_charge_k10;
   Float_t         ttb_ca8_ungroomed_pt;
   Float_t         ttb_ca8_ungroomed_eta;
   Float_t         ttb_ca8_ungroomed_phi;
   Float_t         ttb_ca8_ungroomed_e;
   Float_t         ttb_ca8_ungroomed_gen_pt;
   Float_t         ttb_ca8_ungroomed_gen_eta;
   Float_t         ttb_ca8_ungroomed_gen_phi;
   Float_t         ttb_ca8_ungroomed_gen_e;
   Float_t         ttb_ca8_tau2tau1;
   Float_t         ttb_ca8_tau2tau1_exkT;
   Float_t         ttb_ca8_tau2tau1_pr;
   Float_t         ttb_ca8_GeneralizedECF;
   Float_t         ttb_ca8_mu;
   Float_t         ttb_ca8_mlvj_type0;
   Float_t         ttb_ca8_mlvj_type2;
   Float_t         ttb_ca8_mlvj_type0_met;
   Float_t         ttb_ca8_mlvj_type2_met;
   Float_t         ttb_dR_ca8_bjet_closer;
   Float_t         ttb_dR_ca8_jet_closer;
   Int_t           isttbar;
   Float_t         gen_parton1_px_fromttbar;
   Float_t         gen_parton1_py_fromttbar;
   Float_t         gen_parton1_pz_fromttbar;
   Float_t         gen_parton1_e_fromttbar;
   Float_t         gen_parton1_id_fromttbar;
   Float_t         gen_parton2_px_fromttbar;
   Float_t         gen_parton2_py_fromttbar;
   Float_t         gen_parton2_pz_fromttbar;
   Float_t         gen_parton2_e_fromttbar;
   Float_t         gen_parton2_id_fromttbar;

   // List of branches
   TBranch        *b_nPV;   //!
   TBranch        *b_event;   //!
   TBranch        *b_event_lumi;   //!
   TBranch        *b_event_runNo;   //!
   TBranch        *b_issignal;   //!
   TBranch        *b_numberJetBin;   //!
   TBranch        *b_numberJetBin2;   //!
   TBranch        *b_numberJetBin3;   //!
   TBranch        *b_numberJetBin4;   //!
   TBranch        *b_numberJetBinGen;   //!
   TBranch        *b_numberJetBinGen2;   //!
   TBranch        *b_numberJetBinGen3;   //!
   TBranch        *b_numberJetBinGen4;   //!
   TBranch        *b_totalEventWeight;   //!
   TBranch        *b_eff_and_pu_Weight;   //!
   TBranch        *b_wSampleWeight;   //!
   TBranch        *b_event_weight;   //!
   TBranch        *b_btag_weight;   //!
   TBranch        *b_btag_weight_up;   //!
   TBranch        *b_btag_weight_dn;   //!
   TBranch        *b_btag_weight_dn_up;   //!
   TBranch        *b_btag_weight_up_dn;   //!
   TBranch        *b_interference_Weight_H600;   //!
   TBranch        *b_interference_Weight_H700;   //!
   TBranch        *b_interference_Weight_H800;   //!
   TBranch        *b_interference_Weight_H900;   //!
   TBranch        *b_interference_Weight_H1000;   //!
   TBranch        *b_cps_Weight_H600;   //!
   TBranch        *b_cps_Weight_H700;   //!
   TBranch        *b_cps_Weight_H800;   //!
   TBranch        *b_cps_Weight_H900;   //!
   TBranch        *b_cps_Weight_H1000;   //!
   TBranch        *b_bsmReweight_cPrime01_brNew00;   //!
   TBranch        *b_bsmReweight_cPrime01_brNew01;   //!
   TBranch        *b_bsmReweight_cPrime01_brNew02;   //!
   TBranch        *b_bsmReweight_cPrime01_brNew03;   //!
   TBranch        *b_bsmReweight_cPrime01_brNew04;   //!
   TBranch        *b_bsmReweight_cPrime01_brNew05;   //!
   TBranch        *b_bsmReweight_cPrime02_brNew00;   //!
   TBranch        *b_bsmReweight_cPrime02_brNew01;   //!
   TBranch        *b_bsmReweight_cPrime02_brNew02;   //!
   TBranch        *b_bsmReweight_cPrime02_brNew03;   //!
   TBranch        *b_bsmReweight_cPrime02_brNew04;   //!
   TBranch        *b_bsmReweight_cPrime02_brNew05;   //!
   TBranch        *b_bsmReweight_cPrime03_brNew00;   //!
   TBranch        *b_bsmReweight_cPrime03_brNew01;   //!
   TBranch        *b_bsmReweight_cPrime03_brNew02;   //!
   TBranch        *b_bsmReweight_cPrime03_brNew03;   //!
   TBranch        *b_bsmReweight_cPrime03_brNew04;   //!
   TBranch        *b_bsmReweight_cPrime03_brNew05;   //!
   TBranch        *b_bsmReweight_cPrime04_brNew00;   //!
   TBranch        *b_bsmReweight_cPrime04_brNew01;   //!
   TBranch        *b_bsmReweight_cPrime04_brNew02;   //!
   TBranch        *b_bsmReweight_cPrime04_brNew03;   //!
   TBranch        *b_bsmReweight_cPrime04_brNew04;   //!
   TBranch        *b_bsmReweight_cPrime04_brNew05;   //!
   TBranch        *b_bsmReweight_cPrime05_brNew00;   //!
   TBranch        *b_bsmReweight_cPrime05_brNew01;   //!
   TBranch        *b_bsmReweight_cPrime05_brNew02;   //!
   TBranch        *b_bsmReweight_cPrime05_brNew03;   //!
   TBranch        *b_bsmReweight_cPrime05_brNew04;   //!
   TBranch        *b_bsmReweight_cPrime05_brNew05;   //!
   TBranch        *b_bsmReweight_cPrime06_brNew00;   //!
   TBranch        *b_bsmReweight_cPrime06_brNew01;   //!
   TBranch        *b_bsmReweight_cPrime06_brNew02;   //!
   TBranch        *b_bsmReweight_cPrime06_brNew03;   //!
   TBranch        *b_bsmReweight_cPrime06_brNew04;   //!
   TBranch        *b_bsmReweight_cPrime06_brNew05;   //!
   TBranch        *b_bsmReweight_cPrime07_brNew00;   //!
   TBranch        *b_bsmReweight_cPrime07_brNew01;   //!
   TBranch        *b_bsmReweight_cPrime07_brNew02;   //!
   TBranch        *b_bsmReweight_cPrime07_brNew03;   //!
   TBranch        *b_bsmReweight_cPrime07_brNew04;   //!
   TBranch        *b_bsmReweight_cPrime07_brNew05;   //!
   TBranch        *b_bsmReweight_cPrime08_brNew00;   //!
   TBranch        *b_bsmReweight_cPrime08_brNew01;   //!
   TBranch        *b_bsmReweight_cPrime08_brNew02;   //!
   TBranch        *b_bsmReweight_cPrime08_brNew03;   //!
   TBranch        *b_bsmReweight_cPrime08_brNew04;   //!
   TBranch        *b_bsmReweight_cPrime08_brNew05;   //!
   TBranch        *b_bsmReweight_cPrime09_brNew00;   //!
   TBranch        *b_bsmReweight_cPrime09_brNew01;   //!
   TBranch        *b_bsmReweight_cPrime09_brNew02;   //!
   TBranch        *b_bsmReweight_cPrime09_brNew03;   //!
   TBranch        *b_bsmReweight_cPrime09_brNew04;   //!
   TBranch        *b_bsmReweight_cPrime09_brNew05;   //!
   TBranch        *b_bsmReweight_cPrime10_brNew00;   //!
   TBranch        *b_bsmReweight_cPrime10_brNew01;   //!
   TBranch        *b_bsmReweight_cPrime10_brNew02;   //!
   TBranch        *b_bsmReweight_cPrime10_brNew03;   //!
   TBranch        *b_bsmReweight_cPrime10_brNew04;   //!
   TBranch        *b_bsmReweight_cPrime10_brNew05;   //!
   TBranch        *b_mass_lvj_type0_met;   //!
   TBranch        *b_mass_lvj_type2_met;   //!
   TBranch        *b_mass_lvj_type0;   //!
   TBranch        *b_mass_lvj_type2;   //!
   TBranch        *b_mass_lvj_type0_met_jes_up;   //!
   TBranch        *b_mass_lvj_type2_met_jes_up;   //!
   TBranch        *b_mass_lvj_type0_met_jes_dn;   //!
   TBranch        *b_mass_lvj_type2_met_jes_dn;   //!
   TBranch        *b_mass_lvj_type0_met_jer;   //!
   TBranch        *b_mass_lvj_type2_met_jer;   //!
   TBranch        *b_mass_lvj_type0_met_jer_up;   //!
   TBranch        *b_mass_lvj_type2_met_jer_up;   //!
   TBranch        *b_mass_lvj_type0_met_jer_dn;   //!
   TBranch        *b_mass_lvj_type2_met_jer_dn;   //!
   TBranch        *b_mass_lvj_type0_met_lep_scale_up;   //!
   TBranch        *b_mass_lvj_type2_met_lep_scale_up;   //!
   TBranch        *b_mass_lvj_type0_met_lep_scale_dn;   //!
   TBranch        *b_mass_lvj_type2_met_lep_scale_dn;   //!
   TBranch        *b_mass_lvj_type0_met_lep_res;   //!
   TBranch        *b_mass_lvj_type2_met_lep_res;   //!
   TBranch        *b_mass_lv_subj_type0_met;   //!
   TBranch        *b_mass_lv_subj_type2_met;   //!
   TBranch        *b_mass_lv_subj_type0;   //!
   TBranch        *b_mass_lv_subj_type2;   //!
   TBranch        *b_mass_ungroomedjet_vbf_j1;   //!
   TBranch        *b_mass_ungroomedjet_vbf_j2;   //!
   TBranch        *b_mass_ungroomedjet_vbf_j1_pr;   //!
   TBranch        *b_mass_ungroomedjet_vbf_j2_pr;   //!
   TBranch        *b_mass_ungroomedjet_closerjet;   //!
   TBranch        *b_mass_ungroomedjet_closerjet_pr;   //!
   TBranch        *b_mass_leptonic_closerjet;   //!
   TBranch        *b_l_pt;   //!
   TBranch        *b_l_eta;   //!
   TBranch        *b_l_charge;   //!
   TBranch        *b_l_phi;   //!
   TBranch        *b_l_pt_scale_up;   //!
   TBranch        *b_l_eta_scale_up;   //!
   TBranch        *b_l_phi_scale_up;   //!
   TBranch        *b_l_pt_scale_dn;   //!
   TBranch        *b_l_eta_scale_dn;   //!
   TBranch        *b_l_phi_scale_dn;   //!
   TBranch        *b_l_pt_res;   //!
   TBranch        *b_l_eta_res;   //!
   TBranch        *b_l_phi_res;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMET_Phi;   //!
   TBranch        *b_pfMET_jes_up;   //!
   TBranch        *b_pfMET_Phi_jes_up;   //!
   TBranch        *b_pfMET_jes_dn;   //!
   TBranch        *b_pfMET_Phi_jes_dn;   //!
   TBranch        *b_pfMET_jer;   //!
   TBranch        *b_pfMET_Phi_jer;   //!
   TBranch        *b_pfMET_jer_up;   //!
   TBranch        *b_pfMET_Phi_jer_up;   //!
   TBranch        *b_pfMET_jer_dn;   //!
   TBranch        *b_pfMET_Phi_jer_dn;   //!
   TBranch        *b_pfMET_lep_scale_up;   //!
   TBranch        *b_pfMET_Phi_lep_scale_up;   //!
   TBranch        *b_pfMET_lep_scale_dn;   //!
   TBranch        *b_pfMET_Phi_lep_scale_dn;   //!
   TBranch        *b_pfMET_lep_res;   //!
   TBranch        *b_pfMET_Phi_lep_res;   //!
   TBranch        *b_nu_pz_type0;   //!
   TBranch        *b_nu_pz_type2;   //!
   TBranch        *b_nu_pz_type0_met;   //!
   TBranch        *b_nu_pz_type2_met;   //!
   TBranch        *b_W_pz_type0;   //!
   TBranch        *b_W_pz_type2;   //!
   TBranch        *b_W_pz_type0_met;   //!
   TBranch        *b_W_pz_type2_met;   //!
   TBranch        *b_nu_pz_gen;   //!
   TBranch        *b_W_pz_gen;   //!
   TBranch        *b_W_pt_gen;   //!
   TBranch        *b_v_pt;   //!
   TBranch        *b_v_mt;   //!
   TBranch        *b_v_eta;   //!
   TBranch        *b_v_phi;   //!
   TBranch        *b_ungroomed_jet_eta;   //!
   TBranch        *b_ungroomed_jet_phi;   //!
   TBranch        *b_ungroomed_jet_pt;   //!
   TBranch        *b_ungroomed_jet_e;   //!
   TBranch        *b_ungroomed_gen_jet_eta;   //!
   TBranch        *b_ungroomed_gen_jet_phi;   //!
   TBranch        *b_ungroomed_gen_jet_pt;   //!
   TBranch        *b_ungroomed_gen_jet_e;   //!
   TBranch        *b_jet_mass_pr;   //!
   TBranch        *b_jet_pt_pr;   //!
   TBranch        *b_jet_charge;   //!
   TBranch        *b_jet_charge_k05;   //!
   TBranch        *b_jet_charge_k07;   //!
   TBranch        *b_jet_charge_k10;   //!
   TBranch        *b_gen_jet_mass_pr;   //!
   TBranch        *b_gen_jet_pt_pr;   //!
   TBranch        *b_jet_grsens_ft;   //!
   TBranch        *b_jet_grsens_tr;   //!
   TBranch        *b_jet_massdrop_pr;   //!
   TBranch        *b_jet_qjetvol;   //!
   TBranch        *b_gen_jet_grsens_ft;   //!
   TBranch        *b_gen_jet_grsens_tr;   //!
   TBranch        *b_gen_jet_massdrop_pr;   //!
   TBranch        *b_gen_jet_qjetvol;   //!
   TBranch        *b_jet_tau2tau1;   //!
   TBranch        *b_jet_tau2tau1_exkT;   //!
   TBranch        *b_jet_tau2tau1_pr;   //!
   TBranch        *b_jet_GeneralizedECF;   //!
   TBranch        *b_gen_jet_tau2tau1;   //!
   TBranch        *b_gen_jet_tau2tau1_exkT;   //!
   TBranch        *b_gen_jet_tau2tau1_pr;   //!
   TBranch        *b_gen_jet_GeneralizedECF;   //!
   TBranch        *b_jet_jetconstituents;   //!
   TBranch        *b_gen_jet_jetconstituents;   //!
   TBranch        *b_jet_rcore4;   //!
   TBranch        *b_jet_rcore5;   //!
   TBranch        *b_jet_rcore6;   //!
   TBranch        *b_jet_rcore7;   //!
   TBranch        *b_gen_jet_rcore4;   //!
   TBranch        *b_gen_jet_rcore5;   //!
   TBranch        *b_gen_jet_rcore6;   //!
   TBranch        *b_gen_jet_rcore7;   //!
   TBranch        *b_jet_pt1frac;   //!
   TBranch        *b_jet_pt2frac;   //!
   TBranch        *b_jet_sjdr;   //!
   TBranch        *b_j_jecfactor_up;   //!
   TBranch        *b_j_jecfactor_dn;   //!
   TBranch        *b_jet_mass_pr_jes_up;   //!
   TBranch        *b_jet_mass_pr_jes_dn;   //!
   TBranch        *b_jet_mass_pr_jer;   //!
   TBranch        *b_jet_mass_pr_jer_up;   //!
   TBranch        *b_jet_mass_pr_jer_dn;   //!
   TBranch        *b_ungroomed_jet_pt_jes_dn;   //!
   TBranch        *b_ungroomed_jet_pt_jes_up;   //!
   TBranch        *b_ungroomed_jet_pt_jer;   //!
   TBranch        *b_ungroomed_jet_pt_jer_dn;   //!
   TBranch        *b_ungroomed_jet_pt_jer_up;   //!
   TBranch        *b_jet_planarlow04;   //!
   TBranch        *b_jet_planarlow05;   //!
   TBranch        *b_jet_planarlow06;   //!
   TBranch        *b_jet_planarlow07;   //!
   TBranch        *b_vbf_maxpt_jj_m;   //!
   TBranch        *b_vbf_maxpt_jj_pt;   //!
   TBranch        *b_vbf_maxpt_jj_eta;   //!
   TBranch        *b_vbf_maxpt_jj_phi;   //!
   TBranch        *b_vbf_maxpt_j1_m;   //!
   TBranch        *b_vbf_maxpt_j1_pt;   //!
   TBranch        *b_vbf_maxpt_j1_eta;   //!
   TBranch        *b_vbf_maxpt_j1_phi;   //!
   TBranch        *b_vbf_maxpt_j2_m;   //!
   TBranch        *b_vbf_maxpt_j2_pt;   //!
   TBranch        *b_vbf_maxpt_j2_eta;   //!
   TBranch        *b_vbf_maxpt_j2_phi;   //!
   TBranch        *b_vbf_maxpt_jj_m_gen;   //!
   TBranch        *b_vbf_maxpt_jj_pt_gen;   //!
   TBranch        *b_vbf_maxpt_jj_eta_gen;   //!
   TBranch        *b_vbf_maxpt_jj_phi_gen;   //!
   TBranch        *b_vbf_maxpt_j1_m_gen;   //!
   TBranch        *b_vbf_maxpt_j1_pt_gen;   //!
   TBranch        *b_vbf_maxpt_j1_eta_gen;   //!
   TBranch        *b_vbf_maxpt_j1_phi_gen;   //!
   TBranch        *b_vbf_maxpt_j2_m_gen;   //!
   TBranch        *b_vbf_maxpt_j2_pt_gen;   //!
   TBranch        *b_vbf_maxpt_j2_eta_gen;   //!
   TBranch        *b_vbf_maxpt_j2_phi_gen;   //!
   TBranch        *b_vbf_maxpt_j1_m_jes_up;   //!
   TBranch        *b_vbf_maxpt_j1_pt_jes_up;   //!
   TBranch        *b_vbf_maxpt_j1_eta_jes_up;   //!
   TBranch        *b_vbf_maxpt_j1_phi_jes_up;   //!
   TBranch        *b_vbf_maxpt_j1_m_jes_dn;   //!
   TBranch        *b_vbf_maxpt_j1_pt_jes_dn;   //!
   TBranch        *b_vbf_maxpt_j1_eta_jes_dn;   //!
   TBranch        *b_vbf_maxpt_j1_phi_jes_dn;   //!
   TBranch        *b_vbf_maxpt_j1_m_jer;   //!
   TBranch        *b_vbf_maxpt_j1_pt_jer;   //!
   TBranch        *b_vbf_maxpt_j1_eta_jer;   //!
   TBranch        *b_vbf_maxpt_j1_phi_jer;   //!
   TBranch        *b_vbf_maxpt_j1_m_jer_up;   //!
   TBranch        *b_vbf_maxpt_j1_pt_jer_up;   //!
   TBranch        *b_vbf_maxpt_j1_eta_jer_up;   //!
   TBranch        *b_vbf_maxpt_j1_phi_jer_up;   //!
   TBranch        *b_vbf_maxpt_j1_m_jer_dn;   //!
   TBranch        *b_vbf_maxpt_j1_pt_jer_dn;   //!
   TBranch        *b_vbf_maxpt_j1_eta_jer_dn;   //!
   TBranch        *b_vbf_maxpt_j1_phi_jer_dn;   //!
   TBranch        *b_vbf_maxpt_j2_m_jes_up;   //!
   TBranch        *b_vbf_maxpt_j2_pt_jes_up;   //!
   TBranch        *b_vbf_maxpt_j2_eta_jes_up;   //!
   TBranch        *b_vbf_maxpt_j2_phi_jes_up;   //!
   TBranch        *b_vbf_maxpt_j2_m_jes_dn;   //!
   TBranch        *b_vbf_maxpt_j2_pt_jes_dn;   //!
   TBranch        *b_vbf_maxpt_j2_eta_jes_dn;   //!
   TBranch        *b_vbf_maxpt_j2_phi_jes_dn;   //!
   TBranch        *b_vbf_maxpt_j2_m_jer;   //!
   TBranch        *b_vbf_maxpt_j2_pt_jer;   //!
   TBranch        *b_vbf_maxpt_j2_eta_jer;   //!
   TBranch        *b_vbf_maxpt_j2_phi_jer;   //!
   TBranch        *b_vbf_maxpt_j2_m_jer_up;   //!
   TBranch        *b_vbf_maxpt_j2_pt_jer_up;   //!
   TBranch        *b_vbf_maxpt_j2_eta_jer_up;   //!
   TBranch        *b_vbf_maxpt_j2_phi_jer_up;   //!
   TBranch        *b_vbf_maxpt_j2_m_jer_dn;   //!
   TBranch        *b_vbf_maxpt_j2_pt_jer_dn;   //!
   TBranch        *b_vbf_maxpt_j2_eta_jer_dn;   //!
   TBranch        *b_vbf_maxpt_j2_phi_jer_dn;   //!
   TBranch        *b_vbf_maxpt_j1_QGLikelihood;   //!
   TBranch        *b_vbf_maxpt_j2_QGLikelihood;   //!
   TBranch        *b_vbf_maxpt_j1_isPileUpMedium;   //!
   TBranch        *b_vbf_maxpt_j2_isPileUpMedium;   //!
   TBranch        *b_vbf_maxpt_j1_isPileUpTight;   //!
   TBranch        *b_vbf_maxpt_j2_isPileUpTight;   //!
   TBranch        *b_vbf_maxpt_j1_bDiscriminatorCSV;   //!
   TBranch        *b_vbf_maxpt_j2_bDiscriminatorCSV;   //!
   TBranch        *b_vbf_maxpt_j1_bDiscriminatorCSV_gen;   //!
   TBranch        *b_vbf_maxpt_j2_bDiscriminatorCSV_gen;   //!
   TBranch        *b_nbjets_csvl_veto;   //!
   TBranch        *b_nbjets_csvm_veto;   //!
   TBranch        *b_nbjets_csvt_veto;   //!
   TBranch        *b_nbjets_ssvhem_veto;   //!
   TBranch        *b_nbjets_csvl_veto_cleaned;   //!
   TBranch        *b_nbjets_csvm_veto_cleaned;   //!
   TBranch        *b_nbjets_csvt_veto_cleaned;   //!
   TBranch        *b_nbjets_ssvhem_veto_cleaned;   //!
   TBranch        *b_njets;   //!
   TBranch        *b_deltaR_lca8jet;   //!
   TBranch        *b_deltaphi_METca8jet;   //!
   TBranch        *b_deltaphi_Vca8jet;   //!
   TBranch        *b_deltaphi_METca8jet_met;   //!
   TBranch        *b_deltaphi_Vca8jet_met;   //!
   TBranch        *b_genHMass;   //!
   TBranch        *b_genHphi;   //!
   TBranch        *b_genHeta;   //!
   TBranch        *b_genHpt;   //!
   TBranch        *b_genTagQuark1E;   //!
   TBranch        *b_genTagQuark1phi;   //!
   TBranch        *b_genTagQuark1eta;   //!
   TBranch        *b_genTagQuark1pt;   //!
   TBranch        *b_genTagQuark2E;   //!
   TBranch        *b_genTagQuark2phi;   //!
   TBranch        *b_genTagQuark2eta;   //!
   TBranch        *b_genTagQuark2pt;   //!
   TBranch        *b_ttb_nak5_same;   //!
   TBranch        *b_ttb_nak5_same_csvl;   //!
   TBranch        *b_ttb_nak5_same_csvm;   //!
   TBranch        *b_ttb_nak5_same_csvt;   //!
   TBranch        *b_ttb_nak5_oppo;   //!
   TBranch        *b_ttb_nak5_oppo_csvl;   //!
   TBranch        *b_ttb_nak5_oppo_csvm;   //!
   TBranch        *b_ttb_nak5_oppo_csvt;   //!
   TBranch        *b_ttb_nak5_oppoveto;   //!
   TBranch        *b_ttb_nak5_oppoveto_csvl;   //!
   TBranch        *b_ttb_nak5_oppoveto_csvm;   //!
   TBranch        *b_ttb_nak5_oppoveto_csvt;   //!
   TBranch        *b_ttb_nak5_sameveto;   //!
   TBranch        *b_ttb_nak5_sameveto_csvl;   //!
   TBranch        *b_ttb_nak5_sameveto_csvm;   //!
   TBranch        *b_ttb_nak5_sameveto_csvt;   //!
   TBranch        *b_ttb_ht;   //!
   TBranch        *b_ttb_ca8_mass_pr;   //!
   TBranch        *b_ttb_ca8_charge;   //!
   TBranch        *b_ttb_ca8_charge_k05;   //!
   TBranch        *b_ttb_ca8_charge_k07;   //!
   TBranch        *b_ttb_ca8_charge_k10;   //!
   TBranch        *b_ttb_ca8_ungroomed_pt;   //!
   TBranch        *b_ttb_ca8_ungroomed_eta;   //!
   TBranch        *b_ttb_ca8_ungroomed_phi;   //!
   TBranch        *b_ttb_ca8_ungroomed_e;   //!
   TBranch        *b_ttb_ca8_ungroomed_gen_pt;   //!
   TBranch        *b_ttb_ca8_ungroomed_gen_eta;   //!
   TBranch        *b_ttb_ca8_ungroomed_gen_phi;   //!
   TBranch        *b_ttb_ca8_ungroomed_gen_e;   //!
   TBranch        *b_ttb_ca8_tau2tau1;   //!
   TBranch        *b_ttb_ca8_tau2tau1_exkT;   //!
   TBranch        *b_ttb_ca8_tau2tau1_pr;   //!
   TBranch        *b_ttb_ca8_GeneralizedECF;   //!
   TBranch        *b_ttb_ca8_mu;   //!
   TBranch        *b_ttb_ca8_mlvj_type0;   //!
   TBranch        *b_ttb_ca8_mlvj_type2;   //!
   TBranch        *b_ttb_ca8_mlvj_type0_met;   //!
   TBranch        *b_ttb_ca8_mlvj_type2_met;   //!
   TBranch        *b_ttb_dR_ca8_bjet_closer;   //!
   TBranch        *b_ttb_dR_ca8_jet_closer;   //!
   TBranch        *b_isttbar;   //!
   TBranch        *b_gen_parton1_px_fromttbar;   //!
   TBranch        *b_gen_parton1_py_fromttbar;   //!
   TBranch        *b_gen_parton1_pz_fromttbar;   //!
   TBranch        *b_gen_parton1_e_fromttbar;   //!
   TBranch        *b_gen_parton1_id_fromttbar;   //!
   TBranch        *b_gen_parton2_px_fromttbar;   //!
   TBranch        *b_gen_parton2_py_fromttbar;   //!
   TBranch        *b_gen_parton2_pz_fromttbar;   //!
   TBranch        *b_gen_parton2_e_fromttbar;   //!
   TBranch        *b_gen_parton2_id_fromttbar;   //!

   data_tree(TTree *tree=0);
   virtual ~data_tree();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual int      GetEntries();
   virtual void     Init(TTree *tree);
};

#endif

data_tree::data_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("light_ntuples/ofile_data.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("light_ntuples/ofile_data.root");
      }
      f->GetObject("data_tree",tree);

   }
   Init(tree);
}

data_tree::~data_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t data_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

int data_tree::GetEntries()
{
    if (!fChain) return 0;
    return fChain->GetEntriesFast();
}

void data_tree::Init(TTree *tree)
{
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("event_lumi", &event_lumi, &b_event_lumi);
   fChain->SetBranchAddress("event_runNo", &event_runNo, &b_event_runNo);
   fChain->SetBranchAddress("issignal", &issignal, &b_issignal);
   fChain->SetBranchAddress("numberJetBin", &numberJetBin, &b_numberJetBin);
   fChain->SetBranchAddress("numberJetBin2", &numberJetBin2, &b_numberJetBin2);
   fChain->SetBranchAddress("numberJetBin3", &numberJetBin3, &b_numberJetBin3);
   fChain->SetBranchAddress("numberJetBin4", &numberJetBin4, &b_numberJetBin4);
   fChain->SetBranchAddress("numberJetBinGen", &numberJetBinGen, &b_numberJetBinGen);
   fChain->SetBranchAddress("numberJetBinGen2", &numberJetBinGen2, &b_numberJetBinGen2);
   fChain->SetBranchAddress("numberJetBinGen3", &numberJetBinGen3, &b_numberJetBinGen3);
   fChain->SetBranchAddress("numberJetBinGen4", &numberJetBinGen4, &b_numberJetBinGen4);
   fChain->SetBranchAddress("totalEventWeight", &totalEventWeight, &b_totalEventWeight);
   fChain->SetBranchAddress("eff_and_pu_Weight", &eff_and_pu_Weight, &b_eff_and_pu_Weight);
   fChain->SetBranchAddress("wSampleWeight", &wSampleWeight, &b_wSampleWeight);
   fChain->SetBranchAddress("event_weight", &event_weight, &b_event_weight);
   fChain->SetBranchAddress("btag_weight", &btag_weight, &b_btag_weight);
   fChain->SetBranchAddress("btag_weight_up", &btag_weight_up, &b_btag_weight_up);
   fChain->SetBranchAddress("btag_weight_dn", &btag_weight_dn, &b_btag_weight_dn);
   fChain->SetBranchAddress("btag_weight_dn_up", &btag_weight_dn_up, &b_btag_weight_dn_up);
   fChain->SetBranchAddress("btag_weight_up_dn", &btag_weight_up_dn, &b_btag_weight_up_dn);
   fChain->SetBranchAddress("interference_Weight_H600", &interference_Weight_H600, &b_interference_Weight_H600);
   fChain->SetBranchAddress("interference_Weight_H700", &interference_Weight_H700, &b_interference_Weight_H700);
   fChain->SetBranchAddress("interference_Weight_H800", &interference_Weight_H800, &b_interference_Weight_H800);
   fChain->SetBranchAddress("interference_Weight_H900", &interference_Weight_H900, &b_interference_Weight_H900);
   fChain->SetBranchAddress("interference_Weight_H1000", &interference_Weight_H1000, &b_interference_Weight_H1000);
   fChain->SetBranchAddress("cps_Weight_H600", &cps_Weight_H600, &b_cps_Weight_H600);
   fChain->SetBranchAddress("cps_Weight_H700", &cps_Weight_H700, &b_cps_Weight_H700);
   fChain->SetBranchAddress("cps_Weight_H800", &cps_Weight_H800, &b_cps_Weight_H800);
   fChain->SetBranchAddress("cps_Weight_H900", &cps_Weight_H900, &b_cps_Weight_H900);
   fChain->SetBranchAddress("cps_Weight_H1000", &cps_Weight_H1000, &b_cps_Weight_H1000);
   fChain->SetBranchAddress("bsmReweight_cPrime01_brNew00", &bsmReweight_cPrime01_brNew00, &b_bsmReweight_cPrime01_brNew00);
   fChain->SetBranchAddress("bsmReweight_cPrime01_brNew01", &bsmReweight_cPrime01_brNew01, &b_bsmReweight_cPrime01_brNew01);
   fChain->SetBranchAddress("bsmReweight_cPrime01_brNew02", &bsmReweight_cPrime01_brNew02, &b_bsmReweight_cPrime01_brNew02);
   fChain->SetBranchAddress("bsmReweight_cPrime01_brNew03", &bsmReweight_cPrime01_brNew03, &b_bsmReweight_cPrime01_brNew03);
   fChain->SetBranchAddress("bsmReweight_cPrime01_brNew04", &bsmReweight_cPrime01_brNew04, &b_bsmReweight_cPrime01_brNew04);
   fChain->SetBranchAddress("bsmReweight_cPrime01_brNew05", &bsmReweight_cPrime01_brNew05, &b_bsmReweight_cPrime01_brNew05);
   fChain->SetBranchAddress("bsmReweight_cPrime02_brNew00", &bsmReweight_cPrime02_brNew00, &b_bsmReweight_cPrime02_brNew00);
   fChain->SetBranchAddress("bsmReweight_cPrime02_brNew01", &bsmReweight_cPrime02_brNew01, &b_bsmReweight_cPrime02_brNew01);
   fChain->SetBranchAddress("bsmReweight_cPrime02_brNew02", &bsmReweight_cPrime02_brNew02, &b_bsmReweight_cPrime02_brNew02);
   fChain->SetBranchAddress("bsmReweight_cPrime02_brNew03", &bsmReweight_cPrime02_brNew03, &b_bsmReweight_cPrime02_brNew03);
   fChain->SetBranchAddress("bsmReweight_cPrime02_brNew04", &bsmReweight_cPrime02_brNew04, &b_bsmReweight_cPrime02_brNew04);
   fChain->SetBranchAddress("bsmReweight_cPrime02_brNew05", &bsmReweight_cPrime02_brNew05, &b_bsmReweight_cPrime02_brNew05);
   fChain->SetBranchAddress("bsmReweight_cPrime03_brNew00", &bsmReweight_cPrime03_brNew00, &b_bsmReweight_cPrime03_brNew00);
   fChain->SetBranchAddress("bsmReweight_cPrime03_brNew01", &bsmReweight_cPrime03_brNew01, &b_bsmReweight_cPrime03_brNew01);
   fChain->SetBranchAddress("bsmReweight_cPrime03_brNew02", &bsmReweight_cPrime03_brNew02, &b_bsmReweight_cPrime03_brNew02);
   fChain->SetBranchAddress("bsmReweight_cPrime03_brNew03", &bsmReweight_cPrime03_brNew03, &b_bsmReweight_cPrime03_brNew03);
   fChain->SetBranchAddress("bsmReweight_cPrime03_brNew04", &bsmReweight_cPrime03_brNew04, &b_bsmReweight_cPrime03_brNew04);
   fChain->SetBranchAddress("bsmReweight_cPrime03_brNew05", &bsmReweight_cPrime03_brNew05, &b_bsmReweight_cPrime03_brNew05);
   fChain->SetBranchAddress("bsmReweight_cPrime04_brNew00", &bsmReweight_cPrime04_brNew00, &b_bsmReweight_cPrime04_brNew00);
   fChain->SetBranchAddress("bsmReweight_cPrime04_brNew01", &bsmReweight_cPrime04_brNew01, &b_bsmReweight_cPrime04_brNew01);
   fChain->SetBranchAddress("bsmReweight_cPrime04_brNew02", &bsmReweight_cPrime04_brNew02, &b_bsmReweight_cPrime04_brNew02);
   fChain->SetBranchAddress("bsmReweight_cPrime04_brNew03", &bsmReweight_cPrime04_brNew03, &b_bsmReweight_cPrime04_brNew03);
   fChain->SetBranchAddress("bsmReweight_cPrime04_brNew04", &bsmReweight_cPrime04_brNew04, &b_bsmReweight_cPrime04_brNew04);
   fChain->SetBranchAddress("bsmReweight_cPrime04_brNew05", &bsmReweight_cPrime04_brNew05, &b_bsmReweight_cPrime04_brNew05);
   fChain->SetBranchAddress("bsmReweight_cPrime05_brNew00", &bsmReweight_cPrime05_brNew00, &b_bsmReweight_cPrime05_brNew00);
   fChain->SetBranchAddress("bsmReweight_cPrime05_brNew01", &bsmReweight_cPrime05_brNew01, &b_bsmReweight_cPrime05_brNew01);
   fChain->SetBranchAddress("bsmReweight_cPrime05_brNew02", &bsmReweight_cPrime05_brNew02, &b_bsmReweight_cPrime05_brNew02);
   fChain->SetBranchAddress("bsmReweight_cPrime05_brNew03", &bsmReweight_cPrime05_brNew03, &b_bsmReweight_cPrime05_brNew03);
   fChain->SetBranchAddress("bsmReweight_cPrime05_brNew04", &bsmReweight_cPrime05_brNew04, &b_bsmReweight_cPrime05_brNew04);
   fChain->SetBranchAddress("bsmReweight_cPrime05_brNew05", &bsmReweight_cPrime05_brNew05, &b_bsmReweight_cPrime05_brNew05);
   fChain->SetBranchAddress("bsmReweight_cPrime06_brNew00", &bsmReweight_cPrime06_brNew00, &b_bsmReweight_cPrime06_brNew00);
   fChain->SetBranchAddress("bsmReweight_cPrime06_brNew01", &bsmReweight_cPrime06_brNew01, &b_bsmReweight_cPrime06_brNew01);
   fChain->SetBranchAddress("bsmReweight_cPrime06_brNew02", &bsmReweight_cPrime06_brNew02, &b_bsmReweight_cPrime06_brNew02);
   fChain->SetBranchAddress("bsmReweight_cPrime06_brNew03", &bsmReweight_cPrime06_brNew03, &b_bsmReweight_cPrime06_brNew03);
   fChain->SetBranchAddress("bsmReweight_cPrime06_brNew04", &bsmReweight_cPrime06_brNew04, &b_bsmReweight_cPrime06_brNew04);
   fChain->SetBranchAddress("bsmReweight_cPrime06_brNew05", &bsmReweight_cPrime06_brNew05, &b_bsmReweight_cPrime06_brNew05);
   fChain->SetBranchAddress("bsmReweight_cPrime07_brNew00", &bsmReweight_cPrime07_brNew00, &b_bsmReweight_cPrime07_brNew00);
   fChain->SetBranchAddress("bsmReweight_cPrime07_brNew01", &bsmReweight_cPrime07_brNew01, &b_bsmReweight_cPrime07_brNew01);
   fChain->SetBranchAddress("bsmReweight_cPrime07_brNew02", &bsmReweight_cPrime07_brNew02, &b_bsmReweight_cPrime07_brNew02);
   fChain->SetBranchAddress("bsmReweight_cPrime07_brNew03", &bsmReweight_cPrime07_brNew03, &b_bsmReweight_cPrime07_brNew03);
   fChain->SetBranchAddress("bsmReweight_cPrime07_brNew04", &bsmReweight_cPrime07_brNew04, &b_bsmReweight_cPrime07_brNew04);
   fChain->SetBranchAddress("bsmReweight_cPrime07_brNew05", &bsmReweight_cPrime07_brNew05, &b_bsmReweight_cPrime07_brNew05);
   fChain->SetBranchAddress("bsmReweight_cPrime08_brNew00", &bsmReweight_cPrime08_brNew00, &b_bsmReweight_cPrime08_brNew00);
   fChain->SetBranchAddress("bsmReweight_cPrime08_brNew01", &bsmReweight_cPrime08_brNew01, &b_bsmReweight_cPrime08_brNew01);
   fChain->SetBranchAddress("bsmReweight_cPrime08_brNew02", &bsmReweight_cPrime08_brNew02, &b_bsmReweight_cPrime08_brNew02);
   fChain->SetBranchAddress("bsmReweight_cPrime08_brNew03", &bsmReweight_cPrime08_brNew03, &b_bsmReweight_cPrime08_brNew03);
   fChain->SetBranchAddress("bsmReweight_cPrime08_brNew04", &bsmReweight_cPrime08_brNew04, &b_bsmReweight_cPrime08_brNew04);
   fChain->SetBranchAddress("bsmReweight_cPrime08_brNew05", &bsmReweight_cPrime08_brNew05, &b_bsmReweight_cPrime08_brNew05);
   fChain->SetBranchAddress("bsmReweight_cPrime09_brNew00", &bsmReweight_cPrime09_brNew00, &b_bsmReweight_cPrime09_brNew00);
   fChain->SetBranchAddress("bsmReweight_cPrime09_brNew01", &bsmReweight_cPrime09_brNew01, &b_bsmReweight_cPrime09_brNew01);
   fChain->SetBranchAddress("bsmReweight_cPrime09_brNew02", &bsmReweight_cPrime09_brNew02, &b_bsmReweight_cPrime09_brNew02);
   fChain->SetBranchAddress("bsmReweight_cPrime09_brNew03", &bsmReweight_cPrime09_brNew03, &b_bsmReweight_cPrime09_brNew03);
   fChain->SetBranchAddress("bsmReweight_cPrime09_brNew04", &bsmReweight_cPrime09_brNew04, &b_bsmReweight_cPrime09_brNew04);
   fChain->SetBranchAddress("bsmReweight_cPrime09_brNew05", &bsmReweight_cPrime09_brNew05, &b_bsmReweight_cPrime09_brNew05);
   fChain->SetBranchAddress("bsmReweight_cPrime10_brNew00", &bsmReweight_cPrime10_brNew00, &b_bsmReweight_cPrime10_brNew00);
   fChain->SetBranchAddress("bsmReweight_cPrime10_brNew01", &bsmReweight_cPrime10_brNew01, &b_bsmReweight_cPrime10_brNew01);
   fChain->SetBranchAddress("bsmReweight_cPrime10_brNew02", &bsmReweight_cPrime10_brNew02, &b_bsmReweight_cPrime10_brNew02);
   fChain->SetBranchAddress("bsmReweight_cPrime10_brNew03", &bsmReweight_cPrime10_brNew03, &b_bsmReweight_cPrime10_brNew03);
   fChain->SetBranchAddress("bsmReweight_cPrime10_brNew04", &bsmReweight_cPrime10_brNew04, &b_bsmReweight_cPrime10_brNew04);
   fChain->SetBranchAddress("bsmReweight_cPrime10_brNew05", &bsmReweight_cPrime10_brNew05, &b_bsmReweight_cPrime10_brNew05);
   fChain->SetBranchAddress("mass_lvj_type0_met", &mass_lvj_type0_met, &b_mass_lvj_type0_met);
   fChain->SetBranchAddress("mass_lvj_type2_met", &mass_lvj_type2_met, &b_mass_lvj_type2_met);
   fChain->SetBranchAddress("mass_lvj_type0", &mass_lvj_type0, &b_mass_lvj_type0);
   fChain->SetBranchAddress("mass_lvj_type2", &mass_lvj_type2, &b_mass_lvj_type2);
   fChain->SetBranchAddress("mass_lvj_type0_met_jes_up", &mass_lvj_type0_met_jes_up, &b_mass_lvj_type0_met_jes_up);
   fChain->SetBranchAddress("mass_lvj_type2_met_jes_up", &mass_lvj_type2_met_jes_up, &b_mass_lvj_type2_met_jes_up);
   fChain->SetBranchAddress("mass_lvj_type0_met_jes_dn", &mass_lvj_type0_met_jes_dn, &b_mass_lvj_type0_met_jes_dn);
   fChain->SetBranchAddress("mass_lvj_type2_met_jes_dn", &mass_lvj_type2_met_jes_dn, &b_mass_lvj_type2_met_jes_dn);
   fChain->SetBranchAddress("mass_lvj_type0_met_jer", &mass_lvj_type0_met_jer, &b_mass_lvj_type0_met_jer);
   fChain->SetBranchAddress("mass_lvj_type2_met_jer", &mass_lvj_type2_met_jer, &b_mass_lvj_type2_met_jer);
   fChain->SetBranchAddress("mass_lvj_type0_met_jer_up", &mass_lvj_type0_met_jer_up, &b_mass_lvj_type0_met_jer_up);
   fChain->SetBranchAddress("mass_lvj_type2_met_jer_up", &mass_lvj_type2_met_jer_up, &b_mass_lvj_type2_met_jer_up);
   fChain->SetBranchAddress("mass_lvj_type0_met_jer_dn", &mass_lvj_type0_met_jer_dn, &b_mass_lvj_type0_met_jer_dn);
   fChain->SetBranchAddress("mass_lvj_type2_met_jer_dn", &mass_lvj_type2_met_jer_dn, &b_mass_lvj_type2_met_jer_dn);
   fChain->SetBranchAddress("mass_lvj_type0_met_lep_scale_up", &mass_lvj_type0_met_lep_scale_up, &b_mass_lvj_type0_met_lep_scale_up);
   fChain->SetBranchAddress("mass_lvj_type2_met_lep_scale_up", &mass_lvj_type2_met_lep_scale_up, &b_mass_lvj_type2_met_lep_scale_up);
   fChain->SetBranchAddress("mass_lvj_type0_met_lep_scale_dn", &mass_lvj_type0_met_lep_scale_dn, &b_mass_lvj_type0_met_lep_scale_dn);
   fChain->SetBranchAddress("mass_lvj_type2_met_lep_scale_dn", &mass_lvj_type2_met_lep_scale_dn, &b_mass_lvj_type2_met_lep_scale_dn);
   fChain->SetBranchAddress("mass_lvj_type0_met_lep_res", &mass_lvj_type0_met_lep_res, &b_mass_lvj_type0_met_lep_res);
   fChain->SetBranchAddress("mass_lvj_type2_met_lep_res", &mass_lvj_type2_met_lep_res, &b_mass_lvj_type2_met_lep_res);
   fChain->SetBranchAddress("mass_lv_subj_type0_met", &mass_lv_subj_type0_met, &b_mass_lv_subj_type0_met);
   fChain->SetBranchAddress("mass_lv_subj_type2_met", &mass_lv_subj_type2_met, &b_mass_lv_subj_type2_met);
   fChain->SetBranchAddress("mass_lv_subj_type0", &mass_lv_subj_type0, &b_mass_lv_subj_type0);
   fChain->SetBranchAddress("mass_lv_subj_type2", &mass_lv_subj_type2, &b_mass_lv_subj_type2);
   fChain->SetBranchAddress("mass_ungroomedjet_vbf_j1", &mass_ungroomedjet_vbf_j1, &b_mass_ungroomedjet_vbf_j1);
   fChain->SetBranchAddress("mass_ungroomedjet_vbf_j2", &mass_ungroomedjet_vbf_j2, &b_mass_ungroomedjet_vbf_j2);
   fChain->SetBranchAddress("mass_ungroomedjet_vbf_j1_pr", &mass_ungroomedjet_vbf_j1_pr, &b_mass_ungroomedjet_vbf_j1_pr);
   fChain->SetBranchAddress("mass_ungroomedjet_vbf_j2_pr", &mass_ungroomedjet_vbf_j2_pr, &b_mass_ungroomedjet_vbf_j2_pr);
   fChain->SetBranchAddress("mass_ungroomedjet_closerjet", &mass_ungroomedjet_closerjet, &b_mass_ungroomedjet_closerjet);
   fChain->SetBranchAddress("mass_ungroomedjet_closerjet_pr", &mass_ungroomedjet_closerjet_pr, &b_mass_ungroomedjet_closerjet_pr);
   fChain->SetBranchAddress("mass_leptonic_closerjet", &mass_leptonic_closerjet, &b_mass_leptonic_closerjet);
   fChain->SetBranchAddress("l_pt", &l_pt, &b_l_pt);
   fChain->SetBranchAddress("l_eta", &l_eta, &b_l_eta);
   fChain->SetBranchAddress("l_charge", &l_charge, &b_l_charge);
   fChain->SetBranchAddress("l_phi", &l_phi, &b_l_phi);
   fChain->SetBranchAddress("l_pt_scale_up", &l_pt_scale_up, &b_l_pt_scale_up);
   fChain->SetBranchAddress("l_eta_scale_up", &l_eta_scale_up, &b_l_eta_scale_up);
   fChain->SetBranchAddress("l_phi_scale_up", &l_phi_scale_up, &b_l_phi_scale_up);
   fChain->SetBranchAddress("l_pt_scale_dn", &l_pt_scale_dn, &b_l_pt_scale_dn);
   fChain->SetBranchAddress("l_eta_scale_dn", &l_eta_scale_dn, &b_l_eta_scale_dn);
   fChain->SetBranchAddress("l_phi_scale_dn", &l_phi_scale_dn, &b_l_phi_scale_dn);
   fChain->SetBranchAddress("l_pt_res", &l_pt_res, &b_l_pt_res);
   fChain->SetBranchAddress("l_eta_res", &l_eta_res, &b_l_eta_res);
   fChain->SetBranchAddress("l_phi_res", &l_phi_res, &b_l_phi_res);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMET_Phi", &pfMET_Phi, &b_pfMET_Phi);
   fChain->SetBranchAddress("pfMET_jes_up", &pfMET_jes_up, &b_pfMET_jes_up);
   fChain->SetBranchAddress("pfMET_Phi_jes_up", &pfMET_Phi_jes_up, &b_pfMET_Phi_jes_up);
   fChain->SetBranchAddress("pfMET_jes_dn", &pfMET_jes_dn, &b_pfMET_jes_dn);
   fChain->SetBranchAddress("pfMET_Phi_jes_dn", &pfMET_Phi_jes_dn, &b_pfMET_Phi_jes_dn);
   fChain->SetBranchAddress("pfMET_jer_", &pfMET_jer_, &b_pfMET_jer);
   fChain->SetBranchAddress("pfMET_Phi_jer_", &pfMET_Phi_jer_, &b_pfMET_Phi_jer);
   fChain->SetBranchAddress("pfMET_jer_up", &pfMET_jer_up, &b_pfMET_jer_up);
   fChain->SetBranchAddress("pfMET_Phi_jer_up", &pfMET_Phi_jer_up, &b_pfMET_Phi_jer_up);
   fChain->SetBranchAddress("pfMET_jer_dn", &pfMET_jer_dn, &b_pfMET_jer_dn);
   fChain->SetBranchAddress("pfMET_Phi_jer_dn", &pfMET_Phi_jer_dn, &b_pfMET_Phi_jer_dn);
   fChain->SetBranchAddress("pfMET_lep_scale_up", &pfMET_lep_scale_up, &b_pfMET_lep_scale_up);
   fChain->SetBranchAddress("pfMET_Phi_lep_scale_up", &pfMET_Phi_lep_scale_up, &b_pfMET_Phi_lep_scale_up);
   fChain->SetBranchAddress("pfMET_lep_scale_dn", &pfMET_lep_scale_dn, &b_pfMET_lep_scale_dn);
   fChain->SetBranchAddress("pfMET_Phi_lep_scale_dn", &pfMET_Phi_lep_scale_dn, &b_pfMET_Phi_lep_scale_dn);
   fChain->SetBranchAddress("pfMET_lep_res", &pfMET_lep_res, &b_pfMET_lep_res);
   fChain->SetBranchAddress("pfMET_Phi_lep_res", &pfMET_Phi_lep_res, &b_pfMET_Phi_lep_res);
   fChain->SetBranchAddress("nu_pz_type0", &nu_pz_type0, &b_nu_pz_type0);
   fChain->SetBranchAddress("nu_pz_type2", &nu_pz_type2, &b_nu_pz_type2);
   fChain->SetBranchAddress("nu_pz_type0_met", &nu_pz_type0_met, &b_nu_pz_type0_met);
   fChain->SetBranchAddress("nu_pz_type2_met", &nu_pz_type2_met, &b_nu_pz_type2_met);
   fChain->SetBranchAddress("W_pz_type0", &W_pz_type0, &b_W_pz_type0);
   fChain->SetBranchAddress("W_pz_type2", &W_pz_type2, &b_W_pz_type2);
   fChain->SetBranchAddress("W_pz_type0_met", &W_pz_type0_met, &b_W_pz_type0_met);
   fChain->SetBranchAddress("W_pz_type2_met", &W_pz_type2_met, &b_W_pz_type2_met);
   fChain->SetBranchAddress("nu_pz_gen", &nu_pz_gen, &b_nu_pz_gen);
   fChain->SetBranchAddress("W_pz_gen", &W_pz_gen, &b_W_pz_gen);
   fChain->SetBranchAddress("W_pt_gen", &W_pt_gen, &b_W_pt_gen);
   fChain->SetBranchAddress("v_pt", &v_pt, &b_v_pt);
   fChain->SetBranchAddress("v_mt", &v_mt, &b_v_mt);
   fChain->SetBranchAddress("v_eta", &v_eta, &b_v_eta);
   fChain->SetBranchAddress("v_phi", &v_phi, &b_v_phi);
   fChain->SetBranchAddress("ungroomed_jet_eta", &ungroomed_jet_eta, &b_ungroomed_jet_eta);
   fChain->SetBranchAddress("ungroomed_jet_phi", &ungroomed_jet_phi, &b_ungroomed_jet_phi);
   fChain->SetBranchAddress("ungroomed_jet_pt", &ungroomed_jet_pt, &b_ungroomed_jet_pt);
   fChain->SetBranchAddress("ungroomed_jet_e", &ungroomed_jet_e, &b_ungroomed_jet_e);
   fChain->SetBranchAddress("ungroomed_gen_jet_eta", &ungroomed_gen_jet_eta, &b_ungroomed_gen_jet_eta);
   fChain->SetBranchAddress("ungroomed_gen_jet_phi", &ungroomed_gen_jet_phi, &b_ungroomed_gen_jet_phi);
   fChain->SetBranchAddress("ungroomed_gen_jet_pt", &ungroomed_gen_jet_pt, &b_ungroomed_gen_jet_pt);
   fChain->SetBranchAddress("ungroomed_gen_jet_e", &ungroomed_gen_jet_e, &b_ungroomed_gen_jet_e);
   fChain->SetBranchAddress("jet_mass_pr", &jet_mass_pr, &b_jet_mass_pr);
   fChain->SetBranchAddress("jet_pt_pr", &jet_pt_pr, &b_jet_pt_pr);
   fChain->SetBranchAddress("jet_charge", &jet_charge, &b_jet_charge);
   fChain->SetBranchAddress("jet_charge_k05", &jet_charge_k05, &b_jet_charge_k05);
   fChain->SetBranchAddress("jet_charge_k07", &jet_charge_k07, &b_jet_charge_k07);
   fChain->SetBranchAddress("jet_charge_k10", &jet_charge_k10, &b_jet_charge_k10);
   fChain->SetBranchAddress("gen_jet_mass_pr", &gen_jet_mass_pr, &b_gen_jet_mass_pr);
   fChain->SetBranchAddress("gen_jet_pt_pr", &gen_jet_pt_pr, &b_gen_jet_pt_pr);
   fChain->SetBranchAddress("jet_grsens_ft", &jet_grsens_ft, &b_jet_grsens_ft);
   fChain->SetBranchAddress("jet_grsens_tr", &jet_grsens_tr, &b_jet_grsens_tr);
   fChain->SetBranchAddress("jet_massdrop_pr", &jet_massdrop_pr, &b_jet_massdrop_pr);
   fChain->SetBranchAddress("jet_qjetvol", &jet_qjetvol, &b_jet_qjetvol);
   fChain->SetBranchAddress("gen_jet_grsens_ft", &gen_jet_grsens_ft, &b_gen_jet_grsens_ft);
   fChain->SetBranchAddress("gen_jet_grsens_tr", &gen_jet_grsens_tr, &b_gen_jet_grsens_tr);
   fChain->SetBranchAddress("gen_jet_massdrop_pr", &gen_jet_massdrop_pr, &b_gen_jet_massdrop_pr);
   fChain->SetBranchAddress("gen_jet_qjetvol", &gen_jet_qjetvol, &b_gen_jet_qjetvol);
   fChain->SetBranchAddress("jet_tau2tau1", &jet_tau2tau1, &b_jet_tau2tau1);
   fChain->SetBranchAddress("jet_tau2tau1_exkT", &jet_tau2tau1_exkT, &b_jet_tau2tau1_exkT);
   fChain->SetBranchAddress("jet_tau2tau1_pr", &jet_tau2tau1_pr, &b_jet_tau2tau1_pr);
   fChain->SetBranchAddress("jet_GeneralizedECF", &jet_GeneralizedECF, &b_jet_GeneralizedECF);
   fChain->SetBranchAddress("gen_jet_tau2tau1", &gen_jet_tau2tau1, &b_gen_jet_tau2tau1);
   fChain->SetBranchAddress("gen_jet_tau2tau1_exkT", &gen_jet_tau2tau1_exkT, &b_gen_jet_tau2tau1_exkT);
   fChain->SetBranchAddress("gen_jet_tau2tau1_pr", &gen_jet_tau2tau1_pr, &b_gen_jet_tau2tau1_pr);
   fChain->SetBranchAddress("gen_jet_GeneralizedECF", &gen_jet_GeneralizedECF, &b_gen_jet_GeneralizedECF);
   fChain->SetBranchAddress("jet_jetconstituents", &jet_jetconstituents, &b_jet_jetconstituents);
   fChain->SetBranchAddress("gen_jet_jetconstituents", &gen_jet_jetconstituents, &b_gen_jet_jetconstituents);
   fChain->SetBranchAddress("jet_rcore4", &jet_rcore4, &b_jet_rcore4);
   fChain->SetBranchAddress("jet_rcore5", &jet_rcore5, &b_jet_rcore5);
   fChain->SetBranchAddress("jet_rcore6", &jet_rcore6, &b_jet_rcore6);
   fChain->SetBranchAddress("jet_rcore7", &jet_rcore7, &b_jet_rcore7);
   fChain->SetBranchAddress("gen_jet_rcore4", &gen_jet_rcore4, &b_gen_jet_rcore4);
   fChain->SetBranchAddress("gen_jet_rcore5", &gen_jet_rcore5, &b_gen_jet_rcore5);
   fChain->SetBranchAddress("gen_jet_rcore6", &gen_jet_rcore6, &b_gen_jet_rcore6);
   fChain->SetBranchAddress("gen_jet_rcore7", &gen_jet_rcore7, &b_gen_jet_rcore7);
   fChain->SetBranchAddress("jet_pt1frac", &jet_pt1frac, &b_jet_pt1frac);
   fChain->SetBranchAddress("jet_pt2frac", &jet_pt2frac, &b_jet_pt2frac);
   fChain->SetBranchAddress("jet_sjdr", &jet_sjdr, &b_jet_sjdr);
   fChain->SetBranchAddress("j_jecfactor_up", &j_jecfactor_up, &b_j_jecfactor_up);
   fChain->SetBranchAddress("j_jecfactor_dn", &j_jecfactor_dn, &b_j_jecfactor_dn);
   fChain->SetBranchAddress("jet_mass_pr_jes_up", &jet_mass_pr_jes_up, &b_jet_mass_pr_jes_up);
   fChain->SetBranchAddress("jet_mass_pr_jes_dn", &jet_mass_pr_jes_dn, &b_jet_mass_pr_jes_dn);
   fChain->SetBranchAddress("jet_mass_pr_jer", &jet_mass_pr_jer, &b_jet_mass_pr_jer);
   fChain->SetBranchAddress("jet_mass_pr_jer_up", &jet_mass_pr_jer_up, &b_jet_mass_pr_jer_up);
   fChain->SetBranchAddress("jet_mass_pr_jer_dn", &jet_mass_pr_jer_dn, &b_jet_mass_pr_jer_dn);
   fChain->SetBranchAddress("ungroomed_jet_pt_jes_dn", &ungroomed_jet_pt_jes_dn, &b_ungroomed_jet_pt_jes_dn);
   fChain->SetBranchAddress("ungroomed_jet_pt_jes_up", &ungroomed_jet_pt_jes_up, &b_ungroomed_jet_pt_jes_up);
   fChain->SetBranchAddress("ungroomed_jet_pt_jer", &ungroomed_jet_pt_jer, &b_ungroomed_jet_pt_jer);
   fChain->SetBranchAddress("ungroomed_jet_pt_jer_dn", &ungroomed_jet_pt_jer_dn, &b_ungroomed_jet_pt_jer_dn);
   fChain->SetBranchAddress("ungroomed_jet_pt_jer_up", &ungroomed_jet_pt_jer_up, &b_ungroomed_jet_pt_jer_up);
   fChain->SetBranchAddress("jet_planarlow04", &jet_planarlow04, &b_jet_planarlow04);
   fChain->SetBranchAddress("jet_planarlow05", &jet_planarlow05, &b_jet_planarlow05);
   fChain->SetBranchAddress("jet_planarlow06", &jet_planarlow06, &b_jet_planarlow06);
   fChain->SetBranchAddress("jet_planarlow07", &jet_planarlow07, &b_jet_planarlow07);
   fChain->SetBranchAddress("vbf_maxpt_jj_m", &vbf_maxpt_jj_m, &b_vbf_maxpt_jj_m);
   fChain->SetBranchAddress("vbf_maxpt_jj_pt", &vbf_maxpt_jj_pt, &b_vbf_maxpt_jj_pt);
   fChain->SetBranchAddress("vbf_maxpt_jj_eta", &vbf_maxpt_jj_eta, &b_vbf_maxpt_jj_eta);
   fChain->SetBranchAddress("vbf_maxpt_jj_phi", &vbf_maxpt_jj_phi, &b_vbf_maxpt_jj_phi);
   fChain->SetBranchAddress("vbf_maxpt_j1_m", &vbf_maxpt_j1_m, &b_vbf_maxpt_j1_m);
   fChain->SetBranchAddress("vbf_maxpt_j1_pt", &vbf_maxpt_j1_pt, &b_vbf_maxpt_j1_pt);
   fChain->SetBranchAddress("vbf_maxpt_j1_eta", &vbf_maxpt_j1_eta, &b_vbf_maxpt_j1_eta);
   fChain->SetBranchAddress("vbf_maxpt_j1_phi", &vbf_maxpt_j1_phi, &b_vbf_maxpt_j1_phi);
   fChain->SetBranchAddress("vbf_maxpt_j2_m", &vbf_maxpt_j2_m, &b_vbf_maxpt_j2_m);
   fChain->SetBranchAddress("vbf_maxpt_j2_pt", &vbf_maxpt_j2_pt, &b_vbf_maxpt_j2_pt);
   fChain->SetBranchAddress("vbf_maxpt_j2_eta", &vbf_maxpt_j2_eta, &b_vbf_maxpt_j2_eta);
   fChain->SetBranchAddress("vbf_maxpt_j2_phi", &vbf_maxpt_j2_phi, &b_vbf_maxpt_j2_phi);
   fChain->SetBranchAddress("vbf_maxpt_jj_m_gen", &vbf_maxpt_jj_m_gen, &b_vbf_maxpt_jj_m_gen);
   fChain->SetBranchAddress("vbf_maxpt_jj_pt_gen", &vbf_maxpt_jj_pt_gen, &b_vbf_maxpt_jj_pt_gen);
   fChain->SetBranchAddress("vbf_maxpt_jj_eta_gen", &vbf_maxpt_jj_eta_gen, &b_vbf_maxpt_jj_eta_gen);
   fChain->SetBranchAddress("vbf_maxpt_jj_phi_gen", &vbf_maxpt_jj_phi_gen, &b_vbf_maxpt_jj_phi_gen);
   fChain->SetBranchAddress("vbf_maxpt_j1_m_gen", &vbf_maxpt_j1_m_gen, &b_vbf_maxpt_j1_m_gen);
   fChain->SetBranchAddress("vbf_maxpt_j1_pt_gen", &vbf_maxpt_j1_pt_gen, &b_vbf_maxpt_j1_pt_gen);
   fChain->SetBranchAddress("vbf_maxpt_j1_eta_gen", &vbf_maxpt_j1_eta_gen, &b_vbf_maxpt_j1_eta_gen);
   fChain->SetBranchAddress("vbf_maxpt_j1_phi_gen", &vbf_maxpt_j1_phi_gen, &b_vbf_maxpt_j1_phi_gen);
   fChain->SetBranchAddress("vbf_maxpt_j2_m_gen", &vbf_maxpt_j2_m_gen, &b_vbf_maxpt_j2_m_gen);
   fChain->SetBranchAddress("vbf_maxpt_j2_pt_gen", &vbf_maxpt_j2_pt_gen, &b_vbf_maxpt_j2_pt_gen);
   fChain->SetBranchAddress("vbf_maxpt_j2_eta_gen", &vbf_maxpt_j2_eta_gen, &b_vbf_maxpt_j2_eta_gen);
   fChain->SetBranchAddress("vbf_maxpt_j2_phi_gen", &vbf_maxpt_j2_phi_gen, &b_vbf_maxpt_j2_phi_gen);
   fChain->SetBranchAddress("vbf_maxpt_j1_m_jes_up", &vbf_maxpt_j1_m_jes_up, &b_vbf_maxpt_j1_m_jes_up);
   fChain->SetBranchAddress("vbf_maxpt_j1_pt_jes_up", &vbf_maxpt_j1_pt_jes_up, &b_vbf_maxpt_j1_pt_jes_up);
   fChain->SetBranchAddress("vbf_maxpt_j1_eta_jes_up", &vbf_maxpt_j1_eta_jes_up, &b_vbf_maxpt_j1_eta_jes_up);
   fChain->SetBranchAddress("vbf_maxpt_j1_phi_jes_up", &vbf_maxpt_j1_phi_jes_up, &b_vbf_maxpt_j1_phi_jes_up);
   fChain->SetBranchAddress("vbf_maxpt_j1_m_jes_dn", &vbf_maxpt_j1_m_jes_dn, &b_vbf_maxpt_j1_m_jes_dn);
   fChain->SetBranchAddress("vbf_maxpt_j1_pt_jes_dn", &vbf_maxpt_j1_pt_jes_dn, &b_vbf_maxpt_j1_pt_jes_dn);
   fChain->SetBranchAddress("vbf_maxpt_j1_eta_jes_dn", &vbf_maxpt_j1_eta_jes_dn, &b_vbf_maxpt_j1_eta_jes_dn);
   fChain->SetBranchAddress("vbf_maxpt_j1_phi_jes_dn", &vbf_maxpt_j1_phi_jes_dn, &b_vbf_maxpt_j1_phi_jes_dn);
   fChain->SetBranchAddress("vbf_maxpt_j1_m_jer", &vbf_maxpt_j1_m_jer, &b_vbf_maxpt_j1_m_jer);
   fChain->SetBranchAddress("vbf_maxpt_j1_pt_jer", &vbf_maxpt_j1_pt_jer, &b_vbf_maxpt_j1_pt_jer);
   fChain->SetBranchAddress("vbf_maxpt_j1_eta_jer", &vbf_maxpt_j1_eta_jer, &b_vbf_maxpt_j1_eta_jer);
   fChain->SetBranchAddress("vbf_maxpt_j1_phi_jer", &vbf_maxpt_j1_phi_jer, &b_vbf_maxpt_j1_phi_jer);
   fChain->SetBranchAddress("vbf_maxpt_j1_m_jer_up", &vbf_maxpt_j1_m_jer_up, &b_vbf_maxpt_j1_m_jer_up);
   fChain->SetBranchAddress("vbf_maxpt_j1_pt_jer_up", &vbf_maxpt_j1_pt_jer_up, &b_vbf_maxpt_j1_pt_jer_up);
   fChain->SetBranchAddress("vbf_maxpt_j1_eta_jer_up", &vbf_maxpt_j1_eta_jer_up, &b_vbf_maxpt_j1_eta_jer_up);
   fChain->SetBranchAddress("vbf_maxpt_j1_phi_jer_up", &vbf_maxpt_j1_phi_jer_up, &b_vbf_maxpt_j1_phi_jer_up);
   fChain->SetBranchAddress("vbf_maxpt_j1_m_jer_dn", &vbf_maxpt_j1_m_jer_dn, &b_vbf_maxpt_j1_m_jer_dn);
   fChain->SetBranchAddress("vbf_maxpt_j1_pt_jer_dn", &vbf_maxpt_j1_pt_jer_dn, &b_vbf_maxpt_j1_pt_jer_dn);
   fChain->SetBranchAddress("vbf_maxpt_j1_eta_jer_dn", &vbf_maxpt_j1_eta_jer_dn, &b_vbf_maxpt_j1_eta_jer_dn);
   fChain->SetBranchAddress("vbf_maxpt_j1_phi_jer_dn", &vbf_maxpt_j1_phi_jer_dn, &b_vbf_maxpt_j1_phi_jer_dn);
   fChain->SetBranchAddress("vbf_maxpt_j2_m_jes_up", &vbf_maxpt_j2_m_jes_up, &b_vbf_maxpt_j2_m_jes_up);
   fChain->SetBranchAddress("vbf_maxpt_j2_pt_jes_up", &vbf_maxpt_j2_pt_jes_up, &b_vbf_maxpt_j2_pt_jes_up);
   fChain->SetBranchAddress("vbf_maxpt_j2_eta_jes_up", &vbf_maxpt_j2_eta_jes_up, &b_vbf_maxpt_j2_eta_jes_up);
   fChain->SetBranchAddress("vbf_maxpt_j2_phi_jes_up", &vbf_maxpt_j2_phi_jes_up, &b_vbf_maxpt_j2_phi_jes_up);
   fChain->SetBranchAddress("vbf_maxpt_j2_m_jes_dn", &vbf_maxpt_j2_m_jes_dn, &b_vbf_maxpt_j2_m_jes_dn);
   fChain->SetBranchAddress("vbf_maxpt_j2_pt_jes_dn", &vbf_maxpt_j2_pt_jes_dn, &b_vbf_maxpt_j2_pt_jes_dn);
   fChain->SetBranchAddress("vbf_maxpt_j2_eta_jes_dn", &vbf_maxpt_j2_eta_jes_dn, &b_vbf_maxpt_j2_eta_jes_dn);
   fChain->SetBranchAddress("vbf_maxpt_j2_phi_jes_dn", &vbf_maxpt_j2_phi_jes_dn, &b_vbf_maxpt_j2_phi_jes_dn);
   fChain->SetBranchAddress("vbf_maxpt_j2_m_jer", &vbf_maxpt_j2_m_jer, &b_vbf_maxpt_j2_m_jer);
   fChain->SetBranchAddress("vbf_maxpt_j2_pt_jer", &vbf_maxpt_j2_pt_jer, &b_vbf_maxpt_j2_pt_jer);
   fChain->SetBranchAddress("vbf_maxpt_j2_eta_jer", &vbf_maxpt_j2_eta_jer, &b_vbf_maxpt_j2_eta_jer);
   fChain->SetBranchAddress("vbf_maxpt_j2_phi_jer", &vbf_maxpt_j2_phi_jer, &b_vbf_maxpt_j2_phi_jer);
   fChain->SetBranchAddress("vbf_maxpt_j2_m_jer_up", &vbf_maxpt_j2_m_jer_up, &b_vbf_maxpt_j2_m_jer_up);
   fChain->SetBranchAddress("vbf_maxpt_j2_pt_jer_up", &vbf_maxpt_j2_pt_jer_up, &b_vbf_maxpt_j2_pt_jer_up);
   fChain->SetBranchAddress("vbf_maxpt_j2_eta_jer_up", &vbf_maxpt_j2_eta_jer_up, &b_vbf_maxpt_j2_eta_jer_up);
   fChain->SetBranchAddress("vbf_maxpt_j2_phi_jer_up", &vbf_maxpt_j2_phi_jer_up, &b_vbf_maxpt_j2_phi_jer_up);
   fChain->SetBranchAddress("vbf_maxpt_j2_m_jer_dn", &vbf_maxpt_j2_m_jer_dn, &b_vbf_maxpt_j2_m_jer_dn);
   fChain->SetBranchAddress("vbf_maxpt_j2_pt_jer_dn", &vbf_maxpt_j2_pt_jer_dn, &b_vbf_maxpt_j2_pt_jer_dn);
   fChain->SetBranchAddress("vbf_maxpt_j2_eta_jer_dn", &vbf_maxpt_j2_eta_jer_dn, &b_vbf_maxpt_j2_eta_jer_dn);
   fChain->SetBranchAddress("vbf_maxpt_j2_phi_jer_dn", &vbf_maxpt_j2_phi_jer_dn, &b_vbf_maxpt_j2_phi_jer_dn);
   fChain->SetBranchAddress("vbf_maxpt_j1_QGLikelihood", &vbf_maxpt_j1_QGLikelihood, &b_vbf_maxpt_j1_QGLikelihood);
   fChain->SetBranchAddress("vbf_maxpt_j2_QGLikelihood", &vbf_maxpt_j2_QGLikelihood, &b_vbf_maxpt_j2_QGLikelihood);
   fChain->SetBranchAddress("vbf_maxpt_j1_isPileUpMedium", &vbf_maxpt_j1_isPileUpMedium, &b_vbf_maxpt_j1_isPileUpMedium);
   fChain->SetBranchAddress("vbf_maxpt_j2_isPileUpMedium", &vbf_maxpt_j2_isPileUpMedium, &b_vbf_maxpt_j2_isPileUpMedium);
   fChain->SetBranchAddress("vbf_maxpt_j1_isPileUpTight", &vbf_maxpt_j1_isPileUpTight, &b_vbf_maxpt_j1_isPileUpTight);
   fChain->SetBranchAddress("vbf_maxpt_j2_isPileUpTight", &vbf_maxpt_j2_isPileUpTight, &b_vbf_maxpt_j2_isPileUpTight);
   fChain->SetBranchAddress("vbf_maxpt_j1_bDiscriminatorCSV", &vbf_maxpt_j1_bDiscriminatorCSV, &b_vbf_maxpt_j1_bDiscriminatorCSV);
   fChain->SetBranchAddress("vbf_maxpt_j2_bDiscriminatorCSV", &vbf_maxpt_j2_bDiscriminatorCSV, &b_vbf_maxpt_j2_bDiscriminatorCSV);
   fChain->SetBranchAddress("vbf_maxpt_j1_bDiscriminatorCSV_gen", &vbf_maxpt_j1_bDiscriminatorCSV_gen, &b_vbf_maxpt_j1_bDiscriminatorCSV_gen);
   fChain->SetBranchAddress("vbf_maxpt_j2_bDiscriminatorCSV_gen", &vbf_maxpt_j2_bDiscriminatorCSV_gen, &b_vbf_maxpt_j2_bDiscriminatorCSV_gen);
   fChain->SetBranchAddress("nbjets_csvl_veto", &nbjets_csvl_veto, &b_nbjets_csvl_veto);
   fChain->SetBranchAddress("nbjets_csvm_veto", &nbjets_csvm_veto, &b_nbjets_csvm_veto);
   fChain->SetBranchAddress("nbjets_csvt_veto", &nbjets_csvt_veto, &b_nbjets_csvt_veto);
   fChain->SetBranchAddress("nbjets_ssvhem_veto", &nbjets_ssvhem_veto, &b_nbjets_ssvhem_veto);
   fChain->SetBranchAddress("nbjets_csvl_veto_cleaned", &nbjets_csvl_veto_cleaned, &b_nbjets_csvl_veto_cleaned);
   fChain->SetBranchAddress("nbjets_csvm_veto_cleaned", &nbjets_csvm_veto_cleaned, &b_nbjets_csvm_veto_cleaned);
   fChain->SetBranchAddress("nbjets_csvt_veto_cleaned", &nbjets_csvt_veto_cleaned, &b_nbjets_csvt_veto_cleaned);
   fChain->SetBranchAddress("nbjets_ssvhem_veto_cleaned", &nbjets_ssvhem_veto_cleaned, &b_nbjets_ssvhem_veto_cleaned);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("deltaR_lca8jet", &deltaR_lca8jet, &b_deltaR_lca8jet);
   fChain->SetBranchAddress("deltaphi_METca8jet", &deltaphi_METca8jet, &b_deltaphi_METca8jet);
   fChain->SetBranchAddress("deltaphi_Vca8jet", &deltaphi_Vca8jet, &b_deltaphi_Vca8jet);
   fChain->SetBranchAddress("deltaphi_METca8jet_met", &deltaphi_METca8jet_met, &b_deltaphi_METca8jet_met);
   fChain->SetBranchAddress("deltaphi_Vca8jet_met", &deltaphi_Vca8jet_met, &b_deltaphi_Vca8jet_met);
   fChain->SetBranchAddress("genHMass", &genHMass, &b_genHMass);
   fChain->SetBranchAddress("genHphi", &genHphi, &b_genHphi);
   fChain->SetBranchAddress("genHeta", &genHeta, &b_genHeta);
   fChain->SetBranchAddress("genHpt", &genHpt, &b_genHpt);
   fChain->SetBranchAddress("genTagQuark1W", &genTagQuark1W, &b_genTagQuark1E);
   fChain->SetBranchAddress("genTagQuark1phi", &genTagQuark1phi, &b_genTagQuark1phi);
   fChain->SetBranchAddress("genTagQuark1eta", &genTagQuark1eta, &b_genTagQuark1eta);
   fChain->SetBranchAddress("genTagQuark1pt", &genTagQuark1pt, &b_genTagQuark1pt);
   fChain->SetBranchAddress("genTagQuarkE", &genTagQuarkE, &b_genTagQuark2E);
   fChain->SetBranchAddress("genTagQuark2phi", &genTagQuark2phi, &b_genTagQuark2phi);
   fChain->SetBranchAddress("genTagQuark2eta", &genTagQuark2eta, &b_genTagQuark2eta);
   fChain->SetBranchAddress("genTagQuark2pt", &genTagQuark2pt, &b_genTagQuark2pt);
   fChain->SetBranchAddress("ttb_nak5_same", &ttb_nak5_same, &b_ttb_nak5_same);
   fChain->SetBranchAddress("ttb_nak5_same_csvl", &ttb_nak5_same_csvl, &b_ttb_nak5_same_csvl);
   fChain->SetBranchAddress("ttb_nak5_same_csvm", &ttb_nak5_same_csvm, &b_ttb_nak5_same_csvm);
   fChain->SetBranchAddress("ttb_nak5_same_csvt", &ttb_nak5_same_csvt, &b_ttb_nak5_same_csvt);
   fChain->SetBranchAddress("ttb_nak5_oppo", &ttb_nak5_oppo, &b_ttb_nak5_oppo);
   fChain->SetBranchAddress("ttb_nak5_oppo_csvl", &ttb_nak5_oppo_csvl, &b_ttb_nak5_oppo_csvl);
   fChain->SetBranchAddress("ttb_nak5_oppo_csvm", &ttb_nak5_oppo_csvm, &b_ttb_nak5_oppo_csvm);
   fChain->SetBranchAddress("ttb_nak5_oppo_csvt", &ttb_nak5_oppo_csvt, &b_ttb_nak5_oppo_csvt);
   fChain->SetBranchAddress("ttb_nak5_oppoveto", &ttb_nak5_oppoveto, &b_ttb_nak5_oppoveto);
   fChain->SetBranchAddress("ttb_nak5_oppoveto_csvl", &ttb_nak5_oppoveto_csvl, &b_ttb_nak5_oppoveto_csvl);
   fChain->SetBranchAddress("ttb_nak5_oppoveto_csvm", &ttb_nak5_oppoveto_csvm, &b_ttb_nak5_oppoveto_csvm);
   fChain->SetBranchAddress("ttb_nak5_oppoveto_csvt", &ttb_nak5_oppoveto_csvt, &b_ttb_nak5_oppoveto_csvt);
   fChain->SetBranchAddress("ttb_nak5_sameveto", &ttb_nak5_sameveto, &b_ttb_nak5_sameveto);
   fChain->SetBranchAddress("ttb_nak5_sameveto_csvl", &ttb_nak5_sameveto_csvl, &b_ttb_nak5_sameveto_csvl);
   fChain->SetBranchAddress("ttb_nak5_sameveto_csvm", &ttb_nak5_sameveto_csvm, &b_ttb_nak5_sameveto_csvm);
   fChain->SetBranchAddress("ttb_nak5_sameveto_csvt", &ttb_nak5_sameveto_csvt, &b_ttb_nak5_sameveto_csvt);
   fChain->SetBranchAddress("ttb_ht", &ttb_ht, &b_ttb_ht);
   fChain->SetBranchAddress("ttb_ca8_mass_pr", &ttb_ca8_mass_pr, &b_ttb_ca8_mass_pr);
   fChain->SetBranchAddress("ttb_ca8_charge", &ttb_ca8_charge, &b_ttb_ca8_charge);
   fChain->SetBranchAddress("ttb_ca8_charge_k05", &ttb_ca8_charge_k05, &b_ttb_ca8_charge_k05);
   fChain->SetBranchAddress("ttb_ca8_charge_k07", &ttb_ca8_charge_k07, &b_ttb_ca8_charge_k07);
   fChain->SetBranchAddress("ttb_ca8_charge_k10", &ttb_ca8_charge_k10, &b_ttb_ca8_charge_k10);
   fChain->SetBranchAddress("ttb_ca8_ungroomed_pt", &ttb_ca8_ungroomed_pt, &b_ttb_ca8_ungroomed_pt);
   fChain->SetBranchAddress("ttb_ca8_ungroomed_eta", &ttb_ca8_ungroomed_eta, &b_ttb_ca8_ungroomed_eta);
   fChain->SetBranchAddress("ttb_ca8_ungroomed_phi", &ttb_ca8_ungroomed_phi, &b_ttb_ca8_ungroomed_phi);
   fChain->SetBranchAddress("ttb_ca8_ungroomed_e", &ttb_ca8_ungroomed_e, &b_ttb_ca8_ungroomed_e);
   fChain->SetBranchAddress("ttb_ca8_ungroomed_gen_pt", &ttb_ca8_ungroomed_gen_pt, &b_ttb_ca8_ungroomed_gen_pt);
   fChain->SetBranchAddress("ttb_ca8_ungroomed_gen_eta", &ttb_ca8_ungroomed_gen_eta, &b_ttb_ca8_ungroomed_gen_eta);
   fChain->SetBranchAddress("ttb_ca8_ungroomed_gen_phi", &ttb_ca8_ungroomed_gen_phi, &b_ttb_ca8_ungroomed_gen_phi);
   fChain->SetBranchAddress("ttb_ca8_ungroomed_gen_e", &ttb_ca8_ungroomed_gen_e, &b_ttb_ca8_ungroomed_gen_e);
   fChain->SetBranchAddress("ttb_ca8_tau2tau1", &ttb_ca8_tau2tau1, &b_ttb_ca8_tau2tau1);
   fChain->SetBranchAddress("ttb_ca8_tau2tau1_exkT", &ttb_ca8_tau2tau1_exkT, &b_ttb_ca8_tau2tau1_exkT);
   fChain->SetBranchAddress("ttb_ca8_tau2tau1_pr", &ttb_ca8_tau2tau1_pr, &b_ttb_ca8_tau2tau1_pr);
   fChain->SetBranchAddress("ttb_ca8_GeneralizedECF", &ttb_ca8_GeneralizedECF, &b_ttb_ca8_GeneralizedECF);
   fChain->SetBranchAddress("ttb_ca8_mu", &ttb_ca8_mu, &b_ttb_ca8_mu);
   fChain->SetBranchAddress("ttb_ca8_mlvj_type0", &ttb_ca8_mlvj_type0, &b_ttb_ca8_mlvj_type0);
   fChain->SetBranchAddress("ttb_ca8_mlvj_type2", &ttb_ca8_mlvj_type2, &b_ttb_ca8_mlvj_type2);
   fChain->SetBranchAddress("ttb_ca8_mlvj_type0_met", &ttb_ca8_mlvj_type0_met, &b_ttb_ca8_mlvj_type0_met);
   fChain->SetBranchAddress("ttb_ca8_mlvj_type2_met", &ttb_ca8_mlvj_type2_met, &b_ttb_ca8_mlvj_type2_met);
   fChain->SetBranchAddress("ttb_dR_ca8_bjet_closer", &ttb_dR_ca8_bjet_closer, &b_ttb_dR_ca8_bjet_closer);
   fChain->SetBranchAddress("ttb_dR_ca8_jet_closer", &ttb_dR_ca8_jet_closer, &b_ttb_dR_ca8_jet_closer);
   fChain->SetBranchAddress("isttbar", &isttbar, &b_isttbar);
   fChain->SetBranchAddress("gen_parton1_px_fromttbar", &gen_parton1_px_fromttbar, &b_gen_parton1_px_fromttbar);
   fChain->SetBranchAddress("gen_parton1_py_fromttbar", &gen_parton1_py_fromttbar, &b_gen_parton1_py_fromttbar);
   fChain->SetBranchAddress("gen_parton1_pz_fromttbar", &gen_parton1_pz_fromttbar, &b_gen_parton1_pz_fromttbar);
   fChain->SetBranchAddress("gen_parton1_e_fromttbar", &gen_parton1_e_fromttbar, &b_gen_parton1_e_fromttbar);
   fChain->SetBranchAddress("gen_parton1_id_fromttbar", &gen_parton1_id_fromttbar, &b_gen_parton1_id_fromttbar);
   fChain->SetBranchAddress("gen_parton2_px_fromttbar", &gen_parton2_px_fromttbar, &b_gen_parton2_px_fromttbar);
   fChain->SetBranchAddress("gen_parton2_py_fromttbar", &gen_parton2_py_fromttbar, &b_gen_parton2_py_fromttbar);
   fChain->SetBranchAddress("gen_parton2_pz_fromttbar", &gen_parton2_pz_fromttbar, &b_gen_parton2_pz_fromttbar);
   fChain->SetBranchAddress("gen_parton2_e_fromttbar", &gen_parton2_e_fromttbar, &b_gen_parton2_e_fromttbar);
   fChain->SetBranchAddress("gen_parton2_id_fromttbar", &gen_parton2_id_fromttbar, &b_gen_parton2_id_fromttbar);
}


