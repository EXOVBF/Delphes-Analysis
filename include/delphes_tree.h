//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 28 10:38:53 2014 by ROOT version 5.34/07
// from TTree Delphes/Analysis tree
// found on file: EWK_0.root
//////////////////////////////////////////////////////////

#ifndef delphes_tree_h
#define delphes_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

using namespace std ;

// Fixed size dimensions of array or collections stored in the TTree if any.

class delphes_tree {
public :
   TTree          *fChain ;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   vector<float>   *lhe_lep_number;
   vector<float>   *lhe_lep_pt;
   vector<float>   *lhe_lep_eta;
   vector<float>   *lhe_lep_phi;
   vector<float>   *lhe_lep_flv;
   vector<float>   *lhe_nu_pt;
   vector<float>   *lhe_nu_eta;
   vector<float>   *lhe_nu_phi;
   vector<float>   *lhe_nu_flv;
   vector<float>   *lhe_p_pt;
   vector<float>   *lhe_p_eta;
   vector<float>   *lhe_p_phi;
   vector<float>   *lhe_p_flv;
   vector<float>   *lhe_p_from_W;
   vector<float>   *lhe_W_pt;
   vector<float>   *lhe_W_eta;
   vector<float>   *lhe_W_phi;
   vector<float>   *lhe_W_mass;
   vector<float>   *lhe_W_pid;
   vector<float>   *lhe_X_pt;
   vector<float>   *lhe_X_eta;
   vector<float>   *lhe_X_phi;
   vector<float>   *lhe_X_mass;
   vector<float>   *lep_pt;
   vector<float>   *lep_eta;
   vector<float>   *lep_phi;
   vector<float>   *lep_flv;
   vector<float>   *lep_isolation;
   vector<float>   *lep_E_no_smearing;
   vector<float>   *lep_E_with_smearing;
   vector<float>   *lep_n_particle_cone;
   vector<float>   *lep_number;
   vector<float>   *lep_X_vertex;
   vector<float>   *lep_Y_vertex;
   vector<float>   *lep_Z_vertex;
   vector<float>   *MET;
   vector<float>   *MET_phi;
   vector<float>   *number_gen_jet_ak5;
   vector<float>   *gen_jet_ak5_pt;
   vector<float>   *gen_jet_ak5_eta;
   vector<float>   *gen_jet_ak5_phi;
   vector<float>   *gen_jet_ak5_mass;
   vector<float>   *gen_jet_ak5_mass_pruned;
   vector<float>   *gen_jet_ak5_btag;
   vector<float>   *gen_jet_ak5_Hcal_over_Ecal;
   vector<float>   *number_jet_ak5;
   vector<float>   *jet_ak5_pt;
   vector<float>   *jet_ak5_eta;
   vector<float>   *jet_ak5_phi;
   vector<float>   *jet_ak5_mass;
   vector<float>   *jet_ak5_mass_pruned;
   vector<float>   *jet_ak5_btag;
   vector<float>   *jet_ak5_Hcal_over_Ecal;
   vector<float>   *number_gen_jet_CA8;
   vector<float>   *gen_jet_CA8_pt;
   vector<float>   *gen_jet_CA8_eta;
   vector<float>   *gen_jet_CA8_phi;
   vector<float>   *gen_jet_CA8_mass;
   vector<float>   *gen_jet_CA8_mass_pruned;
   vector<float>   *gen_jet_CA8_btag;
   vector<float>   *gen_jet_CA8_Hcal_over_Ecal;
   vector<float>   *gen_jet_CA8_tau1;
   vector<float>   *gen_jet_CA8_tau2;
   vector<float>   *gen_jet_CA8_tau3;
   vector<float>   *number_jet_CA8;
   vector<float>   *jet_CA8_pt;
   vector<float>   *jet_CA8_eta;
   vector<float>   *jet_CA8_phi;
   vector<float>   *jet_CA8_mass;
   vector<float>   *jet_CA8_mass_pruned;
   vector<float>   *jet_CA8_btag;
   vector<float>   *jet_CA8_Hcal_over_Ecal;
   vector<float>   *jet_CA8_tau1;
   vector<float>   *jet_CA8_tau2;
   vector<float>   *jet_CA8_tau3;
   vector<float>   *Rho_PU_barrel;
   vector<float>   *Rho_PU_endcap;
   vector<float>   *Rho_PU_foward;
   vector<float>   *nPV;
   vector<float>   *vertex_X;
   vector<float>   *vertex_Y;
   vector<float>   *vertex_Z;
   vector<float>   *sum_pt_square;
   vector<float>   *gen_nPV;
   vector<float>   *gen_is_PU;
   vector<float>   *photon_pt;
   vector<float>   *photon_eta;
   vector<float>   *photon_phi;
   vector<float>   *photon_E;
   vector<float>   *photon_isolation;
   vector<float>   *photon_n_particle_cone;
   vector<float>   *photon_number;
   vector<float>   *photon_Hcal_over_Ecal;

   // List of branches
   TBranch        *b_lhe_lep_number;   //!
   TBranch        *b_lhe_lep_pt;   //!
   TBranch        *b_lhe_lep_eta;   //!
   TBranch        *b_lhe_lep_phi;   //!
   TBranch        *b_lhe_lep_flv;   //!
   TBranch        *b_lhe_nu_pt;   //!
   TBranch        *b_lhe_nu_eta;   //!
   TBranch        *b_lhe_nu_phi;   //!
   TBranch        *b_lhe_nu_flv;   //!
   TBranch        *b_lhe_p_pt;   //!
   TBranch        *b_lhe_p_eta;   //!
   TBranch        *b_lhe_p_phi;   //!
   TBranch        *b_lhe_p_flv;   //!
   TBranch        *b_lhe_p_from_W;   //!
   TBranch        *b_lhe_W_pt;   //!
   TBranch        *b_lhe_W_eta;   //!
   TBranch        *b_lhe_W_phi;   //!
   TBranch        *b_lhe_W_mass;   //!
   TBranch        *b_lhe_W_pid;   //!
   TBranch        *b_lhe_X_pt;   //!
   TBranch        *b_lhe_X_eta;   //!
   TBranch        *b_lhe_X_phi;   //!
   TBranch        *b_lhe_X_mass;   //!
   TBranch        *b_lep_pt;   //!
   TBranch        *b_lep_eta;   //!
   TBranch        *b_lep_phi;   //!
   TBranch        *b_lep_flv;   //!
   TBranch        *b_lep_isolation;   //!
   TBranch        *b_lep_E_no_smearing;   //!
   TBranch        *b_lep_E_with_smearing;   //!
   TBranch        *b_lep_n_particle_cone;   //!
   TBranch        *b_lep_number;   //!
   TBranch        *b_lep_X_vertex;   //!
   TBranch        *b_lep_Y_vertex;   //!
   TBranch        *b_lep_Z_vertex;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_number_gen_jet_ak5;   //!
   TBranch        *b_gen_jet_ak5_pt;   //!
   TBranch        *b_gen_jet_ak5_eta;   //!
   TBranch        *b_gen_jet_ak5_phi;   //!
   TBranch        *b_gen_jet_ak5_mass;   //!
   TBranch        *b_gen_jet_ak5_mass_pruned;   //!
   TBranch        *b_gen_jet_ak5_btag;   //!
   TBranch        *b_gen_jet_ak5_Hcal_over_Ecal;   //!
   TBranch        *b_number_jet_ak5;   //!
   TBranch        *b_jet_ak5_pt;   //!
   TBranch        *b_jet_ak5_eta;   //!
   TBranch        *b_jet_ak5_phi;   //!
   TBranch        *b_jet_ak5_mass;   //!
   TBranch        *b_jet_ak5_mass_pruned;   //!
   TBranch        *b_jet_ak5_btag;   //!
   TBranch        *b_jet_ak5_Hcal_over_Ecal;   //!
   TBranch        *b_number_gen_jet_CA8;   //!
   TBranch        *b_gen_jet_CA8_pt;   //!
   TBranch        *b_gen_jet_CA8_eta;   //!
   TBranch        *b_gen_jet_CA8_phi;   //!
   TBranch        *b_gen_jet_CA8_mass;   //!
   TBranch        *b_gen_jet_CA8_mass_pruned;   //!
   TBranch        *b_gen_jet_CA8_btag;   //!
   TBranch        *b_gen_jet_CA8_Hcal_over_Ecal;   //!
   TBranch        *b_gen_jet_CA8_tau1;   //!
   TBranch        *b_gen_jet_CA8_tau2;   //!
   TBranch        *b_gen_jet_CA8_tau3;   //!
   TBranch        *b_number_jet_CA8;   //!
   TBranch        *b_jet_CA8_pt;   //!
   TBranch        *b_jet_CA8_eta;   //!
   TBranch        *b_jet_CA8_phi;   //!
   TBranch        *b_jet_CA8_mass;   //!
   TBranch        *b_jet_CA8_mass_pruned;   //!
   TBranch        *b_jet_CA8_btag;   //!
   TBranch        *b_jet_CA8_Hcal_over_Ecal;   //!
   TBranch        *b_jet_CA8_tau1;   //!
   TBranch        *b_jet_CA8_tau2;   //!
   TBranch        *b_jet_CA8_tau3;   //!
   TBranch        *b_Rho_PU_barrel;   //!
   TBranch        *b_Rho_PU_endcap;   //!
   TBranch        *b_Rho_PU_foward;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_vertex_X;   //!
   TBranch        *b_vertex_Y;   //!
   TBranch        *b_vertex_Z;   //!
   TBranch        *b_sum_pt_square;   //!
   TBranch        *b_gen_nPV;   //!
   TBranch        *b_gen_is_PU;   //!
   TBranch        *b_photon_pt;   //!
   TBranch        *b_photon_eta;   //!
   TBranch        *b_photon_phi;   //!
   TBranch        *b_photon_E;   //!
   TBranch        *b_photon_isolation;   //!
   TBranch        *b_photon_n_particle_cone;   //!
   TBranch        *b_photon_number;   //!
   TBranch        *b_photon_Hcal_over_Ecal;   //!

   delphes_tree (TTree * /*tree*/ = 0, bool signal = false) : fChain (0) { }
   virtual ~delphes_tree () { }
   virtual void    Init (TTree *tree, bool signal = false) ;
   virtual Int_t   GetEntry (Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree ()->GetEntry (entry, getall) : 0 ; }

} ;

void delphes_tree::Init (TTree * tree, bool signal)
{

   // Set object pointer
   // Set object pointer
   lhe_lep_number = 0;
   lhe_lep_pt = 0;
   lhe_lep_eta = 0;
   lhe_lep_phi = 0;
   lhe_lep_flv = 0;
   lhe_nu_pt = 0;
   lhe_nu_eta = 0;
   lhe_nu_phi = 0;
   lhe_nu_flv = 0;
   lhe_p_pt = 0;
   lhe_p_eta = 0;
   lhe_p_phi = 0;
   lhe_p_flv = 0;
   lhe_p_from_W = 0;
   lhe_W_pt = 0;
   lhe_W_eta = 0;
   lhe_W_phi = 0;
   lhe_W_mass = 0;
   lhe_W_pid = 0;
   if( signal )
   {
      lhe_X_pt = 0;
      lhe_X_eta = 0;
      lhe_X_phi = 0;
      lhe_X_mass = 0;
   }
   lep_pt = 0;
   lep_eta = 0;
   lep_phi = 0;
   lep_flv = 0;
   lep_isolation = 0;
   lep_E_no_smearing = 0;
   lep_E_with_smearing = 0;
   lep_n_particle_cone = 0;
   lep_number = 0;
   lep_X_vertex = 0;
   lep_Y_vertex = 0;
   lep_Z_vertex = 0;
   MET = 0;
   MET_phi = 0;
   number_gen_jet_ak5 = 0;
   gen_jet_ak5_pt = 0;
   gen_jet_ak5_eta = 0;
   gen_jet_ak5_phi = 0;
   gen_jet_ak5_mass = 0;
   gen_jet_ak5_mass_pruned = 0;
   gen_jet_ak5_btag = 0;
   gen_jet_ak5_Hcal_over_Ecal = 0;
   number_jet_ak5 = 0;
   jet_ak5_pt = 0;
   jet_ak5_eta = 0;
   jet_ak5_phi = 0;
   jet_ak5_mass = 0;
   jet_ak5_mass_pruned = 0;
   jet_ak5_btag = 0;
   jet_ak5_Hcal_over_Ecal = 0;
   number_gen_jet_CA8 = 0;
   gen_jet_CA8_pt = 0;
   gen_jet_CA8_eta = 0;
   gen_jet_CA8_phi = 0;
   gen_jet_CA8_mass = 0;
   gen_jet_CA8_mass_pruned = 0;
   gen_jet_CA8_btag = 0;
   gen_jet_CA8_Hcal_over_Ecal = 0;
   gen_jet_CA8_tau1 = 0;
   gen_jet_CA8_tau2 = 0;
   gen_jet_CA8_tau3 = 0;
   number_jet_CA8 = 0;
   jet_CA8_pt = 0;
   jet_CA8_eta = 0;
   jet_CA8_phi = 0;
   jet_CA8_mass = 0;
   jet_CA8_mass_pruned = 0;
   jet_CA8_btag = 0;
   jet_CA8_Hcal_over_Ecal = 0;
   jet_CA8_tau1 = 0;
   jet_CA8_tau2 = 0;
   jet_CA8_tau3 = 0;
   Rho_PU_barrel = 0;
   Rho_PU_endcap = 0;
   Rho_PU_foward = 0;
   nPV = 0;
   vertex_X = 0;
   vertex_Y = 0;
   vertex_Z = 0;
   sum_pt_square = 0;
   gen_nPV = 0;
   gen_is_PU = 0;
   photon_pt = 0;
   photon_eta = 0;
   photon_phi = 0;
   photon_E = 0;
   photon_isolation = 0;
   photon_n_particle_cone = 0;
   photon_number = 0;
   photon_Hcal_over_Ecal = 0;

   // Set branch addresses and branch pointers
   if (!tree) return ;
   fChain = tree ;
   fChain->SetMakeClass (1) ;

   fChain->SetBranchAddress("lhe_lep_number", &lhe_lep_number, &b_lhe_lep_number);
   fChain->SetBranchAddress("lhe_lep_pt", &lhe_lep_pt, &b_lhe_lep_pt);
   fChain->SetBranchAddress("lhe_lep_eta", &lhe_lep_eta, &b_lhe_lep_eta);
   fChain->SetBranchAddress("lhe_lep_phi", &lhe_lep_phi, &b_lhe_lep_phi);
   fChain->SetBranchAddress("lhe_lep_flv", &lhe_lep_flv, &b_lhe_lep_flv);
   fChain->SetBranchAddress("lhe_nu_pt", &lhe_nu_pt, &b_lhe_nu_pt);
   fChain->SetBranchAddress("lhe_nu_eta", &lhe_nu_eta, &b_lhe_nu_eta);
   fChain->SetBranchAddress("lhe_nu_phi", &lhe_nu_phi, &b_lhe_nu_phi);
   fChain->SetBranchAddress("lhe_nu_flv", &lhe_nu_flv, &b_lhe_nu_flv);
   fChain->SetBranchAddress("lhe_p_pt", &lhe_p_pt, &b_lhe_p_pt);
   fChain->SetBranchAddress("lhe_p_eta", &lhe_p_eta, &b_lhe_p_eta);
   fChain->SetBranchAddress("lhe_p_phi", &lhe_p_phi, &b_lhe_p_phi);
   fChain->SetBranchAddress("lhe_p_flv", &lhe_p_flv, &b_lhe_p_flv);
   fChain->SetBranchAddress("lhe_p_from_W", &lhe_p_from_W, &b_lhe_p_from_W);
   fChain->SetBranchAddress("lhe_W_pt", &lhe_W_pt, &b_lhe_W_pt);
   fChain->SetBranchAddress("lhe_W_eta", &lhe_W_eta, &b_lhe_W_eta);
   fChain->SetBranchAddress("lhe_W_phi", &lhe_W_phi, &b_lhe_W_phi);
   fChain->SetBranchAddress("lhe_W_mass", &lhe_W_mass, &b_lhe_W_mass);
   fChain->SetBranchAddress("lhe_W_pid", &lhe_W_pid, &b_lhe_W_pid);
   if( signal )
   {
      fChain->SetBranchAddress("lhe_X_pt", &lhe_X_pt, &b_lhe_X_pt);
      fChain->SetBranchAddress("lhe_X_eta", &lhe_X_eta, &b_lhe_X_eta);
      fChain->SetBranchAddress("lhe_X_phi", &lhe_X_phi, &b_lhe_X_phi);
      fChain->SetBranchAddress("lhe_X_mass", &lhe_X_mass, &b_lhe_X_mass);
   }
   fChain->SetBranchAddress("lep_pt", &lep_pt, &b_lep_pt);
   fChain->SetBranchAddress("lep_eta", &lep_eta, &b_lep_eta);
   fChain->SetBranchAddress("lep_phi", &lep_phi, &b_lep_phi);
   fChain->SetBranchAddress("lep_flv", &lep_flv, &b_lep_flv);
   fChain->SetBranchAddress("lep_isolation", &lep_isolation, &b_lep_isolation);
   fChain->SetBranchAddress("lep_E_no_smearing", &lep_E_no_smearing, &b_lep_E_no_smearing);
   fChain->SetBranchAddress("lep_E_with_smearing", &lep_E_with_smearing, &b_lep_E_with_smearing);
   fChain->SetBranchAddress("lep_n_particle_cone", &lep_n_particle_cone, &b_lep_n_particle_cone);
   fChain->SetBranchAddress("lep_number", &lep_number, &b_lep_number);
   fChain->SetBranchAddress("lep_X_vertex", &lep_X_vertex, &b_lep_X_vertex);
   fChain->SetBranchAddress("lep_Y_vertex", &lep_Y_vertex, &b_lep_Y_vertex);
   fChain->SetBranchAddress("lep_Z_vertex", &lep_Z_vertex, &b_lep_Z_vertex);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("number_gen_jet_ak5", &number_gen_jet_ak5, &b_number_gen_jet_ak5);
   fChain->SetBranchAddress("gen_jet_ak5_pt", &gen_jet_ak5_pt, &b_gen_jet_ak5_pt);
   fChain->SetBranchAddress("gen_jet_ak5_eta", &gen_jet_ak5_eta, &b_gen_jet_ak5_eta);
   fChain->SetBranchAddress("gen_jet_ak5_phi", &gen_jet_ak5_phi, &b_gen_jet_ak5_phi);
   fChain->SetBranchAddress("gen_jet_ak5_mass", &gen_jet_ak5_mass, &b_gen_jet_ak5_mass);
   fChain->SetBranchAddress("gen_jet_ak5_mass_pruned", &gen_jet_ak5_mass_pruned, &b_gen_jet_ak5_mass_pruned);
   fChain->SetBranchAddress("gen_jet_ak5_btag", &gen_jet_ak5_btag, &b_gen_jet_ak5_btag);
   fChain->SetBranchAddress("gen_jet_ak5_Hcal_over_Ecal", &gen_jet_ak5_Hcal_over_Ecal, &b_gen_jet_ak5_Hcal_over_Ecal);
   fChain->SetBranchAddress("number_jet_ak5", &number_jet_ak5, &b_number_jet_ak5);
   fChain->SetBranchAddress("jet_ak5_pt", &jet_ak5_pt, &b_jet_ak5_pt);
   fChain->SetBranchAddress("jet_ak5_eta", &jet_ak5_eta, &b_jet_ak5_eta);
   fChain->SetBranchAddress("jet_ak5_phi", &jet_ak5_phi, &b_jet_ak5_phi);
   fChain->SetBranchAddress("jet_ak5_mass", &jet_ak5_mass, &b_jet_ak5_mass);
   fChain->SetBranchAddress("jet_ak5_mass_pruned", &jet_ak5_mass_pruned, &b_jet_ak5_mass_pruned);
   fChain->SetBranchAddress("jet_ak5_btag", &jet_ak5_btag, &b_jet_ak5_btag);
   fChain->SetBranchAddress("jet_ak5_Hcal_over_Ecal", &jet_ak5_Hcal_over_Ecal, &b_jet_ak5_Hcal_over_Ecal);
   fChain->SetBranchAddress("number_gen_jet_CA8", &number_gen_jet_CA8, &b_number_gen_jet_CA8);
   fChain->SetBranchAddress("gen_jet_CA8_pt", &gen_jet_CA8_pt, &b_gen_jet_CA8_pt);
   fChain->SetBranchAddress("gen_jet_CA8_eta", &gen_jet_CA8_eta, &b_gen_jet_CA8_eta);
   fChain->SetBranchAddress("gen_jet_CA8_phi", &gen_jet_CA8_phi, &b_gen_jet_CA8_phi);
   fChain->SetBranchAddress("gen_jet_CA8_mass", &gen_jet_CA8_mass, &b_gen_jet_CA8_mass);
   fChain->SetBranchAddress("gen_jet_CA8_mass_pruned", &gen_jet_CA8_mass_pruned, &b_gen_jet_CA8_mass_pruned);
   fChain->SetBranchAddress("gen_jet_CA8_btag", &gen_jet_CA8_btag, &b_gen_jet_CA8_btag);
   fChain->SetBranchAddress("gen_jet_CA8_Hcal_over_Ecal", &gen_jet_CA8_Hcal_over_Ecal, &b_gen_jet_CA8_Hcal_over_Ecal);
   fChain->SetBranchAddress("gen_jet_CA8_tau1", &gen_jet_CA8_tau1, &b_gen_jet_CA8_tau1);
   fChain->SetBranchAddress("gen_jet_CA8_tau2", &gen_jet_CA8_tau2, &b_gen_jet_CA8_tau2);
   fChain->SetBranchAddress("gen_jet_CA8_tau3", &gen_jet_CA8_tau3, &b_gen_jet_CA8_tau3);
   fChain->SetBranchAddress("number_jet_CA8", &number_jet_CA8, &b_number_jet_CA8);
   fChain->SetBranchAddress("jet_CA8_pt", &jet_CA8_pt, &b_jet_CA8_pt);
   fChain->SetBranchAddress("jet_CA8_eta", &jet_CA8_eta, &b_jet_CA8_eta);
   fChain->SetBranchAddress("jet_CA8_phi", &jet_CA8_phi, &b_jet_CA8_phi);
   fChain->SetBranchAddress("jet_CA8_mass", &jet_CA8_mass, &b_jet_CA8_mass);
   fChain->SetBranchAddress("jet_CA8_mass_pruned", &jet_CA8_mass_pruned, &b_jet_CA8_mass_pruned);
   fChain->SetBranchAddress("jet_CA8_btag", &jet_CA8_btag, &b_jet_CA8_btag);
   fChain->SetBranchAddress("jet_CA8_Hcal_over_Ecal", &jet_CA8_Hcal_over_Ecal, &b_jet_CA8_Hcal_over_Ecal);
   fChain->SetBranchAddress("jet_CA8_tau1", &jet_CA8_tau1, &b_jet_CA8_tau1);
   fChain->SetBranchAddress("jet_CA8_tau2", &jet_CA8_tau2, &b_jet_CA8_tau2);
   fChain->SetBranchAddress("jet_CA8_tau3", &jet_CA8_tau3, &b_jet_CA8_tau3);
   fChain->SetBranchAddress("Rho_PU_barrel", &Rho_PU_barrel, &b_Rho_PU_barrel);
   fChain->SetBranchAddress("Rho_PU_endcap", &Rho_PU_endcap, &b_Rho_PU_endcap);
   fChain->SetBranchAddress("Rho_PU_foward", &Rho_PU_foward, &b_Rho_PU_foward);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("vertex_X", &vertex_X, &b_vertex_X);
   fChain->SetBranchAddress("vertex_Y", &vertex_Y, &b_vertex_Y);
   fChain->SetBranchAddress("vertex_Z", &vertex_Z, &b_vertex_Z);
   fChain->SetBranchAddress("sum_pt_square", &sum_pt_square, &b_sum_pt_square);
   fChain->SetBranchAddress("gen_nPV", &gen_nPV, &b_gen_nPV);
   fChain->SetBranchAddress("gen_is_PU", &gen_is_PU, &b_gen_is_PU);
   fChain->SetBranchAddress("photon_pt", &photon_pt, &b_photon_pt);
   fChain->SetBranchAddress("photon_eta", &photon_eta, &b_photon_eta);
   fChain->SetBranchAddress("photon_phi", &photon_phi, &b_photon_phi);
   fChain->SetBranchAddress("photon_E", &photon_E, &b_photon_E);
   fChain->SetBranchAddress("photon_isolation", &photon_isolation, &b_photon_isolation);
   fChain->SetBranchAddress("photon_n_particle_cone", &photon_n_particle_cone, &b_photon_n_particle_cone);
   fChain->SetBranchAddress("photon_number", &photon_number, &b_photon_number);
   fChain->SetBranchAddress("photon_Hcal_over_Ecal", &photon_Hcal_over_Ecal, &b_photon_Hcal_over_Ecal);
}


#endif 
