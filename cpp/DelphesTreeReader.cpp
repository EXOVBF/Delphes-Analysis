/****************************************************************************************
    - this program reads the output from delphes and store the interesting parameters in 
	  a new tree, it also apllies some preselection cuts
 
    - before compiling ---> source ../Decay/setup_slc6.sh

    - compile with ---> c++ -O2 -lm `root-config --cflags --glibs` -L /afs/cern.ch/user/s/spigazzi/work/DelphesStuff/delphes_code/ -I /afs/cern.ch/user/s/spigazzi/work/DelphesStuff/delphes_code/ -lDelphes -o DelphesTreeReader DelphesTreeReader.cpp

    - instead of lbrianza/work/SIMONE/delphes_code/ put your Delphes folder

*****************************************************************************************/

#include <iostream>
#include <fstream>
#include <utility>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"

#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include "classes/DelphesClasses.h"

#include "external/ExRootAnalysis/ExRootTreeReader.h"

using namespace std;

//****************************************************************************************

struct eventType
{
    vector<long int> eventID;
    vector<int> lepFlavor;
};	

//****************************************************************************************
// runs over all entries in the delphes tree and applies some preselection cuts
// -> semileptonic final state.

eventType event_preselector(ExRootTreeReader *delphesTree, TClonesArray* branchEl ,TClonesArray* branchMu, TClonesArray* branchMET)
{	
    cout << endl << "################## EVENTS PRESELECTION ##################" << endl;
	
    eventType goodEvent;

    int iEvent=0;
    for(iEvent = 0; iEvent < delphesTree->GetEntries(); iEvent++)
    {
	if (iEvent % 10000 == 0)
	{
	    cout << "iEvent =   " << iEvent << endl;
	}
	delphesTree -> ReadEntry(iEvent);
	if(branchEl->GetEntriesFast() == 1 && branchMu->GetEntriesFast() == 0)
	{
	    Electron* el = (Electron*) branchEl->At(0);
	    //MissingET* met = (MissingET*) branchMET->At(0); 
	    if(el->PT > 27)// && met->MET > 65)
	    {
		goodEvent.eventID.push_back(iEvent);
		goodEvent.lepFlavor.push_back(11);
	    }
	}
	if(branchMu->GetEntriesFast() == 1 && branchEl->GetEntriesFast() == 0)
	{
	    Muon* mu = (Muon*) branchMu->At(0);
	    //MissingET* met = (MissingET*) branchMET->At(0); 
	    if(mu->PT > 24)// &&  met->MET > 50)
	    {
		goodEvent.eventID.push_back(iEvent);
		goodEvent.lepFlavor.push_back(13);
	    }
	}
    }
    cout << "######### events from delphes:   " << iEvent << endl;
    cout << "######### events after preselection cuts:   " << goodEvent.eventID.size() << endl;
    return goodEvent;

}

//****************************************************************************************
// main

int main (int argc, char *argv[])
{
//----------------------------------------------------------------------------------------
//importing delphes libraries

    gSystem->Load("libDelphes");

//----------------------------------------------------------------------------------------
//complex object definitions	
		
    vector<string> inputFiles;
    TChain* delphesNtuples = new TChain("Delphes");
    ExRootTreeReader *delphesTree = new ExRootTreeReader(delphesNtuples);	
    TFile* outputFile = TFile::Open(argv[2],"recreate");
    TTree* easyTree = new TTree("easyDelphes","easyDelphes");

//----------------------------------------------------------------------------------------
// reading input files
	
    if(argc < 3)
    {
	cout << "ERROR: not enough info provided" << endl;
	return 0;
    }
    
    ifstream inputList (argv[1], ios::in);
    string buffer;
    while(inputList >> buffer)
    {		
	inputFiles.push_back(buffer);
	cout << "####### Input file #" << inputFiles.size() << ":   " << inputFiles.back() << endl; 
    }
    //--------- filling the TChain      
    for (int iFiles = 0; iFiles < (int)inputFiles.size(); iFiles++)
    {
	delphesNtuples -> Add((inputFiles.at(iFiles)).c_str());
    }
    delphesNtuples -> BranchRef();

//----------------------------------------------------------------------------------------
//variable management
//
// -> jets: jets of both type in the output tree are stored in decreasing pt order.
// -> leptons: the isolation variable is defined as pt_sum(al particle in a cone)/pt_lep 
	
    //--------- getting objects from the delphes tree
    TClonesArray* branchJet_ak5 = delphesTree->UseBranch("Jet_ak5");
    TClonesArray* branchJet_CA8 = delphesTree->UseBranch("Jet_CA8");
    TClonesArray* branchGenJet_ak5 = delphesTree->UseBranch("GenJet_ak5");
    TClonesArray* branchGenJet_CA8 = delphesTree->UseBranch("GenJet_CA8");	
    TClonesArray* branchEl = delphesTree->UseBranch("Electron");
    TClonesArray* branchMu = delphesTree->UseBranch("Muon");
    TClonesArray* branchMET = delphesTree->UseBranch("MissingET");
    TClonesArray* branchGenPart = delphesTree->UseBranch("Particle");
    TClonesArray* branchRho = delphesTree->UseBranch("Rho");
    TClonesArray* branchVertex = delphesTree->UseBranch("Vertex");
    //--------- creating branches for the new (light) tree
    //--------- lepton
    float lep_pt_tmp=0,lep_eta_tmp=0,lep_phi_tmp=0;
    float lep_flv_tmp=0,lep_isolation_tmp=0,lep_Pin_tmp=0,lep_Pout_tmp=0;
    easyTree -> Branch("lep_pt",&lep_pt_tmp,"lep_pt/F");
    easyTree -> Branch("lep_eta",&lep_eta_tmp,"lep_eta/F");
    easyTree -> Branch("lep_phi",&lep_phi_tmp,"lep_phi/F");
    easyTree -> Branch("lep_flavor",&lep_flv_tmp,"lep_flavor/F");
    easyTree -> Branch("lep_isolation",&lep_isolation_tmp,"lep_isolation/F");
    easyTree -> Branch("lep_p_no_smearing",&lep_Pin_tmp,"lep_p_no_smearing/F");
    easyTree -> Branch("lep_p_with_smearing",&lep_Pout_tmp,"lep_p_with_smearing/F");
    //--------- gen lepton
    float gen_lep_pt_tmp=0,gen_lep_eta_tmp=0,gen_lep_phi_tmp=0,gen_lep_flv_tmp=0;
    int gen_wasTau_tmp = 0;
    easyTree -> Branch("gen_lep_pt",&gen_lep_pt_tmp,"gen_lep_pt/F");
    easyTree -> Branch("gen_lep_eta",&gen_lep_eta_tmp,"gen_lep_eta/F");
    easyTree -> Branch("gen_lep_phi",&gen_lep_phi_tmp,"gen_lep_phi/F");
    easyTree -> Branch("gen_lep_flavor",&gen_lep_flv_tmp,"gen_lep_flavor/F");
    easyTree -> Branch("gen_wasTau",&gen_wasTau_tmp,"gen_wasTau/I");
    //--------- MET
    float MET_tmp=0,MET_phi_tmp=0;	
    easyTree -> Branch("MET",&MET_tmp,"MET/F");
    easyTree -> Branch("MET_phi",&MET_phi_tmp,"MET/F");
    //--------- gen MET (v)
    float gen_MET_tmp=0,gen_MET_eta_tmp=0,gen_MET_phi_tmp=0;	
    easyTree -> Branch("gen_MET",&gen_MET_tmp,"gen_MET/F");
    easyTree -> Branch("gen_MET_eta",&gen_MET_eta_tmp,"gen_MET/F");
    easyTree -> Branch("gen_MET_phi",&gen_MET_phi_tmp,"gen_MET/F");
    //--------- PileUp
    float Rho_tmp=0,Rho_fw_tmp=0,nPV_tmp=0,gen_nPV_tmp=0;
    easyTree -> Branch("RhoPU",&Rho_tmp,"RhoPU/F");
    easyTree -> Branch("RhoPU_foward",&Rho_fw_tmp,"RhoPU_foward/F");
    easyTree -> Branch("nPV",&nPV_tmp,"nPV/F");
    easyTree -> Branch("gen_nPV",&gen_nPV_tmp,"gen_nPV/F");
    //--------- ak5 jets
    vector<float> ak5_pt_tmp, ak5_eta_tmp, ak5_phi_tmp, ak5_m_tmp;
    vector<float> ak5_HC_EC_tmp, ak5_m_pruned_tmp;
    vector<bool> ak5_btag_tmp;
    int numberJetak5_tmp=0; 
    easyTree -> Branch("ak5_jet_pt","vector<float>",&ak5_pt_tmp);
    easyTree -> Branch("ak5_jet_eta","vector<float>",&ak5_eta_tmp);
    easyTree -> Branch("ak5_jet_phi","vector<float>",&ak5_phi_tmp);
    easyTree -> Branch("ak5_jet_m","vector<float>",&ak5_m_tmp);
    easyTree -> Branch("ak5_jet_ratio_hcal_ecal","vector<float>",&ak5_HC_EC_tmp);
    easyTree -> Branch("ak5_jet_m_pruned","vector<float>",&ak5_m_pruned_tmp);
    easyTree -> Branch("ak5_jet_btag","vector<bool>",&ak5_btag_tmp);
    easyTree -> Branch("numberJet_ak5",&numberJetak5_tmp,"numberJet_ak5/I");
    //--------- gen ak5 jets
    vector<float> gen_ak5_pt_tmp, gen_ak5_eta_tmp, gen_ak5_phi_tmp, gen_ak5_m_tmp;
    vector<float> gen_ak5_HC_EC_tmp, gen_ak5_m_pruned_tmp;
    vector<bool> gen_ak5_btag_tmp;
    int gen_numberJetak5_tmp=0; 
    easyTree -> Branch("gen_ak5_jet_pt","vector<float>",&gen_ak5_pt_tmp);
    easyTree -> Branch("gen_ak5_jet_eta","vector<float>",&gen_ak5_eta_tmp);
    easyTree -> Branch("gen_ak5_jet_phi","vector<float>",&gen_ak5_phi_tmp);
    easyTree -> Branch("gen_ak5_jet_m","vector<float>",&gen_ak5_m_tmp);
    easyTree -> Branch("gen_ak5_jet_ratio_hcal_ecal","vector<float>",&gen_ak5_HC_EC_tmp);
    easyTree -> Branch("gen_ak5_jet_m_pruned","vector<float>",&gen_ak5_m_pruned_tmp);
    easyTree -> Branch("gen_ak5_jet_btag","vector<bool>",&gen_ak5_btag_tmp);
    easyTree -> Branch("gen_numberJet_ak5",&gen_numberJetak5_tmp,"gen_numberJet_ak5/I");
    //--------- CA8 jets	
    vector<float> CA8_pt_tmp, CA8_eta_tmp, CA8_phi_tmp, CA8_m_tmp;
    vector<float> CA8_HC_EC_tmp, CA8_m_pruned_tmp;
    vector<float> CA8_tau1_tmp,CA8_tau2_tmp,CA8_tau3_tmp;
    vector<bool> CA8_btag_tmp;
    int numberJetCA8_tmp=0; 
    easyTree -> Branch("CA8_jet_pt","vector<float>",&CA8_pt_tmp);
    easyTree -> Branch("CA8_jet_eta","vector<float>",&CA8_eta_tmp);
    easyTree -> Branch("CA8_jet_phi","vector<float>",&CA8_phi_tmp);
    easyTree -> Branch("CA8_jet_m","vector<float>",&CA8_m_tmp);
    easyTree -> Branch("CA8_jet_ratio_hcal_ecal","vector<float>",&CA8_HC_EC_tmp);
    easyTree -> Branch("CA8_jet_m_pruned","vector<float>",&CA8_m_pruned_tmp);
    easyTree -> Branch("CA8_jet_tau1","vector<float>",&CA8_tau1_tmp);
    easyTree -> Branch("CA8_jet_tau2","vector<float>",&CA8_tau2_tmp);
    easyTree -> Branch("CA8_jet_tau3","vector<float>",&CA8_tau3_tmp);
    easyTree -> Branch("CA8_jet_btag","vector<bool>",&CA8_btag_tmp);
    easyTree -> Branch("numberJet_CA8",&numberJetCA8_tmp,"numberJet_CA8/I");
    //--------- gen CA8 jets
    vector<float> gen_CA8_pt_tmp, gen_CA8_eta_tmp, gen_CA8_phi_tmp, gen_CA8_m_tmp;
    vector<float> gen_CA8_HC_EC_tmp, gen_CA8_m_pruned_tmp;
    vector<float> gen_CA8_tau1_tmp,gen_CA8_tau2_tmp,gen_CA8_tau3_tmp;
    vector<bool> gen_CA8_btag_tmp;
    int gen_numberJetCA8_tmp=0; 
    easyTree -> Branch("gen_CA8_jet_pt","vector<float>",&gen_CA8_pt_tmp);
    easyTree -> Branch("gen_CA8_jet_eta","vector<float>",&gen_CA8_eta_tmp);
    easyTree -> Branch("gen_CA8_jet_phi","vector<float>",&gen_CA8_phi_tmp);
    easyTree -> Branch("gen_CA8_jet_m","vector<float>",&gen_CA8_m_tmp);
    easyTree -> Branch("gen_CA8_jet_ratio_hcal_ecal","vector<float>",&gen_CA8_HC_EC_tmp);
    easyTree -> Branch("gen_CA8_jet_m_pruned","vector<float>",&gen_CA8_m_pruned_tmp);
    easyTree -> Branch("gen_CA8_jet_tau1","vector<float>",&gen_CA8_tau1_tmp);
    easyTree -> Branch("gen_CA8_jet_tau2","vector<float>",&gen_CA8_tau2_tmp);
    easyTree -> Branch("gen_CA8_jet_tau3","vector<float>",&gen_CA8_tau3_tmp);
    easyTree -> Branch("gen_CA8_jet_btag","vector<bool>",&gen_CA8_btag_tmp);
    easyTree -> Branch("gen_numberJet_CA8",&gen_numberJetCA8_tmp,"gen_numberJet_CA8/I");

//----------------------------------------------------------------------------------------
//filling the new (plain) tree
			
    eventType goodEvent = event_preselector(delphesTree,branchEl,branchMu,branchMET);
	
    cout << endl << "################# TREE CREATION STARTED #################" << endl;		

    for(int iEvent = 0; iEvent < (long int)goodEvent.eventID.size(); iEvent++)
    {
	if (iEvent % 1000 == 0)
	{
	    cout << "iEvent =   " << iEvent << endl;
	}
	delphesTree -> ReadEntry(goodEvent.eventID.at(iEvent));
	//--------- lepton branches
	if(goodEvent.lepFlavor.at(iEvent) == 11)
	{
	    //--------- lep = el
	    Electron* el = (Electron*) branchEl->At(0);	
	    lep_pt_tmp = el->PT;
	    lep_eta_tmp = el->Eta;
	    lep_phi_tmp = el->Phi;
	    lep_isolation_tmp = el->Isolation;
	    lep_Pin_tmp = el->P_in;
	    lep_Pout_tmp = el->P_out;
	    lep_flv_tmp = 11;
	    //--------- gen lep
	    GenParticle* genPart = (GenParticle*)(el->Particle).GetObject();	
	    gen_lep_pt_tmp = genPart->PT;
	    gen_lep_eta_tmp = genPart->Eta;
	    gen_lep_phi_tmp = genPart->Phi;
	    gen_lep_flv_tmp = TMath::Abs(genPart->PID); 
	}
	else
	{
	    //--------- lep = mu
	    Muon* mu = (Muon*) branchMu->At(0);
	    lep_pt_tmp = mu->PT;
	    lep_eta_tmp = mu->Eta;
	    lep_phi_tmp = mu->Phi;
	    lep_isolation_tmp = mu->Isolation;
	    lep_Pin_tmp = mu->P_in;
	    lep_Pout_tmp = mu->P_out;
	    lep_flv_tmp = 13;
	    //--------- gen lep
	    GenParticle* genPart = (GenParticle*)(mu->Particle).GetObject();	
	    gen_lep_pt_tmp = genPart->PT;
	    gen_lep_eta_tmp = genPart->Eta;
	    gen_lep_phi_tmp = genPart->Phi;
	    gen_lep_flv_tmp = TMath::Abs(genPart->PID);
	}
	//--------- ak5 jets branches
	numberJetak5_tmp = branchJet_ak5->GetEntriesFast();
	int iak5=0;
	while(iak5 < numberJetak5_tmp)
	{
	    Jet* jet = (Jet*) branchJet_ak5->At(iak5);
	    ak5_pt_tmp.push_back(jet->PT);
	    ak5_eta_tmp.push_back(jet->Eta);
	    ak5_phi_tmp.push_back(jet->Phi);
	    ak5_m_tmp.push_back(jet->Mass);
	    ak5_HC_EC_tmp.push_back(jet->EhadOverEem);
	    ak5_m_pruned_tmp.push_back(jet->PrunedMass);
	    ak5_btag_tmp.push_back(jet->BTag);
	    iak5++;
	}
	//--------- ak5 GenJets branches
	gen_numberJetak5_tmp = branchGenJet_ak5->GetEntriesFast();
	int gen_iak5=0;
	while(gen_iak5 < gen_numberJetak5_tmp)
	{
	    Jet* jet = (Jet*) branchGenJet_ak5->At(gen_iak5);
	    gen_ak5_pt_tmp.push_back(jet->PT);
	    gen_ak5_eta_tmp.push_back(jet->Eta);
	    gen_ak5_phi_tmp.push_back(jet->Phi);
	    gen_ak5_m_tmp.push_back(jet->Mass);
	    gen_ak5_HC_EC_tmp.push_back(jet->EhadOverEem);
	    gen_ak5_m_pruned_tmp.push_back(jet->PrunedMass);
	    gen_ak5_btag_tmp.push_back(jet->BTag);
	    gen_iak5++;
	}
	//--------- CA8 jets branches
	numberJetCA8_tmp = branchJet_CA8->GetEntriesFast();
	int iCA8=0;
	while(iCA8 < numberJetCA8_tmp)
	{
	    Jet* jet = (Jet*) branchJet_CA8->At(iCA8);
	    CA8_pt_tmp.push_back(jet->PT);
	    CA8_eta_tmp.push_back(jet->Eta);
	    CA8_phi_tmp.push_back(jet->Phi);
	    CA8_m_tmp.push_back(jet->Mass);
	    CA8_HC_EC_tmp.push_back(jet->EhadOverEem);
	    CA8_m_pruned_tmp.push_back(jet->PrunedMass);
	    CA8_tau1_tmp.push_back(jet->tau1);
	    CA8_tau2_tmp.push_back(jet->tau2);
	    CA8_tau3_tmp.push_back(jet->tau3);
	    CA8_btag_tmp.push_back(jet->BTag);
	    iCA8++;
	}
	//--------- CA8 GenJets branches
	gen_numberJetCA8_tmp = branchGenJet_CA8->GetEntriesFast();
	int gen_iCA8=0;
	while(gen_iCA8 < gen_numberJetCA8_tmp)
	{
	    Jet* jet = (Jet*) branchGenJet_CA8->At(gen_iCA8);
	    gen_CA8_pt_tmp.push_back(jet->PT);
	    gen_CA8_eta_tmp.push_back(jet->Eta);
	    gen_CA8_phi_tmp.push_back(jet->Phi);
	    gen_CA8_m_tmp.push_back(jet->Mass);
	    gen_CA8_HC_EC_tmp.push_back(jet->EhadOverEem);
	    gen_CA8_m_pruned_tmp.push_back(jet->PrunedMass);
	    gen_CA8_tau1_tmp.push_back(jet->tau1);
	    gen_CA8_tau2_tmp.push_back(jet->tau2);
	    gen_CA8_tau3_tmp.push_back(jet->tau3);
	    gen_CA8_btag_tmp.push_back(jet->BTag);
	    gen_iCA8++;
	}
	//--------- MET
	MissingET* met = (MissingET*) branchMET->At(0);		
	MET_tmp = met->MET;
	MET_phi_tmp	= met->Phi;
	//--------- gen MET (v)
	TLorentzVector* gen_met_4vect = new TLorentzVector(0,0,0,0);
	for(int iPart = 0; iPart < branchGenPart->GetEntriesFast(); iPart++)
	{
	    GenParticle* genPart = (GenParticle*) branchGenPart->At(iPart);
	    int gen_PID = TMath::Abs(genPart->PID);
	    if(gen_PID == 12 || gen_PID == 14 || gen_PID == 16)
	    {
		gen_met_4vect -> SetPx(gen_met_4vect->Px() + genPart->Px);
		gen_met_4vect -> SetPy(gen_met_4vect->Py() + genPart->Py);
		gen_met_4vect -> SetPz(gen_met_4vect->Pz() + genPart->Pz);
		gen_met_4vect -> SetE(gen_met_4vect->E() + genPart->E);
	    }
	    if(gen_PID == 16)
	    {
		gen_wasTau_tmp = 1;
	    }
	}
	gen_MET_tmp = gen_met_4vect -> Et();
	gen_MET_eta_tmp = gen_met_4vect -> Eta();
	gen_MET_phi_tmp = gen_met_4vect -> Phi();
	//--------- Pile Up
	Rho* rho = (Rho*) branchRho->At(0);
	Rho* rho_fw = (Rho*) branchRho->At(1);
	Rho_tmp = rho->Rho;
	Rho_fw_tmp = rho_fw->Rho;
	gen_nPV_tmp = branchVertex->GetEntriesFast();
	vector<float> Z_vertex;
	int iVertex = 0;
	while(iVertex < gen_nPV_tmp)
	{
	    Vertex* vertex = (Vertex*) branchVertex->At(iVertex);
            Z_vertex.push_back(vertex->Z);
	    iVertex++;
	}
	sort(Z_vertex.begin(), Z_vertex.end());
	//--------- as in Delphes vertices closer then 100 micron are considered unresolved while otherwise perfect reconstruction efficency is assumed 
        for(iVertex = 1; iVertex < Z_vertex.size(); iVertex++)
        {
            if(TMath::Abs(Z_vertex.at(iVertex-1)-Z_vertex.at(iVertex)) > 0.1)
            {
		nPV_tmp++;
            }
        }
	//--------- Let's put this entry in the output tree
	easyTree -> Fill();
	//--------- Cleaning all tmp variables for the next entry
	ak5_pt_tmp.erase(ak5_pt_tmp.begin(), ak5_pt_tmp.end());
	ak5_eta_tmp.erase(ak5_eta_tmp.begin(), ak5_eta_tmp.end());
	ak5_phi_tmp.erase(ak5_phi_tmp.begin(), ak5_phi_tmp.end());
	ak5_m_tmp.erase(ak5_m_tmp.begin(), ak5_m_tmp.end());
	ak5_HC_EC_tmp.erase(ak5_HC_EC_tmp.begin(), ak5_HC_EC_tmp.end());
	ak5_m_pruned_tmp.erase(ak5_m_pruned_tmp.begin(), ak5_m_pruned_tmp.end());
	ak5_btag_tmp.erase(ak5_btag_tmp.begin(), ak5_btag_tmp.end());
	gen_ak5_pt_tmp.erase(gen_ak5_pt_tmp.begin(), gen_ak5_pt_tmp.end());
	gen_ak5_eta_tmp.erase(gen_ak5_eta_tmp.begin(), gen_ak5_eta_tmp.end());
	gen_ak5_phi_tmp.erase(gen_ak5_phi_tmp.begin(), gen_ak5_phi_tmp.end());
	gen_ak5_m_tmp.erase(gen_ak5_m_tmp.begin(), gen_ak5_m_tmp.end());
	gen_ak5_HC_EC_tmp.erase(gen_ak5_HC_EC_tmp.begin(), gen_ak5_HC_EC_tmp.end());
	gen_ak5_m_pruned_tmp.erase(gen_ak5_m_pruned_tmp.begin(), gen_ak5_m_pruned_tmp.end());
	gen_ak5_btag_tmp.erase(gen_ak5_btag_tmp.begin(), gen_ak5_btag_tmp.end());
	CA8_pt_tmp.erase(CA8_pt_tmp.begin(), CA8_pt_tmp.end());
	CA8_eta_tmp.erase(CA8_eta_tmp.begin(), CA8_eta_tmp.end());
	CA8_phi_tmp.erase(CA8_phi_tmp.begin(), CA8_phi_tmp.end());
	CA8_m_tmp.erase(CA8_m_tmp.begin(), CA8_m_tmp.end());
	CA8_HC_EC_tmp.erase(CA8_HC_EC_tmp.begin(), CA8_HC_EC_tmp.end());
	CA8_m_pruned_tmp.erase(CA8_m_pruned_tmp.begin(), CA8_m_pruned_tmp.end());
	CA8_btag_tmp.erase(CA8_btag_tmp.begin(), CA8_btag_tmp.end());
	CA8_tau1_tmp.erase(CA8_tau1_tmp.begin(), CA8_tau1_tmp.end());
	CA8_tau2_tmp.erase(CA8_tau2_tmp.begin(), CA8_tau2_tmp.end());
	CA8_tau3_tmp.erase(CA8_tau3_tmp.begin(), CA8_tau3_tmp.end());
	gen_CA8_pt_tmp.erase(gen_CA8_pt_tmp.begin(), gen_CA8_pt_tmp.end());
	gen_CA8_eta_tmp.erase(gen_CA8_eta_tmp.begin(), gen_CA8_eta_tmp.end());
	gen_CA8_phi_tmp.erase(gen_CA8_phi_tmp.begin(), gen_CA8_phi_tmp.end());
	gen_CA8_m_tmp.erase(gen_CA8_m_tmp.begin(), gen_CA8_m_tmp.end());
	gen_CA8_HC_EC_tmp.erase(gen_CA8_HC_EC_tmp.begin(), gen_CA8_HC_EC_tmp.end());
	gen_CA8_m_pruned_tmp.erase(gen_CA8_m_pruned_tmp.begin(), gen_CA8_m_pruned_tmp.end());
	gen_CA8_btag_tmp.erase(gen_CA8_btag_tmp.begin(), gen_CA8_btag_tmp.end());
	gen_CA8_tau1_tmp.erase(gen_CA8_tau1_tmp.begin(), gen_CA8_tau1_tmp.end());
	gen_CA8_tau2_tmp.erase(gen_CA8_tau2_tmp.begin(), gen_CA8_tau2_tmp.end());
	gen_CA8_tau3_tmp.erase(gen_CA8_tau3_tmp.begin(), gen_CA8_tau3_tmp.end());
	nPV_tmp = 0;		
	gen_wasTau_tmp = 0;	
    }	
    easyTree -> Print("easyDelphes");
    outputFile -> Write();		
    delete outputFile;				
}