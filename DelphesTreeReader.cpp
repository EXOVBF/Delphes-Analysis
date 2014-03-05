/****************************************************************************************

	before compiling ---> source ../Decay/setup_slc6.sh

    compile with ---> c++ -O2 -W -Wall `root-config --cflags --glibs` -L /afs/cern.ch/user/l/lbrianza/work/SIMONE/delphes_code/ -I /afs/cern.ch/user/l/lbrianza/work/SIMONE/delphes_code/ -lDelphes -o DelphesTreeReader DelphesTreeReader.cpp

    instead of lbrianza/work/SIMONE/delphes_code/ put your Delphes folder

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
#include "THStack.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

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

eventType event_preselector(ExRootTreeReader *delphesTree, TClonesArray* branchJet, TClonesArray* branchEl ,TClonesArray* branchMu, TClonesArray* branchMET)
{
	eventType goodEvent;

	for(int iEvent = 0; iEvent < delphesTree->GetEntries(); iEvent++)
	{
		delphesTree -> ReadEntry(iEvent);
		if(branchJet->GetEntries() > 1 && branchEl->GetEntries() == 1 && branchMu->GetEntries() == 0)
		{
			Electron* el = (Electron*) branchEl->At(0);
		 	Jet* jet_1 = (Jet*) branchJet->At(0);
			Jet* jet_2 = (Jet*) branchJet->At(1);
			MissingET* met = (MissingET*) branchMET->At(0); 
			if(el->PT > 20 && jet_1->PT > 20 && jet_2->PT > 20 && met->MET > 30)
			{
				goodEvent.eventID.push_back(iEvent);
				goodEvent.lepFlavor.push_back(0);
			}
		}
		if(branchJet->GetEntries() > 1 && branchMu->GetEntries() == 1 && branchEl->GetEntries() == 0)
		{
			Muon* mu = (Muon*) branchMu->At(0);
		 	Jet* jet_1 = (Jet*) branchJet->At(0);
			Jet* jet_2 = (Jet*) branchJet->At(1);
			MissingET* met = (MissingET*) branchMET->At(0); 
			if(mu->PT > 20 && jet_1->PT > 20 && jet_2->PT > 20 && met->MET > 30)
			{
				goodEvent.eventID.push_back(iEvent);
				goodEvent.lepFlavor.push_back(1);
			}
		}
	}
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
//definitions	
		
	vector<string> inputFiles;
	TChain* delphesNtuples = new TChain("Delphes");
	ExRootTreeReader *delphesTree = new ExRootTreeReader(delphesNtuples);	
	TFile* outputFile = TFile::Open(argv[2],"recreate");
	TTree* easyTree = new TTree("Delphes_easy","Delphes_easy");
	float lep_pt_tmp=0,lep_eta_tmp=0,lep_phi_tmp=0;
	float MET_tmp=0,m_jj_tmp=0;

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

//----------------------------------------------------------------------------------------
//filling the output tree 
	
	//--------- getting objects from the delphes tree
	TClonesArray* branchJet = delphesTree->UseBranch("Jet");
	TClonesArray* branchEl = delphesTree->UseBranch("Electron");
	TClonesArray* branchMu = delphesTree->UseBranch("Muon");
	TClonesArray* branchMET = delphesTree->UseBranch("MissingET");
	//--------- creating branches for the new (light) tree
	//--------- lep
	TBranch* lep_pt = easyTree->Branch("lep_pt",&lep_pt_tmp,"lep_pt/F");
	TBranch* lep_eta = easyTree->Branch("lep_eta",&lep_eta_tmp,"lep_eta/F");
	TBranch* lep_phi = easyTree->Branch("lep_phi",&lep_phi_tmp,"lep_phi/F");
	//------- jets
	TBranch* MET = easyTree->Branch("MET",&MET_tmp,"MET/F");
	TBranch* m_jj = easyTree->Branch("max_pt_m_jj",&m_jj_tmp,"max_pt_m_jj/F");
		
	eventType goodEvent = event_preselector(delphesTree,branchJet,branchEl,branchMu,branchMET);		
		
	for(int iEvent = 0; iEvent < (long int)goodEvent.eventID.size(); iEvent++)
	{
		delphesTree -> ReadEntry(goodEvent.eventID.at(iEvent));
		if(goodEvent.lepFlavor.at(iEvent) == 0)
		{
			Electron* el = (Electron*) branchEl->At(0);	
			lep_pt_tmp = el->PT;
			lep_eta_tmp = el->Eta;
			lep_phi_tmp = el->Phi;
		}
		else
		{
			Muon* mu = (Muon*) branchMu->At(0);
			lep_pt_tmp = mu->PT;
			lep_eta_tmp = mu->Eta;
			lep_phi_tmp = mu->Phi;
		}
		MissingET* met = (MissingET*) branchMET->At(0);
		MET_tmp = met->MET;	
	 	Jet* jet_1 = (Jet*) branchJet->At(0);
		Jet* jet_2 = (Jet*) branchJet->At(1);
		TLorentzVector* jj4vector = new TLorentzVector(0,0,0,0);
		jj4vector -> SetPtEtaPhiM(jet_1->PT+jet_2->PT,jet_1->Eta+jet_2->Eta,jet_1->Phi+jet_2->Phi,jet_1->Mass+jet_2->Mass);
		m_jj_tmp = jj4vector->M(); 
		easyTree -> Fill();	
	}	
	
	easyTree -> Print("Delphes_easy");
	outputFile -> Write();		
	delete outputFile;				
}
