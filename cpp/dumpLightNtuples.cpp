/***************************************************************************************** 
    
    this program apply the preselection cut on MC sample from Delphes and store
    the good events in a tree with only the variables used in the analysis

    c++ -O2 -lm `root-config --cflags --glibs` -o dumpLightNtuples dumpLightNtuples.cpp 
    
    or exexute compile.sh cpp/dumpLightNtuples.cpp

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

//****************************************************************************************

TLorentzVector BuildNu4Vector (TLorentzVector lep_4vect, TLorentzVector MET_4vect, int type = 0, float W_mass = 80.385)
{
    double M_W = W_mass;
    double M_lep = lep_4vect.M();
    double e_lep = lep_4vect.E();
    double px_lep = lep_4vect.Px(); 
    double py_lep = lep_4vect.Py();
    double pz_lep = lep_4vect.Pz();
    double px_nu = MET_4vect.Px();
    double py_nu = MET_4vect.Py();
    double pz_nu = 0.;
    double otherSol = 0.;
    double tmpsol1 = 0.;
    double tmpsol2 = 0.;
    double newPtneutrino1 = 0.;
    double newPtneutrino2 = 0.;
    bool isComplex = false;

    double a = M_W*M_W - M_lep*M_lep + 2.0*px_lep*px_nu + 2.0*py_lep*py_nu;
    double A = 4.0*(e_lep*e_lep - pz_lep*pz_lep);
    double B = -4.0*a*pz_lep;
    double C = 4.0*e_lep*e_lep*(px_nu*px_nu + py_nu*py_nu) - a*a;

    double tmproot = B*B - 4.0*A*C;

    if (tmproot<0) 
    {
        isComplex = true;
        pz_nu = - B/(2*A); // take real part of complex roots
        
        otherSol = pz_nu;
        double p_nu = TMath::Sqrt(px_nu*px_nu + py_nu*py_nu + pz_nu*pz_nu);
        double Delta = (M_W*M_W - M_lep*M_lep);
        double alpha = (px_lep*px_nu/p_nu + py_lep*py_nu/p_nu);
        double pt_nu = TMath::Sqrt( px_nu*px_nu + py_nu*py_nu); // old
        double AA = 4.*pz_lep*pz_lep - 4*e_lep*e_lep + 4*alpha*alpha;
        double BB = 4.*alpha*Delta;
        double CC = Delta*Delta;

        double tmpdisc = BB*BB - 4.0*AA*CC;
        double tmpsolpt1 = (-BB + TMath::Sqrt(tmpdisc))/(2.0*AA);
        double tmpsolpt2 = (-BB - TMath::Sqrt(tmpdisc))/(2.0*AA);

        if ( fabs( tmpsolpt1 - pt_nu ) < fabs( tmpsolpt2 - pt_nu) ) 
        { 
            newPtneutrino1 = tmpsolpt1; 
            newPtneutrino2 = tmpsolpt2;
        }
        else 
        { 
            newPtneutrino1 = tmpsolpt2; 
            newPtneutrino2 = tmpsolpt1; 
        }

    }
    else 
    {
        isComplex = false;
        double tmpsol1 = (-B + TMath::Sqrt(tmproot))/(2.0*A);
        double tmpsol2 = (-B - TMath::Sqrt(tmproot))/(2.0*A);

        if (type == 0 ) 
        {
            // two real roots, pick the lower one
            if (tmpsol2 <= tmpsol1) 
            {
                pz_nu = tmpsol2; 
                otherSol = tmpsol1;
            }
            else 
            {
                pz_nu = tmpsol1; 
                otherSol = tmpsol2; 
            }
        }
        if (type == 1 ) 
        {
            // if pz_nu is > 300 pick the most central root
            if ( pz_nu > 300. ) 
            {
                if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2)) 
                {
                    pz_nu = tmpsol1; 
                    otherSol = tmpsol2; 
                }
                else 
                {
                    pz_nu = tmpsol2; 
                    otherSol = tmpsol1; 
                }
            }
            else
            {
                // two real roots, pick the one closest to pz of lepton
                if (TMath::Abs(tmpsol2-pz_lep) < TMath::Abs(tmpsol1-pz_lep)) 
                {
                    pz_nu = tmpsol2; 
                    otherSol = tmpsol1; 
                }
                else
                {
                    pz_nu = tmpsol1; 
                    otherSol = tmpsol2; 
                }
            }
        }
        if (type == 2 ) 
        {
            // pick the most central root
            if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2)) 
            {
                pz_nu = tmpsol1; 
		otherSol = tmpsol2; 
            }
            else 
            {
                pz_nu = tmpsol2; 
                otherSol = tmpsol1;
            }
        }
        if (type == 3 ) 
        {
            // pick the largest value of the cosine
            TVector3 p3w, p3mu;
            p3w.SetXYZ(px_lep+px_nu, py_lep+py_nu, pz_lep+ tmpsol1);
            p3mu.SetXYZ(px_lep, py_lep, pz_lep );

            double sinthcm1 = 2.*(p3mu.Perp(p3w))/M_W;
            p3w.SetXYZ(px_lep+px_nu, py_lep+py_nu, pz_lep+ tmpsol2);
            double sinthcm2 = 2.*(p3mu.Perp(p3w))/M_W;
            double costhcm1 = TMath::Sqrt(1. - sinthcm1*sinthcm1);
            double costhcm2 = TMath::Sqrt(1. - sinthcm2*sinthcm2);

            if ( costhcm1 > costhcm2 ) 
            {
                pz_nu = tmpsol1; 
                otherSol = tmpsol2; 
            }
            else 
            {
                pz_nu = tmpsol2;
                otherSol = tmpsol1; 
            }
        }
    }
    if(isComplex)
    {
        px_nu = newPtneutrino1*TMath::Cos(MET_4vect.Phi());
        py_nu = newPtneutrino1*TMath::Sin(MET_4vect.Phi());
    }
    MET_4vect.SetPxPyPzE(px_nu,py_nu,pz_nu,TMath::Sqrt(px_nu*px_nu + py_nu*py_nu + pz_nu*pz_nu)); 
    return MET_4vect;
}

//*****************************************************************************************

int readDataset (TString sampleName, vector<TString> datasetBaseName)
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
    vector<float> ord_quark_tmp;
    vector<float>::iterator max_tmp;
    bool isFirst;
    //---jets
    float deltaR_jets_tmp=0;
    int n_CA8_tmp=0, n_ak5_tmp=0;
    int good_CA8_tmp=0;
    vector<int> CA8_tmp;
    vector<int> ak5_tmp;
    //---Ws
    float W_lep_tmp=0, W_had_tmp=0;
    float Wl_closerjet_tmp=0, Wh_closerjet_tmp=0; 
    float dR_Wl_j=100, dR_Wh_j=100; 
    //---counters
    int countPSEvents=0, countBSEvents=0;

    //---Load data---
    TChain * ch = new TChain ("Delphes") ;
    for(int i=0; i<datasetBaseName.size(); i++)
    {
        ch->Add(datasetBaseName.at(i)+"*root") ;
    }
    cout << "read " << ch->GetEntries() << " events in " << ch->GetNtrees() 
         << " files of " << sampleName.Data() << " sample\n" ;
    //---Init trees 
    delphes_tree* DT = new delphes_tree();
    if(sampleName.Contains("Signal")) 
	DT->Init(ch, true);
    else 
	DT->Init(ch, false);
    TFile* outFile = TFile::Open("/afs/cern.ch/user/s/spigazzi/work/EXOVBF/Delphes-Analysis/light_ntuples/"+sampleName+".root", "recreate");
    outFile->cd();
    TTree* ET = new TTree("light_tree", "light_tree");
    InitLightTree(ET);
    ET->SetDirectory(0);

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
	isFirst = 0;
        ord_quark_tmp.clear();
        good_CA8_tmp = 0;
	CA8_tmp.clear();
        ak5_tmp.clear();
	dR_Wl_j = 100;
	dR_Wh_j = 100; 
        //---Reco object
        TLorentzVector lep_4vect_tmp, nu_4vect_tmp, J_4vect_tmp, Wl_reco_4vect_tmp, G_reco_4vect_tmp;
        TLorentzVector vbf_j1_4vect_tmp, vbf_j2_4vect_tmp;
	TLorentzVector vbf_q1_4vect_tmp, vbf_q2_4vect_tmp;
        TLorentzVector Wl_closerj_4vect_tmp, Wh_closerj_4vect_tmp;
        //------------------Preselection--------------------------------------------------
	//---only semileptonic events---
	if(DT->lhe_lep_pt->size() < 1)
	    continue;
        //---PV selection---
        for(int iVertex=0; iVertex<DT->nPV->at(0); iVertex++)
        {
            //---choose highets pt^2 vertex
            if( DT->sum_pt_square->at(iVertex) > S_ver_tmp )
            {
                X_ver_tmp = DT->vertex_X->at(iVertex);
                Y_ver_tmp = DT->vertex_Y->at(iVertex);
                Z_ver_tmp = DT->vertex_Z->at(iVertex);
		S_ver_tmp = DT->sum_pt_square->at(iVertex);
            }
        }
	nPV = DT->nPV->at(0);
        //---Selec leptons from PV---
        for(int iLep=0; iLep<DT->lep_number->at(0); iLep++)
        {    
            d_xy_tmp = TMath::Sqrt(pow(X_ver_tmp - DT->lep_X_vertex->at(iLep),2) + pow(X_ver_tmp - DT->lep_X_vertex->at(iLep),2));
            //---tight electron
            if( abs(DT->lep_flv->at(iLep)) == 11 && abs(DT->lep_Z_vertex->at(iLep)-Z_ver_tmp) < 0.1 && DT->lep_pt->at(iLep) > 35 && 
                (DT->lep_eta->at(iLep) < 1.4442 || DT->lep_eta->at(iLep) > 1.566) && DT->lep_isolation->at(iLep)*DT->lep_pt->at(iLep) < 5 &&
                d_xy_tmp < 0.2 && DT->MET->at(0) > 65) 
            {
                good_lep_tmp = iLep;
                n_lept_tmp++;
            }
            //---loose electron
            else if( abs(DT->lep_flv->at(iLep)) == 11 && abs(DT->lep_Z_vertex->at(iLep)-Z_ver_tmp) < 0.1 && DT->lep_pt->at(iLep) > 20 && 
		     DT->lep_eta->at(iLep) < 2.5 && DT->lep_isolation->at(iLep)*DT->lep_pt->at(iLep) < 5 && d_xy_tmp < 0.2) 
            {
                n_lepl_tmp++;
            }
            //---tight muon
            if( abs(DT->lep_flv->at(iLep)) == 13 && abs(DT->lep_Z_vertex->at(iLep)-Z_ver_tmp) < 5 && DT->lep_pt->at(iLep) > 50 && 
		DT->lep_eta->at(iLep) < 2.1 && DT->lep_isolation->at(iLep) < 0.1 && d_xy_tmp < 2 && DT->MET->at(0) > 50)
            {
                good_lep_tmp = iLep;
                n_lept_tmp++;
            }
            //---loose muon
            else if( abs(DT->lep_flv->at(iLep)) == 13 && abs(DT->lep_Z_vertex->at(iLep)-Z_ver_tmp) < 5 && DT->lep_pt->at(iLep) > 10 && 
                     DT->lep_eta->at(iLep) < 2.5 && DT->lep_isolation->at(iLep) < 0.1 && d_xy_tmp < 2)
            {
                n_lepl_tmp++;
            }
        }
	//---apply basic selection on leptons---
	if(n_lepl_tmp > 0 || n_lept_tmp != 1)
	    continue;
	//---select reco_CA8---
        for(int iCA8=0; iCA8<DT->number_jet_CA8->at(0); iCA8++)
        {
	    //---remove hard lepton from CA8 collection
	    float dR_l_j_tmp = DeltaR(DT->lep_eta->at(good_lep_tmp), DT->jet_CA8_eta->at(iCA8),
				      DT->lep_phi->at(good_lep_tmp), DT->jet_CA8_phi->at(iCA8));
            if(DT->jet_CA8_pt->at(iCA8) > 150 && abs(DT->jet_CA8_eta->at(iCA8)) < 2.4 &&
	       dR_l_j_tmp > 0.8)
	    {
		CA8_tmp.push_back(iCA8);
		n_CA8_tmp++;
	    }
        }
	//---apply basic selection on the CA8 jet
	if(n_CA8_tmp < 1)
	    continue;
	good_CA8_tmp = CA8_tmp.at(0);
	//---build reco object
	lep_4vect_tmp.SetPtEtaPhiE(DT->lep_pt->at(good_lep_tmp),
				   DT->lep_eta->at(good_lep_tmp),
				   DT->lep_phi->at(good_lep_tmp),
				   DT->lep_E_with_smearing->at(good_lep_tmp));
	nu_4vect_tmp.SetPtEtaPhiM(DT->MET->at(0),0,DT->MET_phi->at(0),0);
	nu_4vect_tmp = BuildNu4Vector (lep_4vect_tmp, nu_4vect_tmp, nu_pz_comp_type);
	Wl_reco_4vect_tmp = lep_4vect_tmp + nu_4vect_tmp;
	J_4vect_tmp.SetPtEtaPhiM(DT->jet_CA8_pt->at(good_CA8_tmp),
				 DT->jet_CA8_eta->at(good_CA8_tmp),
				 DT->jet_CA8_phi->at(good_CA8_tmp),
				 DT->jet_CA8_mass_pruned->at(good_CA8_tmp));
	G_reco_4vect_tmp = Wl_reco_4vect_tmp + J_4vect_tmp;
	//---select reco_ak5---
	for(int iak5=0; iak5<DT->number_jet_ak5->at(0); iak5++)
	{           
	    //---remove W CA8_jet from ak5 collection 
	    float dR_Wh_j_tmp = DeltaR(J_4vect_tmp.Eta(), DT->jet_ak5_eta->at(iak5), 
				       J_4vect_tmp.Phi(), DT->jet_ak5_phi->at(iak5));
	    //---remove hard lepton from ak5 collection
	    float dR_l_j_tmp = DeltaR(lep_4vect_tmp.Eta(), DT->jet_ak5_eta->at(iak5),
				      lep_4vect_tmp.Phi(), DT->jet_ak5_phi->at(iak5));
	    if(DT->jet_ak5_pt->at(iak5) > 30 && DT->jet_ak5_eta->at(iak5) < 4.7 &&
	       dR_Wh_j_tmp > 0.8 && dR_l_j_tmp > 0.5)
	    {            
		ak5_tmp.push_back(iak5);
		n_ak5_tmp++;
		float dR_Wl_j_tmp = DeltaR(Wl_reco_4vect_tmp.Eta(), 
					   DT->jet_ak5_eta->at(iak5),
					   Wl_reco_4vect_tmp.Phi(), 
					   DT->jet_ak5_phi->at(iak5));
		if(dR_Wl_j_tmp < dR_Wl_j)
		{
		    Wl_closerjet_tmp = iak5;
		    dR_Wl_j = dR_Wl_j_tmp;
		}
		if(dR_Wh_j_tmp < dR_Wh_j)
		{
		    Wh_closerjet_tmp = iak5;
		    dR_Wh_j = dR_Wh_j_tmp;
		}
	    }
	}
	//---apply basic selection ak5 jets
	if(n_ak5_tmp < 2)
            continue;
        countPSEvents++;
        //-----Store reco variables-----		
        //---leptons
	lep_pt = lep_4vect_tmp.Pt();
	lep_eta = lep_4vect_tmp.Eta();
	lep_phi = lep_4vect_tmp.Phi();
        //---MET
	MET = DT->MET->at(0);
	MET_phi = DT->MET_phi->at(0);
	//---leptonic W reco
	lv_mass = Wl_reco_4vect_tmp.M();
	lv_pt = Wl_reco_4vect_tmp.Pt();
	lv_eta = Wl_reco_4vect_tmp.Eta();
	lv_phi = Wl_reco_4vect_tmp.Phi();
	lv_delta_R = DeltaR(lep_4vect_tmp.Eta(), nu_4vect_tmp.Eta(),
			    lep_4vect_tmp.Phi(), nu_4vect_tmp.Phi());
	lv_Mt = TMath::Sqrt(lep_pt*MET*(1-TMath::Cos(DeltaPhi(lep_phi, MET_phi))));
	//---Wl closerjet
	Wl_closerj_4vect_tmp.SetPtEtaPhiM(DT->jet_ak5_pt->at(Wl_closerjet_tmp),
					  DT->jet_ak5_eta->at(Wl_closerjet_tmp),
					  DT->jet_ak5_phi->at(Wl_closerjet_tmp),
					  DT->jet_ak5_mass_pruned->at(Wl_closerjet_tmp));
	lv_closerjet_mass = (Wl_reco_4vect_tmp + Wl_closerj_4vect_tmp).M();
	//---CA8 jet (hadronic W reco)
	CA8_jet_pt = J_4vect_tmp.Pt();
	CA8_jet_eta = J_4vect_tmp.Eta();
	CA8_jet_phi = J_4vect_tmp.Phi();
	CA8_jet_mass = J_4vect_tmp.M();
	CA8_jet_t2t1 = DT->jet_CA8_tau2->at(good_CA8_tmp)/DT->jet_CA8_tau1->at(good_CA8_tmp);
	CA8_jet_t3t2 = DT->jet_CA8_tau3->at(good_CA8_tmp)/DT->jet_CA8_tau2->at(good_CA8_tmp);
	//---CA8 closerjet
	Wh_closerj_4vect_tmp.SetPtEtaPhiM(DT->jet_ak5_pt->at(Wh_closerjet_tmp),
					  DT->jet_ak5_eta->at(Wh_closerjet_tmp),
					  DT->jet_ak5_phi->at(Wh_closerjet_tmp),
					  DT->jet_ak5_mass_pruned->at(Wh_closerjet_tmp));
	CA8_closerjet_mass = (J_4vect_tmp + Wh_closerj_4vect_tmp).M();
	//---Object separation
	lv_J_delta_phi = DeltaPhi(lv_phi, CA8_jet_phi);
	MET_J_delta_phi = DeltaPhi(MET_phi, CA8_jet_phi);
	l_J_delta_R = DeltaR(lep_eta, CA8_jet_eta, lep_phi, CA8_jet_phi);
        //---mlvJ
	lvJ_mass = G_reco_4vect_tmp.M();
	//---vbf jets
	vbf_j1_4vect_tmp.SetPtEtaPhiM(DT->jet_ak5_pt->at(ak5_tmp.at(0)), 
				      DT->jet_ak5_eta->at(ak5_tmp.at(0)),
				      DT->jet_ak5_phi->at(ak5_tmp.at(0)), 
				      DT->jet_ak5_mass_pruned->at(ak5_tmp.at(0)));
	vbf_j2_4vect_tmp.SetPtEtaPhiM(DT->jet_ak5_pt->at(ak5_tmp.at(1)), 
				      DT->jet_ak5_eta->at(ak5_tmp.at(1)),
				      DT->jet_ak5_phi->at(ak5_tmp.at(1)), 
				      DT->jet_ak5_mass_pruned->at(ak5_tmp.at(1)));
	vbf_jet1_pt = vbf_j1_4vect_tmp.Pt();
	vbf_jet1_eta = vbf_j1_4vect_tmp.Eta();
	vbf_jet1_phi = vbf_j1_4vect_tmp.Phi();
	vbf_jet1_mass = vbf_j1_4vect_tmp.M();
	vbf_jet1_btag = DT->jet_ak5_btag->at(ak5_tmp.at(0));
	vbf_jet2_pt = vbf_j2_4vect_tmp.Pt();
	vbf_jet2_eta = vbf_j2_4vect_tmp.Eta();
	vbf_jet2_phi = vbf_j2_4vect_tmp.Phi();
	vbf_jet2_mass = vbf_j2_4vect_tmp.M();
	vbf_jet2_btag = DT->jet_ak5_btag->at(ak5_tmp.at(1));
	vbf_jj_mass = (vbf_j1_4vect_tmp + vbf_j2_4vect_tmp).M();
	vbf_jj_delta_phi = DeltaPhi(vbf_jet1_phi, vbf_jet2_phi); 
	vbf_jj_delta_R = DeltaR(vbf_jet1_eta, vbf_jet2_eta, 
				vbf_jet1_phi, vbf_jet2_phi);
	//-----Store gen variables-----
        //---Select W's
        if(sampleName.Contains("Z") == 0)
        {
            if(TMath::Sign(DT->lhe_W_pid->front(),DT->lhe_lep_flv->at(0)) == DT->lhe_W_pid->front()) 
            {
                W_lep_tmp=1;
                W_had_tmp=0;
            }
            else 
            {
                W_lep_tmp=0;
                W_had_tmp=1;
            }
	    gen_Wl_pt = DT->lhe_W_pt->at(W_lep_tmp);
	    gen_Wl_eta = DT->lhe_W_eta->at(W_lep_tmp);
	    gen_Wl_phi = DT->lhe_W_phi->at(W_lep_tmp);
	    gen_Wl_mass = DT->lhe_W_mass->at(W_lep_tmp);
	    gen_Wh_pt = DT->lhe_W_pt->at(W_lep_tmp);
	    gen_Wh_eta = DT->lhe_W_eta->at(W_lep_tmp);
	    gen_Wh_phi = DT->lhe_W_phi->at(W_lep_tmp);
	    gen_Wh_mass = DT->lhe_W_mass->at(W_lep_tmp); 
	    if(sampleName.Contains("Signal") == 1)
		gen_X_mass = DT->lhe_X_mass->at(0);
	}
	//---leptons
	gen_lep_pt = DT->lhe_lep_pt->at(0);
	gen_lep_eta = DT->lhe_lep_eta->at(0);
	gen_lep_phi = DT->lhe_lep_phi->at(0);
	if(sampleName.Contains("Z") == 0)
	{
	    gen_nu_pt = DT->lhe_nu_pt->at(0);
	    gen_nu_eta = DT->lhe_nu_eta->at(0);
	    gen_nu_phi = DT->lhe_nu_phi->at(0);
	}
	//---order quark wrt pt
	ord_quark_tmp.push_back(0);
	for(int iPart=1; iPart<DT->lhe_p_pt->size(); iPart++)
	{
	    bool insert=0;
	    for(int iq=0; iq<ord_quark_tmp.size(); iq++)
	    {
	if(DT->lhe_p_pt->at(iPart) > DT->lhe_p_pt->at(ord_quark_tmp.at(iq)))
		{
		    ord_quark_tmp.insert(ord_quark_tmp.begin()+iq, iPart);
		    insert = 1;
		    break;
		}
	    }
	    if(insert == 0)
		ord_quark_tmp.push_back(iPart);
	}
	//---quark from W
	if(sampleName.Contains("Z") == 0)
	{
	    isFirst = 1;
	    for(int iPart=0; iPart<ord_quark_tmp.size(); iPart++)
	    {
		if(DT->lhe_p_from_W->at(iPart) == 1 && isFirst == 1)
		{
		    gen_W_q1_pt = DT->lhe_p_pt->at(iPart);
		    gen_W_q1_eta = DT->lhe_p_eta->at(iPart);
		    gen_W_q1_phi = DT->lhe_p_phi->at(iPart);
		    isFirst = 0;
		}
		else if(DT->lhe_p_from_W->at(iPart) == 1 && isFirst == 0)
		{
		    gen_W_q2_pt = DT->lhe_p_pt->at(iPart);
		    gen_W_q2_eta = DT->lhe_p_eta->at(iPart);
		    gen_W_q2_phi = DT->lhe_p_phi->at(iPart);
		    isFirst = -1;
		}
	    }
	}
	//---tag quark
	isFirst = 1;
        for(int iPart=0; iPart<ord_quark_tmp.size(); iPart++)
	{
	    if(DT->lhe_p_from_W->at(iPart) == 0 && isFirst == 1)
	    {
		gen_vbf_q1_pt = DT->lhe_p_pt->at(iPart);
		gen_vbf_q1_eta = DT->lhe_p_eta->at(iPart);
		gen_vbf_q1_phi = DT->lhe_p_phi->at(iPart);
		gen_vbf_q1_btag = 0;
		if(TMath::Abs(DT->lhe_p_flv->at(iPart)) == 5)
		    gen_vbf_q1_btag = 1;
		isFirst = 0;
	    }
	    else if(DT->lhe_p_from_W->at(iPart) == 0 && isFirst == 0)
	    {
		gen_vbf_q2_pt = DT->lhe_p_pt->at(iPart);
		gen_vbf_q2_eta = DT->lhe_p_eta->at(iPart);
		gen_vbf_q2_phi = DT->lhe_p_phi->at(iPart);
		gen_vbf_q2_btag = 0;
		if(TMath::Abs(DT->lhe_p_flv->at(iPart)) == 5)
		    gen_vbf_q2_btag = 1;
		isFirst = -1;
	    }
	}
	vbf_q1_4vect_tmp.SetPtEtaPhiM(gen_vbf_q1_pt, gen_vbf_q1_eta,
				      gen_vbf_q1_phi, 0);
	vbf_q2_4vect_tmp.SetPtEtaPhiM(gen_vbf_q2_pt, gen_vbf_q2_eta,
				      gen_vbf_q2_phi, 0);
	gen_vbf_qq_mass = (vbf_q1_4vect_tmp + vbf_q2_4vect_tmp).M();
	gen_vbf_qq_delta_phi = DeltaPhi(gen_vbf_q1_phi, gen_vbf_q2_phi);
	gen_vbf_qq_delta_R = DeltaR(gen_vbf_q1_eta, gen_vbf_q2_eta,
				    gen_vbf_q1_phi, gen_vbf_q2_phi);
	ET->Fill();
    }    
    cout << "number of events that passed the preselection:  " << countPSEvents << endl;    

    ET->Write();
    outFile->Close();
	
    return 0 ;
}


//****************************************************************************************
//main

int main (int argc, char* argv[]) 
{
//-----------------Config files check-----------------------------------------------------
    TString configFile;
    TString outFile = "output.root";
    if(argc < 2)
    {
        cout << "ERROR: config files !!! " << endl;
        return 0;
    }
    configFile = argv[1];
    cout << "CONFIG samples file:  " << configFile.Data() << endl;
    if(argc > 2) outFile = argv[2];
    cout << "OUTPUT file: " << outFile.Data() << endl;
    
//-----------------Definitions------------------------------------------------------------
    //---samples
    string sample_path;
    string buffer_sample_dir, buffer_sample_name, buffer_type;
    vector<string> sample_dir, sample_name;
//-----------------Read Config file-------------------------------------------------------
    //---samples
    ifstream config_file (configFile.Data(), ios::in);
    while(config_file >> buffer_type)
    {   
        if(strcmp(buffer_type.c_str(),"#") == 0)
        {
            getline(config_file, buffer_type);
        }
        else if(strcmp(buffer_type.c_str(),"P") == 0)
        {
            config_file >> sample_path;
        }
        else if(strcmp(buffer_type.c_str(),"D") == 0)
        {
            config_file >> buffer_sample_name >> buffer_sample_dir;
            sample_dir.push_back(buffer_sample_dir);
            sample_name.push_back(buffer_sample_name);
        }
    }
    config_file.close();
//-----------------Make Plots-------------------------------------------------------------
    for(int iSample=0; iSample<sample_dir.size(); iSample++)
    {
        TString sampleName(sample_name.at(iSample));
        TString samplePath(sample_path);
        vector<TString> sample_set;
        samplePath += sample_dir.at(iSample);
        sample_set.push_back(samplePath);
        //---call to the analyzer function
        cout << "Analyzing:  " << sample_set.at(0).Data() << endl;
        readDataset(sampleName, sample_set);
    }
    return 0 ;
}
