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
#include "../include/easy_tree.h"

//#define gSIGMAMG800 = 
#define gSIGMAW4JETS = 211
#define gSIGMAZgJETS = 3053.71

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
    vector<float> std_vect_tmp;
    vector<float>::iterator max_tmp;
    //---jets
    float deltaR_jets_tmp=0;
    int n_CA8_tmp=0, n_ak5_tmp=0;
    int good_CA8_tmp=0;
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
    delphes_tree* DT ;
    TTree* ET;
    if(sampleName.Contains("Signal")) 
	DT->Init(ch, true);
    else 
	DT->Init(ch, false);
    InitEasyTree(ET);

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
        good_CA8_tmp = 0;
        ak5_tmp.clear();
	dR_Wl_j = 100;
	dR_Wh_j = 100;
        //---Reco leptonic object
        TLorentzVector lep_4vect_tmp, nu_4vect_tmp, J_4vect_tmp, Wl_reco_4vect_tmp, G_reco_4vect_tmp;
        TLorentzVector vbf_jj_4vect_tmp, vbf_j1_4vect_tmp, vbf_j2_4vect_tmp;
        TLorentzVector Wl_closerj_4vect_tmp, Wh_closerj_4vect_tmp;
        //------------------Preselection--------------------------------------------------
        //---PV selection---
        for(int iVertex=0; iVertex<DT->nPV->size(); iVertex++)
        {
            //---choose highets pt^2 vertex
            if( DT->sum_pt_square->at(iVertex) > S_ver_tmp )
            {
                X_ver_tmp = DT->vertex_X->at(iVertex);
                Y_ver_tmp = DT->vertex_Y->at(iVertex);
                Z_ver_tmp = DT->vertex_Z->at(iVertex);
            }
        }
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
	lep_4vect_tmp.SetPtEtaPhiE(DT->lep_pt->at(good_lep_tmp),DT->lep_eta->at(good_lep_tmp),
                                   DT->lep_phi->at(good_lep_tmp),DT->lep_E_with_smearing->at(good_lep_tmp));
        nu_4vect_tmp.SetPtEtaPhiM(DT->MET->at(0),0,DT->MET_phi->at(0),0);
        nu_4vect_tmp = BuildNu4Vector (lep_4vect_tmp, nu_4vect_tmp, nu_pz_comp_type);
        Wl_reco_4vect_tmp = lep_4vect_tmp + nu_4vect_tmp;
	//---select reco_CA8---
        int good_CA8=-1;
        for(int iCA8=0; iCA8<DT->number_jet_CA8->at(0); iCA8++)
        {
            if( DT->jet_CA8_pt->at(iCA8) > 200 && abs(DT->jet_CA8_eta->at(iCA8)) < 2.4 )
            {
                good_CA8_tmp = iCA8;
                n_CA8_tmp++;
            }
        }
	J_4vect_tmp.SetPtEtaPhiM(DT->jet_CA8_pt->at(good_CA8),DT->jet_CA8_eta->at(good_CA8),
                                 DT->jet_CA8_phi->at(good_CA8),DT->jet_CA8_mass_pruned->at(good_CA8));
        G_reco_4vect_tmp = Wl_reco_4vect_tmp + J_4vect_tmp;	
        if(n_CA8_tmp > 0)
        {
            for(int iak5=0; iak5<DT->number_jet_ak5->at(0); iak5++)
            {           
                //---remove W's CA8_jet from ak5 collection 
                deltaR_jets_tmp = DeltaR(DT->jet_CA8_eta->at(good_CA8_tmp), DT->jet_ak5_eta->at(iak5),
					 DT->jet_CA8_phi->at(good_CA8_tmp), DT->jet_ak5_phi->at(iak5));
                if(DT->jet_ak5_pt->at(iak5) > 30 && deltaR_jets_tmp > 0.8 && DT->jet_ak5_eta->at(iak5) < 4.5)
                {            
                    ak5_tmp.push_back(iak5);
                    n_ak5_tmp++;
		    for(int iak5=0; iak5<ak5_tmp.size(); iak5++)
		    {
			float dR_Wl_j_tmp = DeltaR(Wl_reco_4vect_tmp.Eta(), 
						   DT->jet_ak5_eta->at(iak5),
						   Wl_reco_4vect_tmp.Phi(), 
						   DT->jet_ak5_phi->at(iak5));
			float dR_Wh_j_tmp = DeltaR(J_4vect_tmp.Eta(), 
						   DT->jet_ak5_eta->at(iak5),
						   J_4vect_tmp.Phi(), 
						   DT->jet_ak5_phi->at(iak5));
			if(dR_Wl_j_tmp < dR_Wh_j_tmp && dR_Wl_j_tmp < dR_Wl_j)
			{
			    Wl_closerjet_tmp = iak5;
			    dR_Wl_j = dR_Wl_j_tmp;
			}
			if(dR_Wh_j_tmp < dR_Wl_j_tmp && dR_Wh_j_tmp < dR_Wh_j)
			{
			    Wh_closerjet_tmp = iak5;
			    dR_Wh_j = dR_Wh_j_tmp;
			}
		    }
                }
	    }
	}
	vbf_j1_4vect_tmp.SetPtEtaPhiM(DT->jet_ak5_pt->at(ak5_tmp.at(0)), 
				      DT->jet_ak5_eta->at(ak5_tmp.at(0)),
				      DT->jet_ak5_phi->at(ak5_tmp.at(0)), 
				      DT->jet_ak5_mass_pruned->at(ak5_tmp.at(0)));
	vbf_j2_4vect_tmp.SetPtEtaPhiM(DT->jet_ak5_pt->at(ak5_tmp.at(1)), 
				      DT->jet_ak5_eta->at(ak5_tmp.at(1)),
				      DT->jet_ak5_phi->at(ak5_tmp.at(1)), 
				      DT->jet_ak5_mass_pruned->at(ak5_tmp.at(1)));
	vbf_jj_4vect_tmp = vbf_j1_4vect_tmp + vbf_j2_4vect_tmp;
	Wl_closerj_4vect_tmp.SetPtEtaPhiM(DT->jet_ak5_pt->at(Wl_closerjet_tmp),
					  DT->jet_ak5_eta->at(Wl_closerjet_tmp),
					  DT->jet_ak5_phi->at(Wl_closerjet_tmp),
					  DT->jet_ak5_mass_pruned->at(Wl_closerjet_tmp));	
	Wh_closerj_4vect_tmp.SetPtEtaPhiM(DT->jet_ak5_pt->at(Wh_closerjet_tmp),
					  DT->jet_ak5_eta->at(Wh_closerjet_tmp),
					  DT->jet_ak5_phi->at(Wh_closerjet_tmp),
					  DT->jet_ak5_mass_pruned->at(Wh_closerjet_tmp));
        //-----apply preselection cuts-----
        if( n_lepl_tmp > 0 || n_lept_tmp != 1 || n_CA8_tmp < 1 || n_ak5_tmp < 2)
        {
            continue;
        }
        countPSEvents++;
        //-----Store reco variables-----		
        //---leptons
	lep_pt = lep_4vect_tmp.Pt();
	lep_pt = lep_4vect_tmp.Eta();
	lep_pt = lep_4vect_tmp.Phi();
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
	lv_closerjet_mass = (Wl_reco_4vect_tmp + Wl_closerj_4vect_tmp).M();
	//---CA8 jet (hadronic W reco)
	CA8_jet_pt = J_4vect_tmp.Pt();
	CA8_jet_eta = J_4vect_tmp.Eta();
	CA8_jet_phi = J_4vect_tmp.Phi();
	CA8_jet_mass = J_4vect_tmp.M();
	CA8_jet_t2t1 = DT->jet_CA8_tau2->at(good_CA8_tmp)/DT->jet_CA8_tau1->at(good_CA8_tmp);
	CA8_jet_t3t2 = DT->jet_CA8_tau3->at(good_CA8_tmp)/DT->jet_CA8_tau2->at(good_CA8_tmp);
	CA8_closerjet_mass = (J_4vect_tmp + Wh_closerj_4vect_tmp).M();
        //---mlvJ
	lvJ_mass = G_reco_4vect_tmp.M();
	
	//-----Store gen variables-----
/*        //---Select W's
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
	    }*/

        //---tag quark (for W+jets they may also be glouns...)
	q1_tmp = -1;
	q2_tmp = -1;
	q1_pt_tmp = -1;
	q2_pt_tmp = -1;
        for(int iPart=0; iPart<DT->lhe_p_pt->size(); iPart++)
        {
            if(DT->lhe_p_from_W->at(iPart) == 0) 
                std_vect_tmp.push_back(DT->lhe_p_pt->at(iPart));
        }
        while(q2_pt_tmp == -1 && std_vect_tmp.size()>0)
        {
            max_tmp = max_element(std_vect_tmp.begin(), std_vect_tmp.end());
            if (q1_pt_tmp == -1)
            {
                q1_pt_tmp = *max_tmp;
                std_vect_tmp.erase(max_tmp);
            }
            else 
		q2_pt_tmp = *max_tmp;
        }
        for(int iPart=0; iPart<DT->lhe_p_pt->size(); iPart++)
        {
            if(DT->lhe_p_pt->at(iPart) == q1_pt_tmp) q1_tmp = iPart;
            else if(DT->lhe_p_pt->at(iPart) == q2_pt_tmp) q2_tmp = iPart;
        }
	if(q1_tmp == -1)
	    q1_tmp = 0;
	if(q2_tmp == -1)
	    q2_tmp = 0;
        TLorentzVector qq_vbf_4vect_tmp, q1_vbf_4vect_tmp, q2_vbf_4vect_tmp;
        q1_vbf_4vect_tmp.SetPtEtaPhiM(DT->lhe_p_pt->at(q1_tmp),DT->lhe_p_eta->at(q1_tmp),
                                      DT->lhe_p_phi->at(q1_tmp),0);
        q2_vbf_4vect_tmp.SetPtEtaPhiM(DT->lhe_p_pt->at(q2_tmp),DT->lhe_p_eta->at(q2_tmp),
                                      DT->lhe_p_phi->at(q2_tmp),0);
        qq_vbf_4vect_tmp = q1_vbf_4vect_tmp + q2_vbf_4vect_tmp;
    }    
    cout << "number of events that passed the preselection:  " << countPSEvents << endl;    
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
