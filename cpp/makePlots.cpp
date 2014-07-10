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

//****************************************************************************************

int CreateHistos1D (TString suffix, map<TString, TH1F *> & histos, 
                    vector<string>* var, vector<int>* nbin, vector<float>* min, vector<float>* max) 
{
    //---set sample colors
    int sample_color=0;
    if( suffix.Contains("tt") == 1 ) sample_color=3;
    if( suffix.Contains("W") == 1 ) sample_color=2;
    if( suffix.Contains("Z") == 1) sample_color=5;
    //---create histogram
    for(int iVar=0; iVar<var->size(); iVar++)
    {
        if( suffix.Contains("Signal")==0 && ((TString)(var->at(iVar))).Contains("MWW")==1) continue;
        TH1F * h_tmp = new TH1F (TString (var->at(iVar)) + "__"+suffix, var->at(iVar) + "__"+suffix, 
				 nbin->at(iVar), min->at(iVar), max->at(iVar)); 
        h_tmp->SetFillColor(sample_color);
        histos[(var->at(iVar)).c_str()] = h_tmp;
    }
    return histos.size();
}

//****************************************************************************************

int CreateHistos2D (TString suffix, map<TString, TH2F *> & histos, 
                    vector<string>* var, vector<int>* nbin_x, vector<float>* min_x, vector<float>* max_x,
		    vector<int>* nbin_y, vector<float>* min_y, vector<float>* max_y)
{
    for(int iVar=0; iVar<var->size(); iVar++)
    {
        if( strcmp(suffix.Data(),"W4jets")==0 && ((TString)(var->at(iVar))).Contains("WJ")==1) continue;
        if( strcmp(suffix.Data(),"Z")==0 && ((TString)(var->at(iVar))).Contains("WJ")==1) continue;
        TH2F * h_tmp = new TH2F (TString (var->at(iVar)) + "__"+suffix, var->at(iVar) + "__"+suffix, 
				 nbin_x->at(iVar), min_x->at(iVar), max_x->at(iVar), 
				 nbin_y->at(iVar), min_y->at(iVar), max_y->at(iVar)); 
        histos[(var->at(iVar)).c_str()] = h_tmp;
    }
    return histos.size();
}

//****************************************************************************************

void SaveHistos1D (TFile * outfile, map<TString, TH1F *> histos, float cross_section = 1.)
{
    outfile->cd () ;
    for (map<TString, TH1F *>::iterator iMap = histos.begin () ;
         iMap != histos.end () ; ++iMap)
    {
        //iMap->second->Scale (cross_section / iMap->second->Integral ()) ;
        iMap->second->Write () ;
    }
    return ;
}

//****************************************************************************************

void SaveHistos2D (TFile * outfile, map<TString, TH2F *> histos, float cross_section = 1.)
{
    outfile->cd () ;
    for (map<TString, TH2F *>::iterator iMap = histos.begin () ;
         iMap != histos.end () ; ++iMap)
    {
        iMap->second->Write () ;
    }
    return ;
}

//*****************************************************************************************

int readDataset (TString sampleName, vector<TString> datasetBaseName, map<TString, TH1F *> histos_1D, map<TString, TH2F *> histos_2D) 
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
    vector<int> CA8_tmp, ak5_tmp;
    //---Ws
    float deltaR_Wl=0, deltaR_Wlv=0, deltaR_Wj=0, deltaR_Wjgen=0, W_lep_tmp=0, W_had_tmp=0;
    //---counters
    int countPSEvents=0, countBSEvents=0;
    
    TChain * ch = new TChain ("Delphes") ;
    for(int i=0; i<datasetBaseName.size(); i++)
    {
        ch->Add(datasetBaseName.at(i)+"*root") ;
    }
    cout << "read " << ch->GetEntries() << " events in " << ch->GetNtrees() 
         << " files of " << sampleName.Data() << " sample\n" ;
    delphes_tree DT ;
    if(sampleName.Contains("Signal")) DT.Init(ch, true);
    else DT.Init(ch, false);

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
        CA8_tmp.clear();
        ak5_tmp.clear();
        //------------------Preselection--------------------------------------------------
        //---PV selection---
        for(int iVertex=0; iVertex<DT.nPV->size(); iVertex++)
        {
            //---choose highets pt^2 vertex
            if( DT.sum_pt_square->at(iVertex) > S_ver_tmp )
            {
                X_ver_tmp = DT.vertex_X->at(iVertex);
                Y_ver_tmp = DT.vertex_Y->at(iVertex);
                Z_ver_tmp = DT.vertex_Z->at(iVertex);
            }
        }
        //---Selec leptons from PV---        
        for(int iLep=0; iLep<DT.lep_number->at(0); iLep++)
        {        
            //---no cuts plot
            if(abs(DT.lep_flv->at(iLep)) == 11 && DT.lep_eta->at(iLep) < 1.5)
            {
                if(histos_1D["nc_el_pt_BR"]) histos_1D["nc_el_pt_BR"]->Fill(DT.lep_pt->at(iLep));
                if(histos_1D["nc_el_phi_BR"]) histos_1D["nc_el_phi_BR"]->Fill(DT.lep_phi->at(iLep));
            }
            else if(abs(DT.lep_flv->at(iLep)) == 11 && DT.lep_eta->at(iLep) > 1.5)
            {
                if(histos_1D["nc_el_pt_EC"]) histos_1D["nc_el_pt_EC"]->Fill(DT.lep_pt->at(iLep));
                if(histos_1D["nc_el_phi_EC"]) histos_1D["nc_el_phi_EC"]->Fill(DT.lep_phi->at(iLep));
            }
            else if (abs(DT.lep_flv->at(iLep)) == 13  && DT.lep_eta->at(iLep) < 1.5)
            {
                if(histos_1D["nc_mu_pt_BR"]) histos_1D["nc_mu_pt_BR"]->Fill(DT.lep_pt->at(iLep));
                if(histos_1D["nc_mu_phi_BR"]) histos_1D["nc_mu_phi_BR"]->Fill(DT.lep_phi->at(iLep));                
            }
            else if (abs(DT.lep_flv->at(iLep)) == 13 && DT.lep_eta->at(iLep) > 1.5)
            {
                if(histos_1D["nc_mu_pt_EC"]) histos_1D["nc_mu_pt_EC"]->Fill(DT.lep_pt->at(iLep));
                if(histos_1D["nc_mu_phi_EC"]) histos_1D["nc_mu_phi_EC"]->Fill(DT.lep_phi->at(iLep));                
            }
            if(abs(DT.lep_flv->at(iLep)) == 11)
            {
                if(histos_1D["nc_el_pt"]) histos_1D["nc_el_pt"]->Fill(DT.lep_pt->at(iLep));
                if(histos_1D["nc_el_eta"]) histos_1D["nc_el_eta"]->Fill(DT.lep_eta->at(iLep));
            }
            if(abs(DT.lep_flv->at(iLep)) == 13)
            {
                if(histos_1D["nc_mu_pt"]) histos_1D["nc_mu_pt"]->Fill(DT.lep_pt->at(iLep));
                if(histos_1D["nc_mu_eta"]) histos_1D["nc_mu_eta"]->Fill(DT.lep_eta->at(iLep));
            }
            d_xy_tmp = TMath::Sqrt(pow(X_ver_tmp - DT.lep_X_vertex->at(iLep),2) + pow(X_ver_tmp - DT.lep_X_vertex->at(iLep),2));
            //---tight electron
            if( abs(DT.lep_flv->at(iLep)) == 11 && abs(DT.lep_Z_vertex->at(iLep)-Z_ver_tmp) < 0.1 && DT.lep_pt->at(iLep) > 35 && 
                (DT.lep_eta->at(iLep) < 1.4442 || DT.lep_eta->at(iLep) > 1.566) && DT.lep_isolation->at(iLep)*DT.lep_pt->at(iLep) < 5 &&
                d_xy_tmp < 0.2 && DT.MET->at(0) > 65) 
            {
                good_lep_tmp = iLep;
                n_lept_tmp++;
            }
            //---loose electron
            else if( abs(DT.lep_flv->at(iLep)) == 11 && abs(DT.lep_Z_vertex->at(iLep)-Z_ver_tmp) < 0.1 && DT.lep_pt->at(iLep) > 20 && 
		     DT.lep_eta->at(iLep) < 2.5 && DT.lep_isolation->at(iLep)*DT.lep_pt->at(iLep) < 5 && d_xy_tmp < 0.2) 
            {
                n_lepl_tmp++;
            }
            //---tight muon
            if( abs(DT.lep_flv->at(iLep)) == 13 && abs(DT.lep_Z_vertex->at(iLep)-Z_ver_tmp) < 5 && DT.lep_pt->at(iLep) > 50 && 
		DT.lep_eta->at(iLep) < 2.1 && DT.lep_isolation->at(iLep) < 0.1 && d_xy_tmp < 2 && DT.MET->at(0) > 50)
            {
                good_lep_tmp = iLep;
                n_lept_tmp++;
            }
            //---loose muon
            else if( abs(DT.lep_flv->at(iLep)) == 13 && abs(DT.lep_Z_vertex->at(iLep)-Z_ver_tmp) < 5 && DT.lep_pt->at(iLep) > 10 && 
                     DT.lep_eta->at(iLep) < 2.5 && DT.lep_isolation->at(iLep) < 0.1 && d_xy_tmp < 2)
            {
                n_lepl_tmp++;
            }
        }
        //---select reco_CA8---
        int good_CA8=-1;
        for(int iCA8=0; iCA8<DT.number_jet_CA8->at(0); iCA8++)
        {
//            if( DT.jet_CA8_pt->at(iCA8) > 200 && abs(DT.jet_CA8_eta->at(iCA8)) < 2.4 )
	    if( DT.jet_CA8_mass_pruned->at(iCA8) > 65 && DT.jet_CA8_mass_pruned->at(iCA8) < 105 && 
		DT.jet_CA8_pt->at(iCA8) > 200 && abs(DT.jet_CA8_eta->at(iCA8)) < 2.4 )
            {
                CA8_tmp.push_back(iCA8);
                n_CA8_tmp++;
            }
        }/*
        if(n_CA8_tmp>1)
        {
            for(int iCA8=0; iCA8<CA8_tmp.size(); iCA8++)
            {
                int Ws_ak5[2]={-1,-1};
                for(int iak5=0; iak5<DT.number_jet_ak5->at(0); iak5++)
                {           
                    deltaR_jets_tmp = TMath::Sqrt(pow(DT.jet_CA8_eta->at(CA8_tmp.at(iCA8))-DT.jet_ak5_eta->at(iak5),2) +
                                                  pow(DeltaPhi(DT.jet_CA8_phi->at(CA8_tmp.at(iCA8)),DT.jet_ak5_phi->at(iak5)),2));
                    if( DT.jet_ak5_pt->at(iak5) > 30 && deltaR_jets_tmp < 0.8)
                    {                 
                        if(Ws_ak5[0] == -1) Ws_ak5[0]=iak5;
                        else if (Ws_ak5[1] == -1) Ws_ak5[1]=iak5;
                        else break;
                    }
                }
                if(Ws_ak5[0] != -1 && Ws_ak5[1] == -1 && DT.jet_ak5_mass_pruned->at(Ws_ak5[0]) > 65 &&
                   DT.jet_ak5_mass_pruned->at(Ws_ak5[0]) < 105 && good_CA8 == -1)
                {
                    good_CA8 = CA8_tmp.at(iCA8);
                }
                else if(Ws_ak5[0] != -1 && Ws_ak5[1] != -1)
                {
                    TLorentzVector jj_vbf_4vect_tmp, j1_vbf_4vect_tmp, j2_vbf_4vect_tmp;
                    j1_vbf_4vect_tmp.SetPtEtaPhiM(DT.jet_ak5_pt->at(Ws_ak5[0]), DT.jet_ak5_eta->at(Ws_ak5[0]),
                                                  DT.jet_ak5_phi->at(Ws_ak5[0]), DT.jet_ak5_mass_pruned->at(Ws_ak5[0]));
                    j2_vbf_4vect_tmp.SetPtEtaPhiM(DT.jet_ak5_pt->at(Ws_ak5[1]), DT.jet_ak5_eta->at(Ws_ak5[1]),
                                                  DT.jet_ak5_phi->at(Ws_ak5[1]), DT.jet_ak5_mass_pruned->at(Ws_ak5[1]));
                    jj_vbf_4vect_tmp = j1_vbf_4vect_tmp + j2_vbf_4vect_tmp;
                    if(jj_vbf_4vect_tmp.M() > 65 && jj_vbf_4vect_tmp.M() < 105 && good_CA8 == -1) 
                        good_CA8 = CA8_tmp.at(iCA8);
                }
            }
	    }*/
        if(good_CA8 == -1 && n_CA8_tmp > 0) good_CA8 = CA8_tmp.at(0);
        if(n_CA8_tmp > 0)
        {
            for(int iak5=0; iak5<DT.number_jet_ak5->at(0); iak5++)
            {           
                //---remove W's CA8_jet from ak5 collection 
                deltaR_jets_tmp = TMath::Sqrt(pow(DT.jet_CA8_eta->at(good_CA8)-DT.jet_ak5_eta->at(iak5),2) +
                                              pow(DeltaPhi(DT.jet_CA8_phi->at(good_CA8),DT.jet_ak5_phi->at(iak5)),2));
                if(DT.jet_ak5_pt->at(iak5) > 30 && deltaR_jets_tmp > 0.8 && 
		   (DT.jet_ak5_mass_pruned->at(iak5) < 65 || DT.jet_ak5_mass_pruned->at(iak5) > 105)
		    && DT.jet_ak5_eta->at(iak5) < 4.5)
                {            
                    ak5_tmp.push_back(iak5);
                    n_ak5_tmp++;
                }
                if(histos_1D["nc_jet_ak5_eta"]) histos_1D["nc_jet_ak5_eta"]->Fill(DT.jet_ak5_eta->at(iak5));        
            }
	}
        //---select good ak5---
	/*
        for(int iak5=0; iak5<ak5_tmp.size(); iak5++)
        {
	    bool fromW=false;
	    TLorentzVector jj_vbf_4vect_tmp, j1_vbf_4vect_tmp, j2_vbf_4vect_tmp;
            j1_vbf_4vect_tmp.SetPtEtaPhiM(DT.jet_ak5_pt->at(ak5_tmp.at(iak5)), DT.jet_ak5_eta->at(ak5_tmp.at(iak5)),
                                          DT.jet_ak5_phi->at(ak5_tmp.at(iak5)), DT.jet_ak5_mass_pruned->at(ak5_tmp.at(iak5)));
	    for(int i=ak5_tmp.at(iak5)+1; i<DT.number_jet_ak5->at(0); i++)
            {
	        j2_vbf_4vect_tmp.SetPtEtaPhiM(DT.jet_ak5_pt->at(i), DT.jet_ak5_eta->at(i),
					      DT.jet_ak5_phi->at(i), DT.jet_ak5_mass_pruned->at(i));
                jj_vbf_4vect_tmp = j1_vbf_4vect_tmp + j2_vbf_4vect_tmp;
		if(jj_vbf_4vect_tmp.M() > 65 && jj_vbf_4vect_tmp.M() < 105)
		    fromW = true;
	    }
	    if((DT.jet_ak5_mass_pruned->at(iak5) > 65 && DT.jet_ak5_mass_pruned->at(iak5) < 105) || fromW)
	    {
		n_ak5_tmp--;
		ak5_tmp.erase(ak5_tmp.begin()+iak5);
	    }
	    }*/
        //---gen_ak5 nc plots---
        for(int iak5=0; iak5<DT.number_gen_jet_ak5->at(0); iak5++)
        {           
            if( DT.gen_jet_ak5_pt->at(iak5) > 30 && abs(DT.gen_jet_ak5_eta->at(iak5)) <= 1.556 )
            {                
                if(histos_1D["nc_gen_jet_ak5_pt_BR"]) histos_1D["nc_gen_jet_ak5_pt_BR"]->Fill(DT.gen_jet_ak5_pt->at(iak5));
                if(histos_1D["nc_gen_jet_ak5_phi_BR"]) histos_1D["nc_gen_jet_ak5_phi_BR"]->Fill(DT.gen_jet_ak5_phi->at(iak5));
            }
            if( DT.gen_jet_ak5_pt->at(iak5) > 30 && abs(DT.gen_jet_ak5_eta->at(iak5)) > 1.556 )
            {                
                if(histos_1D["nc_gen_jet_ak5_pt_EC"]) histos_1D["nc_gen_jet_ak5_pt_EC"]->Fill(DT.gen_jet_ak5_pt->at(iak5));
                if(histos_1D["nc_gen_jet_ak5_phi_EC"]) histos_1D["nc_gen_jet_ak5_phi_EC"]->Fill(DT.gen_jet_ak5_phi->at(iak5));
            }
            if(histos_1D["nc_gen_jet_ak5_eta"]) histos_1D["nc_gen_jet_ak5_eta"]->Fill(DT.gen_jet_ak5_eta->at(iak5));        
        }
        //---tag quarks nc plots---
        for (int iQuark=0; iQuark<DT.lhe_p_pt->size(); iQuark++)
        {
            if(DT.lhe_p_from_W->at(iQuark) == 0 && histos_1D["nc_q_vbf_phi"])
                histos_1D["nc_q_vbf_phi"]->Fill(DT.lhe_p_phi->at(iQuark));
        }
        //-----apply preselection cuts-----
        if( n_lepl_tmp > 0 || n_lept_tmp != 1 || n_CA8_tmp < 1 || n_ak5_tmp < 2)
        {
            continue;
        }
        countPSEvents++;
        //---preselection plots---
        //---Select W's (remember that l- has +code while l+ has -code...)
        if(sampleName.Contains("Z") == 0)
        {
            if(TMath::Sign(DT.lhe_W_pid->front(),DT.lhe_lep_flv->at(0)) == DT.lhe_W_pid->front()) 
            {
                W_lep_tmp=1;
                W_had_tmp=0;
            }
            else 
            {
                W_lep_tmp=0;
                W_had_tmp=1;
            }
        }
        //---leptonic W matching (don't do this for Z+jets)
        if(DT.lhe_W_pt->size() > 0)
        {
            deltaR_Wl = TMath::Sqrt(pow(DT.lhe_W_eta->at(W_lep_tmp)-DT.lep_eta->at(good_lep_tmp),2)+
                                    pow(DeltaPhi(DT.lhe_W_phi->at(W_lep_tmp),DT.lep_phi->at(good_lep_tmp)),2));
        }
        //---hadronic W matching (don't do this for W+jets)
        if(DT.lhe_W_pt->size() == 2)
        {
            deltaR_Wj = TMath::Sqrt(pow(DT.lhe_W_eta->at(W_had_tmp)-DT.jet_CA8_eta->at(good_CA8),2)+
				    pow(DeltaPhi(DT.lhe_W_phi->at(W_had_tmp),DT.jet_CA8_phi->at(good_CA8)),2));
            deltaR_Wjgen = TMath::Sqrt(pow(DT.lhe_W_eta->at(W_had_tmp)-DT.gen_jet_CA8_eta->at(good_CA8),2)+
				       pow(DeltaPhi(DT.lhe_W_phi->at(W_had_tmp),DT.gen_jet_CA8_phi->at(good_CA8)),2));
            if(histos_2D["ps_Wl_deltaR_vs_ptRatio"])
                histos_2D["ps_Wl_deltaR_vs_ptRatio"]->Fill(DT.lep_pt->at(good_lep_tmp)/
							   DT.lhe_W_pt->at(W_lep_tmp),deltaR_Wl);
            if(histos_2D["ps_WJ_deltaR_vs_ptRatio"])
                histos_2D["ps_WJ_deltaR_vs_ptRatio"]->Fill(DT.jet_CA8_pt->at(good_CA8)/DT.lhe_W_pt->at(W_had_tmp),
							   deltaR_Wj);
            if(histos_2D["ps_gen_WJ_deltaR_vs_ptRatio"])
		histos_2D["ps_gen_WJ_deltaR_vs_ptRatio"]->Fill(DT.gen_jet_CA8_pt->at(good_CA8)/DT.lhe_W_pt->at(W_had_tmp),
							       deltaR_Wjgen);
        }
	//---quark from W
	for(int iPart=0; iPart<DT.lhe_p_pt->size(); iPart++)
        {
            if(DT.lhe_p_from_W->at(iPart) == 1) 
                std_vect_tmp.push_back(DT.lhe_p_pt->at(iPart));
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
	std_vect_tmp.clear();
        for(int iPart=0; iPart<DT.lhe_p_pt->size(); iPart++)
        {
            if(DT.lhe_p_pt->at(iPart) == q1_pt_tmp) q1_tmp = iPart;
            else if(DT.lhe_p_pt->at(iPart) == q2_pt_tmp) q2_tmp = iPart;
        }
	if(q1_tmp == -1)
	    q1_tmp = 0;
	if(q2_tmp == -1)
	    q2_tmp = 0;
	if(histos_1D["ps_W_qq_deltaR"])
	    histos_1D["ps_W_qq_deltaR"]->Fill(TMath::Sqrt(pow(TMath::Abs(DT.lhe_p_eta->at(q2_tmp) - 
									    DT.lhe_p_eta->at(q1_tmp)),2)+
							     pow(DeltaPhi(DT.lhe_p_phi->at(q2_tmp), 
							     DT.lhe_p_phi->at(q1_tmp)),2)));
	//---ak5 from W
	TLorentzVector W_jj;
	W_jj.SetPtEtaPhiM(DT.jet_ak5_pt->at(0)+DT.jet_ak5_pt->at(1),
			  DT.jet_ak5_eta->at(0)+DT.jet_ak5_eta->at(1),
			  DT.jet_ak5_phi->at(0)+DT.jet_ak5_phi->at(1),
			  DT.jet_ak5_mass_pruned->at(0)+DT.jet_ak5_mass_pruned->at(1));
	if(histos_1D["ps_W_jj_M_pruned"])
            histos_1D["ps_W_jj_M_pruned"]->Fill(W_jj.M());
	if(histos_1D["ps_W_jj_deltaR"] && W_jj.M() > 65 && W_jj.M() < 105)
	    histos_1D["ps_W_jj_deltaR"]->Fill(TMath::Sqrt(pow(TMath::Abs(DT.jet_ak5_eta->at(0) - 
									    DT.jet_ak5_eta->at(1)),2)+
							     pow(DeltaPhi(DT.jet_ak5_phi->at(0), 
							     DT.jet_ak5_phi->at(1)),2)));
	if(histos_1D["ps_W_q2j2_delta_phi"] && W_jj.M() > 65 && W_jj.M() < 105)
            histos_1D["ps_W_q2j2_delta_phi"]->Fill(DeltaPhi(DT.lhe_p_phi->at(q2_tmp), 
							    DT.jet_ak5_phi->at(1)));
	if(histos_1D["ps_W_q2j2_delta_eta"] && W_jj.M() > 65 && W_jj.M() < 105)
            histos_1D["ps_W_q2j2_delta_eta"]->Fill(TMath::Abs(DT.lhe_p_eta->at(q2_tmp) - 
							      DT.jet_ak5_eta->at(1)));
	if(histos_1D["ps_W_q1j1_delta_phi"] && W_jj.M() > 65 && W_jj.M() < 105)
            histos_1D["ps_W_q1j1_delta_phi"]->Fill(DeltaPhi(DT.lhe_p_phi->at(q1_tmp), 
							    DT.jet_ak5_phi->at(0)));
	if(histos_1D["ps_W_q1j1_delta_eta"] && W_jj.M() > 65 && W_jj.M() < 105)
            histos_1D["ps_W_q1j1_delta_eta"]->Fill(TMath::Abs(DT.lhe_p_eta->at(q1_tmp) - 
							      DT.jet_ak5_eta->at(0)));
        //---fat jet
        if(histos_1D["ps_jet_CA8_M_pruned"])
            histos_1D["ps_jet_CA8_M_pruned"]->Fill(DT.jet_CA8_mass_pruned->at(good_CA8));
        if(histos_1D["ps_jet_CA8_tau2tau1"])
            histos_1D["ps_jet_CA8_tau2tau1"]->Fill(DT.jet_CA8_tau2->at(good_CA8) /
                                                   DT.jet_CA8_tau1->at(good_CA8));
        //---leptons
        if(histos_1D["ps_l_pt"])
            histos_1D["ps_l_pt"]->Fill(DT.lep_pt->at(good_lep_tmp));
        if(histos_1D["ps_l_eta"])
            histos_1D["ps_l_eta"]->Fill(DT.lep_eta->at(good_lep_tmp));
        if(histos_1D["ps_l_phi"])
            histos_1D["ps_l_phi"]->Fill(DT.lep_phi->at(good_lep_tmp));
        //---MET
        if(histos_1D["ps_MET"])
            histos_1D["ps_MET"]->Fill(DT.MET->at(0));
        if(histos_1D["ps_MET_phi"])
            histos_1D["ps_MET_phi"]->Fill(DT.MET_phi->at(0));
        //---leptonic W reco
        TLorentzVector lep_4vect_tmp, nu_4vect_tmp, J_4vect_tmp, W_reco_4vect_tmp, G_reco_4vect_tmp;
        lep_4vect_tmp.SetPtEtaPhiE(DT.lep_pt->at(good_lep_tmp),DT.lep_eta->at(good_lep_tmp),
                                   DT.lep_phi->at(good_lep_tmp),DT.lep_E_with_smearing->at(good_lep_tmp));
        nu_4vect_tmp.SetPtEtaPhiM(DT.MET->at(0),0,DT.MET_phi->at(0),0);
        nu_4vect_tmp = BuildNu4Vector (lep_4vect_tmp, nu_4vect_tmp, nu_pz_comp_type);
        W_reco_4vect_tmp = lep_4vect_tmp + nu_4vect_tmp;
        if(histos_1D["ps_mlv_reco"])
            histos_1D["ps_mlv_reco"]->Fill(W_reco_4vect_tmp.M());
        if(sampleName.Contains("Z") == 0)
        {
            deltaR_Wlv = TMath::Sqrt(pow(DT.lhe_W_eta->at(W_lep_tmp)-W_reco_4vect_tmp.Eta(),2)+
                                     pow(DT.lhe_W_phi->at(W_lep_tmp)-W_reco_4vect_tmp.Phi(),2));    
            if(histos_2D["ps_Wlv_daltaR_vs_ptRatio"])
                histos_2D["ps_Wlv_daltaR_vs_ptRatio"]->Fill(W_reco_4vect_tmp.Pt()/
                                                            DT.lhe_W_pt->at(W_lep_tmp),deltaR_Wlv);
            if(histos_1D["ps_MW_gen"])
                histos_1D["ps_MW_gen"]->Fill(DT.lhe_W_mass->at(W_lep_tmp));
        
        }
        //---tag jets
        TLorentzVector jj_vbf_4vect_tmp, j1_vbf_4vect_tmp, j2_vbf_4vect_tmp;
        j1_vbf_4vect_tmp.SetPtEtaPhiM(DT.jet_ak5_pt->at(ak5_tmp.at(0)), DT.jet_ak5_eta->at(ak5_tmp.at(0)),
                                      DT.jet_ak5_phi->at(ak5_tmp.at(0)), DT.jet_ak5_mass_pruned->at(ak5_tmp.at(0)));
        j2_vbf_4vect_tmp.SetPtEtaPhiM(DT.jet_ak5_pt->at(ak5_tmp.at(1)), DT.jet_ak5_eta->at(ak5_tmp.at(1)),
                                      DT.jet_ak5_phi->at(ak5_tmp.at(1)), DT.jet_ak5_mass_pruned->at(ak5_tmp.at(1)));
        jj_vbf_4vect_tmp = j1_vbf_4vect_tmp + j2_vbf_4vect_tmp;
        if(histos_1D["ps_jet_ak5_vbf1_pt"])
            histos_1D["ps_jet_ak5_vbf1_pt"]->Fill(DT.jet_ak5_pt->at(ak5_tmp.at(0)));
        if(histos_1D["ps_jet_ak5_vbf1_eta"])
            histos_1D["ps_jet_ak5_vbf1_eta"]->Fill(DT.jet_ak5_eta->at(ak5_tmp.at(0)));
        if(histos_1D["ps_jet_ak5_vbf1_phi"])
            histos_1D["ps_jet_ak5_vbf1_phi"]->Fill(DT.jet_ak5_phi->at(ak5_tmp.at(0)));
        if(histos_1D["ps_jet_ak5_vbf2_pt"])
            histos_1D["ps_jet_ak5_vbf2_pt"]->Fill(DT.jet_ak5_pt->at(ak5_tmp.at(1)));
        if(histos_1D["ps_jet_ak5_vbf2_eta"])
            histos_1D["ps_jet_ak5_vbf2_eta"]->Fill(DT.jet_ak5_eta->at(ak5_tmp.at(1)));
        if(histos_1D["ps_jet_ak5_vbf2_phi"])
            histos_1D["ps_jet_ak5_vbf2_phi"]->Fill(DT.jet_ak5_phi->at(ak5_tmp.at(1)));
        if(histos_1D["ps_mjj_vbf"])
            histos_1D["ps_mjj_vbf"]->Fill(jj_vbf_4vect_tmp.M());
        if(histos_1D["ps_delta_eta_jj_vbf"])
            histos_1D["ps_delta_eta_jj_vbf"]->Fill(TMath::Abs(DT.jet_ak5_eta->at(ak5_tmp.at(1)) - 
                                                              DT.jet_ak5_eta->at(ak5_tmp.at(0))));
        if(histos_1D["ps_delta_phi_jj_vbf"])
            histos_1D["ps_delta_phi_jj_vbf"]->Fill(DeltaPhi(DT.jet_ak5_phi->at(ak5_tmp.at(1)),
							    DT.jet_ak5_phi->at(ak5_tmp.at(0))));
        //---tag quark (for W+jets they may also be glouns...)
	q1_tmp = -1;
	q2_tmp = -1;
	q1_pt_tmp = -1;
	q2_pt_tmp = -1;
        for(int iPart=0; iPart<DT.lhe_p_pt->size(); iPart++)
        {
            if(DT.lhe_p_from_W->at(iPart) == 0) 
                std_vect_tmp.push_back(DT.lhe_p_pt->at(iPart));
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
        for(int iPart=0; iPart<DT.lhe_p_pt->size(); iPart++)
        {
            if(DT.lhe_p_pt->at(iPart) == q1_pt_tmp) q1_tmp = iPart;
            else if(DT.lhe_p_pt->at(iPart) == q2_pt_tmp) q2_tmp = iPart;
        }
	if(q1_tmp == -1)
	    q1_tmp = 0;
	if(q2_tmp == -1)
	    q2_tmp = 0;
        TLorentzVector qq_vbf_4vect_tmp, q1_vbf_4vect_tmp, q2_vbf_4vect_tmp;
        q1_vbf_4vect_tmp.SetPtEtaPhiM(DT.lhe_p_pt->at(q1_tmp),DT.lhe_p_eta->at(q1_tmp),
                                      DT.lhe_p_phi->at(q1_tmp),0);
        q2_vbf_4vect_tmp.SetPtEtaPhiM(DT.lhe_p_pt->at(q2_tmp),DT.lhe_p_eta->at(q2_tmp),
                                      DT.lhe_p_phi->at(q2_tmp),0);
        qq_vbf_4vect_tmp = q1_vbf_4vect_tmp + q2_vbf_4vect_tmp;
        if(histos_1D["ps_q_vbf1_pt"])
            histos_1D["ps_q_vbf1_pt"]->Fill(DT.lhe_p_pt->at(q1_tmp));
        if(histos_1D["ps_q_vbf1_eta"])
            histos_1D["ps_q_vbf1_eta"]->Fill(DT.lhe_p_eta->at(q1_tmp));
        if(histos_1D["ps_q_vbf1_phi"])
            histos_1D["ps_q_vbf1_phi"]->Fill(DT.lhe_p_phi->at(q1_tmp));
        if(histos_1D["ps_q_vbf2_pt"])
            histos_1D["ps_q_vbf2_pt"]->Fill(DT.lhe_p_pt->at(q2_tmp));
        if(histos_1D["ps_q_vbf2_eta"])
            histos_1D["ps_q_vbf2_eta"]->Fill(DT.lhe_p_eta->at(q2_tmp));
        if(histos_1D["ps_q_vbf2_phi"])
            histos_1D["ps_q_vbf2_phi"]->Fill(DT.lhe_p_phi->at(q2_tmp));
        if(histos_1D["ps_mqq_vbf"])
            histos_1D["ps_mqq_vbf"]->Fill(qq_vbf_4vect_tmp.M());
        if(histos_1D["ps_delta_eta_qq_vbf"])
            histos_1D["ps_delta_eta_qq_vbf"]->Fill(TMath::Abs(DT.lhe_p_eta->at(q2_tmp) - 
                                                              DT.lhe_p_eta->at(q1_tmp)));
        if(histos_1D["ps_delta_phi_qq_vbf"])
            histos_1D["ps_delta_phi_qq_vbf"]->Fill(DeltaPhi(DT.lhe_p_phi->at(q2_tmp), 
							    DT.lhe_p_phi->at(q1_tmp)));
        if(histos_1D["ps_delta_R_qq_vbf"])
            histos_1D["ps_delta_R_qq_vbf"]->Fill(TMath::Sqrt(pow(TMath::Abs(DT.lhe_p_eta->at(q2_tmp) - 
									    DT.lhe_p_eta->at(q1_tmp)),2)+
							     pow(DeltaPhi(DT.lhe_p_phi->at(q2_tmp), 
									  DT.lhe_p_phi->at(q1_tmp)),2)));
        //---ratio tag jet/tag quark
        if(histos_2D["ps_ratio_tag_j_q_vbf"])
            histos_2D["ps_ratio_tag_j_q_vbf"]->Fill(DT.lhe_p_pt->at(q1_tmp)/DT.jet_ak5_pt->at(ak5_tmp.at(0)),
                                                    DT.lhe_p_pt->at(q2_tmp)/DT.jet_ak5_pt->at(ak5_tmp.at(1)));
        //---mlvJ
        J_4vect_tmp.SetPtEtaPhiM(DT.jet_CA8_pt->at(good_CA8),DT.jet_CA8_eta->at(good_CA8),
                                 DT.jet_CA8_phi->at(good_CA8),DT.jet_CA8_mass_pruned->at(good_CA8));
        G_reco_4vect_tmp = W_reco_4vect_tmp + J_4vect_tmp;
        if(histos_1D["ps_mlvJ_reco"])
            histos_1D["ps_mlvJ_reco"]->Fill(G_reco_4vect_tmp.M());
        if(histos_1D["ps_MWW_gen"])
            histos_1D["ps_MWW_gen"]->Fill(DT.lhe_X_mass->at(0));
        //------------------Basic Selection-----------------------------------------------
        
        
        
    }    
    cout << "number of events that passed the preselection:  " << countPSEvents << endl;    
    return 0 ;
}


//****************************************************************************************
//main

int main (int argc, char* argv[]) 
{
//-----------------Config files check-----------------------------------------------------
    TString configFile1;
    TString configFile2;
    TString outFile = "output.root";
    if(argc < 3)
    {
        cout << "ERROR: config files !!! " << endl;
        return 0;
    }
    configFile1 = argv[1];
    configFile2 = argv[2];
    cout << "CONFIG variables file:  " << configFile1.Data() << endl;
    cout << "CONFIG samples file:  " << configFile2.Data() << endl;
    if(argc > 3) outFile = argv[3];
    cout << "OUTPUT file: " << outFile.Data() << endl;
    
//-----------------Definitions------------------------------------------------------------
    //---variables
    string buffer_type;
    int buffer_nbin=0, buffer_nbin_x=0, buffer_nbin_y=0;
    float buffer_min=0, buffer_max=0, buffer_min_x=0, buffer_max_x=0, buffer_min_y=0, buffer_max_y=0;
    vector<int> nbin, nbin_x, nbin_y; 
    vector<float> min, max, min_x, max_x, min_y, max_y;
    string buffer_var_1d, buffer_var_2d;
    vector<string> var_1d, var_2d;
    //---samples
    string sample_path;
    string buffer_sample_dir, buffer_sample_name;
    vector<string> sample_dir, sample_name;
//-----------------Read Config files------------------------------------------------------
    //---variables
    ifstream config_file_1 (configFile1.Data(), ios::in);
    while(config_file_1 >> buffer_type)
    {   
        if(strcmp(buffer_type.c_str(),"#") == 0)
        {
            getline(config_file_1, buffer_type);
        }
        else if(strcmp(buffer_type.c_str(),"1D") == 0)
        {
            config_file_1 >> buffer_var_1d >> buffer_nbin >> buffer_min >> buffer_max;
            var_1d.push_back(buffer_var_1d);
            nbin.push_back(buffer_nbin);
            min.push_back(buffer_min);
            max.push_back(buffer_max);
        }
        else if(strcmp(buffer_type.c_str(),"2D") == 0)
        {
            config_file_1 >> buffer_var_2d >> buffer_nbin_x >> buffer_min_x >> buffer_max_x 
                          >> buffer_nbin_y >> buffer_min_y >> buffer_max_y;
            var_2d.push_back(buffer_var_2d);
            nbin_x.push_back(buffer_nbin_x);
            min_x.push_back(buffer_min_x);
            max_x.push_back(buffer_max_x);
            nbin_y.push_back(buffer_nbin_y);
            min_y.push_back(buffer_min_y);
            max_y.push_back(buffer_max_y);
        }
    }
    config_file_1.close();
    //---samples
    ifstream config_file_2 (configFile2.Data(), ios::in);
    while(config_file_2 >> buffer_type)
    {   
        if(strcmp(buffer_type.c_str(),"#") == 0)
        {
            getline(config_file_2, buffer_type);
        }
        else if(strcmp(buffer_type.c_str(),"P") == 0)
        {
            config_file_2 >> sample_path;
        }
        else if(strcmp(buffer_type.c_str(),"D") == 0)
        {
            config_file_2 >> buffer_sample_name >> buffer_sample_dir;
            sample_dir.push_back(buffer_sample_dir);
            sample_name.push_back(buffer_sample_name);
        }
    }
    config_file_2.close();

//-----------------Make Plots-------------------------------------------------------------
    TFile * outfile = new TFile (outFile.Data(), "recreate") ;
    for(int iSample=0; iSample<sample_dir.size(); iSample++)
    {
        TString sampleName(sample_name.at(iSample));
        TString samplePath(sample_path);
        map<TString, TH1F*> histos_1D_tmp;
        map<TString, TH2F*> histos_2D_tmp;
        vector<TString> sample_set;
        CreateHistos1D(sampleName, histos_1D_tmp, &var_1d, &nbin, &min, &max);
        CreateHistos2D(sampleName, histos_2D_tmp, &var_2d, &nbin_x, &min_x, &max_x, &nbin_y, &min_y, &max_y);
        samplePath += sample_dir.at(iSample);
        sample_set.push_back(samplePath);
        //---call to the analyzer function
        cout << "Analyzing:  " << sample_set.at(0).Data() << endl;
        readDataset(sampleName, sample_set, histos_1D_tmp, histos_2D_tmp);
        //---call to the plotter
        SaveHistos1D (outfile, histos_1D_tmp);
        SaveHistos2D (outfile, histos_2D_tmp);
    }
    outfile->Close () ;

    return 0 ;
}
