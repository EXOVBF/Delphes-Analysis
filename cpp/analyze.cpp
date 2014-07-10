/***************************************************************************************** 

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

#include "../include/light_tree.h"

const float LUMI = 19250;

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

int makePlots (TString sampleName, vector<TString> datasetBaseName, 
	       map<TString, TH1F *> histos_1D, map<TString, TH2F *> histos_2D) 
{
    
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
        makePlots(sampleName, sample_set, histos_1D_tmp, histos_2D_tmp);
        //---call to the plotter
        SaveHistos1D (outfile, histos_1D_tmp);
        SaveHistos2D (outfile, histos_2D_tmp);
    }
    outfile->Close () ;

    return 0 ;
}


