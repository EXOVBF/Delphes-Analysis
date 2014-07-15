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

#include "../include/delphes_tree.h"
#include "../include/light_tree.h"

//---define global variables
const float gSIGMAttbar = 245.8;
const float ttbarBR = 0.55;
const float gSIGMAW4JETS = 263;
const float gSIGMAZ4JETS = 3503;
const float gSIGMASTOPs = 3.79;
const float gSIGMASTOPt = 56.4;
const float gSIGMASTOPBARs = 1.76;
const float gSIGMASTOPBARt = 30.7;
const float LUMI = 19250;

using namespace std ;


//*****************************************************************************************

int readDataset (TString sampleName, vector<TString> datasetBaseName)
{
//-----------------Definitions------------------------------------------------------------

    float scale;
    //---Load data---
    TChain * ch = new TChain ("Delphes") ;
    for(int i=0; i<datasetBaseName.size(); i++)
    {
        ch->Add(datasetBaseName.at(i)+"*root") ;
    }
    cout << "read " << ch->GetEntries() << " events in " << ch->GetNtrees() 
         << " files of " << sampleName.Data() << " sample\n" ;
    if(sampleName.Contains("ttbar"))
	scale = gSIGMAttbar*LUMI/(ch->GetNtrees()*10000*ttbarBR);
    if(sampleName.Contains("W4jets"))
	scale = gSIGMAW4JETS*LUMI/(ch->GetNtrees()*50000);
    if(sampleName.Contains("Z4jets"))
	scale = gSIGMAZ4JETS*LUMI/(ch->GetNtrees()*10000);
    if(sampleName.Contains("stop_s"))
	scale = gSIGMASTOPs*LUMI/(ch->GetNtrees()*2000);
    if(sampleName.Contains("stop_t"))
	scale = gSIGMASTOPt*LUMI/(ch->GetNtrees()*20000);
    if(sampleName.Contains("stopbar_s"))
	scale = gSIGMASTOPBARs*LUMI/(ch->GetNtrees()*1400);
    if(sampleName.Contains("stopbar_t"))
	scale = gSIGMASTOPBARt*LUMI/(ch->GetNtrees()*20000);

    cout << scale << endl;
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
