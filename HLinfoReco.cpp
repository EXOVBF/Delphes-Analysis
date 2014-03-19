/****************************************************************************************
	- before compiling ---> source ../Decay/setup_slc6.sh

    - compile with ---> c++ -O2 -lm `root-config --cflags --glibs` -o HLinfoReco HLinfoReco.cpp

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

using namespace std;

//****************************************************************************************
// main

int main (int argc, char *argv[])
{
    cout << "HELLO" << endl;


}

