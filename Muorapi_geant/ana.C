#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

void ana(){
gStyle->SetOptFit(1);
bool draw=true;

TChain *ttree;
ttree = new TChain("Hits");
char name_file[500];

sprintf(name_file,"output.root");
ttree->Add(name_file);
cerr<<name_file<<endl;

TTree *T = (TTree*)gROOT->FindObject("Hits");

ofstream qin;
qin.open("data_hits.txt",ios::out);
qin.setf(ios::showpoint|ios::fixed);
qin.precision(3);

Double_t fX;
Double_t fY;
Double_t fZ;
Int_t fDetNum;
Int_t fparentID;

T->SetBranchAddress("fX",    &fX);
T->SetBranchAddress("fY",    &fY);
T->SetBranchAddress("fZ",    &fZ);
T->SetBranchAddress("fDetNum", &fDetNum);
T->SetBranchAddress("fparentID", &fparentID);

for(int ientry=0; ientry<T->GetEntries(); ++ientry)
 {
  T->GetEntry(ientry);
  if (fparentID ==0) {
   cerr<<ientry<<"  "<<fX<<"  "<<fDetNum<<"  "<<fparentID<<endl;
   qin <<fDetNum<<" "<<fX<<" "<<fY<<" "<<fZ<<endl;
   }
 }


}
