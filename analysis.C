#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <time.h>
#include "TROOT.h"
#include "TStyle.h"
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TVector3.h"
#include "TMath.h"
#include "TH1I.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"
#include "TH1F.h"

using namespace std; 

void analysis()
{
 //variable declaration
 //4-Vector:
 TLorentzVector MyProton, MyPip, MyPim, MyProton_Rest, My_Initial, MyPi0, MyPhoton, My3Pion;

 //float and char variable declaration
 char name[500];
 float temp[50];
 float runnum, eventnum, top, beampol, realphi, cosrealtheta, cos_theta_p_cm, beamE, p_px, p_py, p_pz, p_E, pip_px, pip_py, pip_pz, pip_E, pim_px, pim_py, pim_pz, pim_E;
 float cl_nop, cl_pi0, z, kin_beamE, kin_p_px, kin_p_py, kin_p_pz, kin_p_E, kin_pip_px, kin_pip_py, kin_pip_pz, kin_pip_E, kin_pim_px, kin_pim_py, kin_pim_pz, kin_pim_E, q_value, q_error;
 float phistar, costetastar, degreepol; 
 float protonmass = .938272013;
 float three_pi_mass;

 //Histogram declaration
 TH1D *m_3pi = new TH1D("three pion mass","three pion mass",50, 0.65, 0.9);
 m_3pi->Sumw2(); 
 
 //for loadig text file
 ifstream inputfile;
 sprintf(name,"Jlab_data.txt");

 //reading the text file in the while and for loop
 inputfile.open(name,ios::in);
 while (!inputfile.eof()){
  for (int i = 0; i < 38 ; i++)
    { 
     inputfile >> temp[i]; 
    }
  runnum = temp[0];
  eventnum = temp[1];
  top = temp[2];
  beampol = temp[3];
  phistar = temp[4];
  costetastar = temp[5];
  degreepol = temp[6];
  beamE = temp[7]/1000; 
  p_px = temp[8]/1000 ;
  p_py = temp[9]/1000 ;
  p_pz = temp[10]/1000 ;
  p_E = temp[11]/1000 ;
  pip_px = temp[12]/1000 ;
  pip_py = temp[13]/1000 ;
  pip_pz = temp[14]/1000 ;
  pip_E = temp[15]/1000 ;
  pim_px = temp[16]/1000 ;
  pim_py = temp[17]/1000 ;
  pim_pz = temp[18]/1000 ;
  pim_E = temp[19]/1000 ;
  cl_nop = temp[20];
  cl_pi0 = temp[21];
  z = temp[22];
  kin_beamE = temp[23]/1000 ;
  kin_p_px = temp[24]/1000 ;
  kin_p_py = temp[25]/1000 ;
  kin_p_pz = temp[26]/1000 ;
  kin_p_E = temp[27]/1000 ;
  kin_pip_px = temp[28]/1000 ;
  kin_pip_py = temp[29]/1000 ;
  kin_pip_pz = temp[30]/1000 ;
  kin_pip_E = temp[31]/1000 ;
  kin_pim_px = temp[32]/1000 ;
  kin_pim_py = temp[33]/1000 ;
  kin_pim_pz = temp[34]/1000 ;
  kin_pim_E = temp[35]/1000 ;
  q_value = temp[36];
  q_error = temp[37];

  //construct the four vector of the final detected particle
  MyProton.SetPxPyPzE( kin_p_px,kin_p_py, kin_p_pz, kin_p_E );
  MyPip.SetPxPyPzE( kin_pip_px,kin_pip_py, kin_pip_pz, kin_pip_E );
  MyPim.SetPxPyPzE( kin_pim_px,kin_pim_py, kin_pim_pz, kin_pim_E);
 
  //intitial total 4-vector
  My_Initial.SetPxPyPzE(0.0,0.0,kin_beamE, kin_beamE + protonmass);

  //reconstruction of the missing particle 4-vector
  MyPi0 = My_Initial - MyProton - MyPip - MyPim;

  //omega reconstruction
  My3Pion = MyPip + MyPim + MyPi0;
  three_pi_mass = My3Pion.Mag();

  m_3pi->Fill(three_pi_mass);
 
 }

 //Draw and save the histogram
 TCanvas *c = new TCanvas();
 c->cd();
 m_3pi->Draw();
 c->SaveAs("histo_3pion.png");

 //Fitting
 TF1 *fSignal = new TF1("fSignal","gaus",0.65 , 0.9);
 TF1 *fBackground = new TF1("fBackground","pol3",0.65 , 0.9);
 TF1 *fSpectrum = new TF1("fSpectrum","gaus+pol3(3)",0.65, 0.9);

 fSpectrum->SetParameters(50000, 0.78, 0.1, 2000, 0, 0, 0);
 m_3pi->Fit("fSpectrum", "", "", 0.65, 0.9);

 Double_t param[7];
 fSpectrum->GetParameters(param);
 fSignal->SetParameters(&param[0]);
 fBackground->SetParameters(&param[3]);

 //plot in new canvas
 TCanvas *c2 = new TCanvas();
 c2->cd();
 m_3pi->Draw();
 fSpectrum->Draw("SAME");
 fBackground->SetLineColor(kBlack);
 fBackground->Draw("SAME");
 fSignal->SetLineColor(kBlue);
 fSignal->Draw("SAME");
 c2->SaveAs("histo_3pion_fit.png");

}


