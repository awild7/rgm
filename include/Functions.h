#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>
#include <typeinfo>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TDatabasePDG.h>
#include "HipoChain.h"
#include "clas12ana.h"
#include "TRandom3.h"


//Define some masses and energies
auto db=TDatabasePDG::Instance();
const double me = 0.000511;
const double mass_p = db->GetParticle(2212)->Mass();
const double mass_n = db->GetParticle(2112)->Mass();
const double mD = 1.8756;
const double mU = 0.9314941024;
const double m_4He = 4.00260325415 * mU - 2*me;
const double beam_E_6 = 5.98636;
const double beam_E_6_sigma = 0.00299;


double sq(double x){return x*x;}

double G(double x, double N, double mu, double sigma){return (N/(sigma*sqrt(2*M_PI))) * exp(-0.5 * sq((x-mu)/sigma)) ; }

int binX(vector<double> XS, double X){
  for(int i = 0; i < XS.size(); i++){
    if(X<XS[i]){
      return i-1;
    }
  }
  return -1;
}

double getExpProton(double eMom, double eTheta, double ePhi, double pTheta, double pPhi){
  TVector3 v3e;
  v3e.SetMagThetaPhi(eMom,eTheta,ePhi);
  TLorentzVector vLe;
  vLe.SetXYZM(v3e.X(),v3e.Y(),v3e.Z(),me);

  TLorentzVector beam(0,0,beam_E_6,beam_E_6);
  TLorentzVector deut_ptr(0,0,0,mD);
  TLorentzVector balance_ptr = beam + deut_ptr - vLe;

  TVector3 vp;
  vp.SetMagThetaPhi(1,pTheta,pPhi);
  
  double theta_bpar = balance_ptr.Vect().Angle(vp);
  double Eb = balance_ptr.E();
  double Pb = balance_ptr.P();
  double K = ((mass_n*mass_n) - balance_ptr.M2() - mass_p*mass_p) / 2;
  double a = Pb*Pb*cos(theta_bpar)*cos(theta_bpar) - Eb*Eb;
  double b = -2 * K * Pb * cos(theta_bpar);
  double c = K*K - Eb*Eb*mass_p*mass_p;
  double x_min = (-b - sqrt(b*b - 4*a*c))/(2*a);
  double x_max = (-b + sqrt(b*b - 4*a*c))/(2*a);
  return x_min;
}


#endif

