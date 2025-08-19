#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>
#include <typeinfo>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TDatabasePDG.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include "HipoChain.h"
//#include "clas12ana.h"

using namespace std;
using namespace clas12;

//Define some masses and energies
auto db=TDatabasePDG::Instance();
double beam_E = 5.984792;
double beam_E_sigma = 0.00299;
double mN = db->GetParticle(2212)->Mass();

//Define some useful funcitons and phi binning
double sq(double x){return x*x;}
double G(double x, double N, double mu, double sigma){return (N/(sigma*sqrt(2*M_PI))) * exp(-0.5 * sq((x-mu)/sigma)) ; }
vector<double> bE_PhiCD = {-180,-160,-140,-120,-100,-80,-60,-40,-20,0,20,40,60,80,100,120,140,160,180};
int binX(vector<double> XS, double X){
  for(int i = 0; i <= XS.size(); i++){if(X<XS[i]){return i-1;}}
  return -1;
}
void GetLorentzVector_ReconVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());
}


int main(int argc, char ** argv)
{
  if(argc < 2)
    {
      std::cerr << "Usage: ./code outputfile.pdf inputfiles.hipo \n\n\n";
      return -1;
    }
  //Output pdf file
  char * pdfFile = argv[1];
  cout<<"Ouput PDF file "<< pdfFile <<endl;
  //Input hipo files
  clas12root::HipoChain chain;
  for(int k = 2; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12=chain.GetC12Reader();
  auto &c12=chain.C12ref();

  //Create a clas12ana object
  //This handles PID, fiducial cuts, and vertex cuts
  //clas12ana clasAna;
  
  //Lorentz vectors to be used for calculations
  TLorentzVector beam(0,0,beam_E,beam_E);
  TVector3 vbeam(0,0,beam_E);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector proton_ptr(0,0,0,db->GetParticle(2212)->Mass());

  //Define histograms  
  TH1D * h_W=new TH1D("W","W;W;Counts",100,0.5,1.5);
  TH1D * h_Dphi=new TH1D("Dphi","#Delta #phi;#Delta #phi;Counts",100,-20,20);
  TH1D * h_Dtheta=new TH1D("Dtheta","#Delta #theta;#Delta #theta;Counts",100,0,20);
  char temp_name[100];
  char temp_title[100];
  TH1D * h_corr_binPhiCD[18];
  for(int i=0; i<18; i++){
    int min = bE_PhiCD[i];
    int max = bE_PhiCD[i+1];
    sprintf(temp_name,"corr_phi_%d",i);
    sprintf(temp_title,"#Delta p (%d< #phi < %d);#Delta p;Counts",min,max);
    h_corr_binPhiCD[i] = new TH1D(temp_name,temp_title,50,-0.35,0.35);
  }

  int counter = 0;
  while(chain.Next() && (counter<10000000000000))
    {
      //Display completed  
      counter++;
      if((counter%100000) == 0){
	cerr << "\n" <<counter/100000 <<" hundred thousand completed";
      }    

      //Get particles that pass PID, fiducial, and vertex cuts
      //I usually use the RGM repo which containtsall of these functions
      //clasAna.Run(c12);
      //auto electrons = clasAna.getByPid(11);
      //auto protons = clasAna.getByPid(2212);
      //For simplicity, we can just grab CLAS12 electrons and protons
      auto electrons=c12->getByID(11);
      auto protons=c12->getByID(2212);
      if(electrons.size() == 1 && protons.size() > 0)
	{
	  //Fill TLorentzVectors
	  GetLorentzVector_ReconVector(el,electrons[0]);	  
	  GetLorentzVector_ReconVector(proton_ptr,protons[0]);

	  //Get important electron variables
	  TLorentzVector q = beam - el;
          double Q2        = -q.M2();
	  double omega = q.E();
          double xB       = Q2/(2 * mN * (beam.E() - el.E()) );
	  double WSq = (mN*mN) - Q2 + (2*omega*mN);
	  double W = sqrt(WSq);
	  double theta_q = q.Theta()*180/M_PI;

	  //Get the proton variables
	  double mom_p = proton_ptr.P();
	  double theta_p = proton_ptr.Theta()*180/M_PI;
	  double phi_p = proton_ptr.Phi()*180/M_PI;

	  //Calculate the expected elastic electron
	  //and proton momentum from angles.
	  double k_exp = ( mN*beam_E )/( beam_E*(1-cos(el.Theta())) + mN );
	  TVector3 vk_exp;
	  vk_exp.SetMagThetaPhi(k_exp,el.Theta(),el.Phi());
	  TVector3 vp_exp = vbeam-vk_exp;

	  //Calculate the resolution of the proton from scattering angles
	  double pres = (vp_exp.Mag()-mom_p)/vp_exp.Mag();

	  //We will use this to help select QE events.
	  double Delta_theta = (proton_ptr.Theta()-vp_exp.Theta())*180/M_PI;
	  double Delta_phi = (proton_ptr.Phi()-vp_exp.Phi())*180/M_PI;
	  if(Delta_phi<-180){Delta_phi+=360;}
	  else if(Delta_phi>180){Delta_phi-=360;}

	  //Vertex (only necessary since we are not using clasAna)
	  double vtz_e = electrons[0]->par()->getVz();
	  double vtz_p = protons[0]->par()->getVz();

	  //We limit protons to the forward part of the CD
	  //since this is where the QE protons go
	  if(protons[0]->getRegion()!=CD){continue;}
	  if(!((theta_p<62) && (theta_p>38))){continue;}
	  //Vertex Cuts
	  if(vtz_e<-6){continue;}
	  if(vtz_e<0){continue;}
	  if(fabs(vtz_e-vtz_p)>3){continue;}
	  //QE cuts
	  if(W>1.05){continue;}
	  if(fabs(Delta_phi)>3){continue;}	    
	  if(Delta_theta>4){continue;}

	  h_W->Fill(W,1.0);
	  h_Dphi->Fill(Delta_phi,1.0);
	  h_Dtheta->Fill(Delta_theta,1.0);

	  if(binX(bE_PhiCD,phi_p)==-1){continue;}
	  h_corr_binPhiCD[binX(bE_PhiCD,phi_p)]->Fill(pres,1.0);
	}
    }
  
  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////
  TCanvas * myCanvas = new TCanvas("myPage","myPage",1980,1530);  
  char fileName[100];
  sprintf(fileName,"%s[",pdfFile);
  myCanvas->SaveAs(fileName);
  sprintf(fileName,"%s",pdfFile);

  ///////////////////////////////////

  //Draw relavent distributions
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_W->Draw();
  myCanvas->cd(2);
  h_Dphi->Draw();
  myCanvas->cd(3);
  h_Dtheta->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  //Draw resolution distributions
  myCanvas->Divide(5,4);
  for(int i = 0; i < 18; i++){
    myCanvas->cd(i+1);
    h_corr_binPhiCD[i]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  //Draw resolution distributions with fits
  TGraph * g_Phi_Data_sigma_CD = new TGraph;
  myCanvas->Divide(5,4);
  for(int i = 0; i < 18; i++){
    myCanvas->cd(i+1);
    h_corr_binPhiCD[i]->Draw();     
    TF1 * f_phibin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.35,0.35,3);
    f_phibin->SetParameter(0,h_corr_binPhiCD[i]->GetMaximum());
    f_phibin->SetParameter(1,0);
    f_phibin->SetParameter(2,0.1);
    f_phibin->SetParLimits(2,0.0,0.15);
    TFitResultPtr point = h_corr_binPhiCD[i]->Fit(f_phibin,"SrBeqn","",-0.35,0.35);
    f_phibin->Draw("SAME");
    if(h_corr_binPhiCD[i]->GetEntries()<30){continue;}
    double x = (bE_PhiCD[i]+bE_PhiCD[i+1])/2;
    if((point!=-1)){
      g_Phi_Data_sigma_CD->SetPoint(g_Phi_Data_sigma_CD->GetN(),x,100*point->Parameter(2));
    }
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  //Draw momentum resolution as a function of phi
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->GetPad(1)->SetBottomMargin(0.19);
  myCanvas->GetPad(1)->SetLeftMargin(0.19);
  g_Phi_Data_sigma_CD->SetLineWidth(2);
  g_Phi_Data_sigma_CD->SetLineColor(2);
  g_Phi_Data_sigma_CD->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  //Close pdf file
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  return 0;
}
