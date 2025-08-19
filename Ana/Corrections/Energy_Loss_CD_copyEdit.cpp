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
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TDatabasePDG.h>
#include "HipoChain.h"
#include "clas12ana.h"
#include "reweighter.h"
#include "TGraphErrors.h"


using namespace std;
using namespace clas12;

const double c = 29.9792458;

auto db=TDatabasePDG::Instance();
double mass_p = db->GetParticle(2212)->Mass();
double mass_pi = db->GetParticle(-211)->Mass();
double mD = 1.8756;

double beam_E = 5.98;
const double me = 0.000511;
const double mU = 0.9314941024;
const double m_4He = 4.00260325415 * mU - 2*me;

//Parameters
double A[2]= {-0.00555986,-6.06201e-05};
double B[2]= {0.00695634,8.24535e-05};
double C[2]= {0.00155103,1.74283e-05};

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());
}

double G(double x, double N, double mu, double sigma){
  return (N/(sigma*sqrt(2*M_PI))) * exp(-0.5 * sq((x-mu)/sigma)) ; 
}

double func(double x, double a, double b, double c){
  return a + b*x + (c/x); 
}

double funcABC(double x, double a, double b){
  return a + b*x; 
}

double Corr(double mom, double theta){
  double a = funcABC(theta,A[0],A[1]);
  double b = funcABC(theta,B[0],B[1]);
  double c = funcABC(theta,C[0],C[1]);
  double correction = func(mom,a,b,c);
  return correction;
}

double getExp(TLorentzVector balance_ptr, TLorentzVector par){
  double theta_bpar = balance_ptr.Vect().Angle(par.Vect());
  double Eb = balance_ptr.E();
  double Pb = balance_ptr.P();
  double K = ((mass_pi*mass_pi) - balance_ptr.M2() - par.M2()) / 2;
  double a = Pb*Pb*cos(theta_bpar)*cos(theta_bpar) - Eb*Eb;
  double b = -2 * K * Pb * cos(theta_bpar);
  double c = K*K - Eb*Eb*par.M2();
  double x_min = (-b - sqrt(b*b - 4*a*c))/(2*a);
  double x_max = (-b + sqrt(b*b - 4*a*c))/(2*a);
  return x_min;
}

void getGraph(TH2D * h_myhist, TGraphErrors * g_mygraph, TCanvas * myCanvas, char fileName[100]){
  int ctr = 0;
  char temp[100];
  //Now project the histogram    
  for(int j = 0; j < h_myhist->GetXaxis()->GetNbins(); j++){
    //Define x and y(1D histogram)
    double x = h_myhist->GetXaxis()->GetBinCenter(j+1);
    ctr++;
    sprintf(temp,"Proj_num%d",ctr);
    TH1D * proj = h_myhist->ProjectionY(temp,j+1,j+1);
    //proj->Rebin(2);
    //Now preform a guassian fit
    if(proj->GetEntries()<15){continue;}

    TF1 * gFit = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.3,0.3,3);
    double mode = proj->GetBinCenter(proj->GetMaximumBin());
    gFit->SetParameter(0,proj->GetMaximum()/G(0,1,0,0.1));
    gFit->SetParameter(1,0);
    gFit->SetParLimits(1,-0.2,0.2);
    gFit->SetParameter(2,0.01);
    gFit->SetParLimits(2,0.001,0.4);
    
    TFitResultPtr gPoint = proj->Fit(gFit,"SrBeqn","",-0.2,0.2);

    /*
    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    proj->Draw();
    gFit->Draw("SAME");
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
    */

    if((gPoint == 0)){
      if(gPoint->Parameter(1)<0.07){
	g_mygraph->SetPoint(g_mygraph->GetN(),x,gPoint->Parameter(1));
	g_mygraph->SetPointError(g_mygraph->GetN()-1,0,gPoint->Parameter(2));
      }
    }
    proj->Write();
  }
}

void getFunctionMomTFRP(TH2D * h_myhist, TGraphErrors * g_mygraph, TF1 * f_myfunc, TFitResultPtr & p_mypoint, double min, double max,TCanvas * myCanvas, char fileName[100]){
  getGraph(h_myhist,g_mygraph,myCanvas,fileName);

  f_myfunc->SetLineColor(3);
  f_myfunc->SetLineWidth(1);
  f_myfunc->SetParameter(0,0);
  f_myfunc->SetParLimits(0,-0.2,0.2);
  f_myfunc->SetParameter(1,0);
  f_myfunc->SetParLimits(1,-0.2,0.2);
  f_myfunc->SetParameter(2,0);
  f_myfunc->SetParLimits(2,-0.5,+0.05);
  p_mypoint = g_mygraph->Fit(f_myfunc,"SrBeqn","",min,max);
      
}

void getABC(TH2D * h_myhist[16], TGraphErrors * g_mygraph[16], TF1 * f_myfunc[16], TFitResultPtr p_mypoint[16], double min, double max, TGraph * g_Pargraph[3], TF1 * f_Parfunc[3], TFitResultPtr p_Parpoint[3], TF1 * f_Combfunc[16], TCanvas * myCanvas, char fileName[100]){
  
  for(int i=0; i<16; i++){
    double theta = 40.75 + i*5.5;
    getFunctionMomTFRP(h_myhist[i],g_mygraph[i],f_myfunc[i],p_mypoint[i],min,max,myCanvas,fileName);
    for(int j=0; j<3; j++){
      if(p_mypoint[i]!=-1){
	g_Pargraph[j]->SetPoint(g_Pargraph[j]->GetN(),theta,p_mypoint[i]->Parameter(j));
      }
    }
  }
  
  for(int j=0; j<3; j++){
    f_Parfunc[j]->SetLineColor(4);
    f_Parfunc[j]->SetParameter(0,0);
    f_Parfunc[j]->SetParLimits(0,-0.01,0.01);
    f_Parfunc[j]->SetParameter(1,0);
    f_Parfunc[j]->SetParLimits(1,-0.01,0.01);
    p_Parpoint[j] = g_Pargraph[j]->Fit(f_Parfunc[j],"SrBeqn","",30,150);
  }  

  for(int i=0; i<16; i++){
    double theta = 40.75 + i*5.5;
    f_Combfunc[i]->SetLineColor(4);
    f_Combfunc[i]->SetLineWidth(1);
    f_Combfunc[i]->SetParameter(0,f_Parfunc[0]->Eval(theta));
    f_Combfunc[i]->SetParameter(1,f_Parfunc[1]->Eval(theta));
    f_Combfunc[i]->SetParameter(2,f_Parfunc[2]->Eval(theta));
  }
  
  /*
  myCanvas->Divide(4,4);
  for(int i=0; i<16; i++){
    myCanvas->cd(i+1);    
    h_myhist[i]->Draw("colz");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  myCanvas->Divide(4,4);
  for(int i=0; i<16; i++){
    myCanvas->cd(i+1);    
    h_myhist[i]->Draw("colz");
    g_mygraph[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  */
  TStyle *myStyle  = new TStyle("MyStyle","My Root Styles");
  myStyle->SetPalette("kbird",0);
  myStyle->SetTitleSize(0.10, "t");
  myStyle->SetOptStat(0);
  myStyle->cd();

  myCanvas->Divide(4,4);
  for(int i=0; i<16; i++){
    auto pad = myCanvas->cd(i+1);    
    pad->SetBottomMargin(0.19);
    pad->SetLeftMargin(0.19);
    h_myhist[i]->Draw("colz");
    g_mygraph[i]->Draw("SAME");
    f_myfunc[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  char temp_title[100];
  char* names[3] = {"C_{a}","C_{b}","C_{c}"};
  myCanvas->Divide(2,2);
  for(int j=0; j<3; j++){
    sprintf(temp_title,"%s vs. #theta #circ;#theta #circ;%s",names[j],names[j]);
    auto pad = myCanvas->cd(j+1);    
    pad->SetBottomMargin(0.19);
    pad->SetLeftMargin(0.19);
    g_Pargraph[j]->GetXaxis()->CenterTitle();
    g_Pargraph[j]->GetXaxis()->SetTitleSize(0.07);
    g_Pargraph[j]->GetXaxis()->SetTitleOffset(1.2);
    g_Pargraph[j]->GetYaxis()->CenterTitle();
    g_Pargraph[j]->GetYaxis()->SetTitleSize(0.07);
    g_Pargraph[j]->GetYaxis()->SetTitleOffset(1.2);
    g_Pargraph[j]->GetYaxis()->LabelsOption("v");
    g_Pargraph[j]->SetTitle(temp_title);
    g_Pargraph[j]->SetTitle(temp_title);
    g_Pargraph[j]->Draw();
    f_Parfunc[j]->Draw("SAME");    
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(4,4);
  for(int i=0; i<16; i++){
    auto pad = myCanvas->cd(i+1);    
    pad->SetBottomMargin(0.19);
    pad->SetLeftMargin(0.19);
    h_myhist[i]->Draw("colz");
    g_mygraph[i]->Draw("SAME");
    //f_myfunc[i]->Draw("SAME");
    f_Combfunc[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  


    
}


void Usage()
{
  std::cerr << "Usage: ./code outputfile.pdf inputfile.root \n\n\n";
}


int main(int argc, char ** argv)
{

  if(argc < 2)
    {
      Usage();
      return -1;
    }


  char * pdfFile = argv[1];
  cout<<"Ouput PDF file "<< pdfFile <<endl;

  TFile * inFile = new TFile(argv[2]);
  
  double mN = db->GetParticle(2212)->Mass();
  //some particles
  TLorentzVector beam(0,0,beam_E,beam_E);
  TLorentzVector target_ptr(0,0,0,mD);
  //TLorentzVector target_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector proton_ptr_p1(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector proton_ptr_p2(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector pim_ptr(0,0,0,db->GetParticle(-211)->Mass());
  reweighter newWeight(beam_E,6,6,kelly,"AV18");

  //////////////////////////////////
  char temp_name[100];
  char temp_title[100];

  
  vector<TH1*> hist_list;

  TH1D * h_theta = new TH1D("theta","#theta",100,0,180);
  hist_list.push_back(h_theta);
  
  TH2D * h_pCDmom_pCDDeltaP_int_thetaGroup[16];
  TGraphErrors * g_pCDmom_pCDDeltaP_int_thetaGroup[16];
  TF1 * f_pCDmom_pCDDeltaP_int_thetaGroup[16];
  TFitResultPtr p_pCDmom_pCDDeltaP_int_thetaGroup[16];
  TGraph * g_pCDmom_pCDDeltaP_int_Pars[3];
  TF1 * f_pCDmom_pCDDeltaP_int_Pars[3];
  TFitResultPtr p_pCDmom_pCDDeltaP_int_Pars[3];
  TF1 * f_pCDmom_pCDDeltaP_int_combined_thetaGroup[16];
  TH2D * h_pCDmom_pCDDeltaP_int_corrected_thetaGroup[16];
  for(int i=0; i<16; i++){
    double min = (double)i*5.5 + 38.0;
    double max = (double)i*5.5 + 43.5;
    sprintf(temp_name,"h_pCDmom_pCDDeltaP_int_thetaGroup_%d",i+1);
    h_pCDmom_pCDDeltaP_int_thetaGroup[i] = (TH2D*)inFile->Get(temp_name);
    sprintf(temp_title,"(%.1f #circ< #theta < %.1f #circ);p [GeV]; #Delta p_{E.L.} [GeV]",min,max);
    h_pCDmom_pCDDeltaP_int_thetaGroup[i]->SetTitle(temp_title);
    hist_list.push_back(h_pCDmom_pCDDeltaP_int_thetaGroup[i]);
    g_pCDmom_pCDDeltaP_int_thetaGroup[i] = new TGraphErrors();
    sprintf(temp_name,"g_pCDmom_pCDDeltaP_int_thetaGroup_%d",i+1);
    g_pCDmom_pCDDeltaP_int_thetaGroup[i]->SetName(temp_name);
    g_pCDmom_pCDDeltaP_int_thetaGroup[i]->SetLineColor(2);
    sprintf(temp_name,"f_pCDmom_pCDDeltaP_int_thetaGroup_%d",i+1);
    f_pCDmom_pCDDeltaP_int_thetaGroup[i] = new TF1(temp_name,[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);
    sprintf(temp_name,"f_pCDmom_pCDDeltaP_int_combined_thetaGroup_%d",i+1);
    f_pCDmom_pCDDeltaP_int_combined_thetaGroup[i] = new TF1(temp_name,[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);
    sprintf(temp_name,"h_pCDmom_pCDDeltaP_int_corrected_thetaGroup_%d",i+1);
    h_pCDmom_pCDDeltaP_int_corrected_thetaGroup[i] = (TH2D*)inFile->Get(temp_name);
    h_pCDmom_pCDDeltaP_int_corrected_thetaGroup[i]->SetTitle(temp_title);
    hist_list.push_back(h_pCDmom_pCDDeltaP_int_corrected_thetaGroup[i]);
  }
  for(int j=0; j<3; j++){
    g_pCDmom_pCDDeltaP_int_Pars[j] = new TGraph();
    sprintf(temp_name,"g_pCDmom_pCDDeltaP_int_Pars_%d",j+1);
    g_pCDmom_pCDDeltaP_int_Pars[j]->SetName(temp_name);
    g_pCDmom_pCDDeltaP_int_Pars[j]->SetLineColor(3);
    sprintf(temp_name,"f_pCDmom_pCDDeltaP_int_Pars_%d",j+1);
    f_pCDmom_pCDDeltaP_int_Pars[j] = new TF1(temp_name,[&](double *x, double *p){ return funcABC(x[0],p[0],p[1]); },30,150,2);
  }

  

  TH2D * h_pCDtheta_pCDDeltaP[6];
  TGraphErrors * g_pCDtheta_pCDDeltaP[6];
  for(int j=0; j<6; j++){
    sprintf(temp_name,"h_pCDtheta_pCDDeltaP_sector_%d",j+1);
    h_pCDtheta_pCDDeltaP[j] = (TH2D*)inFile->Get(temp_name);
    hist_list.push_back(h_pCDtheta_pCDDeltaP[j]);
    g_pCDtheta_pCDDeltaP[j] = new TGraphErrors();
    sprintf(temp_name,"g_pCDtheta_pCDDeltaP_sector_%d",j+1);
    g_pCDtheta_pCDDeltaP[j]->SetName(temp_name);
    g_pCDtheta_pCDDeltaP[j]->SetLineColor(2);
  }

  TH2D * h_pCDphi_pCDDeltaP[6];
  TGraphErrors * g_pCDphi_pCDDeltaP[6];
  for(int j=0; j<6; j++){
    sprintf(temp_name,"h_pCDphi_pCDDeltaP_sector_%d",j+1);
    h_pCDphi_pCDDeltaP[j] = (TH2D*)inFile->Get(temp_name);
    hist_list.push_back(h_pCDphi_pCDDeltaP[j]);
    g_pCDphi_pCDDeltaP[j] = new TGraphErrors();
    sprintf(temp_name,"g_pCDphi_pCDDeltaP_sector_%d",j+1);
    g_pCDphi_pCDDeltaP[j]->SetName(temp_name);
    g_pCDphi_pCDDeltaP[j]->SetLineColor(2);
  }

  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Sumw2();
    hist_list[i]->SetTitleFont(0.12);
    hist_list[i]->GetXaxis()->CenterTitle();
    hist_list[i]->GetXaxis()->SetTitleSize(0.12);
    hist_list[i]->GetXaxis()->SetLabelSize(0.06);
    hist_list[i]->GetXaxis()->SetTitleOffset(0.6);
    hist_list[i]->GetYaxis()->CenterTitle();
    hist_list[i]->GetYaxis()->SetTitleSize(0.12);
    hist_list[i]->GetYaxis()->SetLabelSize(0.06);
    hist_list[i]->GetYaxis()->SetTitleOffset(0.6);
  }
  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////
  

  int pixelx = 1980;
  int pixely = 1530;
  TCanvas * myCanvas = new TCanvas("myPage","myPage",pixelx,pixely);
  TCanvas * myText = new TCanvas("myText","myText",pixelx,pixely);
  TLatex text;
  text.SetTextSize(0.05);
  
  char fileName[100];
  sprintf(fileName,"%s[",pdfFile);
  myText->SaveAs(fileName);
  sprintf(fileName,"%s",pdfFile);
  
  getABC(h_pCDmom_pCDDeltaP_int_thetaGroup,g_pCDmom_pCDDeltaP_int_thetaGroup,f_pCDmom_pCDDeltaP_int_thetaGroup,p_pCDmom_pCDDeltaP_int_thetaGroup,0.3,5.5,g_pCDmom_pCDDeltaP_int_Pars,f_pCDmom_pCDDeltaP_int_Pars,p_pCDmom_pCDDeltaP_int_Pars,f_pCDmom_pCDDeltaP_int_combined_thetaGroup,myCanvas,fileName);

  TStyle *myStyle  = new TStyle("MyStyle","My Root Styles");
  myStyle->SetPalette("kbird",0);
  myStyle->SetTitleSize(0.10, "t");
  myStyle->SetOptStat(0);
  myStyle->cd();

  
  myCanvas->Divide(4,4);
  for(int i=0; i<16; i++){
    auto pad = myCanvas->cd(i+1);    
    pad->SetBottomMargin(0.19);
    pad->SetLeftMargin(0.19);
    h_pCDmom_pCDDeltaP_int_thetaGroup[i]->Draw("colz");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear(); 

  myCanvas->Divide(4,4);
  for(int i=0; i<16; i++){
    auto pad = myCanvas->cd(i+1);    
    pad->SetBottomMargin(0.19);
    pad->SetLeftMargin(0.19);
    h_pCDmom_pCDDeltaP_int_corrected_thetaGroup[i]->Draw("colz");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear(); 

  
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  cout<<"A=["<< p_pCDmom_pCDDeltaP_int_Pars[0]->Parameter(0)
      <<","<< p_pCDmom_pCDDeltaP_int_Pars[0]->Parameter(1)
      <<"]"<<endl
      <<"B=["<< p_pCDmom_pCDDeltaP_int_Pars[1]->Parameter(0)
      <<","<< p_pCDmom_pCDDeltaP_int_Pars[1]->Parameter(1)
      <<"]"<<endl
      <<"C=["<< p_pCDmom_pCDDeltaP_int_Pars[2]->Parameter(0)
      <<","<< p_pCDmom_pCDDeltaP_int_Pars[2]->Parameter(1)
      <<"]"<<endl;
    
  return 0;
}

