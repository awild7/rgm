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
double A[2]= {-0.000695124,-0.000355869};
double B[2]= {0.00182181,7.77933e-05};
double C[3]= {0.000266541,0.424055,49.068};



void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());
}

double G(double x, double N, double mu, double sigma){
  return (N/(sigma*sqrt(2*M_PI))) * exp(-0.5 * sq((x-mu)/sigma)) ; 
}

double func(double x, double a, double b, double c){
  return a + b*x + (c/x); 
}

double funcAB(double x, double a, double b){
  return a + b*x; 
}

double funcC(double x, double a, double b, double c){
  return a - b/(x-c); 
}

double Corr(double mom, double theta){
  double a = funcAB(theta,A[0],A[1]);
  double b = funcAB(theta,B[0],B[1]);
  double c = funcC(theta,C[0],C[1],C[2]);
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
    gFit->SetParameter(1,mode);
    gFit->SetParLimits(1,mode-0.025,mode+0.025);
    gFit->SetParameter(2,0.01);
    gFit->SetParLimits(2,0.001,0.1);
    
    TFitResultPtr gPoint = proj->Fit(gFit,"SrBeqn","",mode-0.05,0.05);

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
  f_myfunc->SetParLimits(0,-0.2,+0.2);
  f_myfunc->SetParameter(1,0);
  f_myfunc->SetParLimits(1,0,+0.01);
  f_myfunc->SetParameter(2,0);
  f_myfunc->SetParLimits(2,-0.5,+0.05);
  p_mypoint = g_mygraph->Fit(f_myfunc,"SrBeqn","",min,max);
      
}

void getABC(TH2D * h_myhist[16], TGraphErrors * g_mygraph[16], TF1 * f_myfunc[16], TFitResultPtr p_mypoint[16], double min, double max, TGraphErrors * g_Pargraph[3], TF1 * f_Parfunc[3], TFitResultPtr p_Parpoint[3], TF1 * f_Combfunc[16], TCanvas * myCanvas, char fileName[100]){

  TGraphErrors * g_Pargraph_clone[3];
  for(int j=0; j<3; j++){
    g_Pargraph_clone[j]=(TGraphErrors*)g_Pargraph[j]->Clone("clone");
    g_Pargraph_clone[j]->SetLineColor(3);
  }
  
  for(int i=0; i<16; i++){
    double theta = 4.25 + i*2.5;
    getFunctionMomTFRP(h_myhist[i],g_mygraph[i],f_myfunc[i],p_mypoint[i],min,max,myCanvas,fileName);
    for(int j=0; j<3; j++){
      if(p_mypoint[i]!=-1){
	g_Pargraph[j]->SetPoint(g_Pargraph[j]->GetN(),theta,p_mypoint[i]->Parameter(j));
	//g_Pargraph[j]->SetPointError(g_Pargraph[j]->GetN()-1,0,p_mypoint[i]->ParError(j));

	g_Pargraph_clone[j]->SetPoint(g_Pargraph_clone[j]->GetN(),theta,p_mypoint[i]->Parameter(j));
	g_Pargraph_clone[j]->SetPointError(g_Pargraph_clone[j]->GetN()-1,0,p_mypoint[i]->ParError(j));
      }
    }
  }
  
  for(int j=0; j<2; j++){
    f_Parfunc[j]->SetLineColor(4);
    f_Parfunc[j]->SetParameter(0,0);
    f_Parfunc[j]->SetParLimits(0,-0.01,0.01);
    f_Parfunc[j]->SetParameter(1,0);
    f_Parfunc[j]->SetParLimits(1,-0.01,0.01);
    p_Parpoint[j] = g_Pargraph[j]->Fit(f_Parfunc[j],"SrBeqn","",5,41);
  }  
  f_Parfunc[2]->SetLineColor(4);
  f_Parfunc[2]->SetParameter(0,0);
  f_Parfunc[2]->SetParLimits(0,-0.1,0.1);
  f_Parfunc[2]->SetParameter(1,0.01);
  f_Parfunc[2]->SetParLimits(1,0,10);
  f_Parfunc[2]->SetParameter(2,60);
  f_Parfunc[2]->SetParLimits(2,45,100);
  p_Parpoint[2] = g_Pargraph[2]->Fit(f_Parfunc[2],"SrBeqn","",5,41);

  for(int i=0; i<16; i++){
    double theta = 4.25 + i*2.5;
    f_Combfunc[i]->SetLineColor(4);
    f_Combfunc[i]->SetLineWidth(1);
    f_Combfunc[i]->SetParameter(0,f_Parfunc[0]->Eval(theta));
    f_Combfunc[i]->SetParameter(1,f_Parfunc[1]->Eval(theta));
    f_Combfunc[i]->SetParameter(2,f_Parfunc[2]->Eval(theta));
  }
  

  //myStyle->SetPadLeftMargin(1.11);
  //myStyle->SetTitleOffset(1.2);
  //myStyle->SetLabelFont(1182, "xyz"); // size of axis values
  //myStyle->SetNdivisions(4, "xyz");
  //myCanvas->SetBottomMargin(1.1);
  //myCanvas->SetLeftMargin(0.7);

  TStyle *myStyle  = new TStyle("MyStyle","My Root Styles");
  myStyle->SetPalette("kbird",0);
  myStyle->SetTitleSize(0.10, "t");
  myStyle->SetOptStat(0);
  myStyle->cd();

  myCanvas->Divide(4,4,0,0);
  for(int i=1; i<15; i++){
    auto pad = myCanvas->cd(i);    
    pad->SetBottomMargin(0.19);
    pad->SetLeftMargin(0.19);
    h_myhist[i]->Draw("colz");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  myCanvas->Divide(4,4,0,0);
  for(int i=1; i<15; i++){
    auto pad = myCanvas->cd(i);    
    pad->SetBottomMargin(0.19);
    pad->SetLeftMargin(0.19);
    h_myhist[i]->Draw("colz");
    g_mygraph[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  myCanvas->Divide(4,4,0,0);
  for(int i=1; i<15; i++){
    auto pad = myCanvas->cd(i);    
    pad->SetBottomMargin(0.19);
    pad->SetLeftMargin(0.19);
    h_myhist[i]->Draw("colz");
    g_mygraph[i]->Draw("SAME");
    f_myfunc[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  myCanvas->Divide(4,4,0,0);
  for(int i=1; i<15; i++){
    auto pad = myCanvas->cd(i);    
    pad->SetBottomMargin(0.19);
    pad->SetLeftMargin(0.19);
    h_myhist[i]->Draw("colz");
    g_mygraph[i]->Draw("SAME");
    //f_myfunc[i]->Draw("SAME");
    f_Combfunc[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  


  myCanvas->Divide(2,2,0,0);
  char temp_title[100];
  char* names[3] = {"C_{a}","C_{b}","C_{c}"};
  for(int j=0; j<3; j++){
    auto pad = myCanvas->cd(j+1);    
    pad->SetBottomMargin(0.19);
    pad->SetLeftMargin(0.19);
    myCanvas->cd(j+1);
    sprintf(temp_title,"%s vs. #theta #circ;#theta #circ;%s",names[j],names[j]);
    g_Pargraph[j]->GetXaxis()->CenterTitle();
    g_Pargraph[j]->GetXaxis()->SetTitleSize(0.07);
    g_Pargraph[j]->GetXaxis()->SetTitleOffset(1.2);
    g_Pargraph[j]->GetYaxis()->CenterTitle();
    g_Pargraph[j]->GetYaxis()->SetTitleSize(0.07);
    g_Pargraph[j]->GetYaxis()->SetTitleOffset(1.2);
    g_Pargraph[j]->GetYaxis()->LabelsOption("v");
    g_Pargraph[j]->SetTitle(temp_title);
    g_Pargraph[j]->Draw();
    f_Parfunc[j]->Draw("SAME");    
    g_Pargraph_clone[j]->Draw("SAME");
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
  vector<TH1*> hist_list_2;


  //////////////////////////////////
  //Proton FD
  //////////////////////////////////
  TH2D * h_pFDmom_pFDDeltaP_thetaGroup[6][16];
  TGraphErrors * g_pFDmom_pFDDeltaP_thetaGroup[6][16];
  TF1 * f_pFDmom_pFDDeltaP_thetaGroup[6][16];
  TFitResultPtr p_pFDmom_pFDDeltaP_thetaGroup[6][16];
  TGraphErrors * g_pFDmom_pFDDeltaP_Pars[6][3];
  TF1 * f_pFDmom_pFDDeltaP_Pars[6][3];
  TFitResultPtr p_pFDmom_pFDDeltaP_Pars[6][3];
  TF1 * f_pFDmom_pFDDeltaP_combined_thetaGroup[6][16];
  TH2D * h_pFDmom_pFDDeltaP_corrected_thetaGroup[6][16];
  for(int k=0; k<6; k++){
    for(int i=0; i<16; i++){
      double min = (double)i*2.5 + 3.0;
      double max = (double)i*2.5 + 5.5;
      sprintf(temp_name,"h_pFDmom_pFDDeltaP_sector_%d_thetaGroup_%d",k+1,i+1);
      h_pFDmom_pFDDeltaP_thetaGroup[k][i] = (TH2D*)inFile->Get(temp_name);
      //sprintf(temp_title,"#Delta p_{E.L.} vs. p  Sector = %d (%.1f #circ< #theta < %.1f #circ);p [GeV]; #Delta p_{E.L.}",k+1,min,max);
      sprintf(temp_title,"Sector = %d (%.1f #circ< #theta < %.1f #circ);p [GeV]; #Delta p_{E.L.} [GeV]",k+1,min,max);
      h_pFDmom_pFDDeltaP_thetaGroup[k][i]->SetTitle(temp_title);
      hist_list.push_back(h_pFDmom_pFDDeltaP_thetaGroup[k][i]);
      g_pFDmom_pFDDeltaP_thetaGroup[k][i] = new TGraphErrors();
      sprintf(temp_name,"g_pFDmom_pFDDeltaP_sector_%d_thetaGroup_%d",k+1,i+1);
      g_pFDmom_pFDDeltaP_thetaGroup[k][i]->SetName(temp_name);
      g_pFDmom_pFDDeltaP_thetaGroup[k][i]->SetLineColor(2);
      sprintf(temp_name,"f_pFDmom_pFDDeltaP_sector_%d_thetaGroup_%d",k+1,i+1);
      f_pFDmom_pFDDeltaP_thetaGroup[k][i] = new TF1(temp_name,[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);
      sprintf(temp_name,"f_pFDmom_pFDDeltaP_combined_sector_%d_thetaGroup_%d",k+1,i+1);
      f_pFDmom_pFDDeltaP_combined_thetaGroup[k][i] = new TF1(temp_name,[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);
      sprintf(temp_name,"h_pFDmom_pFDDeltaP_corrected_sector_%d_thetaGroup_%d",k+1,i+1);
      h_pFDmom_pFDDeltaP_corrected_thetaGroup[k][i] = (TH2D*)inFile->Get(temp_name);
      //sprintf(temp_title,"Corrected #Delta p_{E.L.} vs. p  Sector = %d (%.1f #circ< #theta < %.1f #circ);p [GeV]; #Delta p_{E.L.}",k+1,min,max);
      sprintf(temp_title,"Sector = %d (%.1f #circ< #theta < %.1f #circ);p [GeV]; #Delta p_{E.L.} [GeV]",k+1,min,max);
      h_pFDmom_pFDDeltaP_corrected_thetaGroup[k][i]->SetTitle(temp_title);
      hist_list.push_back(h_pFDmom_pFDDeltaP_corrected_thetaGroup[k][i]);
    }
    for(int j=0; j<3; j++){
      g_pFDmom_pFDDeltaP_Pars[k][j] = new TGraphErrors();
      sprintf(temp_name,"g_pFDmom_pFDDeltaP_sector_%d_Pars_%d",k+1,j+1);
      g_pFDmom_pFDDeltaP_Pars[k][j]->SetName(temp_name);
      g_pFDmom_pFDDeltaP_Pars[k][j]->SetLineColor(3);
      sprintf(temp_name,"f_pFDmom_pFDDeltaP_sector_%d_Pars_%d",k+1,j+1);
      if(j<2){
	f_pFDmom_pFDDeltaP_Pars[k][j] = new TF1(temp_name,[&](double *x, double *p){ return funcAB(x[0],p[0],p[1]); },5,41,2);
      }
      else{
	f_pFDmom_pFDDeltaP_Pars[k][j] = new TF1(temp_name,[&](double *x, double *p){ return funcC(x[0],p[0],p[1],p[2]); },5,41,3);
      }
    }

    for(int i=0; i<16; i++){
      sprintf(temp_name,"h_pFDmom_pFDDeltaP_sector_%d_thetaGroup_%d",k+1,i+1);
      h_pFDmom_pFDDeltaP_thetaGroup[k][i] = (TH2D*)inFile->Get(temp_name);
      hist_list.push_back(h_pFDmom_pFDDeltaP_thetaGroup[k][i]);
      g_pFDmom_pFDDeltaP_thetaGroup[k][i] = new TGraphErrors();
      sprintf(temp_name,"g_pFDmom_pFDDeltaP_sector_%d_thetaGroup_%d",k+1,i+1);
      g_pFDmom_pFDDeltaP_thetaGroup[k][i]->SetName(temp_name);
      g_pFDmom_pFDDeltaP_thetaGroup[k][i]->SetLineColor(2);
      sprintf(temp_name,"f_pFDmom_pFDDeltaP_sector_%d_thetaGroup_%d",k+1,i+1);
      f_pFDmom_pFDDeltaP_thetaGroup[k][i] = new TF1(temp_name,[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);
      
    }
  }
  
  TH2D * h_pFDmom_pFDDeltaP_int_thetaGroup[16];
  TGraphErrors * g_pFDmom_pFDDeltaP_int_thetaGroup[16];
  TF1 * f_pFDmom_pFDDeltaP_int_thetaGroup[16];
  TFitResultPtr p_pFDmom_pFDDeltaP_int_thetaGroup[16];
  TGraphErrors * g_pFDmom_pFDDeltaP_int_Pars[3];
  TF1 * f_pFDmom_pFDDeltaP_int_Pars[3];
  TFitResultPtr p_pFDmom_pFDDeltaP_int_Pars[3];
  TF1 * f_pFDmom_pFDDeltaP_int_combined_thetaGroup[16];
  TH2D * h_pFDmom_pFDDeltaP_int_corrected_thetaGroup[16];
  for(int i=0; i<16; i++){
    double min = (double)i*2.5 + 3.0;
    double max = (double)i*2.5 + 5.5;
    sprintf(temp_name,"h_pFDmom_pFDDeltaP_int_thetaGroup_%d",i+1);
    h_pFDmom_pFDDeltaP_int_thetaGroup[i] = (TH2D*)inFile->Get(temp_name);
    //sprintf(temp_title,"#Delta p_{E.L.} vs. p (%.1f #circ< #theta < %.1f #circ);p [GeV]; #Delta p_{E.L.}",min,max);
    sprintf(temp_title,"(%.1f #circ< #theta < %.1f #circ);p [GeV]; #Delta p_{E.L.} [GeV]",min,max);
    h_pFDmom_pFDDeltaP_int_thetaGroup[i]->SetTitle(temp_title);
    hist_list.push_back(h_pFDmom_pFDDeltaP_int_thetaGroup[i]);
    g_pFDmom_pFDDeltaP_int_thetaGroup[i] = new TGraphErrors();
    sprintf(temp_name,"g_pFDmom_pFDDeltaP_int_thetaGroup_%d",i+1);
    g_pFDmom_pFDDeltaP_int_thetaGroup[i]->SetName(temp_name);
    g_pFDmom_pFDDeltaP_int_thetaGroup[i]->SetLineColor(2);
    sprintf(temp_name,"f_pFDmom_pFDDeltaP_int_thetaGroup_%d",i+1);
    f_pFDmom_pFDDeltaP_int_thetaGroup[i] = new TF1(temp_name,[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);
    sprintf(temp_name,"f_pFDmom_pFDDeltaP_int_combined_thetaGroup_%d",i+1);
    f_pFDmom_pFDDeltaP_int_combined_thetaGroup[i] = new TF1(temp_name,[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);
    sprintf(temp_name,"h_pFDmom_pFDDeltaP_int_corrected_thetaGroup_%d",i+1);
    h_pFDmom_pFDDeltaP_int_corrected_thetaGroup[i] = (TH2D*)inFile->Get(temp_name);
    //sprintf(temp_title,"Corrected #Delta p_{E.L.} vs. p (%.1f #circ< #theta < %.1f #circ);p [GeV]; #Delta p_{E.L.}",min,max);
    sprintf(temp_title,"(%.1f #circ< #theta < %.1f #circ);p [GeV]; #Delta p_{E.L.} [GeV]",min,max);
    h_pFDmom_pFDDeltaP_int_corrected_thetaGroup[i]->SetTitle(temp_title);
    hist_list.push_back(h_pFDmom_pFDDeltaP_int_corrected_thetaGroup[i]);
  }
  for(int j=0; j<3; j++){
    g_pFDmom_pFDDeltaP_int_Pars[j] = new TGraphErrors();
    sprintf(temp_name,"g_pFDmom_pFDDeltaP_int_Pars_%d",j+1);
    g_pFDmom_pFDDeltaP_int_Pars[j]->SetName(temp_name);
    g_pFDmom_pFDDeltaP_int_Pars[j]->SetLineColor(3);
    sprintf(temp_name,"f_pFDmom_pFDDeltaP_int_Pars_%d",j+1);
    if(j<2){
      f_pFDmom_pFDDeltaP_int_Pars[j] = new TF1(temp_name,[&](double *x, double *p){ return funcAB(x[0],p[0],p[1]); },5,41,2);
    }
    else{
      f_pFDmom_pFDDeltaP_int_Pars[j] = new TF1(temp_name,[&](double *x, double *p){ return funcC(x[0],p[0],p[1],p[2]); },5,41,3);
    }
  }
  

  TH2D * h_pFDtheta_pFDDeltaP[6];
  TGraphErrors * g_pFDtheta_pFDDeltaP[6];
  for(int j=0; j<6; j++){
    sprintf(temp_name,"h_pFDtheta_pFDDeltaP_sector_%d",j+1);
    h_pFDtheta_pFDDeltaP[j] = (TH2D*)inFile->Get(temp_name);
    sprintf(temp_title,"Sector = %d;#theta #circ; #Delta p_{E.L.} [GeV]",j+1);
    h_pFDtheta_pFDDeltaP[j]->SetTitle(temp_title);
    hist_list_2.push_back(h_pFDtheta_pFDDeltaP[j]);
    g_pFDtheta_pFDDeltaP[j] = new TGraphErrors();
    sprintf(temp_name,"g_pFDtheta_pFDDeltaP_sector_%d",j+1);
    g_pFDtheta_pFDDeltaP[j]->SetName(temp_name);
    g_pFDtheta_pFDDeltaP[j]->SetLineColor(2);
  }

  TH2D * h_pFDphi_pFDDeltaP[6];
  TGraphErrors * g_pFDphi_pFDDeltaP[6];
  for(int j=0; j<6; j++){
    sprintf(temp_name,"h_pFDphi_pFDDeltaP_sector_%d",j+1);
    h_pFDphi_pFDDeltaP[j] = (TH2D*)inFile->Get(temp_name);
    sprintf(temp_title,"Sector = %d; #phi #circ ; #Delta p_{E.L.} [GeV]",j+1);
    h_pFDphi_pFDDeltaP[j]->SetTitle(temp_title);
    hist_list_2.push_back(h_pFDphi_pFDDeltaP[j]);
    g_pFDphi_pFDDeltaP[j] = new TGraphErrors();
    sprintf(temp_name,"g_pFDphi_pFDDeltaP_sector_%d",j+1);
    g_pFDphi_pFDDeltaP[j]->SetName(temp_name);
    g_pFDphi_pFDDeltaP[j]->SetLineColor(2);
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

  for(int i=0; i<hist_list_2.size(); i++){
    hist_list_2[i]->Sumw2();
    //hist_list_2[i]->SetTitleFont(0.12);
    hist_list_2[i]->GetXaxis()->CenterTitle();
    hist_list_2[i]->GetXaxis()->SetTitleSize(0.12);
    hist_list_2[i]->GetXaxis()->SetLabelSize(0.06);
    hist_list_2[i]->GetXaxis()->SetTitleOffset(0.6);
    hist_list_2[i]->GetYaxis()->CenterTitle();
    hist_list_2[i]->GetYaxis()->SetTitleSize(0.12);
    hist_list_2[i]->GetYaxis()->SetLabelSize(0.06);
    hist_list_2[i]->GetYaxis()->SetTitleOffset(0.6);
  }

  
  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////
  

  TStyle *myStyle  = new TStyle("MyStyle","My Root Styles");
  myStyle->SetPalette("kbird",0);
  myStyle->SetTitleSize(0.10, "t");
  myStyle->SetOptStat(0);
  myStyle->cd();

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
  
  getABC(h_pFDmom_pFDDeltaP_int_thetaGroup,g_pFDmom_pFDDeltaP_int_thetaGroup,f_pFDmom_pFDDeltaP_int_thetaGroup,p_pFDmom_pFDDeltaP_int_thetaGroup,0.3,5.5,g_pFDmom_pFDDeltaP_int_Pars,f_pFDmom_pFDDeltaP_int_Pars,p_pFDmom_pFDDeltaP_int_Pars,f_pFDmom_pFDDeltaP_int_combined_thetaGroup,myCanvas,fileName);

  for(int k = 0; k < 6; k++){
    getABC(h_pFDmom_pFDDeltaP_thetaGroup[k],g_pFDmom_pFDDeltaP_thetaGroup[k],f_pFDmom_pFDDeltaP_thetaGroup[k],p_pFDmom_pFDDeltaP_thetaGroup[k],0.3,5.5,g_pFDmom_pFDDeltaP_Pars[k],f_pFDmom_pFDDeltaP_Pars[k],p_pFDmom_pFDDeltaP_Pars[k],f_pFDmom_pFDDeltaP_combined_thetaGroup[k],myCanvas,fileName);
  }
  
  myCanvas->Divide(2,2,0,0);
  for(int j=0; j<3; j++){
    auto pad = myCanvas->cd(j+1);    
    pad->SetBottomMargin(0.19);
    pad->SetLeftMargin(0.19);
    g_pFDmom_pFDDeltaP_int_Pars[j]->Draw();
    for(int k=0; k<6; k++){
      g_pFDmom_pFDDeltaP_Pars[k][j]->SetLineColor(k+2);
      g_pFDmom_pFDDeltaP_Pars[k][j]->Draw("SAME");
    }
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  myCanvas->Divide(4,4,0,0);  
  for(int i=0; i<16; i++){
    auto pad = myCanvas->cd(i+1);    
    pad->SetBottomMargin(0.19);
    pad->SetLeftMargin(0.19);
    f_pFDmom_pFDDeltaP_int_combined_thetaGroup[i]->SetLineWidth(3);
    f_pFDmom_pFDDeltaP_int_combined_thetaGroup[i]->Draw();
    for(int k=0; k<6; k++){
      f_pFDmom_pFDDeltaP_combined_thetaGroup[k][i]->SetLineWidth(1);
      f_pFDmom_pFDDeltaP_combined_thetaGroup[k][i]->SetLineColor(k+2);
      f_pFDmom_pFDDeltaP_combined_thetaGroup[k][i]->Draw("SAME");
    }
  }  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  for(int sector = 0; sector<6; sector++){
    myCanvas->Divide(4,4,0,0);
    for(int i=1; i<15; i++){
      auto pad = myCanvas->cd(i);    
      pad->SetBottomMargin(0.19);
      pad->SetLeftMargin(0.19);
      h_pFDmom_pFDDeltaP_corrected_thetaGroup[sector][i]->Draw("colz");
    }
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear(); 
  }

  
  myCanvas->Divide(4,4,0,0);
  for(int i=1; i<15; i++){
    auto pad = myCanvas->cd(i);    
    pad->SetBottomMargin(0.19);
    pad->SetLeftMargin(0.19);
    for(int sector = 0; sector<6; sector++){
      h_pFDmom_pFDDeltaP_int_corrected_thetaGroup[i]->Add(h_pFDmom_pFDDeltaP_corrected_thetaGroup[sector][i],1);
    }
    h_pFDmom_pFDDeltaP_int_corrected_thetaGroup[i]->Draw("colz");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear(); 

  myCanvas->Divide(4,4,0,0);
  for(int i = 1; i < 15; i++){
    double min = (double)i*2.5 + 3.0;
    double max = (double)i*2.5 + 5.5;
    sprintf(temp_name,"h_pFDDeltaP_int_thetaGroup_%d",i+1);
    sprintf(temp_title,"(%.1f #circ< #theta < %.1f #circ); #Delta p_{E.L.} [GeV]; Counts",min,max);
    TH1D * h_proj = h_pFDmom_pFDDeltaP_int_thetaGroup[i]->ProjectionY(temp_name,0,h_pFDmom_pFDDeltaP_int_thetaGroup[i]->GetNbinsX());
    h_proj->SetTitle(temp_title);
    h_proj->Rebin(2);
    sprintf(temp_name,"h_pFDDeltaP_int_corrected_thetaGroup_%d",i+1);
    sprintf(temp_title,"(%.1f #circ< #theta < %.1f #circ); #Delta p_{E.L.} [GeV]; Counts",min,max);
    TH1D * h_proj_corr = h_pFDmom_pFDDeltaP_int_corrected_thetaGroup[i]->ProjectionY(temp_name,0,h_pFDmom_pFDDeltaP_int_corrected_thetaGroup[i]->GetNbinsX());
    h_proj_corr->SetTitle(temp_title);

    myCanvas->cd(i);
    myCanvas->GetPad(i)->SetBottomMargin(0.19);
    myCanvas->GetPad(i)->SetTopMargin(0.19);
    myCanvas->GetPad(i)->SetLeftMargin(0.19);
    h_proj->GetXaxis()->CenterTitle();
    h_proj->GetXaxis()->SetTitleSize(0.12);
    h_proj->GetXaxis()->SetLabelSize(0.06);
    h_proj->GetXaxis()->SetTitleOffset(0.6);
    h_proj->GetYaxis()->CenterTitle();
    h_proj->GetYaxis()->SetTitleSize(0.12);
    h_proj->GetYaxis()->SetLabelSize(0.06);
    h_proj->GetYaxis()->SetTitleOffset(0.6);

    h_proj_corr->SetLineColor(2);
    h_proj_corr->Draw();
    h_proj->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear(); 

  double means_list[14];
  double sigmas_list[14];
  double means_list_corr[14];
  double sigmas_list_corr[14];
  myCanvas->Divide(4,4,0,0);
  for(int i = 1; i < 15; i++){
    double min = (double)i*2.5 + 3.0;
    double max = (double)i*2.5 + 5.5;
    sprintf(temp_name,"h_pFDDeltaP_int_thetaGroup_%d",i+1);
    sprintf(temp_title,"(%.1f #circ< #theta < %.1f #circ); #Delta p_{E.L.} [GeV]; Counts",min,max);
    TH1D * h_proj = h_pFDmom_pFDDeltaP_int_thetaGroup[i]->ProjectionY(temp_name,0,h_pFDmom_pFDDeltaP_int_thetaGroup[i]->GetNbinsX());
    h_proj->SetTitle(temp_title);
    h_proj->Rebin(2);
    sprintf(temp_name,"h_pFDDeltaP_int_corrected_thetaGroup_%d",i+1);
    sprintf(temp_title,"(%.1f #circ< #theta < %.1f #circ); #Delta p_{E.L.} [GeV]; Counts",min,max);
    TH1D * h_proj_corr = h_pFDmom_pFDDeltaP_int_corrected_thetaGroup[i]->ProjectionY(temp_name,0,h_pFDmom_pFDDeltaP_int_corrected_thetaGroup[i]->GetNbinsX());
    h_proj_corr->SetTitle(temp_title);

    TF1 * gFit = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.05,0.05,3);
    double maxbin=h_proj->GetMaximumBin();
    double mode = (i<14)?h_proj->GetBinCenter(maxbin):0.02;
    gFit->SetParameter(0,h_proj->GetMaximum()/G(0,1,0,0.1));
    gFit->SetParameter(1,mode);
    gFit->SetParLimits(1,mode-0.02,mode+0.02);
    gFit->SetParameter(2,0.01);
    gFit->SetParLimits(2,0.001,0.05);
    TFitResultPtr gPoint = h_proj->Fit(gFit,"SrBeqn","",-0.05,0.05);
    if(gPoint==0){
      means_list[i-1]=gPoint->Parameter(1);
      sigmas_list[i-1]=gPoint->Parameter(2);
    }


    TF1 * gFit_corr = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.05,0.05,3);
    double maxbin_corr=h_proj_corr->GetMaximumBin();
    double mode_corr = h_proj_corr->GetBinCenter(maxbin_corr);
    gFit_corr->SetParameter(0,h_proj_corr->GetMaximum()/G(0,1,0,0.1));
    gFit_corr->SetParameter(1,mode_corr);
    gFit_corr->SetParLimits(1,mode_corr-0.02,mode_corr+0.02);
    gFit_corr->SetParameter(2,0.01);
    gFit_corr->SetParLimits(2,0.001,0.05);
    TFitResultPtr gPoint_corr = h_proj_corr->Fit(gFit_corr,"SrBeqn","",-0.05,0.02);
    if(gPoint_corr==0){
      means_list_corr[i-1]=gPoint_corr->Parameter(1);
      sigmas_list_corr[i-1]=gPoint_corr->Parameter(2);
    }



    myCanvas->cd(i);
    myCanvas->GetPad(i)->SetBottomMargin(0.19);
    myCanvas->GetPad(i)->SetTopMargin(0.19);
    myCanvas->GetPad(i)->SetLeftMargin(0.19);
    h_proj->GetXaxis()->CenterTitle();
    h_proj->GetXaxis()->SetTitleSize(0.12);
    h_proj->GetXaxis()->SetLabelSize(0.06);
    h_proj->GetXaxis()->SetTitleOffset(0.6);
    h_proj->GetYaxis()->CenterTitle();
    h_proj->GetYaxis()->SetTitleSize(0.12);
    h_proj->GetYaxis()->SetLabelSize(0.06);
    h_proj->GetYaxis()->SetTitleOffset(0.6);

    h_proj_corr->SetLineColor(2);
    h_proj_corr->Draw();
    h_proj->Draw("SAME");

    gFit->SetLineColor(4);
    gFit->Draw("SAME");

    gFit_corr->SetLineColor(2);
    gFit_corr->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear(); 


  myCanvas->Divide(2,3,0,0);
  for(int i = 0; i < 6; i++){
    auto pad = myCanvas->cd(i+1);    
    pad->SetBottomMargin(0.19);
    pad->SetLeftMargin(0.19);
    h_pFDphi_pFDDeltaP[i]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear(); 

  myCanvas->Divide(2,3,0,0);
  for(int i = 0; i < 6; i++){
    auto pad = myCanvas->cd(i+1);    
    pad->SetBottomMargin(0.19);
    pad->SetLeftMargin(0.19);
    h_pFDtheta_pFDDeltaP[i]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear(); 

  
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  cout<<"params_EnergyLoss_FD= \n {{"
      << p_pFDmom_pFDDeltaP_int_Pars[0]->Parameter(0)<<","
      << p_pFDmom_pFDDeltaP_int_Pars[0]->Parameter(1)<<","
      <<"0},\n"
      <<"{"
      <<p_pFDmom_pFDDeltaP_int_Pars[1]->Parameter(0)<<","
      <<p_pFDmom_pFDDeltaP_int_Pars[1]->Parameter(1)<<","
      <<"0},\n"
      <<"{"
      <<p_pFDmom_pFDDeltaP_int_Pars[2]->Parameter(0)<<","
      <<p_pFDmom_pFDDeltaP_int_Pars[2]->Parameter(1)<<","
      <<p_pFDmom_pFDDeltaP_int_Pars[2]->Parameter(2)
      <<"}}"<<endl;

  for(int i = 0; i < 14; i++){
    cout<<means_list[i]*1000<<"&"<<means_list_corr[i]*1000<<"&"<<sigmas_list[i]*1000<<"&"<<sigmas_list_corr[i]*1000<<"\n";
  }
  
  return 0;
}

/*
	  TLorentzVector balance_ptr = beam + target_ptr - proton_ptr_p1 - proton_ptr_p2;
	  double theta_be = balance_ptr.Vect().Angle(el.Vect());
	  double Eb = balance_ptr.E();
	  double Pb = balance_ptr.P();
	  double K = (pim_ptr.M2() - balance_ptr.M2() - el.M2()) / 2;
	  double a = Pb*Pb*cos(theta_be)*cos(theta_be) - Eb*Eb;
	  double b = -2 * K * Pb * cos(theta_be);
	  double c = K*K - Eb*Eb*el.M2();
	  double x_min = (-b - sqrt(b*b - 4*a*c))/(2*a);
o	  double x_max = (-b + sqrt(b*b - 4*a*c))/(2*a);
*/
