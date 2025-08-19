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
#include "reweighter.h"

using namespace std;
using namespace clas12;

const double c = 29.9792458;

vector<double> bE_MomCD = {0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,3.0};//{0.5,1.0,1.3,1.6,2.0,2.5,3.0};
vector<double> bE_Theta = {8,10,12,14,16,18,20,23,26,29,32,35,38,41,45};
vector<double> bE_ThetaCD = {38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89,92};//{40,45,50,55,60,70,90};
vector<double> bE_Phi = {-35,-15,-5,0,5,10,15,25,35};
int binX(vector<double> XS, double X){
  for(int i = 0; i <= XS.size(); i++){if(X<XS[i]){return i-1;}}
  return -1;
}

auto db=TDatabasePDG::Instance();
double mass_p = db->GetParticle(2212)->Mass();
double mD = 1.8756;

//double beam_E = 5.98636;
double beam_E = 5.984792;
double beam_E_sigma = 0.00299;
const double me = 0.000511;
const double mU = 0.9314941024;
const double m_4He = 4.00260325415 * mU - 2*me;


double G(double x, double N, double mu, double sigma){
  return (N/(sigma*sqrt(2*M_PI))) * exp(-0.5 * sq((x-mu)/sigma)) ; 
}

double FuncPhiDependenceFD(double x, double A, double B, double C, double D){
  return A + B*sin((x*2*M_PI/C)+D); 
}

double FuncPhiDependenceCD(double x, double A, double B, double C){
  return A + B*sin((x*2*M_PI/180)+C);
}

double FuncPhiDependenceCD2(double x, double A, double B, double C, double D, double E){
  return A + B*sin((x*1*M_PI/180)+C) + D*sin((x*2*M_PI/180)+E); 
}

double FuncThetaDependenceCD(double x, double A, double B){
  return A + B/(x); 
}

double FuncMomDependenceCD(double x, double A, double B){
  return A + B*x; 
}

double FuncThetaDependenceCD2(double x, double A){
  return A ; 
}

double ErfFunc(double x, double A, double B, double C, double D){
  return A - B*(1+erf(((-x+D)/C))); 
}


double cotan(double x){ return cos(x)/sin(x);}

double ThetaE(double ThetaP){
  double ThetaP_rad = ThetaP*M_PI/180;
  return (180/M_PI) * 2*atan( (1/((beam_E/mass_p)+1)) * cotan(ThetaP_rad) );
}


double ThetaP(double ThetaE){
  double ThetaE_rad = ThetaE*M_PI/180;
  return (180/M_PI) * atan( (1/((beam_E/mass_p)+1)) * cotan(ThetaE_rad/2) );
}

double SQ(double x){ return x*x;}

double ClosestPoint(double x, double ThetaEp, double ThetaPp){
  return sqrt(SQ(x-ThetaEp) + SQ(ThetaP(x)-ThetaPp));
}

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());
}

void getGraph(TH2D * h_myhist, TGraphErrors * g_mygraph, TCanvas * myCanvas, char * fileName, int & ctr){
  bool islarge=true;
  h_myhist->RebinX(2);
  h_myhist->RebinY(2);
  if(h_myhist->GetEntries()<2000){
    islarge=false;
    h_myhist->RebinX(2);
  }
  char temp[100];
  //Now project the histogram    
  for(int j = 0; j < h_myhist->GetXaxis()->GetNbins(); j++){
    //Define x and y(1D histogram)
    double x = h_myhist->GetXaxis()->GetBinCenter(j+1);
    ctr++;
    sprintf(temp,"Proj_num%d",ctr);
    TH1D * proj = h_myhist->ProjectionY(temp,j+1,j+1);
    double mode = proj->GetBinCenter(proj->GetMaximumBin());
    double mean = proj->GetMean();
    double cent = islarge?mode:mean;
    double stddev = proj->GetStdDev();
    if(stddev>0.05){
      proj->Rebin(2);
    }
    //Now preform a guassian fit
    if(proj->GetEntries()<20){continue;}
    //proj->Smooth(1);
    TF1 * gFit = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },cent-stddev,cent+stddev,3);
    gFit->SetParameter(0,proj->GetMaximum()/G(0,1,0,0.005));
    gFit->SetParLimits(0,proj->GetMaximum()/G(0,1,0,0.005)*0.5,proj->GetMaximum()/G(0,1,0,0.05)*1.5);
    gFit->SetParameter(1,cent);
    gFit->SetParLimits(1,cent-0.5*stddev,cent+0.5*stddev);
    gFit->SetParameter(2,stddev);
    gFit->SetParLimits(2,0.001,0.2);
    
    TFitResultPtr gPoint = proj->Fit(gFit,"SrBeqn","",cent-2*stddev,cent+stddev);
    if(gPoint == 0){
      g_mygraph->SetPoint(g_mygraph->GetN(),x,gPoint->Parameter(1));
      //g_mygraph->SetPoint(g_mygraph->GetN(),x,mode);
      g_mygraph->SetPointError(g_mygraph->GetN()-1,0,gPoint->Parameter(2));
    }
    /*
    myCanvas->Divide(1,1);
    proj->Draw();
    if(gPoint==0){
      gFit->Draw("SAME");
    }
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();
    */
  }
}

void getGraphP(TH2D * h_myhist, TGraphErrors * g_mygraph, TCanvas * myCanvas, char * fileName, int & ctr){
  bool islarge=true;
  h_myhist->RebinX(5);
  h_myhist->RebinY(5);
  char temp[100];
  //Now project the histogram    
  for(int j = 0; j < h_myhist->GetXaxis()->GetNbins(); j++){
    //Define x and y(1D histogram)
    double x = h_myhist->GetXaxis()->GetBinCenter(j+1);
    ctr++;
    sprintf(temp,"Proj_num%d",ctr);
    TH1D * proj = h_myhist->ProjectionY(temp,j+1,j+1);
    double mode = proj->GetBinCenter(proj->GetMaximumBin());
    double mean = proj->GetMean();
    double cent = mean;
    double stddev = proj->GetStdDev();
    //Now preform a guassian fit
    if(proj->GetEntries()<20){continue;}
    //proj->Smooth(1);
    TF1 * gFit = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },cent-2*stddev,cent+2*stddev,3);
    gFit->SetParameter(0,proj->GetMaximum()/G(0,1,0,0.005));
    //gFit->SetParLimits(0,proj->GetMaximum()/G(0,1,0,0.005)*0.5,proj->GetMaximum()/G(0,1,0,0.05)*1.5);
    gFit->SetParameter(1,cent);
    gFit->SetParLimits(1,cent-0.5*stddev,cent+0.5*stddev);
    gFit->SetParameter(2,stddev);
    gFit->SetParLimits(2,0.001,0.2);
    
    TFitResultPtr gPoint = proj->Fit(gFit,"SrBeqn","",cent-2*stddev,cent+2*stddev);
    if(gPoint == 0){
      g_mygraph->SetPoint(g_mygraph->GetN(),x,gPoint->Parameter(1));
      //g_mygraph->SetPoint(g_mygraph->GetN(),x,mode);
      g_mygraph->SetPointError(g_mygraph->GetN()-1,0,gPoint->Parameter(2));
    }
    /*
    myCanvas->Divide(1,1);
    proj->Draw();
    if(gPoint==0){
      gFit->Draw("SAME");
    }
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();
    */
  }
}

void getGraphCD(TH2D * h_myhist, TGraphErrors * g_mygraph, TCanvas * myCanvas, char * fileName, int & ctr){
  //h_myhist->RebinY(2);
  char temp[100];
  //Now project the histogram    
  for(int j = 0; j < h_myhist->GetXaxis()->GetNbins(); j++){
    //Define x and y(1D histogram)
    double x = h_myhist->GetXaxis()->GetBinCenter(j+1);
    ctr++;
    sprintf(temp,"Proj_num%d",ctr);
    TH1D * proj = h_myhist->ProjectionY(temp,j+1,j+1);
    double mode = proj->GetBinCenter(proj->GetMaximumBin());
    double mean = proj->GetMean();
    double cent = mean;
    double stddev = proj->GetStdDev();
    if(stddev>0.05){
      proj->Rebin(2);
    }
    //Now preform a guassian fit
    if(proj->GetEntries()<20){continue;}
    //proj->Smooth(1);
    TF1 * gFit = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.5,0.5,3);
    gFit->SetParameter(0,proj->GetMaximum()/G(0,1,0,0.1));
    //gFit->SetParLimits(0,proj->GetMaximum()/G(0,1,0,0.1)*0.5,proj->GetMaximum()/G(0,1,0,0.05)*1.5);
    gFit->SetParameter(1,cent);
    gFit->SetParLimits(1,-0.4,0.4);
    gFit->SetParameter(2,stddev);
    gFit->SetParLimits(2,0.001,0.5);
    
    TFitResultPtr gPoint = proj->Fit(gFit,"SrBeqn","",-0.5,0.5);
    if(gPoint == 0){
      g_mygraph->SetPoint(g_mygraph->GetN(),x,gPoint->Parameter(1));
      //g_mygraph->SetPoint(g_mygraph->GetN(),x,mode);
      g_mygraph->SetPointError(g_mygraph->GetN()-1,0,gPoint->Parameter(2));
    }
    /*
    myCanvas->Divide(1,1);
    proj->Draw();
    if(gPoint==0){
      gFit->Draw("SAME");
    }
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();
    */
  }
}

void getFunction(TGraphErrors * g_mygraph, TF1 * f_myfunc, TFitResultPtr & p_mypoint, double guess[4],int i, int sector){
  double min = g_mygraph->GetMinimum();
  double max = g_mygraph->GetMaximum();
  f_myfunc->SetLineColor(3);
  f_myfunc->SetLineWidth(1);

  
  //f_myfunc->SetParLimits(0,-0.2,+0.2);  
  double min0 = guess[0]-0.04;
  double max0 = guess[0]+0.04;
  double min1 = ((guess[1] * 0.5)>0.005)?(guess[1] * 0.5):0.005;
  double max1 = guess[1] * 1.05;
  double min2 = ((guess[2]*0.6)>10)?guess[2]*0.6:10;
  double max2 = ((guess[2]*1.1)<150)?(guess[2]*1.1):150;
  double min3 = guess[3]+0;
  double max3 = guess[3]+M_PI;
  if(i>6){
    min2 = ((guess[2]*0.9)>10)?guess[2]*0.9:10;
  }
  
  if(i==11){
    min1 = 0.05;
    max1 = 0.2;
    min2 = 85;
    max2 = 150;
    min3 = -M_PI;
    max3 = +M_PI;
  }

  if((sector==3)||(sector==4)||(sector==6)){
    min3 = -M_PI;
    max3 = +M_PI;
  }


  if((i>7)){
    if((sector==4)||(sector==5)||(sector==6))
    min2 = 100;
    max2 = 130;
  }
  
  if((i<3)&&(sector==6)){    
    min3 = guess[3]-0.2;
    //max3 = guess[3]+0.5;
  }
  if((i>6)&&(sector==6)){    
    max3 = guess[3]+0.5;
  }

  if(sector==2){
    //max2 = 90;
    min2 = 40;
    min1=0.02;
    if(i<6){
    max1=0.022;
    }
  }

  if(sector==2){
    min2 = 40;
    min1=0.02;
    if(i<6){
    max1=0.022;
    }
  }

  f_myfunc->SetParameter(0,guess[0]);
  f_myfunc->SetParLimits(0,min0,max0);
  f_myfunc->SetParameter(1,guess[1]);
  f_myfunc->SetParLimits(1,min1,max1);  
  f_myfunc->SetParameter(2,guess[2]);
  f_myfunc->SetParLimits(2,min2,max2);
  f_myfunc->SetParameter(3,guess[3]);
  f_myfunc->SetParLimits(3,min3,max3);


  p_mypoint = g_mygraph->Fit(f_myfunc,"SrBeqn","",-40,40);
      
}

void getFunctionP(TGraphErrors * g_mygraph, TF1 * f_myfunc, TFitResultPtr & p_mypoint, double guess[4],int i, int sector){
  f_myfunc->SetLineColor(3);
  f_myfunc->SetLineWidth(1);

  
  //f_myfunc->SetParLimits(0,-0.2,+0.2);  
  double min0 = guess[0]-0.1;
  double max0 = guess[0]+0.1;
  double min1 = 0.005;
  double max1 = guess[1] * 1.5;
  double min2 = 100;
  double max2 = 180;
  double min3 = -M_PI;
  double max3 = M_PI;
  if(i<12){
    min3=guess[3]-0.1;
    max3=guess[3]+0.3;
    min2=guess[2]-40;
    max2=guess[2]+10;
  }
  if(min2<40){
    min2=40;
  }
  
  f_myfunc->SetParameter(0,guess[0]);
  f_myfunc->SetParameter(1,guess[1]);
  f_myfunc->SetParameter(2,guess[2]);
  f_myfunc->SetParameter(3,guess[3]);

  f_myfunc->SetParLimits(2,min2,max2);
  f_myfunc->SetParLimits(3,min3,max3);
  /*
  f_myfunc->SetParLimits(0,min0,max0);
  f_myfunc->SetParLimits(1,min1,max1);  
  */

  p_mypoint = g_mygraph->Fit(f_myfunc,"SrBeqn","",-50,40);
      
}


void getFunctionTH2D(TH2D * myhist, TF1 * f_myfunc, TFitResultPtr & p_mypoint, double guess[4],int i, int sector){
  f_myfunc->SetLineColor(3);
  f_myfunc->SetLineWidth(1);

  
  //f_myfunc->SetParLimits(0,-0.2,+0.2);  
  double min0 = guess[0]-0.04;
  double max0 = guess[0]+0.04;
  double min1 = ((guess[1] * 0.5)>0.005)?(guess[1] * 0.5):0.005;
  double max1 = guess[1] * 1.05;
  double min2 = ((guess[2]*0.6)>10)?guess[2]*0.6:10;
  double max2 = ((guess[2]*1.1)<150)?(guess[2]*1.1):150;
  double min3 = guess[3]+0;
  double max3 = guess[3]+M_PI;
  if(i>6){
    min2 = ((guess[2]*0.9)>10)?guess[2]*0.9:10;
  }
  
  if(i==11){
    min1 = 0.05;
    max1 = 0.2;
    min2 = 85;
    max2 = 150;
    min3 = -M_PI;
    max3 = +M_PI;
  }

  if((sector==3)||(sector==4)||(sector==6)){
    min3 = -M_PI;
    max3 = +M_PI;
  }


  if((i>7)){
    if((sector==4)||(sector==5)||(sector==6))
    min2 = 100;
    max2 = 130;
  }
  
  if((i<3)&&(sector==6)){    
    min3 = guess[3]-0.2;
    //max3 = guess[3]+0.5;
  }
  if((i>6)&&(sector==6)){    
    max3 = guess[3]+0.5;
  }
  f_myfunc->SetParameter(0,guess[0]);
  f_myfunc->SetParLimits(0,min0,max0);
  f_myfunc->SetParameter(1,guess[1]);
  f_myfunc->SetParLimits(1,min1,max1);  
  f_myfunc->SetParameter(2,guess[2]);
  f_myfunc->SetParLimits(2,min2,max2);
  f_myfunc->SetParameter(3,guess[3]);
  f_myfunc->SetParLimits(3,min3,max3);


  p_mypoint = myhist->Fit(f_myfunc,"SrBeqn","",-40,40);
      
}

void getLimits(TF1 * f_params, TFitResultPtr p_params, double x, double vals[2]){
  TF1 * f_p = new TF1("params",[&](double *x, double *p){ return ErfFunc(x[0],p[0],p[1],p[2],p[3]);},5,50,4);

  double center = f_params->Eval(x);
  double p0 = p_params->Parameter(0);
  double p1 = p_params->Parameter(1);
  double p2 = p_params->Parameter(2);
  double p3 = p_params->Parameter(3);
  double min = center*0.9;
  double max = center*1.1;
  for(int i = -1; i <=1; i+=2){
  for(int j = -1; j <=1; j+=2){
  for(int k = -1; k <=1; k+=2){
  for(int l = -1; l <=1; l+=2){
    double p0_p = p0 + ((double)i * 0.05 * p0);
    double p1_p = p1 + ((double)j * 0.1 * p1);
    double p2_p = p2 + ((double)k * 0.1 * p2);
    double p3_p = p3 + ((double)l * 0.1 * p3);
    f_p->SetParameter(0,p0_p);
    f_p->SetParameter(1,p1_p);
    f_p->SetParameter(2,p2_p);
    f_p->SetParameter(3,p3_p);
    double eval = f_p->Eval(x);
    //cout<<center<<" "<<eval<<endl;
    if(eval<min){min=eval;}
    if(eval>max){max=eval;}
  }
  }
  }
  }
  vals[0] = min;
  vals[1] = max;
}

void getFunctionRedo(TGraphErrors * g_mygraph, TF1 * f_myfunc, TFitResultPtr & p_mypoint, TF1 * f_params[4], TFitResultPtr p_params[4], TGraph * params_mins_Redo[4], TGraph * params_maxs_Redo[4], double x, int sector){
  f_myfunc->SetLineColor(6);
  f_myfunc->SetLineWidth(1);

  for(int i = 0; i < 4; i++){
    f_myfunc->SetParameter(i,f_params[i]->Eval(x));
    double lims[2];
    getLimits(f_params[i],p_params[i],x,lims);
    params_mins_Redo[i]->SetPoint(params_mins_Redo[i]->GetN(),x,lims[0]);
    params_maxs_Redo[i]->SetPoint(params_maxs_Redo[i]->GetN(),x,lims[1]);
    f_myfunc->SetParLimits(i,lims[0],lims[1]);
  }
  p_mypoint = g_mygraph->Fit(f_myfunc,"SrBeqn","",-40,40);
      
}

void getFunctionCD(TGraphErrors * g_mygraph, TF1 * f_myfunc, TFitResultPtr & p_mypoint, double guess[5]){
  double min = g_mygraph->GetMinimum();
  double max = g_mygraph->GetMaximum();
  f_myfunc->SetLineColor(3);
  f_myfunc->SetLineWidth(1);
  f_myfunc->SetParameter(0,guess[0]);
  f_myfunc->SetParLimits(0,-0.5,0.5);
  f_myfunc->SetParameter(1,0.05);
  f_myfunc->SetParLimits(1,0.001,0.5);
  f_myfunc->SetParameter(2,guess[2]);
  f_myfunc->SetParLimits(2,-M_PI,M_PI);
  p_mypoint = g_mygraph->Fit(f_myfunc,"SrBeqn","",-180,180);
      
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
  TLorentzVector nucleus_ptr(0,0,0,m_4He);
  TLorentzVector deut_ptr(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector proton_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector lead_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector lead_pion_ptr(0,0,0,db->GetParticle(211)->Mass());
  TLorentzVector stationary_proton_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector ntr(0,0,0,db->GetParticle(2112)->Mass());
  reweighter newWeight(beam_E,2,2,kelly,"AV18");

  //////////////////////////////////
  char temp_name[100];
  char temp_title[100];

  vector<TH1*> hist_list;
 
  TH2D * h_e_phi_corr_binSector_binTheta[6][14];
  TGraphErrors * g_e_phi_corr_binSector_binTheta[6][14];
  TF1 * f_e_phi_corr_binSector_binTheta[6][14];
  TFitResultPtr p_e_phi_corr_binSector_binTheta[6][14];
  TF1 * f_e_phi_corr_Combined_binSector_binTheta[6][14];
  TF1 * f_e_phi_corr_Redo_binSector_binTheta[6][14];
  TFitResultPtr p_e_phi_corr_Redo_binSector_binTheta[6][14];
  TF1 * f_e_phi_corr_Redo_Combined_binSector_binTheta[6][14];
  for(int j=1; j<=6; j++){
    for(int i=0; i<14; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      sprintf(temp_name,"e_phi_corr_sector_%d_theta_%d",j,i);
      sprintf(temp_title,"Correction vs. #phi Sector %d (%d< #theta < %d);#phi;Correction [GeV];Counts",j,min,max);
      h_e_phi_corr_binSector_binTheta[j-1][i] = (TH2D*)inFile->Get(temp_name);
      hist_list.push_back(h_e_phi_corr_binSector_binTheta[j-1][i]);
      
      g_e_phi_corr_binSector_binTheta[j-1][i] = new TGraphErrors();
      sprintf(temp_name,"g_e_phi_corr_sector_%d_theta_%d",j,i);
      g_e_phi_corr_binSector_binTheta[j-1][i]->SetName(temp_name);

      sprintf(temp_name,"f_e_phi_corr_sector_%d_theta_%d",j,i);
      f_e_phi_corr_binSector_binTheta[j-1][i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceFD(x[0],p[0],p[1],p[2],p[3]); },-40,40,4);
      
      sprintf(temp_name,"f_e_phi_corr_Combined_sector_%d_theta_%d",j,i);
      f_e_phi_corr_Combined_binSector_binTheta[j-1][i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceFD(x[0],p[0],p[1],p[2],p[3]); },-40,40,4);

      sprintf(temp_name,"f_e_phi_corr_Redo_sector_%d_theta_%d",j,i);
      f_e_phi_corr_Redo_binSector_binTheta[j-1][i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceFD(x[0],p[0],p[1],p[2],p[3]); },-40,40,4);
      
      sprintf(temp_name,"f_e_phi_corr_Redo_Combined_sector_%d_theta_%d",j,i);
      f_e_phi_corr_Redo_Combined_binSector_binTheta[j-1][i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceFD(x[0],p[0],p[1],p[2],p[3]); },-40,40,4);

    }
  }

  TH2D * h_p_phi_corr_binSector_binTheta[6][14];
  TGraphErrors * g_p_phi_corr_binSector_binTheta[6][14];
  TF1 * f_p_phi_corr_binSector_binTheta[6][14];
  TFitResultPtr p_p_phi_corr_binSector_binTheta[6][14];
  TF1 * f_p_phi_corr_Combined_binSector_binTheta[6][14];
  TF1 * f_p_phi_corr_Redo_binSector_binTheta[6][14];
  TFitResultPtr p_p_phi_corr_Redo_binSector_binTheta[6][14];
  TF1 * f_p_phi_corr_Redo_Combined_binSector_binTheta[6][14];
  for(int j=1; j<=6; j++){
    for(int i=0; i<14; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      sprintf(temp_name,"p_phi_corr_sector_%d_theta_%d",j,i);
      sprintf(temp_title,"Correction vs. #phi Sector %d (%d< #theta < %d);#phi;Correction [GeV];Counts",j,min,max);
      h_p_phi_corr_binSector_binTheta[j-1][i] = (TH2D*)inFile->Get(temp_name);
      hist_list.push_back(h_p_phi_corr_binSector_binTheta[j-1][i]);

      g_p_phi_corr_binSector_binTheta[j-1][i] = new TGraphErrors();
      sprintf(temp_name,"g_p_phi_corr_sector_%d_theta_%d",j,i);
      g_p_phi_corr_binSector_binTheta[j-1][i]->SetName(temp_name);

      sprintf(temp_name,"f_p_phi_corr_sector_%d_theta_%d",j,i);
      f_p_phi_corr_binSector_binTheta[j-1][i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceFD(x[0],p[0],p[1],p[2],p[3]); },-40,40,4);
      
      sprintf(temp_name,"f_p_phi_corr_Combined_sector_%d_theta_%d",j,i);
      f_p_phi_corr_Combined_binSector_binTheta[j-1][i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceFD(x[0],p[0],p[1],p[2],p[3]); },-40,40,4);

      sprintf(temp_name,"f_p_phi_corr_Redo_sector_%d_theta_%d",j,i);
      f_p_phi_corr_Redo_binSector_binTheta[j-1][i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceFD(x[0],p[0],p[1],p[2],p[3]); },-40,40,4);
      
      sprintf(temp_name,"f_p_phi_corr_Redo_Combined_sector_%d_theta_%d",j,i);
      f_p_phi_corr_Redo_Combined_binSector_binTheta[j-1][i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceFD(x[0],p[0],p[1],p[2],p[3]); },-40,40,4);

    }
  }

  
  TH2D * h_phi_corr_binThetaCD[18];
  TGraphErrors * g_phi_corr_binThetaCD[18];
  TF1 * f_phi_corr_binThetaCD[18];
  TFitResultPtr p_phi_corr_binThetaCD[18];
  TF1 * f_phi_corr_Combined_binThetaCD[18];
  for(int i=0; i<18; i++){
    int min = bE_ThetaCD[i];
    int max = bE_ThetaCD[i+1];
    sprintf(temp_name,"phi_corr_theta_%d",i);
    sprintf(temp_title,"Correction vs. #phi (%d< #theta < %d);#phi;Correction [GeV];Counts",min,max);
    h_phi_corr_binThetaCD[i] = (TH2D*)inFile->Get(temp_name);
    hist_list.push_back(h_phi_corr_binThetaCD[i]);
    
    g_phi_corr_binThetaCD[i] = new TGraphErrors();
    sprintf(temp_name,"g_phi_corr_theta_%d",i);
    g_phi_corr_binThetaCD[i]->SetName(temp_name);
    
    sprintf(temp_name,"f_phi_corr_theta_%d",i);
    f_phi_corr_binThetaCD[i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceCD(x[0],p[0],p[1],p[2]); },-180,180,3);

    sprintf(temp_name,"f_phi_corr_Combined_theta_%d",i);
    f_phi_corr_Combined_binThetaCD[i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceCD(x[0],p[0],p[1],p[2]); },-180,180,3);    

    //f_phi_corr_binThetaCD[i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceCD(x[0],p[0],p[1],p[2]); },-180,180,5);    
    //f_phi_corr_Combined_binThetaCD[i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceCD(x[0],p[0],p[1],p[2]); },-180,180,5);
  }
  
  TH1D * h_Delta_Int_Eprime_beforeRad = (TH1D*)inFile->Get("Delta_Eprime_beforeRad");
  //"#Delta p_{e} (e,e'p);#Delta p_{e} [GeV];Counts"
  hist_list.push_back(h_Delta_Int_Eprime_beforeRad);
  TH1D * h_Delta_Int_Eprime_Corrected = (TH1D*)inFile->Get("Delta_Eprime_Corrected");//,"#Delta p_{e} (e,e'p);#Delta p_{e} [GeV];Counts",100,-0.2,0.2);
  hist_list.push_back(h_Delta_Int_Eprime_Corrected);  
  TH1D * h_Delta_Int_pMomFD = (TH1D*)inFile->Get("Delta_pMomFD");//,"#Delta p_{p} (e,e'p_{FD});#Delta p_{p} [GeV]; Counts",100,-0.4,0.4);
  hist_list.push_back(h_Delta_Int_pMomFD);
  TH1D * h_Delta_Int_pMomFD_Corrected = (TH1D*)inFile->Get("Delta_pMomFD_Corrected");//,"#Delta p_{p} (e,e'p_{FD});#Delta p_{p} [GeV]; Counts",100,-0.4,0.4);
  hist_list.push_back(h_Delta_Int_pMomFD_Corrected);
  TH1D * h_Delta_Int_pMomCD = (TH1D*)inFile->Get("Delta_pMomCD");//,"#Delta p_{p} (e,e'p_{CD});#Delta p_{p} [GeV]; Counts",100,-0.4,0.4);
  hist_list.push_back(h_Delta_Int_pMomCD);
  TH1D * h_Delta_Int_pMomCD_Corrected = (TH1D*)inFile->Get("Delta_pMomCD_Corrected");//,"#Delta p_{p} (e,e'p_{CD});#Delta p_{p} [GeV]; Counts",100,-0.4,0.4);
  hist_list.push_back(h_Delta_Int_pMomCD_Corrected);

  
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Sumw2();
    hist_list[i]->GetXaxis()->CenterTitle();
    hist_list[i]->GetXaxis()->SetTitleSize(0.10);
    hist_list[i]->GetXaxis()->SetLabelSize(0.06);
    hist_list[i]->GetXaxis()->SetTitleOffset(0.8);
    hist_list[i]->GetYaxis()->CenterTitle();
    hist_list[i]->GetYaxis()->SetTitleSize(0.10);
    hist_list[i]->GetYaxis()->SetLabelSize(0.06);
    hist_list[i]->GetYaxis()->SetTitleOffset(0.8);
    hist_list[i]->GetXaxis()->CenterTitle();
    hist_list[i]->GetYaxis()->CenterTitle();
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

  TStyle *myStyle  = new TStyle("MyStyle","My Root Styles");
  myStyle->SetPalette("kbird",0);
  myStyle->SetTitleSize(0.07, "t");
  myStyle->SetOptStat(0);
  myStyle->cd();

  char fileName[100];
  sprintf(fileName,"%s[",pdfFile);
  myText->SaveAs(fileName);
  sprintf(fileName,"%s",pdfFile);

  ///////////////////////////////////
  int ctrX=0;
  TFitResultPtr p_e_params[6][4];
  for(int j = 0; j < 6; j++){
    TGraph * e_params[4];
    TGraphErrors * e_params_Errors[4];
    TF1 * f_e_params[4];
    double guess[4];
    for(int k = 0; k < 4; k++){
      e_params[k] = new TGraph();
      e_params_Errors[k] = new TGraphErrors();
      sprintf(temp_name,"e_params_%d",k);
      f_e_params[k] = new TF1(temp_name,[&](double *x, double *p){ return ErfFunc(x[0],p[0],p[1],p[2],p[3]);},5,50,4);
      f_e_params[k]->SetParameter(2,20);
      f_e_params[k]->SetParLimits(2,4,30);
      f_e_params[k]->SetParameter(3,20);
      f_e_params[k]->SetParLimits(3,10,25);
      if((k==0)||(k==1)){
	f_e_params[k]->SetParameter(0,0);
	f_e_params[k]->SetParLimits(0,-10,10);
	f_e_params[k]->SetParameter(1,0);
	f_e_params[k]->SetParLimits(1,-10,10);
      }
      else if(k==2){
	f_e_params[k]->SetParameter(0,0);
	f_e_params[k]->SetParLimits(0,-1000,1000);
	f_e_params[k]->SetParameter(1,100);
	f_e_params[k]->SetParLimits(1,-1000,1000);
      }
      else if(k==3){
	f_e_params[k]->SetParameter(0,0);
	f_e_params[k]->SetParLimits(0,-30,30);
	f_e_params[k]->SetParameter(1,0);
	f_e_params[k]->SetParLimits(1,-10,10);
      }

      if(j==1){
	f_e_params[k]->SetParameter(2,20);
	f_e_params[k]->SetParLimits(2,1.5,30);
	//f_e_params[k]->SetParameter(3,23);
	//f_e_params[k]->SetParLimits(3,22,26);
      }
      if(j==2){
	if(k==0){
	f_e_params[k]->SetParameter(2,20);
	f_e_params[k]->SetParLimits(2,1.5,10);
	//f_e_params[k]->SetParameter(3,23);
	//f_e_params[k]->SetParLimits(3,18,26);
	}
	if(k==2){
	f_e_params[k]->SetParameter(2,20);
	f_e_params[k]->SetParLimits(2,1.5,30);
	f_e_params[k]->SetParameter(3,23);
	f_e_params[k]->SetParLimits(3,18,26);
	}
      }
      if(j==5){
	if(k==0){
	f_e_params[k]->SetParameter(2,20);
	f_e_params[k]->SetParLimits(2,1.5,11);
	//f_e_params[k]->SetParameter(3,23);
	//f_e_params[k]->SetParLimits(3,18,26);
	}
      }

      
    }
    
    
    guess[0] = 0.04;
    guess[1] = 0.1;
    guess[2] = 100;
    guess[3] = 0;
    for(int i = 8; i >= 0; i--){
      getGraph(h_e_phi_corr_binSector_binTheta[j][i],g_e_phi_corr_binSector_binTheta[j][i],myCanvas,fileName,ctrX);      
      getFunction(g_e_phi_corr_binSector_binTheta[j][i],f_e_phi_corr_binSector_binTheta[j][i],p_e_phi_corr_binSector_binTheta[j][i],guess,i,j+1);
      ///////////////      
      double x = (bE_Theta[i]+bE_Theta[i+1])/2;
      for(int k = 0; k < 4; k++){
	e_params[k]->SetPoint(e_params[k]->GetN(),x,p_e_phi_corr_binSector_binTheta[j][i]->Parameter(k));
	e_params_Errors[k]->SetPoint(e_params_Errors[k]->GetN(),x,p_e_phi_corr_binSector_binTheta[j][i]->Parameter(k));
	e_params_Errors[k]->SetPointError(e_params_Errors[k]->GetN()-1,0,p_e_phi_corr_binSector_binTheta[j][i]->ParError(k));
	guess[k]=p_e_phi_corr_binSector_binTheta[j][i]->Parameter(k);
      }
      
    }
    
    for(int k = 0; k < 4; k++){
      p_e_params[j][k] = e_params[k]->Fit(f_e_params[k],"SrBeqn","",5,45);
    }
    for(int i = 0; i < 13; i++){
      double x = (bE_Theta[i]+bE_Theta[i+1])/2;
      f_e_phi_corr_Combined_binSector_binTheta[j][i]->SetLineColor(4);
      for(int k = 0; k < 4; k++){      
	f_e_phi_corr_Combined_binSector_binTheta[j][i]->SetParameter(k,f_e_params[k]->Eval(x));
      }
    }

    myCanvas->Divide(3,3,0,0);
    for(int i = 0; i < 9; i++){
      myCanvas->cd(i+1);
      myCanvas->GetPad(i+1)->SetBottomMargin(0.19);
      myCanvas->GetPad(i+1)->SetLeftMargin(0.19);
      h_e_phi_corr_binSector_binTheta[j][i]->Draw("colz");
      g_e_phi_corr_binSector_binTheta[j][i]->SetLineColor(2);
      g_e_phi_corr_binSector_binTheta[j][i]->Draw("SAME");
      f_e_phi_corr_binSector_binTheta[j][i]->Draw("SAME");
      //f_e_phi_corr_Combined_binSector_binTheta[j][i]->SetLineColor(4);
      //f_e_phi_corr_Combined_binSector_binTheta[j][i]->Draw("SAME");
    }    
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();

    char* names[4] = {"C_{a}","C_{b}","C_{c}","C_{d}"};
    myCanvas->Divide(2,2,0,0);
    for(int k = 0; k < 4; k++){
      myCanvas->cd(k+1);      
      myCanvas->GetPad(k+1)->SetBottomMargin(0.19);
      myCanvas->GetPad(k+1)->SetLeftMargin(0.19);
      e_params[k]->SetLineColor(3);
      e_params_Errors[k]->SetLineColor(3);
      f_e_params[k]->SetLineColor(4);      
      sprintf(temp_title,"%s vs. #theta #circ;#theta #circ;%s",names[k],names[k]);
      e_params_Errors[k]->GetXaxis()->CenterTitle();
      e_params_Errors[k]->GetXaxis()->SetTitleSize(0.10);
      e_params_Errors[k]->GetXaxis()->SetLabelSize(0.06);
      e_params_Errors[k]->GetXaxis()->SetTitleOffset(0.8);
      e_params_Errors[k]->GetYaxis()->CenterTitle();
      e_params_Errors[k]->GetYaxis()->SetTitleSize(0.10);
      e_params_Errors[k]->GetYaxis()->SetLabelSize(0.06);
      e_params_Errors[k]->GetYaxis()->SetTitleOffset(0.8);
      e_params_Errors[k]->GetXaxis()->CenterTitle();
      e_params_Errors[k]->GetYaxis()->CenterTitle();
      e_params_Errors[k]->SetTitle(temp_title);
      //e_params[k]->Draw();
      e_params_Errors[k]->Draw();
      f_e_params[k]->Draw("SAME");      
    }
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();
    
    myCanvas->Divide(3,3,0,0);
    for(int i = 0; i < 9; i++){
      myCanvas->cd(i+1);
      myCanvas->GetPad(i+1)->SetBottomMargin(0.19);
      myCanvas->GetPad(i+1)->SetLeftMargin(0.19);
      h_e_phi_corr_binSector_binTheta[j][i]->Draw("colz");
      g_e_phi_corr_binSector_binTheta[j][i]->Draw("SAME");
      //f_e_phi_corr_binSector_binTheta[j][i]->Draw("SAME");
      f_e_phi_corr_Combined_binSector_binTheta[j][i]->SetLineColor(4);
      f_e_phi_corr_Combined_binSector_binTheta[j][i]->Draw("SAME");
    }    
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();
  }

  
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  TFitResultPtr p_p_params[6][4];
  for(int j = 0; j < 6; j++){
    TGraph * p_params[4];
    TGraphErrors * p_params_Errors[4];
    TF1 * f_p_params[4];
    double guess[4];
    for(int k = 0; k < 4; k++){
      p_params[k] = new TGraph();
      p_params_Errors[k] = new TGraphErrors();
      sprintf(temp_name,"p_params_%d",k);
      /*
      f_p_params[k] = new TF1(temp_name,[&](double *x, double *p){ return FuncMomDependenceCD(x[0],p[0],p[1]);},5,50,2);
      if((k==0)||(k==1)||(k==3)){
	f_p_params[k]->SetParameter(0,0);
	f_p_params[k]->SetParLimits(0,-10,10);
	f_p_params[k]->SetParameter(1,0);
	f_p_params[k]->SetParLimits(1,-10,10);
      }
      else if(k==2){
	f_p_params[k]->SetParameter(0,0);
	f_p_params[k]->SetParLimits(0,-1000,1000);
	f_p_params[k]->SetParameter(1,100);
	f_p_params[k]->SetParLimits(1,-1000,1000);
      }
      */
      f_p_params[k] = new TF1(temp_name,[&](double *x, double *p){ return ErfFunc(x[0],p[0],p[1],p[2],p[3]);},5,50,4);
      f_p_params[k]->SetParameter(2,20);
      f_p_params[k]->SetParLimits(2,4,50);
      f_p_params[k]->SetParameter(3,30);
      f_p_params[k]->SetParLimits(3,26,38);
      if((k==0)||(k==1)){
	f_p_params[k]->SetParameter(0,0);
	f_p_params[k]->SetParLimits(0,-10,10);
	f_p_params[k]->SetParameter(1,0);
	f_p_params[k]->SetParLimits(1,-10,10);
      }
      else if(k==2){
	f_p_params[k]->SetParameter(0,0);
	f_p_params[k]->SetParLimits(0,-1000,1000);
	f_p_params[k]->SetParameter(1,100);
	f_p_params[k]->SetParLimits(1,-1000,1000);
      }
      else if(k==3){
	f_p_params[k]->SetParameter(0,0);
	f_p_params[k]->SetParLimits(0,-30,30);
	f_p_params[k]->SetParameter(1,0);
	f_p_params[k]->SetParLimits(1,-10,10);
      }

      if((j==2)||(j==3)){
	f_p_params[k]->SetParameter(2,20);
	f_p_params[k]->SetParLimits(2,2,50);
      }
      
    }
    
    
    guess[0] = 0.04;
    guess[1] = 0.1;
    guess[2] = 100;
    guess[3] = 0;
    for(int i = 12; i >= 7; i--){
      getGraphP(h_p_phi_corr_binSector_binTheta[j][i],g_p_phi_corr_binSector_binTheta[j][i],myCanvas,fileName,ctrX);      
      getFunctionP(g_p_phi_corr_binSector_binTheta[j][i],f_p_phi_corr_binSector_binTheta[j][i],p_p_phi_corr_binSector_binTheta[j][i],guess,i,j+1);
      ///////////////      
      double x = (bE_Theta[i]+bE_Theta[i+1])/2;
      for(int k = 0; k < 4; k++){
	p_params[k]->SetPoint(p_params[k]->GetN(),x,p_p_phi_corr_binSector_binTheta[j][i]->Parameter(k));
	p_params_Errors[k]->SetPoint(p_params_Errors[k]->GetN(),x,p_p_phi_corr_binSector_binTheta[j][i]->Parameter(k));
	p_params_Errors[k]->SetPointError(p_params_Errors[k]->GetN()-1,0,p_p_phi_corr_binSector_binTheta[j][i]->ParError(k));
	guess[k]=p_p_phi_corr_binSector_binTheta[j][i]->Parameter(k);
      }
      
    }
    
    for(int k = 0; k < 4; k++){
      p_p_params[j][k] = p_params[k]->Fit(f_p_params[k],"SrBeqn","",5,45);
    }
    for(int i = 0; i < 13; i++){
      double x = (bE_Theta[i]+bE_Theta[i+1])/2;
      f_p_phi_corr_Combined_binSector_binTheta[j][i]->SetLineColor(4);
      for(int k = 0; k < 4; k++){      
	f_p_phi_corr_Combined_binSector_binTheta[j][i]->SetParameter(k,f_p_params[k]->Eval(x));
      }
    }
    
    myCanvas->Divide(3,2,0,0);    
    for(int i = 7; i < 13; i++){
      myCanvas->cd(i-6);
      myCanvas->GetPad(i-6)->SetBottomMargin(0.19);
      myCanvas->GetPad(i-6)->SetLeftMargin(0.19);
      h_p_phi_corr_binSector_binTheta[j][i]->Draw("colz");
      g_p_phi_corr_binSector_binTheta[j][i]->SetLineColor(2);
      g_p_phi_corr_binSector_binTheta[j][i]->Draw("SAME");
      f_p_phi_corr_binSector_binTheta[j][i]->Draw("SAME");
      //f_p_phi_corr_Combined_binSector_binTheta[j][i]->SetLineColor(4);
      //f_p_phi_corr_Combined_binSector_binTheta[j][i]->Draw("SAME");
    }    
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();

    char* names[4] = {"C_{a}","C_{b}","C_{c}","C_{d}"};
    myCanvas->Divide(2,2,0,0);    
    for(int k = 0; k < 4; k++){
      myCanvas->cd(k+1);      
      myCanvas->GetPad(k+1)->SetBottomMargin(0.19);
      myCanvas->GetPad(k+1)->SetLeftMargin(0.19);
      p_params[k]->SetLineColor(3);
      p_params_Errors[k]->SetLineColor(3);
      f_p_params[k]->SetLineColor(4);      
      sprintf(temp_title,"%s vs. #theta #circ;#theta #circ;%s",names[k],names[k]);
      p_params_Errors[k]->GetXaxis()->CenterTitle();
      p_params_Errors[k]->GetXaxis()->SetTitleSize(0.10);
      p_params_Errors[k]->GetXaxis()->SetLabelSize(0.06);
      p_params_Errors[k]->GetXaxis()->SetTitleOffset(0.8);
      p_params_Errors[k]->GetYaxis()->CenterTitle();
      p_params_Errors[k]->GetYaxis()->SetTitleSize(0.10);
      p_params_Errors[k]->GetYaxis()->SetLabelSize(0.06);
      p_params_Errors[k]->GetYaxis()->SetTitleOffset(0.8);
      p_params_Errors[k]->GetXaxis()->CenterTitle();
      p_params_Errors[k]->GetYaxis()->CenterTitle();
      p_params_Errors[k]->SetTitle(temp_title);
      //p_params[k]->Draw();
      p_params_Errors[k]->Draw();
      f_p_params[k]->Draw("SAME");
    }
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();

    myCanvas->Divide(3,2,0,0);    
    for(int i = 7; i < 13; i++){
      myCanvas->cd(i-6);
      myCanvas->GetPad(i-6)->SetBottomMargin(0.19);
      myCanvas->GetPad(i-6)->SetLeftMargin(0.19);
      h_p_phi_corr_binSector_binTheta[j][i]->Draw("colz");
      g_p_phi_corr_binSector_binTheta[j][i]->SetLineColor(2);
      g_p_phi_corr_binSector_binTheta[j][i]->Draw("SAME");
      //f_p_phi_corr_binSector_binTheta[j][i]->Draw("SAME");
      f_p_phi_corr_Combined_binSector_binTheta[j][i]->SetLineColor(4);
      f_p_phi_corr_Combined_binSector_binTheta[j][i]->Draw("SAME");
    }    
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();

  }

  
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////

double params_Momentum_CD[3][2]={{-0.107267,0.00193308},
				 {0.0267032,0.000545869},
				 {-2.31836,0.0177897}};
  
  TGraph * params_CD[3];
  TGraphErrors * params_CD_Errors[3];
  for(int k = 0; k < 3; k++){
    params_CD[k] = new TGraph();
    params_CD_Errors[k] = new TGraphErrors();
  }
  double guess_CD[3];  
  guess_CD[0]=0;
  guess_CD[1]=0.05;
  guess_CD[2]=1.4;

  for(int i = 0; i < 18; i++){
    
    if(h_phi_corr_binThetaCD[i]->GetEntries()<1000){
      h_phi_corr_binThetaCD[i]->Rebin(20);
    }
    else if(h_phi_corr_binThetaCD[i]->GetEntries()<100000){
      h_phi_corr_binThetaCD[i]->Rebin(10);
    }
    else if(h_phi_corr_binThetaCD[i]->GetEntries()<1000000){
      h_phi_corr_binThetaCD[i]->Rebin(5);
    }
    else if(h_phi_corr_binThetaCD[i]->GetEntries()<10000000){
      h_phi_corr_binThetaCD[i]->Rebin(2);
    }
  }

  
  for(int i = 0; i < 14; i++){
    //getGraphCD(h_phi_corr_binThetaCD[i],g_phi_corr_binThetaCD[i],myCanvas,fileName,ctrX);
    //getFunctionCD(g_phi_corr_binThetaCD[i],f_phi_corr_binThetaCD[i],p_phi_corr_binThetaCD[i],guess_CD);    
    f_phi_corr_binThetaCD[i]->SetLineColor(3);
    f_phi_corr_binThetaCD[i]->SetLineWidth(1);
    f_phi_corr_binThetaCD[i]->SetParameter(0,guess_CD[0]);
    f_phi_corr_binThetaCD[i]->SetParLimits(0,-0.5,0.5);
    f_phi_corr_binThetaCD[i]->SetParameter(1,0.05);
    f_phi_corr_binThetaCD[i]->SetParLimits(1,0.001,0.5);
    f_phi_corr_binThetaCD[i]->SetParameter(2,guess_CD[2]);
    f_phi_corr_binThetaCD[i]->SetParLimits(2,-M_PI,M_PI);
    p_phi_corr_binThetaCD[i] = h_phi_corr_binThetaCD[i]->Fit(f_phi_corr_binThetaCD[i],"SrBeqn","",-180,180);


    double x = (bE_ThetaCD[i]+bE_ThetaCD[i+1])/2;
    for(int k = 0; k < 3; k++){
      params_CD[k]->SetPoint(params_CD[k]->GetN(),x,p_phi_corr_binThetaCD[i]->Parameter(k));
      params_CD_Errors[k]->SetPoint(params_CD_Errors[k]->GetN(),x,p_phi_corr_binThetaCD[i]->Parameter(k));
      params_CD_Errors[k]->SetPointError(params_CD_Errors[k]->GetN()-1,0,p_phi_corr_binThetaCD[i]->ParError(k));
      guess_CD[k]=p_phi_corr_binThetaCD[i]->Parameter(k);
    }
  }

  
  TF1 * f_params_CD[3];
  TFitResultPtr p_params_CD[3];
  for(int k = 0; k < 3; k++){
    sprintf(temp_name,"params_CD_%d",k);
    f_params_CD[k] = new TF1(temp_name,[&](double *x, double *p){ return ErfFunc(x[0],p[0],p[1],p[2],p[3]);},38,80,4);
    //f_params_CD[k]->SetParameter(0,params_Momentum_CD[k][0]);
    //f_params_CD[k]->SetParameter(1,params_Momentum_CD[k][1]);
    f_params_CD[k]->SetParameter(0,0);
    f_params_CD[k]->SetParLimits(0,-10,10);
    f_params_CD[k]->SetParameter(1,0);
    f_params_CD[k]->SetParLimits(1,-10,10);
    f_params_CD[k]->SetParameter(2,50);
    f_params_CD[k]->SetParLimits(2,5,20);
    f_params_CD[k]->SetParameter(3,50);
    f_params_CD[k]->SetParLimits(3,38,70);
    double upper = 80;
    p_params_CD[k] = params_CD[k]->Fit(f_params_CD[k],"SrBeqn","",38,upper);
  }
  
  for(int i = 0; i < 18; i++){
    double x = (bE_ThetaCD[i]+bE_ThetaCD[i+1])/2;
    f_phi_corr_Combined_binThetaCD[i]->SetLineColor(4);
    for(int k = 0; k < 3; k++){      
      f_phi_corr_Combined_binThetaCD[i]->SetParameter(k,f_params_CD[k]->Eval(x));
    }
  }
  

  myCanvas->Divide(3,3,0,0);
  for(int i = 0; i < 9; i++){
    myCanvas->cd(i+1);
    myCanvas->GetPad(i+1)->SetBottomMargin(0.19);
    myCanvas->GetPad(i+1)->SetLeftMargin(0.19);    
    h_phi_corr_binThetaCD[i]->Draw("colz");
    g_phi_corr_binThetaCD[i]->SetLineColor(2);
    g_phi_corr_binThetaCD[i]->Draw("SAME");
    f_phi_corr_binThetaCD[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,3,0,0);
  for(int i = 9; i < 14; i++){
    myCanvas->cd(i-8);
    myCanvas->GetPad(i-8)->SetBottomMargin(0.19);
    myCanvas->GetPad(i-8)->SetLeftMargin(0.19);
    h_phi_corr_binThetaCD[i]->Draw("colz");
    g_phi_corr_binThetaCD[i]->SetLineColor(2);
    g_phi_corr_binThetaCD[i]->Draw("SAME");
    f_phi_corr_binThetaCD[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  char* names_CD[3] = {"C_{a}","C_{b}","C_{c}"};
  myCanvas->Divide(2,2,0,0);
  for(int i = 0; i < 3; i++){
    myCanvas->cd(i+1);
    myCanvas->GetPad(i+1)->SetBottomMargin(0.19);
    myCanvas->GetPad(i+1)->SetLeftMargin(0.19);
    sprintf(temp_title,"%s vs. #theta #circ;#theta #circ;%s",names_CD[i],names_CD[i]);
    params_CD_Errors[i]->SetTitle(temp_title);      
    params_CD_Errors[i]->SetLineColor(3);
    params_CD_Errors[i]->GetXaxis()->CenterTitle();
    params_CD_Errors[i]->GetXaxis()->SetTitleSize(0.10);
    params_CD_Errors[i]->GetXaxis()->SetLabelSize(0.06);
    params_CD_Errors[i]->GetXaxis()->SetTitleOffset(0.8);
    params_CD_Errors[i]->GetYaxis()->CenterTitle();
    params_CD_Errors[i]->GetYaxis()->SetTitleSize(0.10);
    params_CD_Errors[i]->GetYaxis()->SetLabelSize(0.06);
    params_CD_Errors[i]->GetYaxis()->SetTitleOffset(0.8);
    params_CD_Errors[i]->GetXaxis()->CenterTitle();
    params_CD_Errors[i]->GetYaxis()->CenterTitle();
    params_CD_Errors[i]->Draw();
    //params_CD[i]->SetTitle(temp_title);      
    //params_CD[i]->SetLineColor(3);
    //params_CD[i]->Draw();
    f_params_CD[i]->Draw("SAME");      
    f_params_CD[i]->SetLineColor(4);      
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,3,0,0);
  for(int i = 0; i < 9; i++){
    myCanvas->cd(i+1);
    myCanvas->GetPad(i+1)->SetBottomMargin(0.19);
    myCanvas->GetPad(i+1)->SetLeftMargin(0.19);
    h_phi_corr_binThetaCD[i]->Draw("colz");
    g_phi_corr_binThetaCD[i]->SetLineColor(2);
    g_phi_corr_binThetaCD[i]->Draw("SAME");
    f_phi_corr_Combined_binThetaCD[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,3,0,0);
  for(int i = 9; i < 18; i++){
    myCanvas->cd(i-8);
    myCanvas->GetPad(i-8)->SetBottomMargin(0.19);
    myCanvas->GetPad(i-8)->SetLeftMargin(0.19);
    h_phi_corr_binThetaCD[i]->Draw("colz");
    g_phi_corr_binThetaCD[i]->SetLineColor(2);
    g_phi_corr_binThetaCD[i]->Draw("SAME");
    f_phi_corr_Combined_binThetaCD[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  TF1 * f_ab1 = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.2,0.2,3);
  f_ab1->SetParameter(0,h_Delta_Int_Eprime_Corrected->GetMaximum()/G(0,1,0,0.1));
  f_ab1->SetParameter(1,h_Delta_Int_Eprime_Corrected->GetBinCenter(h_Delta_Int_Eprime_Corrected->GetMaximumBin()));
  f_ab1->SetParameter(2,0.05);
  TFitResultPtr p_ab1 = h_Delta_Int_Eprime_Corrected->Fit(f_ab1,"SrBeqn","",-0.025,0.05);
  f_ab1->SetLineColor(2);

  TF1 * f_ab2 = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.2,0.2,3);
  f_ab2->SetParameter(0,h_Delta_Int_Eprime_beforeRad->GetMaximum()/G(0,1,0,0.1));
  f_ab2->SetParameter(1,h_Delta_Int_Eprime_beforeRad->GetBinCenter(h_Delta_Int_Eprime_beforeRad->GetMaximumBin()));
  f_ab2->SetParameter(2,0.05);
  TFitResultPtr p_ab2 = h_Delta_Int_Eprime_beforeRad->Fit(f_ab2,"SrBeqn","",-0.025,0.06);
  f_ab2->SetLineColor(4);

  TF1 * f_ab3 = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.4,0.4,3);
  f_ab3->SetParameter(0,h_Delta_Int_pMomFD_Corrected->GetMaximum()/G(0,1,0,0.1));
  f_ab3->SetParameter(1,h_Delta_Int_pMomFD_Corrected->GetBinCenter(h_Delta_Int_pMomFD_Corrected->GetMaximumBin()));
  f_ab3->SetParameter(2,0.05);
  TFitResultPtr p_ab3 = h_Delta_Int_pMomFD_Corrected->Fit(f_ab3,"SrBeqn","",-0.15,0.1);
  f_ab3->SetLineColor(2);

  TF1 * f_ab4 = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.4,0.4,3);
  f_ab4->SetParameter(0,h_Delta_Int_pMomFD->GetMaximum()/G(0,1,0,0.1));
  f_ab4->SetParameter(1,h_Delta_Int_pMomFD->GetBinCenter(h_Delta_Int_pMomFD->GetMaximumBin()));
  f_ab4->SetParameter(2,0.05);
  TFitResultPtr p_ab4 = h_Delta_Int_pMomFD->Fit(f_ab4,"SrBeqn","",-0.15,0.1);
  f_ab4->SetLineColor(4);

  TF1 * f_ab5 = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.4,0.4,3);
  f_ab5->SetParameter(0,h_Delta_Int_pMomCD_Corrected->GetMaximum()/G(0,1,0,0.1));
  f_ab5->SetParameter(1,h_Delta_Int_pMomCD_Corrected->GetBinCenter(h_Delta_Int_pMomCD_Corrected->GetMaximumBin()));
  f_ab5->SetParameter(2,0.05);
  TFitResultPtr p_ab5 = h_Delta_Int_pMomCD_Corrected->Fit(f_ab5,"SrBeqn","",-0.2,0.2);
  f_ab5->SetLineColor(2);

  TF1 * f_ab6 = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.4,0.4,3);
  f_ab6->SetParameter(0,h_Delta_Int_pMomCD->GetMaximum()/G(0,1,0,0.1));
  f_ab6->SetParameter(1,h_Delta_Int_pMomCD->GetBinCenter(h_Delta_Int_pMomCD->GetMaximumBin()));
  f_ab6->SetParameter(2,0.05);
  TFitResultPtr p_ab6 = h_Delta_Int_pMomCD->Fit(f_ab6,"SrBeqn","",-0.2,0.2);
  f_ab6->SetLineColor(4);


  myCanvas->Divide(1,1,0,0);
  myCanvas->cd(1);
  myCanvas->GetPad(1)->SetBottomMargin(0.19);
  myCanvas->GetPad(1)->SetLeftMargin(0.19);
  h_Delta_Int_Eprime_Corrected->SetLineColor(2);
  h_Delta_Int_Eprime_Corrected->Draw();
  h_Delta_Int_Eprime_beforeRad->Draw("SAME");
  f_ab1->Draw("SAME");
  f_ab2->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1,0,0);
  myCanvas->cd(1);
  myCanvas->GetPad(1)->SetBottomMargin(0.19);
  myCanvas->GetPad(1)->SetLeftMargin(0.19);
  h_Delta_Int_pMomFD_Corrected->SetLineColor(2);
  h_Delta_Int_pMomFD_Corrected->Draw();
  h_Delta_Int_pMomFD->Draw("SAME");
  f_ab3->Draw("SAME");
  f_ab4->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1,0,0);
  myCanvas->cd(1);
  myCanvas->GetPad(1)->SetBottomMargin(0.19);
  myCanvas->GetPad(1)->SetLeftMargin(0.19);
  h_Delta_Int_pMomCD_Corrected->SetLineColor(2);
  h_Delta_Int_pMomCD_Corrected->Draw();
  h_Delta_Int_pMomCD->Draw("SAME");
  f_ab5->Draw("SAME");
  f_ab6->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  cout<<"\n\n";
  cout<<"Electrons&"<<p_ab2->Parameter(1)*1000<<"&"<<p_ab1->Parameter(1)*1000<<"&"<<p_ab2->Parameter(2)*1000<<"&"<<p_ab1->Parameter(2)*1000<<"\n";
  cout<<"FD Protons&"<<p_ab4->Parameter(1)*1000<<"&"<<p_ab3->Parameter(1)*1000<<"&"<<p_ab4->Parameter(2)*1000<<"&"<<p_ab3->Parameter(2)*1000<<"\n";
  cout<<"CD Protons&"<<p_ab6->Parameter(1)*1000<<"&"<<p_ab5->Parameter(1)*1000<<"&"<<p_ab6->Parameter(2)*1000<<"&"<<p_ab5->Parameter(2)*1000<<"\n";
  cout<<"\n\n";

  
  cout<<"params_neg_FD="<<endl<<"{";
  for(int j = 0; j < 6; j++){
    cout<<"{";
    for(int k = 0; k < 4; k++){
      cout<<"{";
      for(int l = 0; l < 4; l++){
	if(l!=0){cout<<",";}
	cout<<p_e_params[j][k]->Parameter(l);	
      }
      cout<<"}";
      if(k!=3){cout<<","<<endl;}
    }
    cout<<"}";
    if(j!=5){cout<<","<<endl;}
  }
  cout<<"}"<<endl;
  
  cout<<"params_pos_FD="<<endl<<"{";
  for(int j = 0; j < 6; j++){
    cout<<"{";
    for(int k = 0; k < 4; k++){
      cout<<"{";
      for(int l = 0; l < 4; l++){
	if(l!=0){cout<<",";}
	cout<<p_p_params[j][k]->Parameter(l);	
      }
      cout<<"}";
      if(k!=3){cout<<","<<endl;}
    }
    cout<<"}";
    if(j!=5){cout<<","<<endl;}
  }
  cout<<"}"<<endl;

  
  cout<<"params_CD="<<endl<<"{";
  for(int k = 0; k < 3; k++){
    cout<<"{";
    for(int l = 0; l < 4; l++){
      if(l!=0){cout<<",";}
      cout<<p_params_CD[k]->Parameter(l);      	
    }
    cout<<"}";
    if(k!=2){cout<<","<<endl;}
  }
  cout<<"}"<<endl;  
  
  return 0;
}


