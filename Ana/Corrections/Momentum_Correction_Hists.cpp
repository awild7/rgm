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
#include <TGraph2D.h>
#include <TF2.h>
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
/*
const double pMnF[6][4][4]{{{0.186693, -0.172461, 4.00326, 24.9998 },
{1.84473, 0.770261, 4, 25 },
{123.059, 39.3685, 10, 15 },
{0.243346, -1.30484, 1.00005, 19.0313 }
  },
{{0.899761, 0.609863, 4, 20.5197 },
{1.13272, 0.415673, 4.00009, 23.5633 },
{113.324, 26.7516, 4, 15.5147 },
{-0.857351, -1.71183, 1, 17.6956 }
  },
{{2.07015, 1.265, 7.52494, 23.0353 },
{1.35334, 0.530038, 4.00006, 22.6681 },
{83.8325, -26.7186, 4.00066, 22.4375 },
{-3, -1.89208, 1, 16.0519 }
  },
{{1.77653, 0.803893, 4.85401, 24.9114 },
{0.41588, -0.0783615, 4.00055, 24.9999 },
{159.596, 60.2332, 3.32834, 16.8565 },
{3, 1.04014, 3.38586, 17.9818 }
  },
{{0.888503, 0.957621, 9.99999, 18.2388 },
{0.00256309, -0.333994, 9.99739, 15 },
{105.15, 9.09147, 4, 25 },
{1.80422, 0.917092, 1, 25 }
  },
{{0.935728, 1.20995, 4.00003, 17.8666 },
{0.240479, -0.661936, 4, 16.9237 },
{40.867, -51.7812, 2, 17.1385 },
{-0.569428, -0.906155, 1, 15.9099 }
  }
};
*/
const double pMnF[6][4][4]{{{0.355763, -0.22643, 9.99968, 15.0028 },
{1.74691, 0.691857, 4, 25 },
{111.161, 25.3496, 4, 15 },
{0.107585, -1.45309, 1.00002, 19.0328 }
  },
{{0.916627, 0.713197, 5.23609, 19.9905 },
{0.975205, 0.286991, 4.00002, 24.5373 },
{86.7248, -6.05024, 4, 25 },
{-0.64038, -1.48147, 1, 18.2531 }
  },
{{1.41609, 0.836034, 5.29416, 22.2758 },
{1.2199, 0.741623, 6.38842, 19.9664 },
{71.5074, -9.31258, 4, 25 },
{-3, 0.149729, 5.81782, 25 }
  },
{{1.5198, 0.690357, 6.18876, 24.9967 },
{0.734563, 0.684944, 4.00012, 19.9929 },
{200, 68.9998, 4, 24.8638 },
{-3, -0.572058, 1, 20.2022 }
  },
{{1.39382, 1.03248, 9.99997, 22.4862 },
{0.211263, -0.0713637, 8.76589, 15.0033 },
{36.0234, -16.9183, 4, 23.249 },
{-1.47604, -1.18077, 9.99999, 15 }
  },
{{1.13099, 0.832053, 5.43072, 20.3934 },
{0.31193, -0.145659, 4, 16.3524 },
{37.7433, -16.4775, 4, 19.3609 },
{-0.550815, -0.787611, 1, 15.9985 }
  }
};

const double f_e_guess[6][14][4]={{{0.536517, 0.310507, 61.9058, 2.88675 },
{0.536517, 0.310507, 67.4919, 2.88675 },
{0.536513, 0.310524, 76.1014, 2.88675 },
{0.536442, 0.31082, 86.7898, 2.88675 },
{0.535658, 0.314105, 97.4782, 2.88129 },
{0.530301, 0.33658, 106.088, 1.60538 },
{0.496943, 0.476606, 112.619, 0.240547 },
{0.379026, 0.971812, 115.619, 0.239921 },
{0.239204, 1.55917, 116.283, 0.239921 },
{0.179706, 1.80912, 116.374, 0.239921 },
{0.1707, 1.84695, 116.382, 0.239921 },
{0.170222, 1.84895, 116.382, 0.239921 },
{0.170214, 1.84899, 116.382, 0.239921 },
{0.170213, 1.84899, 116.382, 0.239921 }
  },
{{-0.328594, 0.29939, 56.225, 2.51221 },
{-0.32813, 0.299394, 59.3492, 2.51221 },
{-0.323591, 0.299471, 67.9838, 2.51221 },
{-0.296212, 0.300453, 82.7477, 2.512 },
{-0.194226, 0.308105, 98.3741, 1.98691 },
{0.0406795, 0.344934, 108.613, -0.740557 },
{0.458868, 0.497151, 113.196, -0.857167 },
{0.804564, 0.831355, 113.95, -0.857168 },
{0.891778, 1.07233, 113.981, -0.857168 },
{0.899546, 1.13448, 113.981, -0.857168 },
{0.899786, 1.14014, 113.981, -0.857168 },
{0.899789, 1.14032, 113.981, -0.857168 },
{0.899789, 1.14032, 113.981, -0.857168 },
{0.899789, 1.14032, 113.981, -0.857168 }
  },
{{-0.476871, 0.30078, 135.09, 0.791959 },
{-0.453004, 0.300797, 135.07, 0.791959 },
{-0.401316, 0.301082, 134.951, 0.791922 },
{-0.302726, 0.304053, 134.421, 0.508727 },
{-0.1371, 0.323168, 132.654, -2.68616 },
{0.107971, 0.399119, 128.229, -2.99995 },
{0.516095, 0.650695, 117.209, -3 },
{1.0742, 1.07795, 98.6237, -3 },
{1.56691, 1.31266, 83.7706, -3 },
{1.89505, 1.35866, 77.4131, -3 },
{2.0599, 1.36183, 75.9583, -3 },
{2.12237, 1.36191, 75.7809, -3 },
{2.14022, 1.36191, 75.7694, -3 },
{2.14427, 1.36191, 75.769, -3 }
  },
{{0.139861, 0.585363, 41.5852, 1.06948 },
{0.139927, 0.585363, 43.0041, 1.0695 },
{0.140504, 0.585361, 50.4364, 1.07123 },
{0.144159, 0.585324, 72.0465, 1.11542 },
{0.16114, 0.584923, 106.987, 1.45537 },
{0.218946, 0.582175, 138.432, 2.26013 },
{0.417636, 0.565056, 155.981, 2.91755 },
{0.888167, 0.504521, 159.181, 2.99928 },
{1.40179, 0.432729, 159.3, 2.99986 },
{1.6871, 0.40218, 159.301, 2.99986 },
{1.76761, 0.397557, 159.301, 2.99986 },
{1.77912, 0.397312, 159.301, 2.99986 },
{1.77995, 0.397308, 159.301, 2.99986 },
{1.77998, 0.397308, 159.301, 2.99986 }
  },
{{-0.678465, 0.430871, 82.4811, -0.0362898 },
{-0.616141, 0.388944, 84.2578, -0.0362898 },
{-0.52317, 0.339794, 86.3406, -0.0362898 },
{-0.395074, 0.286577, 88.5957, -0.0362898 },
{-0.232067, 0.23336, 90.8507, -0.0362898 },
{-0.040482, 0.184208, 92.9335, -0.0362898 },
{0.220218, 0.133157, 95.0968, -0.0362887 },
{0.523035, 0.0904142, 96.9079, 0.651333 },
{0.770216, 0.0660367, 97.9408, 2.8312 },
{0.939196, 0.0543931, 98.4342, 2.83178 },
{1.03594, 0.0497358, 98.6315, 2.83178 },
{1.08233, 0.0481758, 98.6976, 2.83178 },
{1.10096, 0.0477382, 98.7161, 2.83178 },
{1.1077, 0.0476294, 98.7208, 2.83178 }
  },
{{-1.48184, 1.56262, 144.495, 1.24634 },
{-1.46094, 1.54718, 144.495, 1.24634 },
{-1.35915, 1.47772, 144.453, 1.2463 },
{-1.05276, 1.28463, 141.71, 1.05676 },
{-0.482317, 0.952518, 112.798, -0.475845 },
{0.175081, 0.59894, 58.7793, -0.582879 },
{0.717563, 0.331374, 39.6039, -0.582889 },
{0.900724, 0.249746, 39.1932, -0.582889 },
{0.917699, 0.243007, 39.1931, -0.582889 },
{0.918245, 0.242814, 39.1931, -0.582889 },
{0.918251, 0.242812, 39.1931, -0.582889 },
{0.918251, 0.242812, 39.1931, -0.582889 },
{0.918251, 0.242812, 39.1931, -0.582889 },
{0.918251, 0.242812, 39.1931, -0.582889 }
  }
};
const double f_p_guess[6][14][4]={{{-1.37882, 2.45281, 136.487, -2.94586 },
{-1.37882, 2.45281, 136.487, -2.94586 },
{-1.37882, 2.45281, 136.487, -2.94586 },
{-1.37882, 2.45281, 136.487, -2.94586 },
{-1.37882, 2.45281, 136.487, -2.94586 },
{-1.37882, 2.45281, 136.487, -2.94586 },
{-1.37882, 2.45281, 136.487, -2.94544 },
{-1.37882, 2.45281, 136.471, -2.91247 },
{-1.37882, 2.45311, 131.842, -2.57826 },
{-1.37881, 2.52143, 98.8453, -2.08195 },
{-1.36649, 2.90946, 89.3182, -1.97081 },
{-0.245933, 2.99933, 89.2435, -1.9674 },
{2.69243, 2.99987, 89.2435, -1.96738 },
{2.99981, 2.99987, 89.2435, -1.96738 }
  },
{{-1.35681, 1.67556, 57.4565, -5.01461 },
{-1.35681, 1.67556, 57.4565, -5.01461 },
{-1.35681, 1.67556, 57.4565, -5.01461 },
{-1.35681, 1.67556, 57.4565, -5.01461 },
{-1.35681, 1.67556, 57.4565, -5.01461 },
{-1.35681, 1.67556, 57.4565, -5.01396 },
{-1.35681, 1.67085, 57.4681, -4.85454 },
{-1.35681, 0.93923, 61.5704, -1.5642 },
{-1.35681, -2.04395, 94.7092, 0.922803 },
{-1.35681, -2.53941, 105.487, 0.985337 },
{-1.34695, -2.5414, 105.585, 0.985365 },
{-0.326153, -2.5414, 105.585, 0.985365 },
{2.65071, -2.5414, 105.585, 0.985365 },
{2.99977, -2.5414, 105.585, 0.985365 }
  },
{{-0.476224, -5.05467, 102.454, 4.81663 },
{-0.476224, -5.05467, 102.454, 4.81663 },
{-0.476224, -5.05467, 102.454, 4.81663 },
{-0.476224, -5.05467, 102.454, 4.81663 },
{-0.476224, -5.05467, 102.454, 4.81663 },
{-0.476224, -5.05412, 102.454, 4.81561 },
{-0.476224, -4.90933, 102.454, 4.62417 },
{-0.476224, -1.70197, 102.443, 1.63046 },
{-0.476224, 0.875561, 99.1091, -0.14287 },
{-0.476223, 0.945294, 74.7943, -0.176134 },
{-0.47226, 0.945328, 67.6181, -0.176145 },
{0.136469, 0.945328, 67.5603, -0.176145 },
{2.59428, 0.945328, 67.5603, -0.176145 },
{2.99957, 0.945328, 67.5603, -0.176145 }
  },
{{-2.04905, 0.841574, 128.109, 1.10468 },
{-2.04905, 0.841574, 128.109, 1.10468 },
{-2.04905, 0.841574, 128.109, 1.10468 },
{-2.04905, 0.841574, 128.109, 1.10468 },
{-2.04905, 0.841574, 128.109, 1.10468 },
{-2.04904, 0.841574, 128.109, 1.10468 },
{-2.049, 0.841574, 128.109, 1.10468 },
{-2.04678, 0.841574, 128.109, 1.10463 },
{-2.00221, 0.841984, 128.109, 1.05799 },
{-1.63239, 0.974729, 128.109, -0.085838 },
{-0.350324, 1.97078, 128.109, -1.28635 },
{1.51941, 2.27283, 128.088, -1.34106 },
{2.66922, 2.27535, 73.113, -1.34113 },
{2.9789, 2.27535, -3697.83, -1.34113 }
  },
{{-0.161849, 2.67097, 142.608, -1.52901 },
{-0.161849, 2.67097, 142.608, -1.52901 },
{-0.161849, 2.67097, 142.608, -1.52901 },
{-0.161849, 2.67097, 142.608, -1.52901 },
{-0.161849, 2.67097, 142.608, -1.52901 },
{-0.161849, 2.67097, 142.608, -1.52901 },
{-0.161849, 2.67097, 142.608, -1.52901 },
{-0.161849, 2.67097, 142.608, -1.52901 },
{-0.161849, 2.67097, 142.608, -1.52898 },
{-0.161849, 2.67097, 142.608, -1.52688 },
{-0.156867, 2.66916, 142.618, -1.48058 },
{0.392779, 2.11503, 144.189, -1.08785 },
{2.09715, -1.84482, 150.777, 0.205405 },
{2.3122, -2.99729, 151.907, 2.08544 }
  },
{{-1.04758, -1.46923, 56.0239, 2.33031 },
{-1.04758, -1.46923, 56.0239, 2.33031 },
{-1.04758, -1.46923, 56.0239, 2.33031 },
{-1.04758, -1.46923, 56.0239, 2.33031 },
{-1.04758, -1.46923, 56.0239, 2.33031 },
{-1.04758, -1.46923, 56.0239, 2.33031 },
{-1.04758, -1.46923, 56.0239, 2.33031 },
{-1.04756, -1.46923, 56.024, 2.33031 },
{-1.04436, -1.46923, 56.342, 2.33031 },
{-0.931308, -1.46923, 81.1767, 2.33031 },
{0.00603735, -1.46858, 138.638, 2.33031 },
{1.91016, -1.1713, 143.876, 2.33031 },
{2.87628, 1.80422, 143.886, 2.33019 },
{2.9983, 2.9956, 143.886, 2.23281 }
  }
};

const double c = 29.9792458;

vector<double> bE_Theta = {8,10,12,14,16,18,20,23,26,29,32,35,38,41,45};
int binX(vector<double> XS, double X){
  for(int i = 0; i < XS.size(); i++){if(X<XS[i]){return i-1;}}
  return -1;
}
double G(double x, double N, double mu, double sigma){ return (N/(sigma*sqrt(2*M_PI))) * exp(-0.5 * sq((x-mu)/sigma)); }
double FuncPhiDependenceFD(double x, double A, double B, double C, double D){ return A + B*sin((x*2*M_PI/C)+D); }
double ErfFunc(double x, double A, double B, double C, double D){ return A - B*(1+erf(((-x+D)/C))); }
double SQ(double x){ return x*x;}
void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){ p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());}


//void centralValue(TH2D * h_myhist, TGraphErrors * g_mygraph, TCanvas * myCanvas, char * fileName, int & ctr){
void centralValue(TH2D * h_myhist, TGraphErrors * g_mygraph, int i , int & ctr){

  bool islarge=(h_myhist->GetEntries()<2000)?false:true;
  if(i>3){
    h_myhist->RebinY(2);
  }
  if(i>5){
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
    if(i==0&&proj->GetEntries()<2000){continue;}
    if(proj->GetEntries()<20){continue;}
    //proj->Smooth(1);
    TF1 * gFit = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },cent-(stddev*2),cent+(stddev*2),3);
    gFit->SetParameter(0,proj->GetMaximum()/G(0,1,0,1));
    gFit->SetParameter(1,cent);
    gFit->SetParameter(2,stddev);
    gFit->SetParLimits(1,cent-0.5*stddev,cent+0.5*stddev);
    gFit->SetParLimits(2,stddev*0.5,stddev*1.2);
    
    TFitResultPtr gPoint = proj->Fit(gFit,"SrBeqn","",cent-1*stddev,cent+0.5*stddev);
    if(gPoint != -1){
      g_mygraph->SetPoint(g_mygraph->GetN(),x,gPoint->Parameter(1));
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


void getFunctionConst(TGraphErrors * g_mygraph, TF1 * f_myfunc, TFitResultPtr & p_mypoint, const double guess[4],int i, int sector){
  double min0 = -3;
  double max0 = 3;

  double guess1 = (guess[1]<0)?-guess[1]:guess[1];
  double min1 = 0;
  double max1 = 2;

  double guess2 = (guess[2]<0)?50:guess[2];
  double min2 = ((guess[2]*0.6)>30)?guess[2]*0.6:30;
  double max2 = (guess[2]<30)?150:((guess[2]*1.5)<150)?guess[2]*1.5:150;
  guess2 = ((guess[2]<min2)||(guess[2]>max2))?(max2+min2)/2:guess[2];
  double min3 = ((guess[3]-1)>-M_PI)?guess[3]-1:-M_PI;
  double max3 = ((guess[3]+1)< M_PI)?guess[3]+1: M_PI;


  if(sector==1){
    min2=guess2*0.9;
    max2=guess2*1.1;
  }
  if(sector==2){
    min2=guess2*0.9;
    max2=guess2*1.1;
  }
  if(sector==3){
    min2=guess2*0.9;
    max2=guess2*1.1;
  }
  if(sector==5){
    min2=guess2*0.9;
    max2=guess2*1.1;
  }

  
  if(sector==4){
    max1=0.7;
    min3 = 0;
    max3 =  2*M_PI;
  }
  if(sector==5){
    min2 = 80;
    guess2 =100;
  }
  
  f_myfunc->SetParameter(0,guess[0]);
  f_myfunc->SetParameter(1,guess1);
  f_myfunc->SetParameter(2,guess2);
  f_myfunc->SetParameter(3,guess[3]);
  f_myfunc->SetParLimits(0,min0,max0);
  f_myfunc->SetParLimits(1,min1,max1);  
  f_myfunc->SetParLimits(2,min2,max2);
  f_myfunc->SetParLimits(3,min3,max3);

  p_mypoint = g_mygraph->Fit(f_myfunc,"SrBeqn","",-40,40);
}

void getFunction_WithGoodGuess(TGraphErrors * g_mygraph, TF1 * f_myfunc, TFitResultPtr & p_mypoint, const double guess[4],int i, int sector){

  double guess0 = guess[0];
  double guess1 = guess[1];
  double guess2 = guess[2];
  double guess3 = guess[3];
  if(sector==4){
    //guess3-=2*M_PI;
  }
  double min0 = guess0-0.5;
  double max0 = guess0+0.5;
  double min1 = guess1-0.5;
  double max1 = guess1+0.5;
  double min2 = guess2*0.75;
  double max2 = guess2*1.25;
  double min3 = guess3-0.2*M_PI;
  double max3 = guess3+0.2*M_PI;
  f_myfunc->SetParameter(0,guess0);
  f_myfunc->SetParameter(1,guess1);
  f_myfunc->SetParameter(2,guess2);
  f_myfunc->SetParameter(3,guess3);
  f_myfunc->SetParLimits(0,min0,max0);
  f_myfunc->SetParLimits(1,min1,max1);  
  f_myfunc->SetParLimits(2,min2,max2);
  f_myfunc->SetParLimits(3,min3,max3);
  p_mypoint = g_mygraph->Fit(f_myfunc,"SrBeqn","",-40,40);
}

void getFunction_ProtonConst(TGraphErrors * g_mygraph, TF1 * f_myfunc, TFitResultPtr & p_mypoint, const double guess[4],int i, int sector){

  double min0 = -3;
  double max0 = 3;
  double min1 = -3;
  double max1 = 3;
  double min2 = 60;
  double max2 = 150;
  double min3 = -M_PI;
  double max3 =  M_PI;
  f_myfunc->SetParameter(0,guess[0]);
  f_myfunc->SetParameter(1,guess[1]);
  f_myfunc->SetParameter(2,guess[2]);
  f_myfunc->SetParameter(3,guess[3]);
  f_myfunc->SetParLimits(0,min0,max0);
  f_myfunc->SetParLimits(1,min1,max1);  
  f_myfunc->SetParLimits(2,min2,max2);
  f_myfunc->SetParLimits(3,min3,max3);
  
  p_mypoint = g_mygraph->Fit(f_myfunc,"SrBeqn","",-40,40);
}

void getSecondFunction(TGraphErrors * g_mygraph, TF1 * f_myfunc, TFitResultPtr & p_mypoint, int k, int sector){

  double fitmin = 6;
  double fitmax = 50;
  
  double guess0 = 0;
  double min0 = -3;
  double max0 = 3;
  double guess1 = 0;
  double min1 = -3;
  double max1 = 3;
  double guess2 = 7;
  double min2 = 4;
  double max2 = 10;
  double guess3 = 20;
  double min3 = 15;
  double max3 = 25;

  if(k==2){
    min0=-200;
    max0=200;
    min1=-200;
    max1=200;
  }
  
  if(k==3){
    min2=1;
    max2=10;
  }


  if((sector==4) && (k==2)){
    //guess0=150;
    //min0=140;
    //guess1=00;
    //max1=40;
    min2=2;
    guess3=12;
    min3=10;
    max3=20;
    fitmax=23;
  }

  if((sector==6) && (k==2)){
    guess2 = 2;
    min2 = 2;
    max2 = 5;
    fitmin=12;
  }
  
  f_myfunc->SetParameter(0,guess0);
  f_myfunc->SetParameter(1,guess1);
  f_myfunc->SetParameter(2,guess2);
  f_myfunc->SetParameter(3,guess3);
  f_myfunc->SetParLimits(0,min0,max0);
  f_myfunc->SetParLimits(1,min1,max1);
  f_myfunc->SetParLimits(2,min2,max2);
  f_myfunc->SetParLimits(3,min3,max3);

  
  p_mypoint = g_mygraph->Fit(f_myfunc,"SrBeqn","",fitmin,fitmax);

}

void getSecondFunction_New(TGraphErrors * g_mygraph, TF1 * f_myfunc, TFitResultPtr & p_mypoint, int k, int sector){

  double fitmin = 6;
  double fitmax = 50;
  
  double guess0 = 0;
  double min0 = -5;
  double max0 = 5;
  double guess1 = 0;
  double min1 = -5;
  double max1 = 5;
  double guess2 = 7;
  double min2 = 4;
  double max2 = 10;
  double guess3 = 20;
  double min3 = 15;
  double max3 = 25;
  if(k==2){
    min0=-200;
    max0=200;
    min1=-200;
    max1=200;
  }
  
  if(k==3){
    min2=1;
    max2=10;
  }
  f_myfunc->SetParameter(0,guess0);
  f_myfunc->SetParameter(1,guess1);
  f_myfunc->SetParameter(2,guess2);
  f_myfunc->SetParameter(3,guess3);
  f_myfunc->SetParLimits(0,min0,max0);
  f_myfunc->SetParLimits(1,min1,max1);
  f_myfunc->SetParLimits(2,min2,max2);
  f_myfunc->SetParLimits(3,min3,max3);

  
  p_mypoint = g_mygraph->Fit(f_myfunc,"SrBeqn","",fitmin,fitmax);

}

void getSecondFunction_Proton(TGraphErrors * g_mygraph, TF1 * f_myfunc, TFitResultPtr & p_mypoint, int k, int sector){


  f_myfunc->SetParameter(0,0);
  f_myfunc->SetParameter(1,0);
  f_myfunc->SetParameter(2,10);
  f_myfunc->SetParameter(3,35);
  /*
  if(k==2){
    f_myfunc->SetParLimits(0,-50,50);
    f_myfunc->SetParLimits(1,20,150);
  }
  */

  if(k!=2){
    f_myfunc->SetParLimits(0,-3,3);
    f_myfunc->SetParLimits(1,-3,3);
  }  
  f_myfunc->SetParLimits(2,2,20);
  f_myfunc->SetParLimits(3,24,45);
  p_mypoint = g_mygraph->Fit(f_myfunc,"SrBeqn","",22,40);

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


  //////////////////////////////////
  char temp_name[100];
  char temp_title[100];
  int ctrX=0;

  vector<TH1*> hist_list;
    
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

  
  
  TH2D * h2_e_cvp_FIT0[6][14];
  TGraphErrors * g1_e_cvp_FIT1[6][14];
  TF1 * f_e_frac0[6][14];
  TFitResultPtr p_e_frac[6][14];
  TF1 * f1_e_cvp_FIT4[6][14];
  TGraph2D * g2_e_cvpvt_FIT1[6];
  for(int j=1; j<=6; j++){

    g2_e_cvpvt_FIT1[j-1] = new TGraph2D();
    sprintf(temp_name,"g2_e_cvpvt_FIT1_%d",j);
    g2_e_cvpvt_FIT1[j-1]->SetName(temp_name);
    
    for(int i=0; i<14; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      double theta = (bE_Theta[i+1]+bE_Theta[i])/2;
      
      sprintf(temp_title,"Sector = %d (%d #circ< #theta < %d #circ);#phi;Correction;Counts",j,min,max);
      sprintf(temp_name,"e_frac_%d_%d",j,i);
      h2_e_cvp_FIT0[j-1][i] = (TH2D*)inFile->Get(temp_name);
      h2_e_cvp_FIT0[j-1][i]->SetName(temp_name);
      hist_list.push_back(h2_e_cvp_FIT0[j-1][i]);

      g1_e_cvp_FIT1[j-1][i] = new TGraphErrors();
      sprintf(temp_name,"g1_e_cvp_FIT1_%d_%d",j,i);
      g1_e_cvp_FIT1[j-1][i]->SetName(temp_name);
      centralValue(h2_e_cvp_FIT0[j-1][i],g1_e_cvp_FIT1[j-1][i],i,ctrX); 
      
      for(int k = 0; k < g1_e_cvp_FIT1[j-1][i]->GetN(); k++){
	double phi,z;
	g1_e_cvp_FIT1[j-1][i]->GetPoint(k,phi,z);
	g2_e_cvpvt_FIT1[j-1]->SetPoint(g2_e_cvpvt_FIT1[j-1]->GetN(),theta,phi,z);
      }
      
      sprintf(temp_name,"f_e_frac0_%d_%d",j,i);
      f_e_frac0[j-1][i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceFD(x[0],p[0],p[1],p[2],p[3]); },-40,40,4);
      getFunctionConst(g1_e_cvp_FIT1[j-1][i],f_e_frac0[j-1][i],p_e_frac[j-1][i],f_e_guess[j-1][i],i,j);      
      
      sprintf(temp_name,"f1_e_cvp_FIT4_%d_%d",j,i);
      f1_e_cvp_FIT4[j-1][i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceFD(x[0],p[0],p[1],p[2],p[3]); },-40,40,4);
    }
  }
  
  TH2D * h_p_frac[6][14];
  TGraphErrors * g_p_frac[6][14];
  TF1 * f_p_frac[6][14];
  TFitResultPtr p_p_frac[6][14];
  TF1 * f_p_frac_comb[6][14];
  TGraph2D * g2_p_frac[6];
  for(int j=1; j<=6; j++){

    g2_p_frac[j-1] = new TGraph2D();
    sprintf(temp_name,"g2_p_frac_%d",j);
    g2_p_frac[j-1]->SetName(temp_name);

    for(int i=0; i<14; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      sprintf(temp_name,"p_frac_%d_%d",j,i);
      sprintf(temp_title,"Sector = %d (%d #circ< #theta < %d #circ);#phi;Correction;Counts",j,min,max);
      h_p_frac[j-1][i] = (TH2D*)inFile->Get(temp_name);
      h_p_frac[j-1][i]->SetName(temp_name);
      hist_list.push_back(h_p_frac[j-1][i]);

      g_p_frac[j-1][i] = new TGraphErrors();
      sprintf(temp_name,"g_p_frac_%d_%d",j,i);
      g_p_frac[j-1][i]->SetName(temp_name);
      centralValue(h_p_frac[j-1][i],g_p_frac[j-1][i],i,ctrX); 

      sprintf(temp_name,"f_p_frac_%d_%d",j,i);
      f_p_frac[j-1][i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceFD(x[0],p[0],p[1],p[2],p[3]); },-40,40,4);
      getFunction_ProtonConst(g_p_frac[j-1][i],f_p_frac[j-1][i],p_p_frac[j-1][i],f_p_guess[j-1][i],i,j);      
      
      sprintf(temp_name,"f_p_frac_comb_%d_%d",j,i);
      f_p_frac_comb[j-1][i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceFD(x[0],p[0],p[1],p[2],p[3]); },-40,40,4);
    }
  }

  ////////////
  //Fit Params
  ////////////
  ////////////
  //Fit Params
  ////////////
  ////////////
  //Fit Params
  ////////////  
  
  TF2 * f2_e_cvpvt_FIT2[6];
  TFitResultPtr p2_e_cvpvt_FIT2[6];
  TF1 * f1_e_cvt_FIT2[6][4];  
  for(int j=1; j<=6; j++){
    ////////////////////////////
    ////////////////////////////
    ////////////////////////////
    sprintf(temp_name,"f2_e_cvpvt_FIT2_%d",j);
    f2_e_cvpvt_FIT2[j-1] = new TF2(temp_name,[&](double *x, double *p){ return FuncPhiDependenceFD(x[1],ErfFunc(x[0],p[0],p[1],p[2],p[3]),ErfFunc(x[0],p[4],p[5],p[6],p[7]),ErfFunc(x[0],p[8],p[9],p[10],p[11]),ErfFunc(x[0],p[12],p[13],p[14],p[15]));},5,50,-40,40,16);
    f2_e_cvpvt_FIT2[j-1]->SetParameter(0,pMnF[j-1][0][0]);
    f2_e_cvpvt_FIT2[j-1]->SetParameter(1,pMnF[j-1][0][1]);
    f2_e_cvpvt_FIT2[j-1]->SetParameter(2,pMnF[j-1][0][2]);
    f2_e_cvpvt_FIT2[j-1]->SetParameter(3,pMnF[j-1][0][3]);
    f2_e_cvpvt_FIT2[j-1]->SetParameter(4,pMnF[j-1][1][0]);
    f2_e_cvpvt_FIT2[j-1]->SetParameter(5,pMnF[j-1][1][1]);
    f2_e_cvpvt_FIT2[j-1]->SetParameter(6,pMnF[j-1][1][2]);
    f2_e_cvpvt_FIT2[j-1]->SetParameter(7,pMnF[j-1][1][3]);
    f2_e_cvpvt_FIT2[j-1]->SetParameter(8,pMnF[j-1][2][0]);
    f2_e_cvpvt_FIT2[j-1]->SetParameter(9,pMnF[j-1][2][1]);
    f2_e_cvpvt_FIT2[j-1]->SetParameter(10,pMnF[j-1][2][2]);
    f2_e_cvpvt_FIT2[j-1]->SetParameter(11,pMnF[j-1][2][3]);
    f2_e_cvpvt_FIT2[j-1]->SetParameter(12,pMnF[j-1][3][0]);
    f2_e_cvpvt_FIT2[j-1]->SetParameter(13,pMnF[j-1][3][1]);
    f2_e_cvpvt_FIT2[j-1]->SetParameter(14,pMnF[j-1][3][2]);
    f2_e_cvpvt_FIT2[j-1]->SetParameter(15,pMnF[j-1][3][3]);
    f2_e_cvpvt_FIT2[j-1]->SetParLimits(0,-3,3);
    f2_e_cvpvt_FIT2[j-1]->SetParLimits(1,-3,3);
    f2_e_cvpvt_FIT2[j-1]->SetParLimits(2,4,10);
    f2_e_cvpvt_FIT2[j-1]->SetParLimits(3,15,25);
    f2_e_cvpvt_FIT2[j-1]->SetParLimits(4,-3,3);
    f2_e_cvpvt_FIT2[j-1]->SetParLimits(5,-3,3);
    f2_e_cvpvt_FIT2[j-1]->SetParLimits(6,4,10);
    f2_e_cvpvt_FIT2[j-1]->SetParLimits(7,15,25);
    switch(j){
    case 1:
      f2_e_cvpvt_FIT2[j-1]->SetParLimits(8,00,150);
      f2_e_cvpvt_FIT2[j-1]->SetParLimits(9,0,50);
      break;
    case 2:
      f2_e_cvpvt_FIT2[j-1]->SetParLimits(8,00,150);
      f2_e_cvpvt_FIT2[j-1]->SetParLimits(9,0,50);
      break;
    case 3:
      f2_e_cvpvt_FIT2[j-1]->SetParLimits(8,80,150);
      f2_e_cvpvt_FIT2[j-1]->SetParLimits(9,-40,40);
      break;
    case 4:
      f2_e_cvpvt_FIT2[j-1]->SetParLimits(8,130,170);
      f2_e_cvpvt_FIT2[j-1]->SetParLimits(9,60,100);
      break;
    case 5:
      f2_e_cvpvt_FIT2[j-1]->SetParLimits(8,80,120);
      f2_e_cvpvt_FIT2[j-1]->SetParLimits(9,2,30);
      break;
    default:
      f2_e_cvpvt_FIT2[j-1]->SetParLimits(8,-50,50);
      f2_e_cvpvt_FIT2[j-1]->SetParLimits(9,-50,50);
      break;
    }
    switch(j){
    case 4:
      f2_e_cvpvt_FIT2[j-1]->SetParLimits(10,4,7);
      break;
    default:
      f2_e_cvpvt_FIT2[j-1]->SetParLimits(10,6,15);      
    }
    f2_e_cvpvt_FIT2[j-1]->SetParLimits(11,15,25);
    f2_e_cvpvt_FIT2[j-1]->SetParLimits(12,-3,3);
    f2_e_cvpvt_FIT2[j-1]->SetParLimits(13,-3,3);
    f2_e_cvpvt_FIT2[j-1]->SetParLimits(14,4,10);
    f2_e_cvpvt_FIT2[j-1]->SetParLimits(15,15,25);
    p2_e_cvpvt_FIT2[j-1] = g2_e_cvpvt_FIT1[j-1]->Fit(f2_e_cvpvt_FIT2[j-1],"SrBeqn");
    ////////////////////////////
    ////////////////////////////
    ////////////////////////////

    for(int k=0; k<4; k++){
      f1_e_cvt_FIT2[j-1][k] = new TF1(temp_name,[&](double *x, double *p){ return ErfFunc(x[0],p[0],p[1],p[2],p[3]);},5,50,4);
      f1_e_cvt_FIT2[j-1][k]->SetParameter(0,p2_e_cvpvt_FIT2[j-1]->Parameter(0+4*k));
      f1_e_cvt_FIT2[j-1][k]->SetParameter(1,p2_e_cvpvt_FIT2[j-1]->Parameter(1+4*k));
      f1_e_cvt_FIT2[j-1][k]->SetParameter(2,p2_e_cvpvt_FIT2[j-1]->Parameter(2+4*k));
      f1_e_cvt_FIT2[j-1][k]->SetParameter(3,p2_e_cvpvt_FIT2[j-1]->Parameter(3+4*k));
      
    }
  } 
   
  TF1 * f1_e_cvp_FIT3[6][14];
  TFitResultPtr p1_e_cvp_FIT3[6][14];
  for(int j=1; j<=6; j++){
    for(int i=0; i<9; i++){
      double theta = (bE_Theta[i+1]+bE_Theta[i])/2;
      double guess[4];

      sprintf(temp_name,"f1_e_cvp_FIT3_%d_%d",j,i);
      f1_e_cvp_FIT3[j-1][i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceFD(x[0],p[0],p[1],p[2],p[3]); },-40,40,4);

      guess[0] = ErfFunc(theta,p2_e_cvpvt_FIT2[j-1]->Parameter(0),p2_e_cvpvt_FIT2[j-1]->Parameter(1),p2_e_cvpvt_FIT2[j-1]->Parameter(2),p2_e_cvpvt_FIT2[j-1]->Parameter(3));
      guess[1] = ErfFunc(theta,p2_e_cvpvt_FIT2[j-1]->Parameter(4),p2_e_cvpvt_FIT2[j-1]->Parameter(5),p2_e_cvpvt_FIT2[j-1]->Parameter(6),p2_e_cvpvt_FIT2[j-1]->Parameter(7));
      guess[2] = ErfFunc(theta,p2_e_cvpvt_FIT2[j-1]->Parameter(8),p2_e_cvpvt_FIT2[j-1]->Parameter(9),p2_e_cvpvt_FIT2[j-1]->Parameter(10),p2_e_cvpvt_FIT2[j-1]->Parameter(11));
      guess[3] = ErfFunc(theta,p2_e_cvpvt_FIT2[j-1]->Parameter(12),p2_e_cvpvt_FIT2[j-1]->Parameter(13),p2_e_cvpvt_FIT2[j-1]->Parameter(14),p2_e_cvpvt_FIT2[j-1]->Parameter(15));      
      getFunction_WithGoodGuess(g1_e_cvp_FIT1[j-1][i],f1_e_cvp_FIT3[j-1][i],p1_e_cvp_FIT3[j-1][i],guess,i,j);      
    }
  }

  TGraphErrors * g1_e_cvt_FIT3[6][4];
  TF1 * f1_e_cvt_FIT4[6][4];
  TFitResultPtr p1_e_cvt_FIT4[6][4];
  for(int j=1; j<=6; j++){
    for(int k=0; k<4; k++){
      g1_e_cvt_FIT3[j-1][k] = new TGraphErrors();
      for(int i=0; i<9; i++){
	int x = (bE_Theta[i]+bE_Theta[i+1])/2.0;
        g1_e_cvt_FIT3[j-1][k]->SetPoint(g1_e_cvt_FIT3[j-1][k]->GetN(),x,p1_e_cvp_FIT3[j-1][i]->Parameter(k));
      }
      sprintf(temp_name,"g_e_par_%d_%d",j,k);
      f1_e_cvt_FIT4[j-1][k] = new TF1(temp_name,[&](double *x, double *p){ return ErfFunc(x[0],p[0],p[1],p[2],p[3]);},5,50,4);
      getSecondFunction_New(g1_e_cvt_FIT3[j-1][k],f1_e_cvt_FIT4[j-1][k],p1_e_cvt_FIT4[j-1][k],k,j);
    }
  }
  /*
  TGraphErrors * g_e_par[6][4];
  for(int j=1; j<=6; j++){
    for(int k=0; k<4; k++){
      g_e_par[j-1][k] = new TGraphErrors();
      for(int i=0; i<9; i++){
	int x = (bE_Theta[i]+bE_Theta[i+1])/2.0;
	g_e_par[j-1][k]->SetPoint(g_e_par[j-1][k]->GetN(),x,p_e_frac[j-1][i]->Parameter(k));
	//g_e_par[j-1][k]->SetPointError(g_e_par[j-1][k]->GetN()-1,0,p_e_frac[j-1][i]->ParError(k));
      }
  
    }
  }
  */
  /*
  for(int j=1; j<=6; j++){
    for(int k=0; k<4; k++){
      g_e_par[j-1][k] = new TGraphErrors();
      for(int i=0; i<9; i++){
	int x = (bE_Theta[i]+bE_Theta[i+1])/2.0;
	g_e_par[j-1][k]->SetPoint(g_e_par[j-1][k]->GetN(),x,p_e_frac[j-1][i]->Parameter(k));
	//g_e_par[j-1][k]->SetPointError(g_e_par[j-1][k]->GetN()-1,0,p_e_frac[j-1][i]->ParError(k));
      }
      sprintf(temp_name,"g_e_par_%d_%d",j,k);
      f1_e_cvt_FIT4[j-1][k] = new TF1(temp_name,[&](double *x, double *p){ return ErfFunc(x[0],p[0],p[1],p[2],p[3]);},5,50,4);
      getSecondFunction(g_e_par[j-1][k],f1_e_cvt_FIT4[j-1][k],p1_e_cvt_FIT4[j-1][k],k,j);
    }
  }
  */
  /*
  for(int j=1; j<=6; j++){
    for(int i=0; i<14; i++){
      double theta = (bE_Theta[i+1]+bE_Theta[i])/2;
      double guess[4];
      guess[0] = ErfFunc(theta,p2_e_cvpvt_FIT2[j-1]->Parameter(0),p2_e_cvpvt_FIT2[j-1]->Parameter(1),p2_e_cvpvt_FIT2[j-1]->Parameter(2),p2_e_cvpvt_FIT2[j-1]->Parameter(3));
      guess[1] = ErfFunc(theta,p2_e_cvpvt_FIT2[j-1]->Parameter(4),p2_e_cvpvt_FIT2[j-1]->Parameter(5),p2_e_cvpvt_FIT2[j-1]->Parameter(6),p2_e_cvpvt_FIT2[j-1]->Parameter(7));
      guess[2] = ErfFunc(theta,p2_e_cvpvt_FIT2[j-1]->Parameter(8),p2_e_cvpvt_FIT2[j-1]->Parameter(9),p2_e_cvpvt_FIT2[j-1]->Parameter(10),p2_e_cvpvt_FIT2[j-1]->Parameter(11));
      guess[3] = ErfFunc(theta,p2_e_cvpvt_FIT2[j-1]->Parameter(12),p2_e_cvpvt_FIT2[j-1]->Parameter(13),p2_e_cvpvt_FIT2[j-1]->Parameter(14),p2_e_cvpvt_FIT2[j-1]->Parameter(15));      
      getFunction_WithGoodGuess(g1_e_cvp_FIT1[j-1][i],f1_e_cvp_FIT3[j-1][i],p1_e_cvp_FIT3[j-1][i],guess,i,j);      
    }
  }

  */
  ////////////////////////////////
  TGraphErrors * g_p_par[6][4];
  TF1 * f_p_par[6][4];
  TFitResultPtr p_p_par[6][4];
  for(int j=1; j<=6; j++){
    for(int k=0; k<4; k++){
      g_p_par[j-1][k] = new TGraphErrors();
      for(int i=7; i<13; i++){
	int x = (bE_Theta[i]+bE_Theta[i+1])/2.0;
	g_p_par[j-1][k]->SetPoint(g_p_par[j-1][k]->GetN(),x,p_p_frac[j-1][i]->Parameter(k));
	//g_p_par[j-1][k]->SetPointError(g_p_par[j-1][k]->GetN()-1,0,p_p_frac[j-1][i]->ParError(k));
      }
      sprintf(temp_name,"g_p_par_%d_%d",j,k);
      f_p_par[j-1][k] = new TF1(temp_name,[&](double *x, double *p){ return ErfFunc(x[0],p[0],p[1],p[2],p[3]);},5,50,4);
      getSecondFunction_Proton(g_p_par[j-1][k],f_p_par[j-1][k],p_p_par[j-1][k],k,j);
    }
  }

  //////////////
  //Final Params
  //////////////
  //////////////
  //Final Params
  //////////////
  //////////////
  //Final Params
  //////////////
  for(int j=1; j<=6; j++){
    for(int i=0; i<14; i++){
      double x = (bE_Theta[i]+bE_Theta[i+1])/2;
      for(int k = 0; k < 4; k++){      
	f1_e_cvp_FIT4[j-1][i]->SetParameter(k,f1_e_cvt_FIT4[j-1][k]->Eval(x));
      }
    }
  }

  for(int j=1; j<=6; j++){
    for(int i=0; i<14; i++){
      double x = (bE_Theta[i]+bE_Theta[i+1])/2;
      for(int k = 0; k < 4; k++){      
	f_p_frac_comb[j-1][i]->SetParameter(k,f_p_par[j-1][k]->Eval(x));
      }
    }
  }

  
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
  TCanvas * myCanvas = new TCanvas("myPage","myPage",1980,1530);
  TStyle *myStyle  = new TStyle("MyStyle","My Root Styles");
  myStyle->SetPalette("kbird",0);
  myStyle->SetTitleSize(0.07, "t");
  myStyle->SetOptStat(0);
  myStyle->cd();

  std::string open = std::string(pdfFile) + "[";
  std::string write = std::string(pdfFile);
  std::string close = std::string(pdfFile) + "]";

  myCanvas->SaveAs(open.c_str());


  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////

  myCanvas->Divide(1,1,0,0);
  myCanvas->cd(1);
  myCanvas->GetPad(1)->SetBottomMargin(0.19);
  myCanvas->GetPad(1)->SetLeftMargin(0.19);
  h_Delta_Int_Eprime_Corrected->SetLineColor(2);
  h_Delta_Int_Eprime_Corrected->Draw();
  h_Delta_Int_Eprime_beforeRad->Draw("SAME");
  f_ab1->Draw("SAME");
  f_ab2->Draw("SAME");
  myCanvas->Print(write.c_str(),"pdf");
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
  myCanvas->Print(write.c_str(),"pdf");
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
  myCanvas->Print(write.c_str(),"pdf");
  myCanvas->Clear();
  
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  double x_ab1[2] = {8,28};
  double x_ab2[2] = {24,39};

  double y_ab1[2] = {-2.5,2.5};  
  double y_ab2[2] = {-2.5,2.5};  
  double y_ab3[2] = {0,250};  
  double y_ab4[2] = {-M_PI,M_PI};  

  TGraph * r_ab1 = new TGraph(2,x_ab1,y_ab1);
  r_ab1->SetTitle("C_{a} vs. #theta #circ;#theta #circ;C_{a}");
  TGraph * r_ab2 = new TGraph(2,x_ab1,y_ab2);
  r_ab2->SetTitle("C_{b} vs. #theta #circ;#theta #circ;C_{b}");
  TGraph * r_ab3 = new TGraph(2,x_ab1,y_ab3);
  r_ab3->SetTitle("C_{c} vs. #theta #circ;#theta #circ;C_{c}");
  TGraph * r_ab4 = new TGraph(2,x_ab1,y_ab4);
  r_ab4->SetTitle("C_{d} vs. #theta #circ;#theta #circ;C_{d}");
  TGraph * r_ab5 = new TGraph(2,x_ab2,y_ab1);
  r_ab5->SetTitle("C_{a} vs. #theta #circ;#theta #circ;C_{a}");
  TGraph * r_ab6 = new TGraph(2,x_ab2,y_ab2);
  r_ab6->SetTitle("C_{b} vs. #theta #circ;#theta #circ;C_{b}");
  TGraph * r_ab7 = new TGraph(2,x_ab2,y_ab3);
  r_ab7->SetTitle("C_{c} vs. #theta #circ;#theta #circ;C_{c}");
  TGraph * r_ab8 = new TGraph(2,x_ab2,y_ab4);
  r_ab8->SetTitle("C_{d} vs. #theta #circ;#theta #circ;C_{d}");

  vector<TGraph*> r_abc_1;
  r_abc_1.push_back(r_ab1);
  r_abc_1.push_back(r_ab2);
  r_abc_1.push_back(r_ab3);
  r_abc_1.push_back(r_ab4);
  for(int i = 0; i < r_abc_1.size(); i++){
    r_abc_1[i]->SetLineColor(0);
    r_abc_1[i]->GetXaxis()->CenterTitle();
    r_abc_1[i]->GetXaxis()->SetTitleSize(0.10);
    r_abc_1[i]->GetXaxis()->SetLabelSize(0.06);
    r_abc_1[i]->GetXaxis()->SetTitleOffset(0.8);
    r_abc_1[i]->GetYaxis()->CenterTitle();
    r_abc_1[i]->GetYaxis()->SetTitleSize(0.10);
    r_abc_1[i]->GetYaxis()->SetLabelSize(0.06);
    r_abc_1[i]->GetYaxis()->SetTitleOffset(0.8);
  }
  
  
  vector<TGraph*> r_abc_2;
  r_abc_2.push_back(r_ab5);
  r_abc_2.push_back(r_ab6);
  r_abc_2.push_back(r_ab7);
  r_abc_2.push_back(r_ab8);
  for(int i = 0; i < r_abc_2.size(); i++){
    r_abc_2[i]->SetLineColor(0);
    r_abc_2[i]->GetXaxis()->CenterTitle();
    r_abc_2[i]->GetXaxis()->SetTitleSize(0.10);
    r_abc_2[i]->GetXaxis()->SetLabelSize(0.06);
    r_abc_2[i]->GetXaxis()->SetTitleOffset(0.8);
    r_abc_2[i]->GetYaxis()->CenterTitle();
    r_abc_2[i]->GetYaxis()->SetTitleSize(0.10);
    r_abc_2[i]->GetYaxis()->SetLabelSize(0.06);
    r_abc_2[i]->GetYaxis()->SetTitleOffset(0.8);
  }
  
  ///////////////////////////////////
  ///////////////////////////////////

  for(int j = 1; j <= 6; j++){
    myCanvas->Divide(3,3,0,0);
    for(int i = 0; i < 9; i++){
      double theta = (bE_Theta[i+1]+bE_Theta[i])/2;

      myCanvas->cd(i+1);
      myCanvas->GetPad(i+1)->SetBottomMargin(0.19);
      myCanvas->GetPad(i+1)->SetLeftMargin(0.19);
      h2_e_cvp_FIT0[j-1][i]->Draw("colz");//The 2D Hist
      g1_e_cvp_FIT1[j-1][i]->SetLineColor(2); 
      g1_e_cvp_FIT1[j-1][i]->Draw("SAME");//The 1D Graph of h2_e_cvp_FIT0
      f_e_frac0[j-1][i]->SetLineColor(6);  
      //f_e_frac0[j-1][i]->Draw("SAME");//Function fitted to g1_e_cvp_FIT1
      f1_e_cvp_FIT4[j-1][i]->SetLineColor(6);
      //f1_e_cvp_FIT4[j-1][i]->Draw("SAME");//Function parameters from f1_e_cvt_FIT4
      
      TGraph* g_slice = new TGraph(100);
      for (int i = 0; i < 100; ++i) {
	double y = -40 + i * (80) / (99);
	double z = f2_e_cvpvt_FIT2[j-1]->Eval(theta, y);
	g_slice->SetPoint(i, y, z);
      }      
      g_slice->SetLineColor(3);
      g_slice->Draw("SAME"); //Slice of f2_e_cvpvt_FIT2 which is a fit to the full g2_e_cvpvt_FIT1 that covers the full TH2D
      
      f1_e_cvp_FIT3[j-1][i]->SetLineColor(4);
      f1_e_cvp_FIT3[j-1][i]->Draw("SAME"); //Function fitted to g1_e_cvp_FIT1 but with a guess from f2_e_cvpvt_FIT2
      
    }    
    myCanvas->Print(write.c_str(),"pdf");
    myCanvas->Clear();

    myCanvas->Divide(2,2,0,0);
    for(int k = 0; k < 4; k++){
      myCanvas->cd(k+1);
      myCanvas->GetPad(k+1)->SetBottomMargin(0.19);
      myCanvas->GetPad(k+1)->SetLeftMargin(0.19);
      r_abc_1[k]->Draw();
      //g_e_par[j-1][k]->SetLineColor(3);
      //g_e_par[j-1][k]->Draw("SAME");//Nothing
      f1_e_cvt_FIT4[j-1][k]->SetLineColor(6);
      f1_e_cvt_FIT4[j-1][k]->Draw("SAME");//Fit to g1_e_cvt_FIT3
      f1_e_cvt_FIT2[j-1][k]->SetLineColor(3);
      f1_e_cvt_FIT2[j-1][k]->Draw("SAME");//From f2_e_cvpvt_FIT2
      g1_e_cvt_FIT3[j-1][k]->SetLineColor(4);
      g1_e_cvt_FIT3[j-1][k]->Draw("SAME");//Params from f1_e_cvp_FIT3
    }    
    myCanvas->Print(write.c_str(),"pdf");
    myCanvas->Clear();    


    myCanvas->Divide(3,3,0,0);
    for(int i = 0; i < 9; i++){
      double theta = (bE_Theta[i+1]+bE_Theta[i])/2;
      myCanvas->cd(i+1);
      myCanvas->GetPad(i+1)->SetBottomMargin(0.19);
      myCanvas->GetPad(i+1)->SetLeftMargin(0.19);
      h2_e_cvp_FIT0[j-1][i]->Draw("colz");//The 2D Hist
      g1_e_cvp_FIT1[j-1][i]->SetLineColor(2); 
      g1_e_cvp_FIT1[j-1][i]->Draw("SAME");//The 1D Graph of h2_e_cvp_FIT0
      f1_e_cvp_FIT4[j-1][i]->SetLineColor(6);
      f1_e_cvp_FIT4[j-1][i]->Draw("SAME");//Function parameters from f1_e_cvt_FIT4           
    }    
    myCanvas->Print(write.c_str(),"pdf");
    myCanvas->Clear();

  }
  
  for(int j = 1; j <= 6; j++){
    myCanvas->Divide(3,2,0,0);    
    for(int i = 7; i < 13; i++){
      myCanvas->cd(i-6);
      myCanvas->GetPad(i-6)->SetBottomMargin(0.19);
      myCanvas->GetPad(i-6)->SetLeftMargin(0.19);
      h_p_frac[j-1][i]->Draw("colz");
      g_p_frac[j-1][i]->SetLineColor(2);
      g_p_frac[j-1][i]->Draw("SAME");
      f_p_frac[j-1][i]->SetLineColor(3);
      f_p_frac[j-1][i]->Draw("SAME");
      f_p_frac_comb[j-1][i]->SetLineColor(4);
      f_p_frac_comb[j-1][i]->Draw("SAME");
    }    
    myCanvas->Print(write.c_str(),"pdf");
    myCanvas->Clear();

    myCanvas->Divide(2,2,0,0);
    for(int k = 0; k < 4; k++){
      myCanvas->cd(k+1);
      myCanvas->GetPad(k+1)->SetBottomMargin(0.19);
      myCanvas->GetPad(k+1)->SetLeftMargin(0.19);
      r_abc_2[k]->Draw();
      g_p_par[j-1][k]->SetLineColor(3);
      g_p_par[j-1][k]->Draw("SAME");
      f_p_par[j-1][k]->SetLineColor(4);
      f_p_par[j-1][k]->Draw("SAME");
    }    
    myCanvas->Print(write.c_str(),"pdf");
    myCanvas->Clear();    
  }
  
  myCanvas->Print(close.c_str(),"pdf");

  ////////////////////////////////////////////
  
  cout<<"const double params_Momentum_negparts_FD[6][4][4]{";
  for(int j=1; j<=6; j++){
    cout<<"{";
    for(int k=0; k<4; k++){
      cout<<"{";
      for(int l = 0; l < 4; l++){      
	cout<<p1_e_cvt_FIT4[j-1][k]->Parameter(l);	
	if (l < 3) cout << ", ";
      }
      cout << " }";
      if (k < 3) cout << ",";
      cout << "\n";
    }
    cout << "  }";
    if (j < 6) cout << ",";
    cout << "\n";
  }
  cout<<"};\n";

  cout<<"const double params_Momentum_posparts_FD[6][4][4]{";
  for(int j=1; j<=6; j++){
    cout<<"{";
    for(int k=0; k<4; k++){
      cout<<"{";
      for(int l = 0; l < 4; l++){      
	cout<<p_p_par[j-1][k]->Parameter(l);	
	if (l < 3) cout << ", ";
      }
      cout << " }";
      if (k < 3) cout << ",";
      cout << "\n";
    }
    cout << "  }";
    if (j < 6) cout << ",";
    cout << "\n";
  }
  cout<<"};\n";

  /*
  cout<<"const double f_e_guess[6][14][4]={";
  for(int j=1; j<=6; j++){
    cout<<"{";
    for(int i=0; i<14; i++){
      double x = (bE_Theta[i]+bE_Theta[i+1])/2;
      cout<<"{";
      for(int k = 0; k < 4; k++){      
	cout<<f1_e_cvt_FIT4[j-1][k]->Eval(x);
	if (k < 3) cout << ", ";
      }
      cout << " }";
      if (i < 13) cout << ",";
      cout << "\n";
    }
    cout << "  }";
    if (j < 6) cout << ",";
    cout << "\n";
  }
  cout<<"};\n";

  cout<<"const double f_p_guess[6][14][4]={";
  for(int j=1; j<=6; j++){
    cout<<"{";
    for(int i=0; i<14; i++){
      double x = (bE_Theta[i]+bE_Theta[i+1])/2;
      cout<<"{";
      for(int k = 0; k < 4; k++){      
	cout<<f_p_par[j-1][k]->Eval(x);
	if (k < 3) cout << ", ";
      }
      cout << " }";
      if (i < 13) cout << ",";
      cout << "\n";
    }
    cout << "  }";
    if (j < 6) cout << ",";
    cout << "\n";
  }
  cout<<"};\n";
  */
  return 0;

}


