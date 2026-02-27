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
#include "Corrections.h"


// For Fitting
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/UnBinData.h"
#include "Fit/Chi2FCN.h"
#include "Fit/FitResult.h"
#include "Fit/DataOptions.h"
#include "Fit/FitConfig.h"

// For defining the functions
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include <vector>

using namespace std;
using namespace clas12;

const double c = 29.9792458;

auto db=TDatabasePDG::Instance();
const double beam_E = 5.98636;
const double mass_n = db->GetParticle(2112)->Mass();
const double mass_p = db->GetParticle(2212)->Mass();
const double mass_pi = db->GetParticle(-211)->Mass();
//const double mD = 1.8756;
const double me = 0.000511;
const double mU = 0.9314941024;
const double m_2H = 2.01410178 * mU - me;
const double m_3H = 3.01604928199 * mU - me;
const double m_3He = 3.0160293 * mU - 2*me;
const double m_4He = 4.00260325415 * mU - 2*me;


void Usage()
{
  std::cerr << "Usage: ./code isMC A outputfile.root outputfile.pdf inputfiles.hipo \n\n\n";
}

int main(int argc, char ** argv)
{
  
  if(argc < 5)
    {
      Usage();
      return -1;
    }
  
  int isMC = atoi(argv[1]);
  int nucleus_A = atoi(argv[2]);
  TString outFile = argv[3];
  char * pdfFile = argv[4];

  cout<<"Ouput file "<< outFile <<endl;
  cout<<"Ouput PDF file "<< pdfFile <<endl;


  clas12ana clasAna;
  clasAna.printParams();
    
  //clas12root::HipoChain chain;
  /*
  */
  clas12root::HipoChain chain;
  for(int k = 5; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12=chain.GetC12Reader();

  
  int counter = 0;
  int cutcounter = 0;

  auto &c12=chain.C12ref();
  
  double mN = db->GetParticle(2212)->Mass();
  //some particles
  TLorentzVector beam(0,0,beam_E,beam_E);
  TLorentzVector nucleus_ptr(0,0,0,m_4He);
  TLorentzVector deut_recoil_ptr(0,0,0,m_2H);

  TLorentzVector deut_ptr(0,0,0,m_2H);
  TLorentzVector hel3_ptr(0,0,0,m_3He);
  TLorentzVector hel4_ptr(0,0,0,m_4He);

  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector lead_ptr(0,0,0,db->GetParticle(2212)->Mass());
  int Z=2;
  int N=2;
  if(isMC){
    Z=nucleus_A/2;
    N=nucleus_A/2;
  }
  reweighter newWeight(beam_E,Z,N,kelly,"AV18");
  char temp_name[100];
  char temp_title[100];


  vector<TH1*> hist_list;

  TH1D * h_Q2 = new TH1D("Q2","Q2;Q2;Counts",100,1.0,5.0);
  hist_list.push_back(h_Q2);
  TH1D * h_xB = new TH1D("xB","xB;xB;Counts",100,1.0,2.5);
  hist_list.push_back(h_xB);
  TH1D * h_pL = new TH1D("pL","p_{Lead};p_{Lead};Counts",100,0.9,3.5);
  hist_list.push_back(h_pL);
  TH1D * h_tL = new TH1D("tL","#theta_{Lead};#theta_{Lead};Counts",100,0,90);
  hist_list.push_back(h_tL);
  TH1D * h_pM = new TH1D("pM","p_{miss};p_{miss};Counts",100,0.5,1.5);
  hist_list.push_back(h_pM);
  TH1D * h_tM = new TH1D("tM","#theta_{miss};#theta;Counts",100,0,180);
  hist_list.push_back(h_tM);
  TH1D * h_mM = new TH1D("mM","m_{miss};m_{miss};Counts",100,0.0,2.0);
  hist_list.push_back(h_mM);
  TH1D * h_tpq = new TH1D("tpq","#theta_{pq};#theta_{pq};Counts",100,0,40);
  hist_list.push_back(h_tpq);
  TH1D * h_poq = new TH1D("poq","p/q;p/q;Counts",100,0.5,1.0);
  hist_list.push_back(h_poq);
  TH2D * h_poq_tpq = new TH2D("poq_tpq","#theta_{pq} vs. p/q;#theta_{pq};p/q;Counts",100,0.5,1.0,100,0,40);
  hist_list.push_back(h_poq);


  
  TH1D * h_mM3 = new TH1D("mM3","mM3;mM3;Counts",100,0.0,3.0);
  hist_list.push_back(h_mM3);

  TH2D * h_p_mass_FD = new TH2D("p_mass_FD","Mass vs. p FD",100,0.0,3.0,100,0.0,4.0);
  hist_list.push_back(h_p_mass_FD);
  TH2D * h_p_mass_CD = new TH2D("p_mass_CD","Mass vs. p CD",100,0.0,3.0,100,0.0,4.0);
  hist_list.push_back(h_p_mass_CD);
  
  TH2D * h_p_beta_FD = new TH2D("p_beta_FD","#beta vs. p FD",100,0.0,3.0,100,0.2,1.2);
  hist_list.push_back(h_p_beta_FD);
  TH2D * h_p_beta_CD = new TH2D("p_beta_CD","#beta vs. p CD",100,0.0,3.0,100,0.2,1.2);
  hist_list.push_back(h_p_beta_CD);


  TH2D * h_vtz_e_p = new TH2D("vtz_e_p","Electron Vertex vs. Proton Vertex;Vertex_{e} [cm];Vertex_{p} [cm]",100,-10,10,100,-10,10);
  hist_list.push_back(h_vtz_e_p);
  TH2D * h_vtz_e_d = new TH2D("vtz_e_d","Electron Vertex vs. Deuteron Vertex;Vertex_{e} [cm];Vertex_{d} [cm]",100,-10,10,100,-10,10);
  hist_list.push_back(h_vtz_e_d);
  TH2D * h_vtz_p_d = new TH2D("vtz_p_d","Proton Vertex vs. Deuteron Vertex;Vertex_{p} [cm];Vertex_{d} [cm]",100,-10,10,100,-10,10);
  hist_list.push_back(h_vtz_p_d);

  TH1D * h_diffvtz_pd_bc = new TH1D("diffvtz_pd_bc","#Delta Vtz_{pd};#Delta Vtz_{pd} [cm];Counts",100,-5,5);
  hist_list.push_back(h_diffvtz_pd_bc);
  TH1D * h_diffvtz_ed_bc = new TH1D("diffvtz_ed_bc","#Delta Vtz_{ed};#Delta Vtz_{ed} [cm];Counts",100,-5,5);
  hist_list.push_back(h_diffvtz_ed_bc);

  TH2D * h_diffvtz_ed_pd_bc = new TH2D("diffvtz_ed_pd_bc","#Delta Vtz_{pd} vs. #Delta Vtz_{ed};#Delta Vtz_{ed} [cm];#Delta Vtz_{pd} [cm]",100,-5,5,100,-5,5);
  hist_list.push_back(h_diffvtz_ed_pd_bc);
  TH2D * h_diffvtz_ed_pd_ac = new TH2D("diffvtz_ed_pd_ac","#Delta Vtz_{pd} vs. #Delta Vtz_{ed};#Delta Vtz_{ed} [cm];#Delta Vtz_{pd} [cm]",100,-5,5,100,-5,5);
  hist_list.push_back(h_diffvtz_ed_pd_ac);

  
  TH1D * h_mM4 = new TH1D("mM4","mM4;mM4;Counts",100,0.0,3.0);
  hist_list.push_back(h_mM4);
  
  TH1D * h_t_pd = new TH1D("t_pd","#theta_{miss,d};#theta_{miss,d};Counts",180,0,180);
  hist_list.push_back(h_t_pd);
  TH1D * h_ct_pd = new TH1D("ct_pd","cos(#theta_{miss,d});cos(#theta_{miss,d});Counts",50,0.4,1);
  hist_list.push_back(h_ct_pd);
  TH1D * h_r_pd = new TH1D("r_pd","p_{miss}-p_{d};p_{miss}-p_{d};Counts",100,-1.0,1.0);
  hist_list.push_back(h_r_pd);
  TH1D * h_t_ld = new TH1D("t_ld","#theta_{lead,d};#theta_{lead,d};Counts",180,0,180);
  hist_list.push_back(h_t_ld);

  TH2D * h_r_t_pd = new TH2D("r_t_pd",";p_{miss}-p_{d};#theta_{miss,d};Counts",100,-1.0,1.0,180,0,180);
  hist_list.push_back(h_r_t_pd);
  TH2D * h_r_ct_pd = new TH2D("r_ct_pd",";p_{miss}-p_{d};cos(#theta_{miss,d});Counts",100,-1.0,1.0,50,0.4,1.0);
  hist_list.push_back(h_r_t_pd);
  
  TH2D * h_tt_pd = new TH2D("tt_pd",";#theta_{miss};#theta_{d}",100,0,180,100,0,180);
  hist_list.push_back(h_tt_pd);
  TH2D * h_pp_pd = new TH2D("pp_pd",";p_{miss};p_{d}",100,0.4,1.5,100,0.4,1.5);
  hist_list.push_back(h_pp_pd);
  TH2D * h_xB_pmiss = new TH2D("xB_pmiss",";x_{B};p_{miss}",100,1.0,2.5,100,0.0,1.5);
  hist_list.push_back(h_xB_pmiss);

  TH1D * h_pcmx = new TH1D("pcmx","p_{cmx};p_{cmx};p_{cmx}",100,-1.0,1.0);
  hist_list.push_back(h_pcmx);
  TH1D * h_pcmy = new TH1D("pcmy","p_{cmy};p_{cmy};p_{cmy}",100,-1.0,1.0);
  hist_list.push_back(h_pcmy);
  TH1D * h_pcmz = new TH1D("pcmz","p_{cmz};p_{cmz};p_{cmz}",100,-1.0,1.0);
  hist_list.push_back(h_pcmz);
  
  
  TH1D * h_xB_ep = new TH1D("xB_ep","xB;xB;Counts",50,1.0,3.0);
  hist_list.push_back(h_xB_ep);
  TH1D * h_xB_epd = new TH1D("xB_epd","xB;xB;Counts",50,1.0,3.0);
  hist_list.push_back(h_xB_epd);

  TH1D * h_pM_ep = new TH1D("pM_ep","pM;pM;Counts",50,0.5,0.85);
  hist_list.push_back(h_pM_ep);
  TH1D * h_pM_epd = new TH1D("pM_epd","pM;pM;Counts",50,0.5,0.85);
  hist_list.push_back(h_pM_epd);

  TH1D * h_alphalead = new TH1D("alphalead","alphalead",100,0.0,2.0);
  hist_list.push_back(h_alphalead);
  TH1D * h_alphamiss = new TH1D("alphamiss","alphamiss",100,0.0,2.0);
  hist_list.push_back(h_alphamiss);
  TH1D * h_alphadeut = new TH1D("alphadeut","alphadeut",100,0.0,2.0);
  hist_list.push_back(h_alphadeut);
  TH1D * h_alphacm = new TH1D("alphacm","alphacm",100,-4.0,4.0);
  hist_list.push_back(h_alphacm);

  TH2D * h_pMissN_mMissN = new TH2D("pMissN_mMissN","m_{Miss} vs. p_{Miss} (e,e'pd);p_{Miss};m_{Miss}",100,0.0,1.0,100,0,1.5);
  hist_list.push_back(h_pMissN_mMissN);
  
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Sumw2();
    hist_list[i]->GetXaxis()->CenterTitle();
    hist_list[i]->GetYaxis()->CenterTitle();
  }

  double num = 0;
  double den = 0;
  
  while(chain.Next() && counter <100000000000)
    {

      double weight = 1;
      if(isMC){
	double original_weight = c12->mcevent()->getWeight(); //used if MC events have a weight
	weight = original_weight * newWeight.get_weight_ep(c12->mcparts());
      }

      //Display completed  
      counter++;
      if((counter%100000) == 0){
	cerr << "\n" <<counter/100000 <<" hundred thousand completed";
      }    
      if((counter%10000) == 0){
	cerr << ".";
      }    

      clasAna.Run(c12);
      auto electrons = clasAna.getByPid(11);
      auto protons = clasAna.getByPid(2212);
      auto pims = clasAna.getByPid(-211);
      auto pips = clasAna.getByPid(211);
      auto deuterons = clasAna.getByPid(45);
      //auto deuterons = clasAna.getByPid(45);
      auto particles = c12->getDetParticles(); 

      if(electrons.size() != 1){continue;}

      GetLorentzVector_Corrected(el,electrons[0],isMC);
      TLorentzVector q = beam - el;
      double Q2        = -q.M2();
      double omega = q.E();
      double xB       = Q2/(2 * mass_p * (beam.E() - el.E()) );
      double WSq = (mN*mN) - Q2 + (2*omega*mN);
      double W = sqrt(WSq);
      double vtz_e = electrons[0]->par()->getVz();
	  
      int sector_e = electrons[0]->getSector()-1;
      double mom_e = el.P();
      double theta_e = el.Theta() * 180 / M_PI;
      double phi_e = el.Vect().Phi() * 180 / M_PI;
      double shift_e = 7.5;


      
      clasAna.getLeadRecoilSRC(beam,deut_ptr,el);
      auto lead    = clasAna.getLeadSRC();
      auto recoil  = clasAna.getRecoilSRC();
      
      if(lead.size()!=1){continue;}
      GetLorentzVector_Corrected(lead_ptr,lead[0],isMC);

      TLorentzVector miss = q + deut_ptr - lead_ptr;
      TLorentzVector miss_fromHelium3_ofDeuteron = q + hel3_ptr - lead_ptr;
      
      double mmiss2 = miss.M2();
      double mmiss= sqrt(mmiss2);
      double alphamiss = (miss.E() - miss.Vect().Dot(q.Vect().Unit()))/mN;
      TVector3 miss_neg = -miss.Vect();
      TVector3 vz = miss_neg.Unit();
      TVector3 vy = miss_neg.Cross(q.Vect()).Unit();
      TVector3 vx = vz.Cross(vy).Unit();
      double mom_miss = miss.P();
      double theta_miss = miss.Theta() * 180 / M_PI;
      
      double E0miss = sqrt(mom_miss*mom_miss + mN*mN)-mN;
      double beta_lead = lead[0]->par()->getBeta();
      double mom_lead = lead_ptr.P();
      double momT_lead = lead_ptr.Vect().Perp();
      double theta_lead = lead_ptr.Theta() * 180 / M_PI;
      double phi_lead = lead_ptr.Phi() * 180 / M_PI;
      double vtz_lead = lead[0]->par()->getVz();
      double EP = lead_ptr.E();
      double EB = omega + nucleus_ptr.M() - EP;
      double TB = EB - sqrt(EB*EB - mom_miss*mom_miss);
      double TP = EP - sqrt(EP*EP - mom_lead*mom_lead);
      double E1miss = omega - TP - TB;
      double thetamissq = miss_neg.Angle(q.Vect())*180/M_PI;
      double thetapq = lead_ptr.Vect().Angle(q.Vect())*180/M_PI;
      double poq = mom_lead/q.P();

      TLorentzVector miss_LC = lead_ptr - q;

      TVector3 u = q.Vect().Unit();
      double pmm = miss_LC.E() - miss_LC.Vect().Dot(u);
      double pmp = miss_LC.Vect().Perp(u);
      double pmiss = miss.P();
      double kmiss = sqrt(mN*mN*((pmp*pmp+mN*mN)/(pmm*(2*mN-pmm))) - mN*mN);
      
      double phidiff = (q.Vect().Phi()*180/M_PI)-phi_lead;
      if(phidiff<-180){phidiff+=360;}
      if(phidiff>180){phidiff-=360;}
            
      bool rec = false;
      if(recoil.size()==1){rec = true;}

      //Other cuts we could do
      //if(theta_lead>37){continue;}
      //if(xB<1.0){continue;}
      //if(mmiss<0.65){continue;}
      //if(mmiss>1.10){continue;}
      //if(kmiss<0.3){continue;}            
      //if(mom_lead<1.0){continue;}
      //if(miss.P()<0.5){continue;}
      //if(miss.P()>1.0){continue;}

      if(Q2<1.5){continue;}
      if(theta_miss<40){continue;}
      if(theta_miss>120){continue;}
      if(lead[0]->getRegion()!=CD){continue;}

      h_Q2->Fill(Q2,weight);
      h_xB->Fill(xB,weight);
      h_pL->Fill(mom_lead,weight);
      h_tL->Fill(theta_lead,weight);
      h_pM->Fill(pmiss,weight);
      h_tM->Fill(theta_miss,weight);
      h_mM->Fill(mmiss,weight);
      h_tpq->Fill(thetapq,weight);
      h_poq->Fill(poq,weight);
      h_poq_tpq->Fill(poq,thetapq,weight);
      
      h_xB_ep->Fill(xB,weight);
      h_pM_ep->Fill(miss.P(),weight);


      if(deuterons.size()!=1){continue;}
      GetLorentzVector_Corrected(deut_recoil_ptr,deuterons[0],isMC);
      //GetLorentzVector_ReconVector(deut_recoil_ptr,deuterons[0]);

      double mom = deut_recoil_ptr.P();
      double theta = deut_recoil_ptr.Theta() * 180 / M_PI;
      double phi = deut_recoil_ptr.Phi() * 180 / M_PI;
      double beta = deuterons[0]->par()->getBeta();
      double gamma = 1/sqrt(1-(beta*beta));
      double mass_pid = mom/(beta*gamma);
      
      TVector3 dr = deut_recoil_ptr.Vect();
      double theta_missd = miss.Vect().Angle(dr)*180/M_PI;
      double theta_ld = lead_ptr.Vect().Angle(dr)*180/M_PI;
      double vtz_d = deuterons[0]->par()->getVz();
      
      TVector3 v_rec = deut_recoil_ptr.Vect();
      TVector3 v_cm   = miss_neg + v_rec;
      
      TLorentzVector miss_fromHelium4_ofHelium3 = q + hel4_ptr - lead_ptr - deut_recoil_ptr;
      
      //if(mass_pid>2.5){continue;}
	
      if(deuterons[0]->getRegion()==FD){
	h_p_beta_FD->Fill(mom,beta,weight);
	h_p_mass_FD->Fill(mom,mass_pid,weight);
      }
      else if(deuterons[0]->getRegion()==CD){
	h_p_beta_CD->Fill(mom,beta,weight);
	h_p_mass_CD->Fill(mom,mass_pid,weight);
      }
      if(deuterons[0]->getRegion()==FD){continue;}

      
      h_vtz_e_p->Fill(vtz_e,vtz_lead,weight);
      h_vtz_e_d->Fill(vtz_e,vtz_d,weight);
      h_vtz_p_d->Fill(vtz_lead,vtz_d,weight);
      
      h_diffvtz_pd_bc->Fill(vtz_lead-vtz_d,weight);
      h_diffvtz_ed_bc->Fill(vtz_e-vtz_d,weight);
      h_diffvtz_ed_pd_bc->Fill(vtz_e-vtz_d,vtz_lead-vtz_d,weight);
      if(fabs(vtz_e-vtz_d)>1.5){continue;}
      if(fabs(vtz_lead-vtz_d)>1.5){continue;}	
      h_diffvtz_ed_pd_ac->Fill(vtz_e-vtz_d,vtz_lead-vtz_d,weight);
	
	
      h_mM3->Fill(miss_fromHelium3_ofDeuteron.M(),weight);
      //if(miss_fromHelium3_ofDeuteron.M()<(m_2H-0.2)){continue;}
      
      h_mM4->Fill(miss_fromHelium4_ofHelium3.M(),weight);
      //if(miss_fromHelium4_ofHelium3.M()<(mass_n-0.2)){continue;}
      
      h_t_pd->Fill(theta_missd,weight);
      h_ct_pd->Fill(cos(theta_missd*M_PI/180),weight);
      h_r_pd->Fill(pmiss-mom,weight);	  
      h_t_ld->Fill(theta_ld,weight);
      
      h_r_t_pd->Fill(pmiss-mom,theta_missd,weight);
      h_r_ct_pd->Fill(pmiss-mom,cos(theta_missd*M_PI/180),weight);
      
      
      h_tt_pd->Fill(theta_miss,theta,weight);
      h_pp_pd->Fill(pmiss,mom,weight);
      h_xB_pmiss->Fill(xB,pmiss,weight);
      
      h_xB_epd->Fill(xB,weight);
      h_pM_epd->Fill(miss.P(),weight);
      
      h_pcmx->Fill(v_cm.Dot(vx),weight);
      h_pcmy->Fill(v_cm.Dot(vy),weight);
      h_pcmz->Fill(v_cm.Dot(vz),weight);
      
      double alpha_lead = (lead_ptr.E()-lead_ptr.Vect().Dot(u))/lead_ptr.M();
      double alpha_deut = (deut_recoil_ptr.E()-deut_recoil_ptr.Vect().Dot(u))/deut_recoil_ptr.M();
      //double alpha_lead = (lead_ptr.E()-lead_ptr.Vect().Dot(u))/lead_ptr.M();
      h_alphalead->Fill(alpha_lead,weight);
      h_alphadeut->Fill(alpha_deut,weight);
      //h_alphacm->Fill(alphacm,weight);
      
      h_pMissN_mMissN->Fill(v_cm.Mag(),miss_fromHelium4_ofHelium3.M(),weight);
      
      

    }
  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////

  TFile *f = new TFile(outFile,"RECREATE");
  f->cd();
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Write();
  }

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

  
  myCanvas->Divide(3,3);
  myCanvas->cd(1);    
  h_Q2->Draw();
  myCanvas->cd(2);    
  myCanvas->cd(2)->SetLogy();
  h_xB->Draw();
  myCanvas->cd(3);    
  h_pL->Draw();
  myCanvas->cd(4);    
  h_tL->Draw();
  myCanvas->cd(5);    
  h_pM->Draw();
  myCanvas->cd(6);    
  h_tM->Draw();
  myCanvas->cd(7);    
  h_mM->Draw();
  myCanvas->cd(8);    
  h_tpq->Draw();
  myCanvas->cd(9);    
  h_poq->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_poq_tpq->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(2,2);
  myCanvas->cd(1);    
  myCanvas->cd(1)->SetLogz();
  h_p_mass_FD->Draw("SAME");
  myCanvas->cd(2);    
  myCanvas->cd(2)->SetLogz();
  h_p_mass_CD->Draw("SAME");
  myCanvas->cd(3);    
  h_p_beta_FD->Draw("SAME");
  myCanvas->cd(4);    
  h_p_beta_CD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(2,2);
  myCanvas->cd(1);    
  h_vtz_e_p->Draw("SAME");
  myCanvas->cd(2);    
  h_vtz_e_d->Draw("SAME");
  myCanvas->cd(3);    
  h_vtz_p_d->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(2,2);
  myCanvas->cd(1);    
  h_diffvtz_ed_pd_bc->Draw("SAME");
  myCanvas->cd(2);    
  h_diffvtz_ed_pd_ac->Draw("SAME");
  myCanvas->cd(3);    
  h_diffvtz_pd_bc->Draw("SAME");
  myCanvas->cd(4);    
  h_diffvtz_ed_bc->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_mM3->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_mM4->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_t_pd->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_ct_pd->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_r_pd->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_t_ld->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_r_t_pd->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_r_ct_pd->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_tt_pd->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pp_pd->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_xB_pmiss->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  myCanvas->cd(1)->SetLogy();
  h_xB_ep->Rebin(4);
  h_xB_epd->Rebin(4);
  h_xB_ep->Draw();
  h_xB_epd->SetLineColor(2);
  h_xB_epd->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  myCanvas->cd(1)->SetLogy();
  h_xB_epd->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pcmx->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pcmy->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pcmz->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_alphalead->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_alphadeut->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pMissN_mMissN->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  double x_ab1[2] = {1.2,2.6};
  double y_ab1[2] = {0,5.0};  
  TGraph * r_ab1 = new TGraph(2,x_ab1,y_ab1);
  r_ab1->SetTitle("epd/ep;xB;epd/ep [%]");
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  r_ab1->SetLineColor(0);
  r_ab1->Draw();
  h_xB_epd->Divide(h_xB_ep);
  h_xB_epd->Scale(100.0);
  h_xB_epd->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pM_ep->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pM_epd->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  double x_ab2[2] = {0.6,1.3};
  double y_ab2[2] = {0,5.0};  
  TGraph * r_ab2 = new TGraph(2,x_ab2,y_ab2);
  r_ab2->SetTitle("epd/ep;pM;epd/ep [%]");
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  r_ab2->SetLineColor(0);
  //r_ab2->Draw();
  h_pM_epd->Divide(h_pM_ep);
  h_pM_epd->Scale(100.0);
  h_pM_epd->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  f->Close();

  return 0;
}
