#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>
#include <typeinfo>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TDatabasePDG.h>
#include "HipoChain.h"
#include "HipoChainWriter.h"
#include "clas12ana.h"
#include "Functions.h"

using namespace std;
using namespace clas12;

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());

}

void Usage()
{
  std::cerr << "Usage: ./code outputfile.hipo inputfiles.hipo  \n\n\n";

}


int main(int argc, char ** argv)
{

  if(argc < 3)
    {
      Usage();
      return -1;
    }

  //This is the clas12ana class that helps us
  //cut on detector level and SRC variables.
  clas12ana clasAna;
  clasAna.printParams();

  //Take the output file name
  char * outName = argv[1];
  cout<<"Ouput file "<< outName <<endl;

  //make c12writer for the output hipo file
  clas12root::HipoChainWriter chain(outName);
  //Now add the input files to the chain
  for(int k = 2; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  //Some necessary tags for the chain
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  //And this is the object we use to get event
  //details.
  auto config_c12=chain.GetC12Reader();
  auto &c12=chain.C12ref();

  ////////////////////////////////////////////////

  //Define some useful variables
  int counter = 0;
  auto db=TDatabasePDG::Instance();
  double mass_p = db->GetParticle(2212)->Mass();
  double mD = 1.8756;
  double beam_E = 5.98636;
  TLorentzVector beam(0,0,beam_E,beam_E);
  TLorentzVector deut_ptr(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  
  //Now loop over all events in the input files.
  while(chain.Next())
    {
      //Display number of completed events
      counter++;
      if((counter%1000000) == 0){
	cerr << "\n" <<counter/1000000 <<" million completed";
      }    
      if((counter%100000) == 0){
	cerr << ".";
      }    

      //Here is where you run the clas12ana class.
      //It does the pid, fiducial, and vertex cuts
      //on all particles and returns particles
      //by PID number.
      clasAna.Run(c12);
      auto electrons = clasAna.getByPid(11);
      auto protons = clasAna.getByPid(2212);
      //auto pip = clasAna.getByPid(211);
      //auto pim = clasAna.getByPid(-211);

      //In this skim we only keep events with
      //exactly 1 proton.
      if(electrons.size() == 1)
	{
          SetLorentzVector(el,electrons[0]);
	  TLorentzVector q = beam - el;
          double Q2 = -q.M2();
          double xB = Q2/(2 * mass_p * (beam.E() - el.E()));
	  if(xB<1.5){continue;}
	  if(Q2<1){continue;}
	  clasAna.getLeadRecoilSRC(beam,deut_ptr,el);
	  auto lead    = clasAna.getLeadSRC();
	  auto recoil  = clasAna.getRecoilSRC();
	  //if(lead.size()!=1){continue;}
	  //Finally you write the event you want to keep
	  chain.WriteEvent();	  
	}
    }
  return 0;
}

