#include "TLorentzVector.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TLegend.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <set>
#include "studyPIDs.h"
#include "studyXQ2Rec.h"
#include "TGraph.h"
#include "TFile.h"

#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h> 

using namespace std;

int main(int argc, char** argv)
{

  srand(time(NULL));
  if(argc<4)
    {
      cout <<"filename and electron+hadron beam energy required " <<endl;
      exit(0);
    }
  Float_t x, y, Q2, W, Nu;

  const int maxParticles=200;
  ULong64_t evt_prt_count;
  Long64_t pdg[maxParticles];
  Bool_t smear_has_e[maxParticles];
  Bool_t smear_has_p[maxParticles];
  Double_t smear_orig_tot_e[maxParticles];
  Double_t tot_e[maxParticles];
  Double_t dir_x[maxParticles];
  Double_t dir_y[maxParticles];
  Double_t dir_z[maxParticles];
  Double_t px[maxParticles];
  Double_t py[maxParticles];
  Double_t pz[maxParticles];
  
  Double_t smear_orig_px[maxParticles];
  Double_t smear_orig_py[maxParticles];
  Double_t smear_orig_pz[maxParticles];
  
  TFile _file0(argv[1]);
  
  int beamEnergyI = atoi(argv[2]);
  int hadronBeamEnergyI = atoi(argv[3]);
  double beamEnergy=0.0+beamEnergyI;
  double hadronBeamEnergy=0.0+hadronBeamEnergyI;
  
  cout <<"using beam energy of " << beamEnergy <<" and hadron beam energy of " << hadronBeamEnergy<<endl;
  TTree* myTree=(TTree*)_file0.Get("events/tree");
  myTree->SetBranchAddress("evt_prt_count",&evt_prt_count); //I
  myTree->SetBranchAddress("pdg",pdg); //Long64_t
  myTree->SetBranchAddress("smear_has_e",smear_has_e); //O
  myTree->SetBranchAddress("smear_has_p",smear_has_p); //O
  
  myTree->SetBranchAddress("smear_orig_tot_e",smear_orig_tot_e); //D
  myTree->SetBranchAddress("tot_e",tot_e); //D
  myTree->SetBranchAddress("dir_x",dir_x); //D
  myTree->SetBranchAddress("dir_y",dir_y);
  myTree->SetBranchAddress("dir_z",dir_z);
  myTree->SetBranchAddress("px",px);
  myTree->SetBranchAddress("py",py);
  myTree->SetBranchAddress("pz",pz);

  myTree->SetBranchAddress("smear_orig_px",smear_orig_px);
  myTree->SetBranchAddress("smear_orig_py",smear_orig_py);
  myTree->SetBranchAddress("smear_orig_pz",smear_orig_pz);

  Int_t nevents=myTree->GetEntries();
  cout <<"Hello, World!" <<endl;

    
  for(long i=0;i<nevents;i++)
    {
      
#ifdef MAX_EVENTS
      if(i>MAX_EVENTS)
	break;
#endif
      if(!(i%10000))
	cout <<"processing event nr " << i << " of " << nevents << "(" << 100*i/(float)nevents<< "% )"<<endl;
      myTree->GetEntry(i);
      //      if(Z[0]>0.2 && Z[1]>0.1)

      if(evt_prt_count>maxParticles)
	{
	  cout <<"too many particles" <<endl;
	  continue;
	}


      
      //      cout <<"looking at " << evt_prt_count<< " particles in event " <<endl;
      //for all particles in the event



      ///find scattered electron
      double maxE=-1;
      int eIndex=-1;
      double ePx,ePy,ePz,eE;
      double ePxOrig,ePyOrig,ePzOrig,eEOrig;
      
      for(int pi=0;pi<evt_prt_count;pi++)
	{
	  if(pdg[pi]==11 && tot_e[pi]>maxE)
	    {
	      maxE=tot_e[pi];
	      eIndex=pi;
	    }
	}
      if(eIndex>-1)
	{
	  //	  cout <<"found scattered electron at index " << maxEIndex <<" maxE: "<< maxE <<endl;
	  //seems like the scattered electron is always at index 0
	}
      if(eIndex<0)
	continue;
      if(!smear_has_e[eIndex] || !smear_has_p[eIndex])
	{
	  //	  not in the detector
	  continue;
	}
      ePx=px[eIndex];
      ePy=py[eIndex];
      ePz=pz[eIndex];
      eE=tot_e[eIndex];
      ePxOrig=smear_orig_px[eIndex];
      ePyOrig=smear_orig_py[eIndex];
      ePzOrig=smear_orig_pz[eIndex];
      eEOrig=smear_orig_tot_e[eIndex];
      cout <<"original e px: " << ePxOrig <<" py: " << ePyOrig<< " pz: "<< ePzOrig<<endl;
      
      //calculate Q2 just from original and scattered electron
      Kins kinSmeared=getKinsFromScatElectron(beamEnergy,hadronBeamEnergy,ePx,ePy,ePz,eE);
      Kins kinOrig=getKinsFromScatElectron(beamEnergy,hadronBeamEnergy,ePxOrig,ePyOrig,ePzOrig,eEOrig);

      float s=63; //sqrtS for 10x100;
      //      float s=45; //sqrtS for 5x100;
      //      float s=141; //sqrtS for 18x275;
      //      float s=28.6; //sqrtS for 5x41;
            Kins kinJB=getKinsJB(evt_prt_count,px,py,pz,tot_e,pdg,smear_has_e, smear_has_p,beamEnergy,s);
      
      
      for(int pi=0;pi<evt_prt_count;pi++)
	{
	  if(pdg[pi]==22)
	    {
	      //no px, py, pz assigned for photon
	      //	      cout <<"found photon wiht e: " << tot_e[pi] << " px: " << px[pi] << " py: " << py[pi] << " pz: " << pz[pi] <<endl;
	    }
	  //look for electron

	    {
	      //	      cout <<"found photon " << endl;
	      //	      if(smear_has_e[pi])
		{
		  TLorentzVector lvphoton(dir_x[pi]*tot_e[pi],dir_y[pi]*tot_e[pi],dir_z[pi]*tot_e[pi],tot_e[pi]);
		  if(lvphoton.Eta()<-4.5|| lvphoton.Eta()>4.5 || tot_e[pi]<0.2)
		    {
		      continue;
		    }

		  //		  cout <<" found photon number " << photonE.size() <<endl;
		}
	      
	    }
	}
    }

    
  
}

template<class T> T* allocateArray(int dim1)
{
  T* ret=new T[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=0;
    }
}


template<class T> T** allocateArray(int dim1, int dim2)
{
  T** ret=new T*[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=new T[dim2];
    }

  for(int i=0;i<dim1;i++)
    {
      for(int j=0;j<dim2;j++)
	{
	  ret[i][j]=0;
	  
	}
    }
  return ret;
};

template<class T> T*** allocateArray(int dim1, int dim2, int dim3)
{
  T*** ret=new T**[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocateArray<T>(dim2,dim3);
    }
  return ret;
};

template<class T> T**** allocateArray(int dim1, int dim2, int dim3, int dim4)
{
  T**** ret=new T***[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocateArray<T>(dim2,dim3,dim4);
    }
  return ret;
};


template<class T> T***** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5)
{
 T***** ret=new T****[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocateArray<T>(dim2,dim3,dim4,dim5);
    }
  return ret;
};


template<class T> T****** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6)
{
 T****** ret=new T*****[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocateArray<T>(dim2,dim3,dim4,dim5,dim6);
    }
  return ret;
};

template<class T> T******* allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6, int dim7)
{
 T******* ret=new T******[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocateArray<T>(dim2,dim3,dim4,dim5,dim6,dim7);
    }
  return ret;
};
