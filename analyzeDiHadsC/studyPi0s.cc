//0x33 piPlus1_piPlus2
//0x34 piPlus_piMinus
//0x37 piPlus_KMinus
//0x44 piMinus1_piMinus2
//0x63 KPlus_piPlus
//0x64 KPlus_piMinus
//0x66 KPlus1_KPlus2
//0x67 KPlus_KMinus
//0x74 KMinus_piMinus
//0x77 KMinus1_KMinus2

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
#include "studyPi0s.h"
#include "TGraph.h"
#include "TFile.h"

#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h> 

using namespace std;

int main(int argc, char** argv)
{

  srand(time(NULL));
  if(argc<2)
    {
      cout <<"filename required " <<endl;
      exit(0);
    }

  vector<float> binningEta;
  vector<float> binningPt;

  binningEta.push_back(-2);
  binningEta.push_back(-1);
  binningEta.push_back(-0);
  binningEta.push_back(1);
  binningEta.push_back(2);
  binningEta.push_back(1100);

  binningPt.push_back(0.6);
  binningPt.push_back(1.5);
  binningPt.push_back(2.5);
  binningPt.push_back(100.0);

  
  Float_t x, y, Q2, W, Nu;

  const int maxParticles=200;
  ULong64_t evt_prt_count;
  Long64_t pdg[maxParticles];
  Bool_t smear_has_e[maxParticles];
  Double_t smear_orig_tot_e[maxParticles];
  Double_t tot_e[maxParticles];
  Double_t dir_x[maxParticles];
  Double_t dir_y[maxParticles];
  Double_t dir_z[maxParticles];

  TFile _file0(argv[1]);
  TTree* myTree=(TTree*)_file0.Get("events/tree");
  myTree->SetBranchAddress("evt_prt_count",&evt_prt_count); //I
  myTree->SetBranchAddress("pdg",pdg); //Long64_t
  myTree->SetBranchAddress("smear_has_e",smear_has_e); //O
  myTree->SetBranchAddress("smear_orig_tot_e",smear_orig_tot_e); //D
    myTree->SetBranchAddress("tot_e",tot_e); //D
  myTree->SetBranchAddress("dir_x",dir_x); //D
  myTree->SetBranchAddress("dir_y",dir_y);
  myTree->SetBranchAddress("dir_z",dir_z);
  
 
  Int_t nevents=myTree->GetEntries();
  cout <<"Hello, World!" <<endl;


  vector<float> photonE;
  vector<float> phX;
  vector<float> phY;
  vector<float> phZ;

  TH1D pi0Mass("pi0Mass","pi0Mass",500,0.12,0.15);
  TH1D etaMass("etaMass","etaMass",500,0.5,0.6);

  //-4.5 -1.5 0.5, 0.5-4.5
  enum rapRegion{forward,barrel, backward, endRapRegion};

  //  float upperRegion={-1.5,0.5,4.5};
  vector<float> etaBinning;
  etaBinning.push_back(-1.5);
  etaBinning.push_back(0.5);
  etaBinning.push_back(4.5);

  char buffer[300];
  
  TH1D** pi0Histos=new TH1D*[4];
  TH1D** etaHistos=new TH1D*[4];
  for(int i=0;i<endRapRegion;i++)
    {
      sprintf(buffer,"pi0HistoBin_%d",i);
      pi0Histos[i]=new TH1D(buffer, buffer, 200,0.05,0.25);
      sprintf(buffer,"etaHistoBin_%d",i);
      etaHistos[i]=new TH1D(buffer, buffer, 200,0.3,0.8);
    }


  

  
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


      photonE.clear();
      phX.clear();
      phY.clear();
      phZ.clear();
      //      cout <<"looking at " << evt_prt_count<< " particles in event " <<endl;
      //for all particles in the event
      for(int pi=0;pi<evt_prt_count;pi++)
	{
	  //look for photons
	  if(pdg[pi]==22)
	    {
	      //	      cout <<"found photon " << endl;
	      //	      if(smear_has_e[pi])
		{
		  TLorentzVector lvphoton(dir_x[pi]*tot_e[pi],dir_y[pi]*tot_e[pi],dir_z[pi]*tot_e[pi],tot_e[pi]);
		  if(lvphoton.Eta()<-4.5|| lvphoton.Eta()>4.5 || tot_e[pi]<0.2)
		    {
		      continue;
		    }

		  photonE.push_back(tot_e[pi]);
		  phX.push_back(dir_x[pi]);
		  phY.push_back(dir_y[pi]);
		  phZ.push_back(dir_z[pi]);
		  //		  cout <<" found photon with e: " << tot_e[pi] <<" dirX: " << dir_x[pi] << " y: " << dir_y[pi] << " dirz: " << dir_z[pi] <<endl;
		  //		  cout <<" found photon number " << photonE.size() <<endl;
		}
	      
	    }
	}
      vector<TLorentzVector*> pi0s=constructPi0s(photonE,phX,phY,phZ);
      for(int i=0;i<pi0s.size();i++)
	{

	int etaBin=getBin(etaBinning,pi0s[i]->Eta());
	  if(etaBin>=0)
	    {
	      pi0Histos[etaBin]->Fill(pi0s[i]->M());
	      etaHistos[etaBin]->Fill(pi0s[i]->M());
	    }
	  
	  //pi0Mass.Fill(pi0s[i]->M());
	  //	  etaMass.Fill(pi0s[i]->M());

	  
	  delete pi0s[i];
	}
    }

  TCanvas c1;
  //  TF1 f("f","pol4(0)+gaus(5)");
  TF1 f("f","gaus(0)+pol4(3)");
  f.SetParameter(5,500);
  f.SetParameter(0,15000);
  f.SetParameter(1,0.13957);
  f.SetParameter(2,005);
  //  f.SetParameter(7,);

  for(int i=0;i<3;i++)
    {
      pi0Histos[i]->Draw();
      //      pi0Histos[i]->Fit(&f);
      sprintf(buffer,"pi0MassBin_%d.png",i);
      c1.SaveAs(buffer);
    }

  for(int i=0;i<3;i++)
    {
      etaHistos[i]->Draw();
      //      etaHistos[i]->Fit(&f);
      sprintf(buffer,"etaMassBin_%d.png",i);
      c1.SaveAs(buffer);
    }
  
  //  pi0Mass.Fit("f");
  //  pi0Mass.Draw();
  c1.SaveAs("pi0Mass.png");
  f.SetParameter(5,0.547862);
  //    f.SetParameter(5,14000);
  //  etaMass.Fit(&f);
  etaMass.Draw();
  c1.SaveAs("etaMass.png");
  
  float** totalRecsEta=allocateArray<float>(3,binningEta.size());
  float** totalRecsPt=allocateArray<float>(3,binningPt.size());
  
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
