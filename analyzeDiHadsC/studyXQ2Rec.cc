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

  TH2D corrConventionalQ2("convQ2","convQ2",100 ,1,200,100 ,1,200);
  TH2D corrConventionalX("convX","convX",100 ,0.0,1.0,100 ,-1.0,1.0);

  TH2D corrJBQ2("jbQ2","jbQ2",100 ,1,200,100 ,1,200);
  TH2D corrJBX("jbX","jbX",100 ,0.0,1.0,100 ,-1.0,1.0);

  TH2D corrDAQ2("daQ2","daQ2",100 ,1,200,100 ,1,200);
  TH2D corrDAX("daX","daX",100 ,0.0,1.0,100 ,-1.0,1.0);

  TH2D corrMixedQ2("mixedQ2","mixedQ2",100 ,1,200,100 ,1,200);
  TH2D corrMixedX("mixedX","mixedX",100 ,0.0,1.0,100 ,-1.0,1.0);


  

  TH2D corrConventionalQ2SmallY("convQ2SmallY","convQ2SmallY",100 ,1,200,100 ,1,200);
  TH2D corrConventionalXSmallY("convXSmallY","convXSmallY",100 ,0.0,1.0,100 ,-1.0,1.0);

  TH2D corrJBQ2SmallY("jbQ2SmallY","jbQ2SmallY",100 ,1,200,100 ,1,200);
  TH2D corrJBXSmallY("jbXSmallY","jbXSmallY",100 ,0.0,1.0,100 ,-1.0,1.0);

  TH2D corrDAQ2SmallY("daQ2SmallY","daQ2SmallY",100 ,1,200,100 ,1,200);
  TH2D corrDAXSmallY("daXSmallY","daXSmallY",100 ,0.0,1.0,100 ,-1.0,1.0);

  TH2D corrMixedQ2SmallY("mixedQ2SmallY","mixedQ2SmallY",100 ,1,200,100 ,1,200);
  TH2D corrMixedXSmallY("mixedXSmallY","mixedXSmallY",100 ,0.0,1.0,100 ,-1.0,1.0);

  
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
      //      cout <<"original e px: " << ePxOrig <<" py: " << ePyOrig<< " pz: "<< ePzOrig<<endl;

      //float sqrtS=63.2527; //sqrtS for 10x100;
      //      float sqrtS=45; //sqrtS for 5x100;
         float sqrtS=140.716; //sqrtS for 18x275;
      //      float sqrtS=28.6; //sqrtS for 5x41;
      
      //calculate Q2 just from original and scattered electron
      Kins kinSmeared=getKinsFromScatElectron(beamEnergy,hadronBeamEnergy,ePx,ePy,ePz,eE);
      Kins kinOrig=getKinsFromScatElectron(beamEnergy,hadronBeamEnergy,ePxOrig,ePyOrig,ePzOrig,eEOrig);
      Kins kinJBSmeared=getKinsJB(evt_prt_count,px,py,pz,dir_x,dir_y,dir_z,tot_e,pdg,smear_has_e, smear_has_p,beamEnergy,sqrtS*sqrtS);
      Kins kinJB=getKinsJB(evt_prt_count,smear_orig_px,smear_orig_py,smear_orig_pz,dir_x,dir_y,dir_z,smear_orig_tot_e,pdg,smear_has_e, smear_has_p,beamEnergy,sqrtS*sqrtS,false);

      HadronicVars hvOrig=getHadronicVars(evt_prt_count,smear_orig_px,smear_orig_py,smear_orig_pz,dir_x,dir_y,dir_z,smear_orig_tot_e,pdg,smear_has_e, smear_has_p,false);
      HadronicVars hv=getHadronicVars(evt_prt_count,px,py,pz,dir_x,dir_y,dir_z,tot_e,pdg,smear_has_e, smear_has_p);
      
      Kins kinJB2=getKinsJB2(hvOrig,beamEnergy, sqrtS*sqrtS);
      
      Kins kinsDAOrig=getKinsDA(ePxOrig,ePyOrig,ePzOrig,eEOrig,hvOrig.theta, sqrtS*sqrtS);
      Kins kinsDASmeared=getKinsDA(ePx,ePy,ePz,eE,hv.theta, sqrtS*sqrtS);
      //      cout <<" kin orig x: " << kinOrig.x <<" Q2 conv: " << kinOrig.Q2 <<" y: " << kinOrig.y<<endl;
      //      cout <<" kin orig 2 x " << kinJB2.x <<" Q2: "<< kinJB2.Q2 <<endl;
      //      cout <<" kin conventional x: " << kinSmeared.x <<" JB: " << kinJB.x <<" Q2 conv: " << kinSmeared.Q2 <<" JB: " << kinJB.Q2 <<endl;
      //      cout <<" conv y: " << kinSmeared.y <<" JB: " << kinJB.y <<endl;
      //      cout <<"kins DA: x: " << kinsDAOrig.x <<" Q2: "<< kinsDAOrig.Q2 <<" y: " << kinsDAOrig.y<<endl;

      float mixedXSmeared=kinSmeared.Q2/(sqrtS*sqrtS*kinJBSmeared.y);
      float mixedXOrig=kinOrig.Q2/(sqrtS*sqrtS*kinJB.y);      
      
      if(kinOrig.Q2>4)
	{
	  //split y in above 0.01 and below 0.01
	  if(kinOrig.y>0.01)
	    {
	      corrConventionalQ2.Fill(kinOrig.Q2,(kinOrig.Q2-kinSmeared.Q2)/kinOrig.Q2);
	      corrConventionalX.Fill(kinOrig.x,(kinOrig.x-kinSmeared.x)/kinOrig.x);
	      
	      corrJBQ2.Fill(kinOrig.Q2,(kinOrig.Q2-kinJBSmeared.Q2)/kinOrig.Q2);
	      corrJBX.Fill(kinOrig.x,(kinOrig.x-kinJBSmeared.x)/kinOrig.x);
	      
	      corrDAQ2.Fill(kinOrig.Q2,(kinOrig.Q2-kinsDASmeared.Q2)/kinOrig.Q2);
	      corrDAX.Fill(kinOrig.x,(kinOrig.x-kinsDASmeared.x)/kinOrig.x);
	      
	      corrMixedX.Fill(kinOrig.x,(kinOrig.x-mixedXSmeared)/kinOrig.x);
	      
	      
	    }
	  else
	    {
	      corrConventionalQ2SmallY.Fill(kinOrig.Q2,(kinOrig.Q2-kinSmeared.Q2)/kinOrig.Q2);
	      //	      corrConventionalXSmallY.Fill(kinOrig.x,kinSmeared.x);
	      corrConventionalXSmallY.Fill(kinOrig.x,(kinOrig.x-kinSmeared.x)/kinOrig.x);
	      
	      corrJBQ2SmallY.Fill(kinOrig.Q2,(kinOrig.Q2-kinJBSmeared.Q2)/kinOrig.Q2);
	      corrJBXSmallY.Fill(kinOrig.x,(kinOrig.x-kinJBSmeared.x)/kinOrig.x);
	      
	      corrDAQ2SmallY.Fill(kinOrig.Q2,(kinOrig.Q2-kinsDASmeared.Q2)/kinOrig.Q2);
	      corrDAXSmallY.Fill(kinOrig.x,(kinOrig.x-kinsDASmeared.x)/kinOrig.x);
	      
	      corrMixedXSmallY.Fill(kinOrig.x,(kinOrig.x-mixedXSmeared)/kinOrig.x);
	    }
	}

      
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


  char buffer[200];
  TCanvas c1;
   c1.SetLogx();
//  c1.SetLogy();

  c1.SetLogz();
  corrConventionalQ2.Draw("colz");
  sprintf(buffer,"correlationsConventionalQ2.png");
  c1.SaveAs(buffer);
  corrConventionalX.Draw("colz");
  sprintf(buffer,"correlationsConventionalX.png");
  c1.SaveAs(buffer);

  
  corrJBQ2.Draw("colz");
  sprintf(buffer,"correlationsJBQ2.png");
  c1.SaveAs(buffer);

  
  corrJBX.Draw("colz");
  sprintf(buffer,"correlationsJBX.png");
  c1.SaveAs(buffer);

  
  corrDAQ2.Draw("colz");
  sprintf(buffer,"correlationsDAQ2.png");
  c1.SaveAs(buffer);

  corrDAX.Draw("colz");
    sprintf(buffer,"correlationsDAX.png");
  c1.SaveAs(buffer);

  corrMixedX.Draw("colz");
  sprintf(buffer,"correlationsMixedX.png");
  c1.SaveAs(buffer);



  
  corrConventionalQ2SmallY.Draw("colz");
  sprintf(buffer,"correlationsConventionalQ2SmallY.png");
  c1.SaveAs(buffer);
  corrConventionalXSmallY.Draw("colz");
  sprintf(buffer,"correlationsConventionalXSmallY.png");
  c1.SaveAs(buffer);

  
  corrJBQ2SmallY.Draw("colz");
  sprintf(buffer,"correlationsJBQ2SmallY.png");
  c1.SaveAs(buffer);

  
  corrJBXSmallY.Draw("colz");
  sprintf(buffer,"correlationsJBXSmallY.png");
  c1.SaveAs(buffer);

  
  corrDAQ2SmallY.Draw("colz");
  sprintf(buffer,"correlationsDAQ2SmallY.png");
  c1.SaveAs(buffer);

  corrDAXSmallY.Draw("colz");
    sprintf(buffer,"correlationsDAXSmallY.png");
  c1.SaveAs(buffer);

  corrMixedXSmallY.Draw("colz");
  sprintf(buffer,"correlationsMixedXSmallY.png");
  c1.SaveAs(buffer);
	  

  

    
  
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
