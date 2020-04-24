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

  Float_t x, y, Q2, W, Nu;
    Float_t x_pythia, Q2_pythia, y_pythia, W_pythia, Nu_pythia;

    Float_t eleE,eleP,elePt,eleEta,elePhi;
    Float_t eleVertex[3];

    Float_t hadE[2];
    Float_t hadP[2];
    Float_t hadPt[2];
    Float_t hadEta[2];
    Float_t hadPhi[2];
    Float_t hadVertex[2][3];

    Float_t PhMag; // dihadron total momentum
    Float_t PhPerpMag; // transverse component of dihadron total momentum (perp frame)
    Float_t PhEta; // pseudorapidity of dihadron pair
    Float_t PhPhi; // azimuth of dihadron pair
    Float_t RMag; // dihadron relative momentum
    Float_t RTMag; // transverse componet of relative momentum (T-frame)
    Float_t RPerpMag; // transverse componet of relative momentum (perp-frame)

  
    Float_t PhiH; // angle[ reaction_plane, Ph^q ]
    Float_t PhiR; // angle[ reaction_plane, RT^q ]
    Float_t PhiS; // angle[ reaction_plane, S^q ]
    Float_t Z[2];

  // carried by the hadron
      Float_t Zpair; // fraction of energy of fragmenting parton

  // carried by the hadron pair
    Float_t Mh; // dihadron invariant mass
    Float_t Mmiss; // missing mass
    Float_t xF; // feynman-x
    Float_t alpha; //opening angle between hadrons
    Float_t hadXF[2]; // feynman-x for each hadron
    
    Float_t zeta; // lab-frame energy sharing
    Float_t theta; // CoM-frame angle between Ph and P1
    Float_t MRterm[2]; // (used for computing theta)
    
    Int_t evnum,evnum_local;
    Int_t spinE,spinP;
    Int_t particleCnt;

    Int_t runnum;
    Int_t eleStatus;
    Float_t eleChi2pid;
    Bool_t eleFidPCAL[3];
    Bool_t eleFidDC[3];
    Int_t pairType,hadOrder;
    Int_t hadIdx[2];
    Int_t hadI[2];
    Int_t hadStatus[2];



    TChain* chAll=new TChain("tree");
    string fileType="pythia";


    chAll->SetBranchAddress("pairType",&pairType);
    chAll->SetBranchAddress("hadOrder",&hadOrder);
    chAll->SetBranchAddress("evnum",&evnum);
    chAll->SetBranchAddress("particleCnt",&particleCnt);

    chAll->SetBranchAddress("spinE",&spinE);
    chAll->SetBranchAddress("spinP",&spinP);

    chAll->SetBranchAddress("Q2",&Q2);
    chAll->SetBranchAddress("x",&x);
    chAll->SetBranchAddress("y",&y);
    chAll->SetBranchAddress("W",&W);
    chAll->SetBranchAddress("Nu",&Nu);
    if(fileType=="pythia") {
      chAll->SetBranchAddress("Q2_pythia",&Q2_pythia);
      chAll->SetBranchAddress("x_pythia",&x_pythia);
      chAll->SetBranchAddress("y_pythia",&y_pythia);
      chAll->SetBranchAddress("W_pythia",&W_pythia);
      chAll->SetBranchAddress("Nu_pythia",&Nu_pythia);
    };

    chAll->SetBranchAddress("eleE",&eleE);
    chAll->SetBranchAddress("eleP",&eleP);
    chAll->SetBranchAddress("elePt",&elePt);
    chAll->SetBranchAddress("eleEta",&eleEta);
    chAll->SetBranchAddress("elePhi",&elePhi);
    chAll->SetBranchAddress("eleVertex",eleVertex);
    chAll->SetBranchAddress("hadE",hadE);
    chAll->SetBranchAddress("hadP",hadP);
    chAll->SetBranchAddress("hadPt",hadPt);
    chAll->SetBranchAddress("hadEta",hadEta);
    chAll->SetBranchAddress("hadPhi",hadPhi);
    chAll->SetBranchAddress("hadVertex",hadVertex);

    chAll->SetBranchAddress("Mh",&Mh);
    chAll->SetBranchAddress("Mmiss",&Mmiss);
    chAll->SetBranchAddress("Z",Z);
    chAll->SetBranchAddress("Zpair",&Zpair);
    chAll->SetBranchAddress("xF",&xF);
    chAll->SetBranchAddress("hadXF",hadXF);
    chAll->SetBranchAddress("alpha",&alpha);
    chAll->SetBranchAddress("theta",&theta);
    chAll->SetBranchAddress("zeta",&zeta);
    chAll->SetBranchAddress("PhiH",&PhiH);
    chAll->SetBranchAddress("PhiR",&PhiR);
    chAll->SetBranchAddress("PhiS",&PhiS);

    chAll->SetBranchAddress("PhiH",&PhiH);
    chAll->SetBranchAddress("PhiR",&PhiR);
    chAll->SetBranchAddress("PhiS",&PhiS);
    
    chAll->SetBranchAddress("Ph",&PhMag);
    chAll->SetBranchAddress("PhPerp",&PhPerpMag);
    chAll->SetBranchAddress("PhEta",&PhEta);
    chAll->SetBranchAddress("PhPhi",&PhPhi);
    chAll->SetBranchAddress("R",&RMag);
    chAll->SetBranchAddress("RT",&RTMag);
    chAll->SetBranchAddress("RPerp",&RPerpMag);
 

    
  char* rootPath=argv[1];
  //  chAll->Add((string(rootPath)+"/*.root").c_str());
  chAll->Add(rootPath);
  Int_t nevents=chAll->GetEntries();
  cout <<"Hello, World!" <<endl;


  TH1D h1("xF","xF",100,0,1.0);
  TH1D h1Q2("Q2","Q2",100,0,5.0);
  TH1D h1Q2Pythia("Q2Pythia","Q2Pythia",100,0,5.0);
  TH1D h1Q2Diff("Q2Diff","Q2Diff",100,-2.0,2.0);


  TH2D xVsQ2("xVsQ2","xVsQ2",200,0.0005,1.0,200,0.1,100);
  
  for(long i=0;i<nevents;i++)
    {
#ifdef MAX_EVENTS
      if(i>MAX_EVENTS)
	break;
#endif


      if(!(i%10000))
	cout <<"processing event nr " << i << " of " << nevents << "(" << 100*i/(float)nevents<< "% )"<<endl;
      chAll->GetEntry(i);
      //      if(Z[0]>0.2 && Z[1]>0.1)
	{
	  h1.Fill(hadXF[0]);
	  h1.Fill(hadXF[1]);
	}




	pair<int,int> pids=getPids(pairType);
	pair<int,int> charges=getCharges(pairType);
	int firstRecPid=getRecPID(hadEta[0],hadPt[0],pids.first);
	int secondRecPid=getRecPID(hadEta[1],hadPt[1],pids.second);
	int newType=getPairCode(firstRecPid,secondRecPid,charges.first,charges.second);


	
      h1Q2.Fill(Q2);
      h1Q2Pythia.Fill(Q2_pythia);
      h1Q2Diff.Fill(Q2-Q2_pythia);
      xVsQ2.Fill(x,Q2);
    }

  TCanvas c1;
  h1.Draw();
  c1.SaveAs("hXf.png");

  h1Q2.Draw();
  c1.SaveAs("Q2.png");

  h1Q2Pythia.Draw();
  c1.SaveAs("Q2Pythia.png");

    h1Q2Diff.Draw();
  c1.SaveAs("Q2Diff.png");

  c1.SetLogx();
  c1.SetLogy();
  c1.SetLogz();
  xVsQ2.Draw("colz");
  c1.SaveAs("xVsQ2.root");
    c1.SaveAs("xVsQ2.png");
  
}
