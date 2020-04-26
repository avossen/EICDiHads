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
#include "TGraph.h"

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


  int totalIndex=3;
  
  //fraction of pi pi that are pik
  float** numPiPiEta=allocateArray<float>(4,binningEta.size());
  float** numPiPiPt=allocateArray<float>(4,binningPt.size());

  float** numPiKEta=allocateArray<float>(4,binningEta.size());
  float** numPiKPt=allocateArray<float>(4,binningPt.size());

  float** numKKEta=allocateArray<float>(4,binningEta.size());
  float** numKKPt=allocateArray<float>(4,binningPt.size());

  float*** numRecEta=allocateArray<float>(3,4,binningEta.size());
  float*** numRecPt=allocateArray<float>(3,4,binningPt.size());




  //  cout<<"eta size: "<< binningEta.size() <<" pt size: " << binningPt.size() <<endl;
  
  //  float* meanEta=allocateArray<float>(binningEta.size());
  //  float* meanPt= allocateArray<float>(binningPt.size());
  float* meanEta=new float[binningEta.size()];
  float* meanPt=new float[binningPt.size()];

  //  float* countEta=allocateArray<float>(binningEta.size());
  //  float* countPt=allocateArray<float>(binningPt.size());
  //  float* countEta=new float[binningEta.size()];
  float* countEta=new float[binningEta.size()];
  float* countPt=new float[binningPt.size()];
  
  //fraction of pi pi that kk
  //fraction of pipi that stay pipi
  //  cout <<"size of binningpt: "<< binningPt.size()<<endl;
  //  cout <<"size of binningEta: "<< binningEta.size()<<endl;
  
  for(long i=0;i<nevents;i++)
    {
#ifdef MAX_EVENTS
      if(i>MAX_EVENTS)
	break;
#endif
      //      cout <<"yf1size of binningpt: "<< binningPt.size()<<endl;
      //      cout <<"yf1size of binningEta: "<< binningEta.size()<<endl;

      
      //      cout <<"getting eta bin "<<endl;
      int etaBin=getBin(binningEta,PhEta);
      //      cout <<"binningEta size: "<< binningEta.size() << " pheta: " << PhEta <<" bin: "<< etaBin <<endl;
      //      cout <<"getting pt bin for "<< PhPerpMag<<"numPt bins: " << binningPt.size()<<endl;
      //            cout <<"1yf1size of binningpt: "<< binningPt.size()<<endl;
      //      cout <<"1yf1size of binningEta: "<< binningEta.size()<<endl;

      int ptBin=getBin(binningPt,PhPerpMag);

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
      //      cout <<"3yf1size of binningpt: "<< binningPt.size()<<endl;
      //      cout <<"3yf1size of binningEta: "<< binningEta.size()<<endl;

      int firstRecPid=getRecPID(hadEta[0],hadPt[0],pids.first);
      int secondRecPid=getRecPID(hadEta[1],hadPt[1],pids.second);
      int newType=getPairCode(firstRecPid,secondRecPid,charges.first,charges.second);
      int recPair=getRecPair(firstRecPid,secondRecPid);



      
      //      cout <<"4yf1size of binningpt: "<< binningPt.size()<<endl;
      //      cout <<"4yf1size of binningEta: "<< binningEta.size()<<endl;
      //      cout <<"inc with etabin: " << etaBin <<endl;
      countEta[etaBin]++;
      meanEta[etaBin]+=PhEta;

      //      cout <<"5yf1size of binningpt: "<< binningPt.size()<<endl;
      //      cout <<"5yf1size of binningEta: "<< binningEta.size()<<endl;
      //      cout <<" ptBin: " << ptBin <<endl;
	
      countPt[ptBin]++;
      meanPt[ptBin]+=PhPerpMag;
	//      cout <<"6yf1size of binningpt: "<< binningPt.size()<<endl;
	//      cout <<"6yf1size of binningEta: "<< binningEta.size()<<endl;

	//	cout <<"etaBin: "<< etaBin << ", ptBin: "<< ptBin <<endl;

	//	cout <<"recPair: "<< recPair <<endl;

      //      if(firstRecPid==Pi&&secondRecPid==Pi)
            if(pids.first==Pi&&pids.second==Pi)
	  {
	    //	    cout <<"xf1size of binningpt: "<< binningPt.size()<<endl;
	    //	    cout <<"xf1size of binningEta: "<< binningEta.size()<<endl;

	    numPiPiEta[3][etaBin]++;
	    //	    cout <<"f1size of binningpt: "<< binningPt.size()<<endl;
	    //	    cout <<"f1size of binningEta: "<< binningEta.size()<<endl;

	    numPiPiEta[recPair][etaBin]++;
	    //	    cout <<"s1size of binningpt: "<< binningPt.size()<<endl;
	    //	    cout <<"s1size of binningEta: "<< binningEta.size()<<endl;

	    numPiPiPt[3][ptBin]++;
	    //	    cout <<"a1size of binningpt: "<< binningPt.size()<<endl;
	    //	    cout <<"a1size of binningEta: "<< binningEta.size()<<endl;

	    numPiPiPt[recPair][ptBin]++;


	  }
	//	  cout <<"1size of binningpt: "<< binningPt.size()<<%endl;
	//  cout <<"1size of binningEta: "<< binningEta.size()<<endl;
	    //	if((firstRecPid==Pi && secondRecPid==K) || (firstRecPid==K && secondRecPid==Pi))
	    if((pids.first==Pi && pids.second==K) || (pids.first==K && pids.second==Pi))
	  {
	    numPiKEta[3][etaBin]++;
	    numPiKEta[recPair][etaBin]++;

	    numPiKPt[3][ptBin]++;
	    numPiKPt[recPair][ptBin]++;

	  }
	//	  cout <<"2size of binningpt: "<< binningPt.size()<<endl;
	//	  cout <<"2size of binningEta: "<< binningEta.size()<<endl;
	//	  cout <<"etabin: "<< etaBin <<" ptBin: "<< ptBin <<" recPair: "<< recPair<< endl;
		//		if(firstRecPid==K && secondRecPid==K)
	    if(pids.first==K && pids.second==K)
	  {
	    //total counts
	    numKKEta[3][etaBin]++;
	    numKKEta[recPair][etaBin]++;

	    numKKPt[3][ptBin]++;
	    numKKPt[recPair][ptBin]++;
	  }

	//	  cout <<"3size of binningpt: "<< binningPt.size()<<endl;
	//  cout <<"3size of binningEta: "<< binningEta.size()<<endl;

      h1Q2.Fill(Q2);
      h1Q2Pythia.Fill(Q2_pythia);
      h1Q2Diff.Fill(Q2-Q2_pythia);
      xVsQ2.Fill(x,Q2);
    }

  float** totalRecsEta=allocateArray<float>(3,binningEta.size());
  float** totalRecsPt=allocateArray<float>(3,binningPt.size());

  for(int i=0;i<binningEta.size();i++)
    {
      meanEta[i]=meanEta[i]/countEta[i];
      for(int pairType=0;pairType<pairTypeEnd;pairType++)
	{
	  totalRecsEta[pairType][i]+=(numPiPiEta[pairType][i]+numPiKEta[pairType][i]+numKKEta[pairType][i]);
	}
      for(int recPType=0;recPType<pairTypeEnd;recPType++)
	{
	  numRecEta[recPType][PiPi][i]=numPiPiEta[recPType][i]/totalRecsEta[recPType][i];
	  numRecEta[recPType][PiK][i]=numPiKEta[recPType][i]/totalRecsEta[recPType][i];
	  numRecEta[recPType][KK][i]=numKKEta[recPType][i]/totalRecsEta[recPType][i];
	}

      for(int j=0;j<3;j++)
	{
	  numPiPiEta[j][i]/=numPiPiEta[3][i];
	  
	}
      for(int j=0;j<3;j++)
	{
	  numPiKEta[j][i]/=numPiKEta[3][i];
	}
      for(int j=0;j<3;j++)
	{
	  numKKEta[j][i]/=numKKEta[3][i];
	}
    }

  for(int i=0;i<binningPt.size();i++)
    {
      meanPt[i]=meanPt[i]/countPt[i];
      for(int pairType=0;pairType<pairTypeEnd;pairType++)
	{
	  totalRecsPt[pairType][i]+=(numPiPiPt[pairType][i]+numPiKPt[pairType][i]+numKKPt[pairType][i]);
	}

       for(int recPType=0;recPType<pairTypeEnd;recPType++)
	{
	  numRecPt[recPType][PiPi][i]=numPiPiPt[recPType][i]/totalRecsPt[recPType][i];
	  numRecPt[recPType][PiK][i]=numPiKPt[recPType][i]/totalRecsPt[recPType][i];
	  numRecPt[recPType][KK][i]=numKKPt[recPType][i]/totalRecsPt[recPType][i];

	}

      for(int j=0;j<3;j++)
	{
	  numPiPiPt[j][i]/=numPiPiPt[3][i];
	}
      for(int j=0;j<3;j++)
	{
	  numPiKPt[j][i]/=numPiKPt[3][i];
	}
      for(int j=0;j<3;j++)
	{
	  numKKPt[j][i]/=numKKPt[3][i];
	}
    }


  TGraph pipi_pipiFractionEta(binningEta.size(),meanEta, numPiPiEta[0]);
  TGraph pipi_piKFractionEta(binningEta.size(),meanEta, numPiPiEta[1]);
  TGraph pipi_KKFractionEta(binningEta.size(),meanEta, numPiPiEta[2]);

  TGraph piK_pipiFractionEta(binningEta.size(),meanEta, numPiKEta[0]);
  TGraph piK_piKFractionEta(binningEta.size(),meanEta, numPiKEta[1]);
  TGraph piK_KKFractionEta(binningEta.size(),meanEta, numPiKEta[2]);


  TGraph KK_pipiFractionEta(binningEta.size(),meanEta, numKKEta[0]);
  TGraph KK_piKFractionEta(binningEta.size(),meanEta, numKKEta[1]);
  TGraph KK_KKFractionEta(binningEta.size(),meanEta, numKKEta[2]);



  TGraph pipi_pipiFractionPt(binningPt.size(),meanPt, numPiPiPt[0]);
  TGraph pipi_piKFractionPt(binningPt.size(),meanPt, numPiPiPt[1]);
  TGraph pipi_KKFractionPt(binningPt.size(),meanPt, numPiPiPt[2]);

  TGraph piK_pipiFractionPt(binningPt.size(),meanPt, numPiKPt[0]);
  TGraph piK_piKFractionPt(binningPt.size(),meanPt, numPiKPt[1]);
  TGraph piK_KKFractionPt(binningPt.size(),meanPt, numPiKPt[2]);


  TGraph KK_pipiFractionPt(binningPt.size(),meanPt, numKKPt[0]);
  TGraph KK_piKFractionPt(binningPt.size(),meanPt, numKKPt[1]);
  TGraph KK_KKFractionPt(binningPt.size(),meanPt, numKKPt[2]);

  TGraph*** recGraphsEta=new TGraph**[pairTypeEnd];
  for(int pt=0;pt<pairTypeEnd;pt++)
    {
      recGraphsEta[pt]=new TGraph*[pairTypeEnd];
      for(int pt2=0;pt2<pairTypeEnd;pt2++)
	{
	  recGraphsEta[pt][pt2]=new TGraph(binningEta.size(),meanEta,numRecEta[pt][pt2]);
	  setGraphProps(recGraphsEta[pt][pt2],pt,pt2);
	}
    }
  TGraph*** recGraphsPt=new TGraph**[pairTypeEnd];
  for(int pt=0;pt<pairTypeEnd;pt++)
    {
      recGraphsPt[pt]=new TGraph*[pairTypeEnd];
      for(int pt2=0;pt2<pairTypeEnd;pt2++)
	{
	  recGraphsPt[pt][pt2]=new TGraph(binningPt.size(),meanPt,numRecPt[pt][pt2]);
	  setGraphProps(recGraphsPt[pt][pt2],pt,pt2);
	}
    }
  TCanvas c1;
  for(int pt=0;pt<pairTypeEnd;pt++)
    {
      for(int pt2=0;pt2<pairTypeEnd;pt2++)
	{
	  if(pt2==0)
	    {
	      recGraphsPt[pt][pt2]->Draw("AP");
	    }
	  if(pt2==1 || pt2==2)
	    {
	      recGraphsPt[pt][pt2]->Draw("P SAME");
	    }
	}
      if(pt==0)
	{
	  c1.SaveAs("recPiPiPt.png");
	}
      if(pt==1)
	{
	  c1.SaveAs("recPiKPt.png");
	}

      if(pt==2)
	{
	  c1.SaveAs("recKKPt.png");
	}

    }
  
  for(int pt=0;pt<pairTypeEnd;pt++)
    {
      for(int pt2=0;pt2<pairTypeEnd;pt2++)
	{
	  if(pt2==0)
	    {
	      recGraphsEta[pt][pt2]->Draw("AP");
	    }
	  if(pt2==1 || pt2==2)
	    {
	      recGraphsEta[pt][pt2]->Draw("P SAME");
	    }
	}
      if(pt==0)
	{
	  c1.SaveAs("recPiPiEta.png");
	}
      if(pt==1)
	{
	  c1.SaveAs("recPiKEta.png");
	}

      if(pt==2)
	{
	  c1.SaveAs("recKKEta.png");
	}

    }
  
  

  pipi_pipiFractionEta.GetYaxis()->SetRangeUser(0.0,1.0);
  pipi_pipiFractionEta.SetMarkerStyle(20);
  //  pipi_pipiFractionEta.SetMarkerSize(2);
  pipi_pipiFractionEta.SetMarkerColor(kBlue);
  pipi_piKFractionEta.SetMarkerStyle(21);
  pipi_piKFractionEta.SetMarkerColor(kRed);
  pipi_KKFractionEta.SetMarkerStyle(22);
  pipi_KKFractionEta.SetMarkerColor(kGreen);
  pipi_pipiFractionEta.Draw("AP");
  pipi_piKFractionEta.Draw("SAME P");
  pipi_KKFractionEta.Draw("SAME P");

  for(int i=0;i<6;i++)
    {
      //      cout <<"eta bin " << i <<" piK: "<< pipi_piKFractionEta
    }
  
  //  c1.Get
  c1.SaveAs("fractionsPiPiEta.png");

  piK_pipiFractionEta.GetYaxis()->SetRangeUser(0.0,1.0);
  piK_pipiFractionEta.SetMarkerStyle(20);
  //  piK_pipiFractionEta.SetMarkerSize(2);
  piK_pipiFractionEta.SetMarkerColor(kBlue);
  piK_piKFractionEta.SetMarkerStyle(21);
  piK_piKFractionEta.SetMarkerColor(kRed);
  piK_KKFractionEta.SetMarkerStyle(22);
  piK_KKFractionEta.SetMarkerColor(kGreen);
  piK_pipiFractionEta.Draw("AP");
  piK_piKFractionEta.Draw("SAME P");
  piK_KKFractionEta.Draw("SAME P");

  c1.SaveAs("fractionsPiKEta.png");

  
  KK_pipiFractionEta.GetYaxis()->SetRangeUser(0.0,1.0);
  KK_pipiFractionEta.SetMarkerStyle(20);
  //  pipi_pipiFractionEta.SetMarkerSize(2);
  KK_pipiFractionEta.SetMarkerColor(kBlue);
  
  KK_piKFractionEta.SetMarkerStyle(21);
  KK_piKFractionEta.SetMarkerColor(kRed);
  KK_KKFractionEta.SetMarkerStyle(22);
  KK_KKFractionEta.SetMarkerColor(kGreen);
  KK_pipiFractionEta.Draw("AP");
  KK_piKFractionEta.Draw("SAME P");
  KK_KKFractionEta.Draw("SAME P");

  c1.SaveAs("fractionsKKEta.png");

  for(int i=0;i<6;i++)
    {
      cout <<" PiPi pipi fraction: " << numPiPiEta[0][i] <<endl;
      cout <<" PiPi pipi fraction: " << numPiPiEta[1][i] <<endl;
      cout <<" PiPi pipi fraction: " << numPiPiEta[2][i] <<endl;
      
      cout <<" PiK pipi fraction: " << numPiKEta[0][i] <<endl;
      cout <<" PiK pipi fraction: " << numPiKEta[1][i] <<endl;
      cout <<" PiK pipi fraction: " << numPiKEta[2][i] <<endl;

      cout <<" KK pipi fraction: " << numKKEta[0][i] <<endl;
      cout <<" KK pipi fraction: " << numKKEta[1][i] <<endl;
      cout <<" KK pipi fraction: " << numKKEta[2][i] <<endl;

    }

  
  pipi_pipiFractionPt.GetYaxis()->SetRangeUser(0.0,1.0);
  pipi_pipiFractionPt.SetMarkerStyle(20);
  //  pipi_pipiFractionPt.SetMarkerSize(2);
  pipi_pipiFractionPt.SetMarkerColor(kBlue);
  pipi_piKFractionPt.SetMarkerStyle(21);
  pipi_piKFractionPt.SetMarkerColor(kRed);
  pipi_KKFractionPt.SetMarkerStyle(22);
  pipi_KKFractionPt.SetMarkerColor(kGreen);
  pipi_pipiFractionPt.Draw("AP");
  pipi_piKFractionPt.Draw("SAME P");
  pipi_KKFractionPt.Draw("SAME P");

  for(int i=0;i<6;i++)
    {
      //      cout <<"eta bin " << i <<" piK: "<< pipi_piKFractionPt
    }
  
  //  c1.Get
  c1.SaveAs("fractionsPiPiPt.png");

  piK_pipiFractionPt.GetYaxis()->SetRangeUser(0.0,1.0);
  piK_pipiFractionPt.SetMarkerStyle(20);
  //  piK_pipiFractionPt.SetMarkerSize(2);
  piK_pipiFractionPt.SetMarkerColor(kBlue);
  piK_piKFractionPt.SetMarkerStyle(21);
  piK_piKFractionPt.SetMarkerColor(kRed);
  piK_KKFractionPt.SetMarkerStyle(22);
  piK_KKFractionPt.SetMarkerColor(kGreen);
  piK_pipiFractionPt.Draw("AP");
  piK_piKFractionPt.Draw("SAME P");
  piK_KKFractionPt.Draw("SAME P");

  c1.SaveAs("fractionsPiKPt.png");

  
  KK_pipiFractionPt.GetYaxis()->SetRangeUser(0.0,1.0);
  KK_pipiFractionPt.SetMarkerStyle(20);
  //  pipi_pipiFractionPt.SetMarkerSize(2);
  KK_pipiFractionPt.SetMarkerColor(kBlue);
  
  KK_piKFractionPt.SetMarkerStyle(21);
  KK_piKFractionPt.SetMarkerColor(kRed);
  KK_KKFractionPt.SetMarkerStyle(22);
  KK_KKFractionPt.SetMarkerColor(kGreen);
  KK_pipiFractionPt.Draw("AP");
  KK_piKFractionPt.Draw("SAME P");
  KK_KKFractionPt.Draw("SAME P");

  c1.SaveAs("fractionsKKPt.png");

  
  
  
  
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
