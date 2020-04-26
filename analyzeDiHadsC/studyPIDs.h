#ifndef STUDYPID_HH
#define STUDYPID_HH
#include <utility>
#include <string>
#include "TGraph.h"

using namespace std;

enum particleType{Pi,K};
enum pairType{PiPi,PiK,KK,pairTypeEnd};


void setGraphProps(TGraph* graph, int pt1, int pt2)
{
  	  if(pt2==0)
           {
	      graph->GetYaxis()->SetRangeUser(0.0,1.05);
	      graph->SetMarkerStyle(20);
	      graph->SetMarkerColor(kBlue);
	    }
	  if(pt2==1)
	    {
	      graph->SetMarkerStyle(21);
	      graph->SetMarkerColor(kRed);
	    }
	  if(pt2==2)
	    {
	      graph->SetMarkerStyle(22);
	      graph->SetMarkerColor(kGreen);
	    }
}

int getRecPID(float eta, float pt, int pid)
{

  //3 sigma: 99.7
  //2 sigma: 95.0
  int k=rand()%1000;

  if(k<50)
    {
      if(pid==Pi)
	return K;
      if(pid==K)
	return Pi;
    }
  
  return pid;  
}

int getRecPair(int pid1, int pid2)
{

  if(pid1==Pi && pid2==Pi)
    return PiPi;
  if(pid1==Pi && pid2==K)
    return PiK;
  if(pid1==K && pid2==Pi)
    return PiK;
  if(pid1==K && pid2==K)
    return KK;
}

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

int getPairCode(int firstPid, int secondPid, int charge1, int charge2)
{
  //++
  if(charge1>0 && charge2>0)
    {
      if(firstPid==Pi && secondPid==Pi)
	return 0x33;
      if(firstPid==K && secondPid==Pi)
	return 0x63;
      if(firstPid==K && secondPid==K)
	return 0x66;
    }
  //+ -
    if(charge1>0 && charge2<0)
      {
	if(firstPid==Pi && secondPid==K)
	  return 0x37;
	if(firstPid==Pi && secondPid==Pi)
	  return 0x34;
	if(firstPid==K && secondPid==Pi)
	  return 0x64;
	if(firstPid==K && secondPid==K)
	  return 0x67;
      }
    // --
    if(charge1<0 && charge2<0)
      {

      if(firstPid==Pi && secondPid==Pi)
	return 0x44;
      if(firstPid==K && secondPid==Pi)
	return 0x74;
      if(firstPid==K && secondPid==K)
	return 0x77;
	
      }
}


pair<int,int> getPids(int hadT){
  switch(hadT)
    {
    case 0x33:
      return pair<int,int>(Pi,Pi);
      break;
    case 0x34:
      return pair<int,int>(Pi,Pi);
      break;
    case 0x37:
      return pair<int,int>(Pi,K);
      break;
    case 0x44:
      return pair<int,int>(Pi,Pi);
      break;
    case 0x63:
      return pair<int,int>(K,Pi);
      break;
    case 0x64:
      return pair<int,int>(K,Pi);
      break;
    case 0x66:
      return pair<int,int>(K,K);
      break;
    case 0x67:
      return pair<int,int>(K,K);
      break;
    case 0x74:
      return pair<int,int>(K,Pi);
      break;
    case 0x77:
      return pair<int,int>(K,K);
      break;
      
    }
};

pair<int,int> getCharges(int hadT){
  switch(hadT)
    {
    case 0x33:
      return pair<int,int>(1,1);
      break;
    case 0x34:
      return pair<int,int>(1,-1);
      break;
    case 0x37:
      return pair<int,int>(1,-1);
      break;
    case 0x44:
      return pair<int,int>(-1,-1);
      break;
    case 0x63:
      return pair<int,int>(1,1);
      break;
    case 0x64:
      return pair<int,int>(1,-1);
      break;
    case 0x66:
      return pair<int,int>(1,1);
      break;
    case 0x67:
      return pair<int,int>(1,-1);
      break;
    case 0x74:
      return pair<int,int>(-1,-1);
      break;
    case 0x77:
      return pair<int,int>(-1,-1);
      break;
    } 
}

template<class T> T* allocateArray(int dim1);
template<class T> T** allocateArray(int dim1, int dim2);
template<class T> T*** allocateArray(int dim1, int dim2, int dim3);
template<class T> T**** allocateArray(int dim1, int dim2, int dim3, int dim4);
template<class T> T***** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5);
template<class T> T****** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6);
template<class T> T******* allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6, int dim7);
int getBin(vector<float>& b1, float value)
{
  int coo1=-1;

  for(unsigned int i=0;i<b1.size();i++)
    {
      if(value<=b1[i])
	{
	coo1=i;
	break;
	}
    }
  /*  if(coo1<0)
    {
        cout <<"wrong coo: val: " << value <<endl;
	}*/
  //  cout <<"value: " << value <<" coo: " << coo1 <<endl;
  return coo1;
}


#endif
