#ifndef STUDYPI0_H
#define STUDYPI0_H

#include "TLorentzVector.h"
#include "TVector3.h"



vector<TLorentzVector*> constructPi0s(vector<float>& e, vector<float>& dirX, vector<float>& dirY, vector<float>& dirZ)
{
  
  vector<TLorentzVector*> ret;
  if(e.size()<2)
    return ret;

  //    cout <<"combining " << e.size()-1 << " photons " <<endl;

  for(int i=0;i<(e.size()-1);i++)
    {
      //  cout <<" loop with i: " << i <<endl;
      for(int j=i+1;j<e.size();j++)
	{
	  //  cout <<" j:  : " << j <<endl;
	  //m^2=2E1E2(1-cos\theta)
	  //(or just from the lorentz vectors)
	  //	  TLorentzVector lv1(
	  TVector3 v1(dirX[i],dirY[i],dirZ[i]);
	  TVector3 v2(dirX[j],dirY[j],dirZ[j]);

	  float theta=v1.Angle(v2);
	  float pi0Mass=2*e[i]*e[j]*(1-cos(theta));
	  if(pi0Mass>0)
	    pi0Mass=sqrt(pi0Mass);
	  else
	    pi0Mass=-1;

	  TLorentzVector lv1(dirX[i]*e[i],dirY[i]*e[i],dirZ[i]*e[i],e[i]);
	  TLorentzVector lv2(dirX[j]*e[j],dirY[j]*e[j],dirZ[j]*e[j],e[j]);

	  TLorentzVector* lvPi0= new TLorentzVector(lv1+lv2);
	  //	  cout <<"comparing pi0masses: " <<pi0Mass <<" and " << lvPi0->M() <<" diff: "<< pi0Mass-lvPi0->M()<<endl;
	  ret.push_back(lvPi0);
	}
    }
  return ret;
}

  
#endif
