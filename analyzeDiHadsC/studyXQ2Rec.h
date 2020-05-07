#ifndef STUDYXQ2REC_HH
#define STUDYXQ2REC_HH


#define PionMass 0.13957
#define KaonMass 0.493677
#define ElectronMass 0.000511
#define ProtonMass 0.938272088
#define MuonMass 0.10566
#define NeutronMass 0.93957

struct Kins
{
  float x;
  float Q2;
  float y;
  float nu;
  float W;
};


struct HadronicVars
{
  double sumPx;
  double sumPy;
  double sumPz;
  double sumEMinusPz;
  double theta;
};

float getParticleMass(int pdg)
{
  if(fabs(pdg)==211)
    return PionMass;
  if(fabs(pdg)==321)
    return KaonMass;
  if(fabs(pdg)==11)
    return ElectronMass;
  if(fabs(pdg)==13)
    return MuonMass;
  if(fabs(pdg)==2212)
    return ProtonMass;
  if(fabs(pdg)==2112)
    return NeutronMass;

  return 0;
  
}

float isTrack(int pdg)
{
  int fpdg=fabs(pdg);
  if(fpdg==211 || fpdg==321|| fpdg==11 || fpdg==13 || fpdg==2212)
    return true;
  return false;
}

HadronicVars getHadronicVars(int evt_prt_count,Double_t* px,Double_t* py,Double_t* pz, Double_t* dir_x, Double_t* dir_y, Double_t* dir_z, Double_t* tot_e,Long64_t* pdg, Bool_t* smear_has_e, Bool_t* smear_has_p, bool smeared=true)
{
  HadronicVars ret;
  double sumPx=0;
  double sumPy=0;
  double sumPz=0;
  double sumEMinusPz=0;
  
  for(int pi=0;pi<evt_prt_count;pi++)
    {
      //exclude scattered electron (should always be first in any case)
      if(pdg[pi]==11)
	continue;
      if(isTrack(pdg[pi]) && (!smeared || smear_has_p[pi]))
	{

	  double mass=getParticleMass(pdg[pi]);
	  //	  cout <<"is track, mass: " << mass <<endl;
	  double particleE=sqrt(px[pi]*px[pi]+py[pi]*py[pi]+pz[pi]*pz[pi]+mass*mass);
	  sumPx+=px[pi];
	  sumPy+=py[pi];
	  sumPz+=pz[pi];
	  sumEMinusPz+=(particleE-pz[pi]);
	}
      if(pdg[pi]==22 &&(!smeared || smear_has_e[pi]))
	{
	  //	  cout <<" is photon " << endl;
	  //	  TLorentzVector lvphoton(dir_x[pi]*tot_e[pi],dir_y[pi]*tot_e[pi],dir_z[pi]*tot_e[pi],tot_e[pi]);
	  sumPx+=dir_x[pi]*tot_e[pi];
	  sumPy+=dir_y[pi]*tot_e[pi];
	  sumPz+=dir_z[pi]*tot_e[pi];
	  sumEMinusPz+=(tot_e[pi]-dir_z[pi]*tot_e[pi]);
	}
    }

  ret.sumPx=sumPx;
  ret.sumPy=sumPy;
  ret.sumPz=sumPz;
  ret.sumEMinusPz=sumEMinusPz;
  TVector3 v(sumPx,sumPy,sumPz);
  ret.theta=v.Theta();
  return ret;
}


Kins getKinsJB2(HadronicVars v, double beamEnergy, double s)
{
  Kins ret;
  double pt2=v.sumPx*v.sumPx+v.sumPy*v.sumPy;
  ret.y=v.sumEMinusPz/(2*beamEnergy);
  if(ret.y>0)
    {
      ret.Q2=pt2/(1-ret.y);
      ret.x=ret.Q2/(s*ret.y);
    }
  return ret;
}

Kins getKinsJB(int evt_prt_count,Double_t* px,Double_t* py,Double_t* pz, Double_t* dir_x, Double_t* dir_y, Double_t* dir_z, Double_t* tot_e,Long64_t* pdg, Bool_t* smear_has_e, Bool_t* smear_has_p, float beamEnergy,float s, bool smeared=true)
{
  Kins ret;
  double sumPx=0;
  double sumPy=0;
  double sumPz=0;
  double sumEMinusPz=0;
  for(int pi=0;pi<evt_prt_count;pi++)
    {
      //      cout <<" particle " << pi << " pdg: "<< pdg[pi] <<endl;
      //exclude scattered electron (should always be first in any case)
      if(pdg[pi]==11)
	continue;
      if(isTrack(pdg[pi]) && (!smeared || smear_has_p[pi]))
	{

	  double mass=getParticleMass(pdg[pi]);
	  //	  cout <<"is track, mass: " << mass <<endl;
	  double particleE=sqrt(px[pi]*px[pi]+py[pi]*py[pi]+pz[pi]*pz[pi]+mass*mass);
	  sumPx+=px[pi];
	  sumPy+=py[pi];
	  sumEMinusPz+=(particleE-pz[pi]);
	}
      if(pdg[pi]==22 &&(!smeared || smear_has_e[pi]))
	{
	  //	  cout <<" is photon " << endl;
	  //	  TLorentzVector lvphoton(dir_x[pi]*tot_e[pi],dir_y[pi]*tot_e[pi],dir_z[pi]*tot_e[pi],tot_e[pi]);
	  sumPx+=dir_x[pi]*tot_e[pi];
	  sumPy+=dir_y[pi]*tot_e[pi];
	  sumEMinusPz+=(tot_e[pi]-dir_z[pi]*tot_e[pi]);
	}
    }
  double pt2=sumPx*sumPx+sumPy*sumPy;
  ret.y=sumEMinusPz/(2*beamEnergy);
  if(ret.y>0)
    {
      ret.Q2=pt2/(1-ret.y);
      ret.x=ret.Q2/(s*ret.y);
    }
  return ret;
}


//electron goes into the negative direction
Kins getKinsFromScatElectron(double eBeamEnergy, double hadronBeamEnergy,double scatElPx,double scatElPy,double scatElPz,double scatElectronEnergy)
{
  Kins ret;
  double m_e = 0.000510998950;
  double m_p = 0.938272088;
  TLorentzVector lv_e;// = new TLorentzVector();
  TLorentzVector lv_beam;// = new TLorentzVector();
  TLorentzVector lv_target;// = new TLorentzVector();

  //e2 =m2+p2 -->p2=e2-m2
  double longBeamMomentum=sqrt(eBeamEnergy*eBeamEnergy-m_e*m_e);
  double longHadBeamMomentum=sqrt(hadronBeamEnergy*hadronBeamEnergy-m_p*m_p);

  cout <<" longbeamMomentum: " << longBeamMomentum << " beam Mom: " << longHadBeamMomentum <<endl;
  lv_e.SetPxPyPzE(scatElPx, scatElPy, scatElPz, scatElectronEnergy);
  lv_beam.SetPxPyPzE(0.0, 0.0, -longBeamMomentum, eBeamEnergy);
  lv_target.SetPxPyPzE(0.0, 0.0, longHadBeamMomentum, hadronBeamEnergy);
  
    TLorentzVector q=lv_beam-lv_e;
  //TLorentzVector q(lv_beam);
  // q=q-lv_e;
  //  q.sub(lv_e);
  // System.out.println("q x " + q.px()+ " y: "+q.py() + " pz: " +q.pz() + " e: "+
  // q.e() );
  // q= new LorentzVector(lv_e);
  // q.sub(lv_beam);
  ret.Q2 = (-1) * q.M2();
  // need to multiply lorentz vectors... doesn't seem to be implemented in the
  // jlab clas
  // since this is lv_target*q /m_p and the momentum of the target is zero, it is
  // just the product of the energy parts, divided by m_p:
  ret.nu = lv_target.E() * q.E() / m_p;
  //only true for fixed target
  ret.y = (lv_target*q)/(lv_target*lv_beam);
  // System.out.println("target x " + lv_target.px()+ " y: "+lv_target.py() + "
  // pz: " +lv_target.pz() + " e: "+ lv_target.e() );
  ret.x = ret.Q2 / (2 * (lv_target*q));
  ret.W = sqrt(m_p+ret.Q2 * (1 - ret.x) / ret.x);
  return ret;
}
//following https://arxiv.org/pdf/1208.6087v2.pdf
Kins getKinsDA(double scatElPx,double scatElPy,double scatElPz, double scatElectronEnergy,double tP, double s)
{
  Kins ret;
  TLorentzVector lv_e;
  lv_e.SetPxPyPzE(scatElPx, scatElPy, scatElPz, scatElectronEnergy);
  double t=lv_e.Theta();
  double nominator=4*lv_e.E()*lv_e.E()*cos(t/2)*cos(t/2);
  double denom=sin(t/2)*sin(t/2)+sin(t/2)*cos(t/2)*tan(tP/2);
  ret.Q2=nominator/denom;

  //this is from the Bluemlein paper, but doesn't seem to work
  //  ret.y=1.0-(sin(t/2))/(sin(t/2)+cos(t/2)*tan(tP/2));
  //Zeus paper:
  ret.y=tan(tP/2);
  ret.y/=(tan(t/2)+tan(tP/2));
  
  if(ret.y>0)
    {
      ret.x=ret.Q2/(s*ret.y);
    }
  
  
  return ret;
  
}



#endif
