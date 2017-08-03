// first test (Joern Putschke)

#include "ktTrackEff.hh"
#include "TFile.h"
#include "TRandom.h"

#include <Riostream.h>
//#include<iostream> // needed for io
//#include<sstream>  // needed for internal io
//#include<vector> 

//REMARK JP: Quick implementation for STAR only,
//           should be setup in a more general way !!!!

ClassImp(ktTrackEff)

ktTrackEff::ktTrackEff(TString mfName)
{

  fName=mfName;

  //centrality bins: 0 is 0-5%, 1 is 5-10%, 2 is 10-20%
  TFile *fEff = new TFile(fName,"read");
  //DEBUG:
  //fEff->ls();

  effY07pteta[0] = (TH2D*)fEff->Get("ptetaScale_0");
  effY07pteta[0]->SetName("effY07pteta_0"); effY07pteta[0]->SetDirectory(0);
  effY07pteta[1] = (TH2D*)fEff->Get("ptetaScale_1");
  effY07pteta[1]->SetName("effY07pteta_1"); effY07pteta[1]->SetDirectory(0);
  effY07pteta[2] = (TH2D*)fEff->Get("ptetaScale_2");
  effY07pteta[2]->SetName("effY07pteta_2"); effY07pteta[2]->SetDirectory(0);
  effY07eta[0] = (TH1D*)fEff->Get("etaScale_0");
  effY07eta[0]->SetName("effY07eta_0"); effY07eta[0]->SetDirectory(0);
  effY07eta[1] = (TH1D*)fEff->Get("etaScale_1");
  effY07eta[1]->SetName("effY07eta_1"); effY07eta[1]->SetDirectory(0);
  effY07eta[2] = (TH1D*)fEff->Get("etaScale_2");
  effY07eta[2]->SetName("effY07eta_2"); effY07eta[2]->SetDirectory(0);
  fEff->Close();
  effY04[0] = GetEffY04(0);
  effY04[0]->SetName("effY04_0");
  effY04[1] = GetEffY04(1);
  effY04[1]->SetName("effY04_1");
  effY04[2] = GetEffY04(2);
  effY04[2]->SetName("effY04_2");
  
  effY06 = GetEffY06();
  effY06->SetName("effY06");

  sysUn=0;

  //DEBUG:
  //std::cout<<std::endl;
  //std::cout<<"Default constructor of ktTrackEff ..."<<std::endl;
}

ktTrackEff::~ktTrackEff()
{
  delete effY06;
 
  delete effY04[0];
  delete effY04[1];
  delete effY04[2];

  delete effY07pteta[0];
  delete effY07pteta[1];
  delete effY07pteta[2];
  
  delete effY07eta[0];
  delete effY07eta[1];
  delete effY07eta[2];
 
  // funny way way not working via delete[] !???
  //delete[] effY04;
  //delete[] effY07pteta;
  //delete[] effY07eta;

  //DEBUG:
  //std::cout<<std::endl;
  //std::cout<<"Destructor of ktTrackEff ..."<<std::endl;
}

//void ktTrackEff::Clear()
//{
//}

void ktTrackEff::PrintInfo()
{
  std::cout<<"STAR Track Eff. Info:"<<std::endl;
  std::cout<<"---------------------"<<std::endl;
  std::cout<<" Run 4 --> Run 7 eff correction file = "<<fName<<std::endl;
  if (sysUn!=0)
    {
      std::cout<<" Sys Uncertainty = "<<sysUn<<std::endl;
      std::cout<<" Diff. eff. uncertainty AuAu = 4%"<<std::endl;
      std::cout<<" Diff. eff. uncertainty AuAu and pp (applied twice) = 3%"<<std::endl; 
    }
}


TF2* ktTrackEff::GetEffY06()
{

  TF2* funcpp = new TF2("ppEfficiency","[0]-0.06-[1]*exp([2]/x)+[3]*exp(-0.5*((x-[4])/[5])**2)/sqrt(2*pi*[5]*[5])-[6]*exp(-0.5*((x-[7])/[8])**2)/sqrt(2*pi*[8]*[8])+([9]-[10]*(y-[11])^2-[12]*(y-[11])^4-[13]*(y-[11])^6-[14]*(y-[11])^8)*exp(-[15]*x)",0.,10.,-1.,1.);

  Double_t parset[] = {0.869233,0.0223402,0.44061,0.558762,0.145162,0.508033,110.008,-4.63659,1.73765,0.0452674,-0.101279,0.0081551,0.945287,-2.00949,1.61746,1.39352};

  ((TF2*)funcpp)->SetParameters(parset);

  return funcpp;
}



TF2* ktTrackEff::GetEffY04(Int_t cb)
{
  Double_t fMaxPtPara = 5;

  char name[30];
  sprintf(name,"EfficiencyFunction%i",cb);

  TF2* func = new TF2(name,"[0]+[1]*x^2+[2]*x^4+[3]*x^6+[4]*x^8+[5]*exp([6]*y)+[7]*y+[8]*y*y +  [9]*exp(-((abs(x)-[10])^2)  /[11] - ((  abs(y)-[12]  )^2) /[13]) ",-1.,1.,0.,fMaxPtPara);

  Double_t parset4[]={ 0.698809, 0.0347652, -0.00819333, 0.112736, -0.356907, -1.62506, -7.26695, 0.0436162, -0.00453185, 0.249514, 0.308879, 0.133046, 0.295414, 0.0019349 };
  Double_t parset5[]={ 0.652242, 0.0760859, -0.0784171, 0.0393619, -0.247293, -2.04786, -8.96039, 0.0603416, -0.006971, -0.0435101, 0.131005, 0.00053132, 0.74369, 0.00576589 };
  Double_t parset6[]={ 0.631911, 0.117639, -0.29002, 0.522928, -0.569609, -1.3921, -6.73044, 0.0588101, -0.00686795, 0.110982, 0.2951, 0.14493, 0.295612, 0.00290843 };

  if(cb == 0)
    ((TF2*)func)->SetParameters(parset6);
  else if(cb == 1)
    ((TF2*)func)->SetParameters(parset5);
  else if(cb == 2)
    ((TF2*)func)->SetParameters(parset4);

  if(cb < 0 || cb > 2)
    std::cout << "Error: Nonsensical Centrality Class!" << std::endl;

  return func;
}

Double_t ktTrackEff::EffAAY07_20(Double_t eta, Double_t mPt)
{
  Double_t e1=EffAAY07(eta,mPt,0);
  Double_t e2=EffAAY07(eta,mPt,1);
  Double_t e3=EffAAY07(eta,mPt,2);
  
  return (e1*5+e2*5+e3*10)/20.;
}

Double_t ktTrackEff::EffAAY07(Double_t eta, Double_t mPt, Int_t centBin)
{
  Double_t effWeight=1.0;

  if(mPt < 5.)
    effWeight = effY04[centBin]->Eval(eta,mPt);
  else
    effWeight = effY04[centBin]->Eval(eta,5.0);
  if(mPt > 1.5)
    effWeight *= effY07eta[centBin]->GetBinContent(effY07eta[centBin]->GetXaxis()->FindBin(eta));
  else
    effWeight *= effY07pteta[centBin]->GetBinContent(effY07pteta[centBin]->GetXaxis()->FindBin(mPt),effY07pteta[centBin]->GetYaxis()->FindBin(eta));

  return effWeight;
}

Double_t ktTrackEff::EffPPY06(Double_t eta, Double_t mPt)
{
  Double_t effWeight=1.0;
  
  effWeight = effY06->Eval(mPt,eta);
  
  return effWeight;
}


Double_t ktTrackEff::EffRatio(Double_t eta, Double_t mPt, Int_t centBin)
{
  return EffAAY07(eta,mPt,centBin)/EffPPY06(eta,mPt);
}

Double_t ktTrackEff::EffRatio_20(Double_t eta, Double_t mPt)
{
  Double_t effRatio=EffAAY07_20(eta,mPt)/EffPPY06(eta,mPt);
  //std::cout<<effRatio<<std::endl;

  if(sysUn==0)
    {
      return effRatio;
    }
  else if (TMath::Abs(sysUn)==1)
    {
      Double_t sysRatio=EffRatio_20_Unc(eta,mPt);
      //std::cout<<sysRatio<<std::endl;
      Double_t effRatioSys=effRatio+sysUn*sysRatio;
      //std::cout<<effRatioSys<<std::endl;

      if (effRatioSys>1)
	{
	  return 1.0;
	}
      else
	{
	  return effRatioSys;
	}
    }
  else
    {
      std::cout<<" Error !!!"<<std::endl;
    }

  return effRatio;
}

Double_t ktTrackEff::EffRatio_20_Unc(Double_t eta, Double_t mPt)
{
  //Fix values according to Jet-Hadron analysis Note
  Double_t unAuAu=0.04;
  Double_t unEff=0.03;
  
  Double_t epp=EffPPY06(eta,mPt);
  Double_t eauau=EffAAY07_20(eta,mPt);

  Double_t factor=TMath::Sqrt(unEff*unEff/(epp*epp)+(unEff*unEff+unAuAu*unAuAu)/(eauau*eauau));
  Double_t effRatio=EffAAY07_20(eta,mPt)/EffPPY06(eta,mPt);
  
  return effRatio*factor;
}



