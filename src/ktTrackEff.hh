// first test (Joern Putschke)

#ifndef ROOT_ktTrackEff
#define ROOT_ktTrackEff

#include "TObject.h"
#include "TBuffer.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TString.h"

class ktTrackEff : public TObject
{

private:

  TString fName;

  //centrality bins: 0 is 0-5%, 1 is 5-10%, 2 is 10-20%
  TF2* effY04[3]; // Run 4 parameterization
  TH2D* effY07pteta[3]; // Run 7 HT / Run 4 MB pt-eta map, used for pt <= 1.5 GeV/c
  TH1D* effY07eta[3]; // Run 7 HT / Run 4 MB eta map, used for pt > 1.5 GeV/c
  TF2* effY06; // Run 6 parameterization


  Int_t sysUn;

  public:

  TF2* GetEffY06();
  TF2* GetEffY04(Int_t cb);
   
  ktTrackEff(TString mfName="./run7eff.root");
  virtual ~ktTrackEff();

  void SetSysUncertainty(Int_t mSys) {sysUn=mSys;}
  void PrintInfo();
  //void Clear();

  Double_t EffAAY07(Double_t eta, Double_t mPt, Int_t centBin=0);
  Double_t EffAAY07_20(Double_t eta, Double_t mPt);
  Double_t EffPPY06(Double_t eta, Double_t mPt);
  Double_t EffRatio(Double_t eta, Double_t mPt, Int_t centBin=0);
  Double_t EffRatio_20(Double_t eta, Double_t mPt);
  Double_t EffRatio_20_Unc(Double_t eta, Double_t mPt);
  //Double_t EffRatio_20_Unc();

  ClassDef(ktTrackEff,1)
};

#endif
