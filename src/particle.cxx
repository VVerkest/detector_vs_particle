//  Veronica Verkest          July, 2017
//  adapted from p6GeantTrackEfficiency.cxx (May 3, 2017)
//  analysis of jet&constituent number and pt

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/Selector.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include <iostream>
#include <sstream>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TH3.h>
#include <TClonesArray.h>
#include <TLatex.h>
#include <TMathText.h>
#include <Trandom.h>
#include <chrono>

#include "TStopwatch.h"
#include "ktTrackEff.hh"

// TStarJetPico
#include "TStarJetPicoReader.h"
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoEventCuts.h"
#include "TStarJetPicoPrimaryTrack.h"
#include "TStarJetPicoTower.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTowerCuts.h"
#include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"
#include "TStarJetPicoTriggerInfo.h"
#include "TStarJetPicoUtils.h"

using namespace fastjet;
using namespace std;

int main() {

  string outFileName = "out/particle_pt.root";
  
  TStopwatch TimeKeeper;  //Start a timer
  TimeKeeper.Start( );
  
  TH1::SetDefaultSumw2( );  // Histograms will calculate gaussian errors
  TH2::SetDefaultSumw2( );
  TH3::SetDefaultSumw2( );

  double pi = 3.14159265358979;
  double efficiency;         // dependent upon track eta and pt

  double R = 0.4;     //  JET DEFINITION
  JetDefinition jet_def(antikt_algorithm, R);

  TFile *fout = new TFile( (outFileName).c_str() ,"RECREATE");  // Create ONE output file
  TTree *P6EfficiencyJets = new TTree("P6EfficiencyJets","P6EfficiencyJets"); // Create a raw jet output tree
  TChain* chain = new TChain( "JetTreeMc" );
  chain->Add( "AddedGeantPythia/picoDst_25_35*" );

  string collisionType = "pp";
  TStarJetPicoReader reader;
  transform(collisionType.begin(), collisionType.end(), collisionType.begin(), ::tolower);
    
  reader.SetInputChain( chain );    // set the chain
  // apply hadronic correction - subtract 100% of charged track energy from towers
  reader.SetApplyFractionHadronicCorrection( false );
  //   reader.SetFractionHadronicCorrection( hadronicCorrectionFraction );
  reader.SetRejectTowerElectrons( kFALSE );
    
  TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();

  evCuts->SetRefMultCut( 0 );
    
  // Tracks cuts
  TStarJetPicoTrackCuts* trackCuts = reader.GetTrackCuts();
  trackCuts->SetDCACut( 100 );
  trackCuts->SetMinNFitPointsCut( -1 );
  trackCuts->SetFitOverMaxPointsCut( -1 );
  trackCuts->SetMaxPtCut ( 9999.0 );
    
  std::cout << "Using these track cuts:" << std::endl;
  std::cout << " dca : " << trackCuts->GetDCACut(  ) << std::endl;
  std::cout << " nfit : " <<   trackCuts->GetMinNFitPointsCut( ) << std::endl;
  std::cout << " nfitratio : " <<   trackCuts->GetFitOverMaxPointsCut( ) << std::endl;
  std::cout << " maxpt : " << trackCuts->GetMaxPtCut (  ) << std::endl;
    
  // Towers
  TStarJetPicoTowerCuts* towerCuts = reader.GetTowerCuts();
  towerCuts->SetMaxEtCut( 9999.0 );
  towerCuts->AddBadTowers( "src/Combined_y7_PP_Nick.txt" );

  std::cout << "Using these tower cuts:" << std::endl;
  std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;
  std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut() << std::endl;
    
  reader.SetProcessV0s(false);
  reader.Init( 400 ); //runs through all events with -1

  vector<PseudoJet> Jets;
  vector<PseudoJet> Constituents;
  
  //  P6EfficiencyJets->Branch("Jets", &Jets);

  // Data classes
  TStarJetVectorContainer<TStarJetVector>* container;
  TStarJetVector* sv;
  TStarJetPicoEventHeader* header;
  TStarJetPicoEvent* event;
  TClonesArray* towers;

  vector<PseudoJet> all_rawJets;

  int nEvents = 0;
  int nHardJets = 0;
  int nMatchedHard = 0;

  while ( reader.NextEvent() ) {       //  BEGIN EVENT LOOP!
    all_rawJets.clear();   //  clear all containers
    if ( nEvents == 0 ) { cout<< "Writing to:  " << outFileName <<endl; }
    
    nEvents++;
    reader.PrintStatus(10);      // Print out reader status every 10 seconds
    event = reader.GetEvent();           // Get the event header and event
    header = event->GetHeader();

    container = reader.GetOutputContainer();      // Get the output container from the reader

    double vertexZ = header->GetPrimaryVertexZ();      // Find vertex Z bin
    if ( vertexZ > 30 || vertexZ <= -30 ) { continue; }      //    |Vz| <= 30.0 cm

    TStarJetVector* sv;
    for ( int i=0; i < container->GetEntries() ; ++i ) {
      sv = container->Get(i);
      PseudoJet current = PseudoJet( *sv );
      current.set_user_index( sv->GetCharge() );

      //  All particles
      if ( sv->IsCharged() == false )  { continue; }  // removes neutral particles
      //      if ( current.pt() < 0.2)                { continue; }  // removes particles below 0.2 GeV
      if ( abs(current.eta()) > 1.0 )      { continue; }  // removes particles with eta>|1|

      
      //  pT-DEPENDENT EFFICIENCY

      ktTrackEff* tEff=new ktTrackEff();
      double eta = current.eta();
      double pT = current.pt();
      //   cout << tEff->EffPPY06(eta, pT) << endl;

      efficiency = tEff->EffPPY06(eta, pT);
      double effic_num = gRandom->Uniform(0.0,1.0);  // generate random number
      if ( effic_num > (1-efficiency) ) {
	all_rawJets.push_back(current);
      }
  
      ClusterSequence c_all_rawJets(all_rawJets, jet_def);          //  CLUSTER ALL JETS
      vector<PseudoJet> allRawJets = sorted_by_pt(c_all_rawJets.inclusive_jets());     // EXTRACT CLUSTERED JETS

      for (int i=0; i<allRawJets.size(); ++i) {
	Jets.push_back(allRawJets[i]);
	vector<PseudoJet> allrawjetcons = allRawJets[i].constituents();     //  Leading Jet Constituents
	for (int j = 0; j < allrawjetcons.size(); ++ j) {
	  if ( abs( allrawjetcons[j].eta() )<=(1.0) ) {
	    Constituents.push_back(allrawjetcons[j]);
	  }
	}
      }
    }
  }    //  END EVENT LOOP

  fout->Write();
  
  return 0;
}
