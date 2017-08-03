//  Veronica Verkest          July 11, 2017
//  Compare p6+efficiency to p6+GEANT
//  Functions in src/functions.cxx

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/Selector.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include <string>
#include <iostream>
#include <sstream>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLatex.h>
#include <Trandom.h>

#include "TStopwatch.h"
#include "ktTrackEff.hh"
#include "functions.hh"

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
using namespace det_vs_part;

int main () {

  string outFileName = "out/detector_vs_particle_analysis.root";
  
  TStopwatch TimeKeeper;  //Start a timer
  TimeKeeper.Start( );
  
  TH1::SetDefaultSumw2( );  // Histograms will calculate gaussian errors
  TH2::SetDefaultSumw2( );
  TH3::SetDefaultSumw2( );

  double numEvents = -1;       //  NUMBER OF EVENTS  (-1 runs all)
  double pi = 3.141592653589793238462643383;
  double efficiency;                 // dependent upon track eta and pt
  double R = 0.4;                     //  JET DEFINITION
  JetDefinition jet_def(antikt_algorithm, R);
  ktTrackEff* tEff=new ktTrackEff();
  
  TH2D* pTMatrix = new TH2D("pTMatrix","Detector p_{T} vs. Particle p_{T};x;y",200,0,100,200,0,100);   //  SET TITLES!
  TH2D* deltaPt_vs_deltaR = new TH2D("deltaPt_vs_deltaR","#Delta p_{T} vs. #Delta R",100,-40,40,50,-0.1,R);
  
  //  readers
  TStarJetPicoReader detReader;    //  Pythia+Geant
  TStarJetPicoReader partReader;   //  Pythia
  
  TChain* detChain = new TChain( "JetTree" );                    //  CORRESPONDING GEANT DATA  (detector)
  // detChain->Add( "AddedGeantPythia/picoDst_25_35*" );
  detChain->Add( "AddedGeantPythia/picoDst*" );
  TChain* partChain = new TChain( "JetTreeMc" );               //  PURE PYTHIA DATA  (particle)
  // partChain->Add( "AddedGeantPythia/picoDst_25_35*" );
  partChain->Add( "AddedGeantPythia/picoDst*" );
  
  containers * trees = new containers();        trees->SetBranches();          // Add branches to trees

  InitReaderGeant( detReader, detChain, numEvents );
  InitReaderPythia( partReader, partChain, numEvents );

  vector<PseudoJet> det_rawJets, part_rawJets, e_Jets, g_Jets, e_Cons, g_Cons;

  // Data classes
  TStarJetVectorContainer<TStarJetVector> * det_container;         TStarJetVector* det_sv;
  TStarJetPicoEventHeader* det_header;         TStarJetPicoEvent* det_event;
  TStarJetVectorContainer<TStarJetVector> * part_container;        TStarJetVector* part_sv;
  TStarJetPicoEventHeader* part_header;        TStarJetPicoEvent* part_event;
  

  int nEvents = 0;       int n_detJets = 0;       int n_partJets = 0;       int eventID = 0;       int detEventID;       int partEventID;

  while ( detReader.NextEvent() ) {       //  BEGIN EVENT LOOP!    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    detEventID = detReader.GetNOfCurrentEvent();

    if ( partReader.ReadEvent( detEventID ) != 1 ) continue;   //  ENSURES PART AND DET ON SAME EVENT!
    
    det_rawJets.clear(), part_rawJets.clear(), e_Jets.clear(), e_Cons.clear(), g_Jets.clear(), g_Cons.clear();   //  clear all containers

    nEvents++;    eventID++;
    detReader.PrintStatus(10);     partReader.PrintStatus(10);      // Print out reader status every 10 seconds
       
    det_event = detReader.GetEvent();    det_header = det_event->GetHeader();           // Get the event header and event
    det_container = detReader.GetOutputContainer();      // Get the output container from the reader
    part_event = partReader.GetEvent();    part_header = part_event->GetHeader();           // Get the event header and event
    part_container = partReader.GetOutputContainer();      // Get the output container from the reader

    double det_vertexZ = det_header->GetPrimaryVertexZ();      double part_vertexZ = part_header->GetPrimaryVertexZ();      // Find vertex Z
    if ( det_vertexZ > 30 || det_vertexZ <= -30 || part_vertexZ > 30 || part_vertexZ <= -30 ) { continue; }      //    |Vz| <= 30.0 cm

    //  GATHER CHARGED PYTHIA PARTICLES
    for ( int i=0; i < part_container->GetEntries() ; ++i ) {
      part_sv = part_container->Get(i);
      PseudoJet part_current = PseudoJet( *part_sv );
      part_current.set_user_index( part_sv->GetCharge() );

      //  All particles
      if ( part_sv->IsCharged() == false )  { continue; }  // removes neutral particles
      if ( abs(part_current.eta()) > 1.0 )      { continue; }  // removes particles with eta>|1|

      part_rawJets.push_back(part_current);
    }

      if ( part_rawJets.size() == 0 ) { continue; }
    
      //  GATHER ALL CHARGED GEANT PARTICLES
      for ( int i=0; i < det_container->GetEntries() ; ++i ) {
	det_sv = det_container->Get(i);
	PseudoJet det_current = PseudoJet( *det_sv );
	det_current.set_user_index( det_sv->GetCharge() );

	//  All particles
	if ( det_sv->IsCharged() == false )  { continue; }  // removes neutral particles
	if ( det_current.pt() < 0.2)                { continue; }  // removes particles below 0.2 GeV (Geant only!)
	if ( abs(det_current.eta()) > 1.0 )      { continue; }  // removes particles with eta>|1|

	det_rawJets.push_back(det_current);
      }

      //   cout << part_rawJets.size() << " particles for Particle Jets" << endl;
      //   cout << det_rawJets.size() << " particles for Detector Jets" << endl;
      
      if ( det_rawJets.size() == 0 ) { continue; }

      Selector etaSelector = SelectorAbsEtaMax( 1.0-R );
      Selector ptSelector = SelectorPtMin(2.0);
      Selector etaPtSelector = etaSelector && ptSelector;

      ClusterSequence geantCluster(det_rawJets, jet_def);           //  CLUSTER GEANT JETS WITHIN |ETA|<1
      vector<PseudoJet> detRawJets = sorted_by_pt( etaPtSelector(geantCluster.inclusive_jets()) );     // EXTRACT GEANT JETS
    
      ClusterSequence efficiencyCluster(part_rawJets, jet_def);    //  CLUSTER EFFICIENCY JETS WITHIN |ETA|<1
      vector<PseudoJet> partRawJets = sorted_by_pt( etaPtSelector(efficiencyCluster.inclusive_jets()) );     // EXTRACT EFFICIENCY JETS

      if ( detRawJets.size() == 0 || partRawJets.size() == 0 ) { continue; }

      Selector jetCircleMatcher = SelectorCircle( R );
      vector<PseudoJet> candidateJets, localJets, matchedJets, unmatchedJets, detectorJet, particleJet;
      double deltaR, deltaPt;
   
      particleJet.push_back( partRawJets[0] );
      
      jetCircleMatcher.set_reference( partRawJets[0] );
      localJets = sorted_by_pt( jetCircleMatcher( detRawJets ) );

      if ( localJets.size() == 0 ) { continue; }
      
      deltaR = partRawJets[0].delta_R( localJets[0] );     //  CALCULATE DELTA R
      detectorJet.push_back( localJets[0] );                         //  TAKE HIGHEST PT JET WITHIN R
    
      deltaPt = partRawJets[0].pt() - localJets[0].pt();

      g_Jets.push_back(detectorJet[0]);               //  Leading Geant (detector) Jet
      g_Cons = detectorJet[0].constituents();     //  Leading Geant Jet Constituents
      e_Jets.push_back( particleJet[0] );               //  Leading Efficiency (particle) Jet
      e_Cons = particleJet[0].constituents();     //  Leading Efficiency Jet Constituents


      
      trees->d_eventID = eventID;
      trees->g_eventID = eventID;                        trees->e_eventID = eventID;
      trees->g_jetPt = detectorJet[0].pt();              trees->e_jetPt = particleJet[0].pt();
      trees->g_jetEta = detectorJet[0].eta();           trees->e_jetEta = particleJet[0].eta();
      trees->g_jetPhi = detectorJet[0].phi_std();    trees->e_jetPhi = particleJet[0].phi_std();
      trees->g_jetE = detectorJet[0].E();                 trees->e_jetE = particleJet[0].E();
      trees->g_nCons = g_Cons.size();                 trees->e_nCons = e_Cons.size();

      trees->pt_diff = ( detectorJet[0].pt() - particleJet[0].pt() );
      trees->eta_diff = ( detectorJet[0].eta() - particleJet[0].eta() );
      trees->phi_diff = ( detectorJet[0].phi_std() - particleJet[0].phi_std() );
      trees->e_diff = ( detectorJet[0].E() - particleJet[0].E() );
      trees->nCons_diff = ( g_Cons.size() - e_Cons.size() );
      trees->delta_R = deltaR;

      pTMatrix->Fill( detectorJet[0].pt(), particleJet[0].pt() );
      deltaPt_vs_deltaR->Fill( deltaPt, deltaR );
    
      n_detJets ++;    n_partJets ++;

      trees->geantJets->Fill();    trees->efficiencyJets->Fill();    trees->detVSpart->Fill();

      //  cout << n_detJets << " detector Jets & " << n_partJets << " particle Jets" << endl << endl;
	
    }        //  END EVENT LOOP!    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    cout << endl << endl;
    cout << " Of " << nEvents << " events, " << n_detJets << " lead jet pairs have been found for the detector-vs.-particle analysis" << endl;
    cout<< "Writing to:  " << outFileName << endl << endl;

    TFile *fout = new TFile( (outFileName).c_str() ,"RECREATE");  // Create ONE output file
    pTMatrix->Write();
    deltaPt_vs_deltaR->Write();
    trees->write();
    fout->Close();

    return 0;
  }
