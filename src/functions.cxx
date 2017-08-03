//  functions.cxx
//  Veronica Verkest July 11, 2017
//  Originally by Isaac Mooney on 5/4/17.

#include "functions.hh"

namespace det_vs_part {

  //  CONVERT PSEUDOJET INTO TLORENTZ VECTOR!
  // TLorentzVector ConvertPJtoTLV( std::vector<fastjet::PseudoJet> & PJ, TLorentzVector &TLV ) {
  //   double pt = PJ.pt();
  //   double eta = PJ.eta();
  //   double phi = PJ.phi_std();
  //   double e = PJ.e();
  //   TLV->SetPtEtaPhiE ( pt, eta, phi, e );
  //   return TLV;
  //   }
  //////////////////////////////////////////////////////////////////////////////////////

  //  INITIATE PYTHIA READER
  void InitReaderPythia( TStarJetPicoReader & reader, TChain* chain, int nEvents ) {

    std::string collisionType = "pp";
      
    // First tolower() on the analysisType
    // shouldnt be necessary....
    std::transform(collisionType.begin(), collisionType.end(), collisionType.begin(), ::tolower);
    
    // set the chain
    reader.SetInputChain( chain );
    // apply hadronic correction - subtract 100% of charged track energy from towers
    reader.SetApplyFractionHadronicCorrection( true );
    reader.SetFractionHadronicCorrection( 0.9999 );
    reader.SetRejectTowerElectrons( kFALSE );
    
    // Event and track selection
    // -------------------------
    
    TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
    evCuts->SetVertexZCut ( 30 );
    evCuts->SetRefMultCut( 0 );
    evCuts->SetMaxEventPtCut( 30 );
    evCuts->SetMaxEventEtCut( 30 );

    // Tracks cuts
    TStarJetPicoTrackCuts* trackCuts = reader.GetTrackCuts();
    trackCuts->SetDCACut( 100 );
    trackCuts->SetMinNFitPointsCut( -1 );
    trackCuts->SetFitOverMaxPointsCut( -1 );    
    
    std::cout << "Using these track cuts:" << std::endl;
    std::cout << " dca : " << trackCuts->GetDCACut(  ) << std::endl;
    std::cout << " nfit : " <<   trackCuts->GetMinNFitPointsCut( ) << std::endl;
    std::cout << " nfitratio : " <<   trackCuts->GetFitOverMaxPointsCut( ) << std::endl;
    
    // Towers
    TStarJetPicoTowerCuts* towerCuts = reader.GetTowerCuts();
    towerCuts->SetMaxEtCut( 9999.0 );
    towerCuts->AddBadTowers( "src/dummy_tower_list.txt" );

    std::cout << "Using these tower cuts:" << std::endl;
    std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;
    std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut() << std::endl;
    
    // V0s: Turn off
    reader.SetProcessV0s(false);
    
    // Initialize the reader
    reader.Init( nEvents ); //runs through all events with -1
  }
  
  //////////////////////////////////////////////////////////////////////////////////////
 //  INITIATE GEANT READER
  void InitReaderGeant( TStarJetPicoReader & reader, TChain* chain, int nEvents ) {

    std::string collisionType = "pp";
      
    // First tolower() on the analysisType
    // shouldnt be necessary....
    std::transform(collisionType.begin(), collisionType.end(), collisionType.begin(), ::tolower);
    
    // set the chain
    reader.SetInputChain( chain );
    // apply hadronic correction - subtract 100% of charged track energy from towers
    reader.SetApplyFractionHadronicCorrection( true );
    reader.SetFractionHadronicCorrection( 0.9999 );
    reader.SetRejectTowerElectrons( kFALSE );
    
    // Event and track selection
    // -------------------------
    
    TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
    evCuts->SetVertexZCut ( 30 );
    evCuts->SetRefMultCut( 0 );
    evCuts->SetMaxEventPtCut( 30 );
    evCuts->SetMaxEventEtCut( 30 );

    // Tracks cuts
    TStarJetPicoTrackCuts* trackCuts = reader.GetTrackCuts();
    trackCuts->SetDCACut( 3.0 );
    trackCuts->SetMinNFitPointsCut( 20 );  // ~15-20
    trackCuts->SetFitOverMaxPointsCut( 0.52 );    //   Fitting to TPC
    
    std::cout << "Using these track cuts:" << std::endl;
    std::cout << " dca : " << trackCuts->GetDCACut(  ) << std::endl;
    std::cout << " nfit : " <<   trackCuts->GetMinNFitPointsCut( ) << std::endl;
    std::cout << " nfitratio : " <<   trackCuts->GetFitOverMaxPointsCut( ) << std::endl;
    
    // Towers
    TStarJetPicoTowerCuts* towerCuts = reader.GetTowerCuts();
    towerCuts->SetMaxEtCut( 9999.0 );
    towerCuts->AddBadTowers( "src/dummy_tower_list.txt" );

    std::cout << "Using these tower cuts:" << std::endl;
    std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;
    std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut() << std::endl;
    
    // V0s: Turn off
    reader.SetProcessV0s(false);
    
    // Initialize the reader
    reader.Init( nEvents ); //runs through all events with -1
  }
  //////////////////////////////////////////////////////////////////////////////////////


  //  IMPLEMENTATION FOR CLASS 'ANALYSIS::CONTAINERS'
  containers::containers() {
    //trees
    efficiencyJets   = new TTree("efficiencyJets","efficiencyJets");
    geantJets  = new TTree("geantJets","geantJets");
    detVSpart = new TTree("detVSpart","detVSpart");
  }
        
  //getters & setters
  void containers::SetBranches() {
    //branches
    efficiencyJets->Branch("p_eventID", &p_eventID);    efficiencyJets->Branch("p_jetPt", &p_jetPt);
    efficiencyJets->Branch("p_jetEta", &p_jetEta);           efficiencyJets->Branch("p_jetPhi", &p_jetPhi);
    efficiencyJets->Branch("p_jetE", &p_jetE);                 efficiencyJets->Branch("p_nCons", &p_nCons);
    geantJets->Branch("g_eventID", &g_eventID);          geantJets->Branch("g_jetPt", &g_jetPt);
    geantJets->Branch("g_jetEta", &g_jetEta);                 geantJets->Branch("g_jetPhi", &g_jetPhi);
    geantJets->Branch("g_jetE", &g_jetE);                       geantJets->Branch("g_nCons", &g_nCons);

    detVSpart->Branch("d_eventID", &d_eventID);
    detVSpart->Branch("pt_diff", &pt_diff);
    detVSpart->Branch("eta_diff", &eta_diff);
    detVSpart->Branch("phi_diff", &phi_diff);
    detVSpart->Branch("e_diff", &e_diff);
    detVSpart->Branch("nCons_diff", &nCons_diff);
    detVSpart->Branch("delta_R", &delta_R);
  }
    
  void containers::write() {
    efficiencyJets->Write();
    geantJets->Write();
    detVSpart->Write();
  }
  //////////////////////////////////////////////////////////////////////////////////////

}
