//  functions.hh
//  Veronica Verkest July 11, 2017
//  Original by Isaac Mooney on 5/4/17.

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include "TROOT.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TClonesArray.h"
#include "TLatex.h"
#include "TMathText.h"
#include "TProfile.h"

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

#include <iostream>
#include <sstream>

#ifndef functions_hh
#define functions_hh

namespace det_vs_part {
  class containers;
  
  void InitReaderPythia( TStarJetPicoReader & reader, TChain* chain, int nEvents );
  void InitReaderGeant( TStarJetPicoReader & reader, TChain* chain, int nEvents );

  // TLorentzVector ConvertPJtoTLV( std::vector<fastjet::PseudoJet> & PJ, TLorentzVector &TLV );
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~CONTAINER CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
  class containers {
  private:
  public:
    //trees
    TTree *efficiencyJets, *geantJets, *detVSpart;
        
    //branches
    double p_eventID, p_jetPt, p_jetEta, p_jetPhi, p_jetE;
    double g_eventID, g_jetPt, g_jetEta, g_jetPhi, g_jetE;
    double d_eventID, pt_diff, eta_diff, phi_diff, e_diff, delta_R;
    int p_nCons, g_nCons, nCons_diff;
    
    containers();               //constructor (initializes the trees, reserves space for vectors)
    virtual ~containers() {}    //destructor (shouldn't need to implement?)
        
    //getters & setters
    void SetBranches();         //setter
    void Clear();               //clears the particle vectors for each event
    void write();               //writes the trees to file
  };

}

#endif
