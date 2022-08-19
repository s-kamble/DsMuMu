// -*- C++ -*-
//
// Package:    DsMuMuRootupler
// Class:      DsMuMuRootupler
//
//====================================================
// Original author:  Samadhan Kamble                 |
//         created:  Sun, July 31, 2022 0300 hrs IST |
//         <samadhan.kamble@cern.ch>                 | 
//====================================================



// system include files
#include <memory>


// user include files
#include "myAnalyzers/DsMuMuPAT/src/DsMuMu.h"


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"


#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/Association.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/CLHEP/interface/Migration.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"


#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"


#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"

#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

//user include files: for kinematic fit
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"

#include "TrackingTools/IPTools/interface/IPTools.h"


#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"


#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"


#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"


#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h"


#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/PxPyPzE4D.h"


#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <vector>
#include <utility>
#include <string>
#include <iostream>
#include <map>



using namespace reco;
using namespace edm;
using namespace std;
using namespace pat;


//Define histograms
struct HistArgs{
  char name[128];
  char title[128];
  int n_bins;
  double x_min;
  double x_max;
};

enum HistName{
  h_dimumass,
  h_dimumass_vtxfit,
  h_dimu_vtxprob,
  h_dimuDca,
  h_Dsmass,
  h_Dspt,
/*  h_trkpt,
  h_trketa,
  h_trkmass,
  h_trkmass_vtxfit,*/
  nHistNameSize
};

HistArgs hist_args[nHistNameSize] = {
  // name, title, n_bins, x_min, x_max
  {"h_dimumass", "DiMuon Invariant Mass; m_{#mu^{+}#mu^{-}} [GeV]", 1000, 0, 100},
  {"h_dimumass_vtxfit", "DiMuon Invariant Mass after vertex fit; m_{#mu^{+}#mu^{-}} [GeV]", 1000, 0, 100},
  {"h_dimu_vtxprob", "h_dimu_vtxprob", 100, 0, 1.0},
  {"h_dimuDca", "Dimuon Distance of closest approach", 1000, 0, 10},
  {"h_Dsmass", "Ds Invariant Mass; m_{#k^{+}#k^{-}#pi^{-}} [GeV]", 1000, 0, 100},
  {"h_Dspt" , "Ds p_{T}; p_{T} [GeV]" , 1000 , 0 , 100}
//  {"h_trkpt" , "Track p_{T}; p_{T} [GeV]" , 1000 , 0 , 100},
//  {"h_trketa" , "Track #eta; #eta" , 1000 , -3. , 3.},
//  {"h_trkmass", "Three track invariant mass; m_{KK#pi} [GeV]", 1000, 0, 100},
//  {"h_trkmass_vtxfit", "Three track invariant mass after vertex fit; m_{KK#pi} [GeV]", 1000, 0, 100},
};

TH1D *histos[nHistNameSize];


//
// constants, enums and typedefs
//
 
typedef math::Error<3>::type CovarianceMatrix;

//
// static data member definitions
//
const double PI = 3.141592653589793;


//
//constructors and destructor
//
DsMuMu::DsMuMu(const edm::ParameterSet& iConfig):
  //labels 
  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  trackCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Track"))),
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),
/*  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),      
  triggerPrescales_ (consumes<pat::PackedTriggerPrescales> (iConfig.getParameter<edm::InputTag>("prescales"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),  
  trigTable_( iConfig.getParameter<std::vector<std::string> >("TriggerNames")),*/
  
  genParticles_(consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"))),
  packedGenToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter <edm::InputTag> ("packedGenParticles"))),


  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
  doMC_ ( iConfig.getUntrackedParameter<bool>("doMC",false) ),


  //particle properties
  MuonMass_(iConfig.getUntrackedParameter<double>("MuonMass")),
  MuonMassErr_(iConfig.getUntrackedParameter<double>("MuonMassErr")),
  KaonMass_(iConfig.getUntrackedParameter<double>("KaonMass")),
  KaonMassErr_(iConfig.getUntrackedParameter<double>("KaonMassErr")),
  PionMass_(iConfig.getUntrackedParameter<double>("PionMass")),
  PionMassErr_(iConfig.getUntrackedParameter<double>("PionMassErr")),

//  BcMass_(iConfig.getUntrackedParameter<double>("BcMass")),  
  
  //pre-selection cuts
  MuonMinPt_(iConfig.getUntrackedParameter<double>("MuonMinPt")),
  MuonMaxEta_(iConfig.getUntrackedParameter<double>("MuonMaxEta")),
  
/*  TrkMinPt_(iConfig.getUntrackedParameter<double>("TrkMinPt")),
  TrkMaxEta_(iConfig.getUntrackedParameter<double>("TrkMaxEta")),

  MuMuMinPt_(iConfig.getUntrackedParameter<double>("MuMuMinPt")),
  MuMuMinInvMass_(iConfig.getUntrackedParameter<double>("MuMuMinInvMass")),
  MuMuMaxInvMass_(iConfig.getUntrackedParameter<double>("MuMuMaxInvMass")),
  
  DsMinMass_(iConfig.getUntrackedParameter<double>("DsMinMass")),
  DsMaxMass_(iConfig.getUntrackedParameter<double>("DsMaxMass")),
  
  BcMinMass_(iConfig.getUntrackedParameter<double>("BcMinMass")),
  BcMaxMass_(iConfig.getUntrackedParameter<double>("BcMaxMass")),  */ 



  tree_(0),

  nBc(0), nMu(0),
  
  Bc_mass(0), Bc_px(0), Bc_py(0), Bc_pz(0), Bc_pt(0), Bc_eta(0), Bc_phi(0), Bc_Prob(0), Bc_charge(0),
  Bc_chi2(0), Bc_DecayVtxCL(0), Bc_DecayVtxX(0), Bc_DecayVtxY(0), Bc_DecayVtxZ(0),
  Bc_DecayVtxXE(0), Bc_DecayVtxYE(0), Bc_DecayVtxZE(0),
  Bc_DecayVtxXYE(0), Bc_DecayVtxXZE(0), Bc_DecayVtxYZE(0),


  Bc_Ds_mass(0), Bc_Ds_px(0), Bc_Ds_py(0), Bc_Ds_pz(0), Bc_Ds_pt(0), Bc_Ds_eta(0), Bc_Ds_phi(0), Bc_Ds_Prob(0), Bc_Ds_charge(0),
  Bc_Ds_chi2(0), Bc_Ds_DecayVtxCL(0), Bc_Ds_DecayVtxX(0), Bc_Ds_DecayVtxY(0), Bc_Ds_DecayVtxZ(0),
  Bc_Ds_DecayVtxXE(0), Bc_Ds_DecayVtxYE(0), Bc_Ds_DecayVtxZE(0),
  Bc_Ds_DecayVtxXYE(0), Bc_Ds_DecayVtxXZE(0), Bc_Ds_DecayVtxYZE(0),


  Bc_mumu_mass(0), Bc_mumu_px(0), Bc_mumu_py(0), Bc_mumu_pz(0), Bc_mumu_pt(0), Bc_mumu_eta(0), Bc_mumu_phi(0), Bc_mumu_Prob(0), Bc_mumu_charge(0),
  Bc_mumu_chi2(0), Bc_mumu_DecayVtxCL(0), Bc_mumu_DecayVtxX(0), Bc_mumu_DecayVtxY(0), Bc_mumu_DecayVtxZ(0),
  Bc_mumu_DecayVtxXE(0), Bc_mumu_DecayVtxYE(0), Bc_mumu_DecayVtxZE(0),
  Bc_mumu_DecayVtxXYE(0), Bc_mumu_DecayVtxXZE(0), Bc_mumu_DecayVtxYZE(0),


  Bc_Ds_pt1(0), Bc_Ds_px1(0), Bc_Ds_py1(0), Bc_Ds_pz1(0), Bc_Ds_charge1(0),
  Bc_Ds_pt2(0), Bc_Ds_px2(0), Bc_Ds_py2(0), Bc_Ds_pz2(0), Bc_Ds_charge2(0),
  Bc_Ds_pt3(0), Bc_Ds_px3(0), Bc_Ds_py3(0), Bc_Ds_pz3(0), Bc_Ds_charge3(0),

  Bc_Ds_px1_track(0), Bc_Ds_py1_track(0), Bc_Ds_pz1_track(0),
  Bc_Ds_px2_track(0), Bc_Ds_py2_track(0), Bc_Ds_pz2_track(0),
  Bc_Ds_px3_track(0), Bc_Ds_py3_track(0), Bc_Ds_pz3_track(0),

  Bc_mumu_pt1(0), Bc_mumu_px1(0), Bc_mumu_py1(0), Bc_mumu_pz1(0), Bc_mumu_charge1(0),
  Bc_mumu_pt2(0), Bc_mumu_px2(0), Bc_mumu_py2(0), Bc_mumu_pz2(0), Bc_mumu_charge2(0),


  nVtx(0), 
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
  priVtxXYE(0), priVtxXZE(0), priVtxYZE(0),

  //event information
  run(0), event(0), lumiblock(0),

  trk1px(0), trk1py(0), trk1pz(0), trk1chg(0),
  trk2px(0), trk2py(0), trk2pz(0), trk2chg(0)


  
{
  //now do what ever initialization is needed
}


DsMuMu::~DsMuMu()
{
  //do anything here that needs to be done at desctruction time
  //(e.g. close files, deallocate resources etc.)
}


//
//member functions
//

// ------------ method called for each event -------------
void DsMuMu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  event_counter_ += 1;
  
  //*********************************
  // Get event content information
  //*********************************
  
  //get magnetic field
  edm::ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  /////////////////////
  //  Kinematic fit
  /////////////////////
  edm::ESHandle<TransientTrackBuilder> theB; //theBuilder
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  edm::Handle<View<pat::PackedCandidate>> thePATTrackHandle;
  iEvent.getByToken(trackCollection_label, thePATTrackHandle);

  edm::Handle<View<pat::Muon>> thePATMuonHandle;
  iEvent.getByToken(dimuon_Label, thePATMuonHandle);
  
  std::cout<<"No. of muons (size of the muon collection): "<<thePATMuonHandle->size()<<endl;
  std::cout<<"No. of tracks (size of the track collection): "<<thePATTrackHandle->size()<<endl;
  
  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByToken(genParticles_, pruned);

  edm::Handle<pat::PackedGenParticleCollection> packed;
  iEvent.getByToken(packedGenToken_,packed);



  //Get the primary vertex

  reco::Vertex bestVtx;

  //getting the first primary vertex of the container
  edm::Handle<std::vector<reco::Vertex> > recVtxs;
  iEvent.getByToken(primaryVertices_Label, recVtxs);
  bestVtx = *(recVtxs->begin());

  priVtxX = bestVtx.x();
  priVtxY = bestVtx.y();
  priVtxZ = bestVtx.z();
  priVtxXE = bestVtx.covariance(0, 0);
  priVtxYE = bestVtx.covariance(1, 1);
  priVtxZE = bestVtx.covariance(2, 2);
  priVtxXYE = bestVtx.covariance(0, 1);
  priVtxXZE = bestVtx.covariance(0, 2);
  priVtxYZE = bestVtx.covariance(1, 2);
  priVtxCL = ChiSquaredProbability((double)(bestVtx.chi2()),(double)(bestVtx.ndof()));

  nVtx = recVtxs->size();  
 
  lumiblock = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event();


  //************************************************
  //let's begin by looking for dimuons (mu+ and mu-)
  //************************************************

  unsigned int nMu_tmp = thePATMuonHandle->size();
  
  for(edm::View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1)  //muon 1
    {
      for(edm::View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2)  //muon 2
        {
         //skip the pairing of muons with itself
         if(iMuon1==iMuon2) continue;
         
         //pair only opposite signed muons 
         if((iMuon1->charge())*(iMuon2->charge()) == 1) continue;
	  
	 //get tracks from the muons 
         TrackRef glbTrackP;	  
         TrackRef glbTrackM;	

         if(iMuon1->charge() == 1) {glbTrackP = iMuon1->track();}
         if(iMuon1->charge() == -1){glbTrackM = iMuon1->track();}

         if(iMuon2->charge() == 1) {glbTrackP = iMuon2->track();}
         if(iMuon2->charge() == -1){glbTrackM = iMuon2->track();}
   
         //check for the track reference
         if(glbTrackP.isNull() || glbTrackM.isNull()) 
           {
	    //std::cout << "continue due to no track reference" << endl;
	    continue;
         }

         //check for muon minimum-pt and muon max-eta
         if(iMuon1->track()->pt()<MuonMinPt_) continue;
         if(iMuon2->track()->pt()<MuonMinPt_) continue;
         if(fabs(iMuon1->eta())>MuonMaxEta_ || fabs(iMuon2->eta())>MuonMaxEta_) continue;
      
         if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
         if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;

         //let's check the vertex and mass
         reco::TransientTrack muon1TT((*theB).build(glbTrackP));
         reco::TransientTrack muon2TT((*theB).build(glbTrackM));

         //Trajectory states to calculate DCA for the 2 muons
         FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
         FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();
         if(!muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid()) continue;

         //measure distance between tracks at their closest approach
         ClosestApproachInRPhi cApp;
         cApp.calculate(mu1State, mu2State);
         if(!cApp.status()) continue; //continue --> bad status of closest approach
         float dca = fabs( cApp.distance() );
         histos[h_dimuDca]->Fill(dca);	  
         if(dca < 0. || dca > 0.5) continue;
         //cout<<" closest approach  "<<dca<<endl;


         TLorentzVector Mu1_4V, Mu2_4V, DiMu_4V;
      
         Mu1_4V.SetXYZM(iMuon1->px(),iMuon1->py(),iMuon1->pz(),MuonMass_);
         Mu2_4V.SetXYZM(iMuon2->px(),iMuon2->py(),iMuon2->pz(),MuonMass_);

         DiMu_4V = Mu1_4V + Mu2_4V;
      
         histos[h_dimumass]->Fill(DiMu_4V.M());

 
 
         // ####################################################
         // # Try to vertex the two muons to get dimuon vertex #
         // ####################################################
         

         //create a kinematic particle factory instance
         KinematicParticleFactoryFromTransientTrack pFactory;
         
         //initial chi2 and ndf before kinematic fits.
         float chi = 0.;
         float ndf = 0.;

         //make dimuon particles
         vector<RefCountedKinematicParticle> muonParticles;
 
         try{
             muonParticles.push_back(pFactory.particle(muon1TT,MuonMass_,chi,ndf,MuonMassErr_));
             muonParticles.push_back(pFactory.particle(muon2TT,MuonMass_,chi,ndf,MuonMassErr_));
         }
         catch(...) { 
             std::cout<<"Exception caught ... continuing 1 "<<std::endl; 
             continue;
         }
	 
         //create the vertex fitter instance
         KinematicParticleVertexFitter fitter;
   
         //create dimuon vertex fit tree
         RefCountedKinematicTree mumuVertexFitTree;
         
         try{
             mumuVertexFitTree = fitter.fit(muonParticles); 
         }
         catch (...) { 
             std::cout<<"Exception caught ... continuing 2 "<<std::endl; 
             continue;
	  }
	  
         //check validity of created fit tree
	 if(!mumuVertexFitTree->isValid()){ 
	    
	     std::cout << "Continue    ------------>      Invalid dimuon vertex fit tree" << std::endl;
	     continue; 
	 }
          
         //accessing the tree components
         mumuVertexFitTree->movePointerToTheTop();
         RefCountedKinematicParticle mumu_KP = mumuVertexFitTree->currentParticle();
         RefCountedKinematicVertex mumu_KV = mumuVertexFitTree->currentDecayVertex();
         
 
         //check validity of fitted vertex
         if(!mumu_KV->vertexIsValid()){
             continue;
         }
         
         //check chi2 validity for fitted vertex
         if(mumu_KV->chiSquared()<0 || mumu_KV->chiSquared()>50 ){
             std::cout << " continue --> negative chi2 for dimuon fit"<< endl;
             continue;
         }

         //calculate vertex fit probability and apply cut 
         double mumu_Prob_tmp  = TMath::Prob(mumu_KV->chiSquared(),(int)mumu_KV->degreesOfFreedom());
         histos[h_dimu_vtxprob]->Fill(mumu_Prob_tmp);
         if(mumu_Prob_tmp<0.1){
	    std::cout << "continue --> bad dimuon vtx fit prob" << endl;
	    continue;
         }
          
         histos[h_dimumass_vtxfit]->Fill(mumu_KP->currentState().mass());
          
          
          
         // #########################################################################
         // # Now we have the dimuon candidates. Let us look for the Ds candidates  #
         // #########################################################################
          

  /*        int it = 0;

          //create trklist with minimum 3 charged tracks 
	  vector <int> trklist;
	  trklist.reserve(2000);

	  for(it = 0; (unsigned)it<thePATTrackHandle->size(); it++)
             {
	      pat::PackedCandidate trackView((*thePATTrackHandle)[it]);
	      if(trackView.charge()==0) continue;  
	      
              //pT cut on each track
              //if(trackView.pt()<1.0) continue;  
	      if(trackView.pt()<0.5) continue; 
              if(!(trackView.hasTrackDetails())) continue;
	      if(!(trackView.bestTrack()->quality(reco::Track::highPurity))) continue; 
	      
              //exclude muon tracks
              if(IsTheSamePFtk(trackView,*iMuon1) || IsTheSamePFtk(trackView,*iMuon2)) continue;
	      trklist.push_back(it);
	     }
	  std::cout<<"no. of tracks (trklist size): "<<trklist.size()<<std::endl;
	 
          if(trklist.size()<3)continue;

          //looping on preselected trackcollection         
 	  for(unsigned int it1 = 0; it1<trklist.size(); it1++)
             {	  
	      for(unsigned int it2 = it1+1; it2<trklist.size(); it2++)
                 {	      
	          for(unsigned int it3 = it1+2; it3<trklist.size(); it3++)
                     {
	              if((it1 == it2) || (it1 == it3) || (it2 == it3)) continue;
	              
                      pat::PackedCandidate iTrack1((*thePATTrackHandle)[trklist.at(it1)]);
	  	      pat::PackedCandidate iTrack2((*thePATTrackHandle)[trklist.at(it2)]);
	  	      pat::PackedCandidate iTrack3((*thePATTrackHandle)[trklist.at(it3)]);
	  	     
//                      if(iTrack1.charge()==iTrack2.charge()) continue;
		     
                      //let's check the vertex and mass 
                      reco::TransientTrack pion1TT((*theB).build(iTrack1.bestTrack()));
		      reco::TransientTrack pion2TT((*theB).build(iTrack2.bestTrack()));
		      reco::TransientTrack pion3TT((*theB).build(iTrack3.bestTrack()));         
*/          

///////////////////////////////////////////////////////////////////////////////////////////////////////


        for(View<pat::PackedCandidate>::const_iterator iTrack1 = thePATTrackHandle->begin();  iTrack1 != thePATTrackHandle->end(); ++iTrack1 )	
           {
	    if(iTrack1->charge()==0) continue;		  
	    if(!(iTrack1->trackHighPurity())) continue;
	    if(iTrack1->pt()<0.5) continue;
            if(!(iTrack1->hasTrackDetails())) continue;
            
            for(View<pat::PackedCandidate>::const_iterator iTrack2 = iTrack1+1;  iTrack2 != thePATTrackHandle->end(); ++iTrack2 ) 
              {
	       if(iTrack1==iTrack2) continue;
	       if(iTrack2->charge()==0) continue;
	       if(!(iTrack2->trackHighPurity())) continue;
               if(iTrack2->pt()<0.5) continue;
               if(!(iTrack2->hasTrackDetails())) continue;
	     //  if(iTrack1->charge() == iTrack2->charge()) continue;
	     
               for(View<pat::PackedCandidate>::const_iterator iTrack3 = iTrack2+1;  iTrack3 != thePATTrackHandle->end(); ++iTrack3 )
                  {
		   if(iTrack3==iTrack2) continue;
		   if(iTrack3==iTrack1) continue; 
		   if(iTrack3->charge()==0) continue;
		   if(!(iTrack3->trackHighPurity())) continue;
		   if(iTrack3->pt()<0.5) continue; 
                   if(!(iTrack3->hasTrackDetails())) continue;
                  

		   //let's checks if our muons do not use the same tracks as we are using now
		   if ( IsTheSame(*iTrack1,*iMuon1) || IsTheSame(*iTrack1,*iMuon2) ) continue;
		   if ( IsTheSame(*iTrack2,*iMuon1) || IsTheSame(*iTrack2,*iMuon2) ) continue;
		   if ( IsTheSame(*iTrack3,*iMuon1) || IsTheSame(*iTrack3,*iMuon2) ) continue;
                  
   
                  //#####################################################
                  //Now let's see if these three tracks make a vertex	
                  //#####################################################

		  reco::TransientTrack pion1TT((*theB).build(iTrack1->pseudoTrack()));
		  reco::TransientTrack pion2TT((*theB).build(iTrack2->pseudoTrack()));
		  reco::TransientTrack pion3TT((*theB).build(iTrack3->pseudoTrack()));


//////////////////////////////////////////////////////////////////////////////////////////////////////////




                      TLorentzVector pion14V,pion24V,pion34V, Ds4V; 
	       
                      pion14V.SetXYZM(pion1TT.track().px(),pion1TT.track().py(),pion1TT.track().pz(),KaonMass_);
                      pion24V.SetXYZM(pion2TT.track().px(),pion2TT.track().py(),pion2TT.track().pz(),KaonMass_);
                      pion34V.SetXYZM(pion3TT.track().px(),pion3TT.track().py(),pion3TT.track().pz(),PionMass_);
                   
                      Ds4V=pion14V+pion24V+pion34V;
 
                      histos[h_Dsmass]->Fill(Ds4V.M());
                      histos[h_Dspt]->Fill(Ds4V.Pt());
                 
                      if(Ds4V.M()<1.40 || Ds4V.M()>2.5) continue;
                      if(Ds4V.Pt()<0)continue;
		
          //            std::cout << "\n" << __LINE__ << " : @@@ I have 3 good charged tracks. I'm trying to vertex them to Ds @@@" << std::endl;                    

                      //initial chi2 and ndf before kinematic fits.
                      float chi = 0.;
                      float ndf = 0.;
		
                      //make Ds particles
                      vector<RefCountedKinematicParticle> DsParticles;
                      
                      try {
                           DsParticles.push_back(pFactory.particle(pion1TT,KaonMass_,chi,ndf,KaonMassErr_));
                           DsParticles.push_back(pFactory.particle(pion2TT,KaonMass_,chi,ndf,KaonMassErr_));
                           DsParticles.push_back(pFactory.particle(pion3TT,PionMass_,chi,ndf,PionMassErr_));
			  }
                      catch(...) {
                           std::cout<<"Exception caught ... continuing 3 "<<std::endl;
                           continue;
                          }

                      RefCountedKinematicTree DsVertexFitTree;
                      KinematicParticleVertexFitter fitter;
                
                      //create vertex fit tree with Ds particles       
                      try {
                           DsVertexFitTree = fitter.fit(DsParticles);
                          }
                      catch(...) {
		           std::cout<<"Exception caught ... continuing 4 "<<std::endl;                   
		           continue;
	                  }
               
                      //check validity of vertex fit tree
                      if(!DsVertexFitTree->isValid()) 
                        {
//		         std::cout << "invalid Ds vertex fit tree" << std::endl;
		         continue; 
                        }
          
                      //access tree contents
                      DsVertexFitTree->movePointerToTheTop();	

                      RefCountedKinematicParticle Ds_KP = DsVertexFitTree->currentParticle();
                      RefCountedKinematicVertex Ds_KV = DsVertexFitTree->currentDecayVertex();

                      //check chi2 validity of fitted vertex
                      //if(Ds_KV->chiSquared()<0 || Ds_KV->chiSquared()>50 ){
                  //    if(Ds_KV->chiSquared()<0){// || Ds_KV->chiSquared()>50 ){
                    //     std::cout << " continue --> invalid chi2 for Ds fit: "<< Ds_KV->chiSquared() << endl;
                      //   continue;
                     // }
                    

                      //check validity of fitted vertex 
                      if (Ds_KV->vertexIsValid() == false)
                        {
                         std::cout << __LINE__ << " : continue --> invalid vertex from the Ds vertex fit" << std::endl;
                         continue;
                        }

                      //Ds mass cut    
                      if(Ds_KP->currentState().mass()<1.90 || Ds_KP->currentState().mass()>2.04) continue;
          
                      double Ds_Prob_tmp  = TMath::Prob((double)Ds_KV->chiSquared(),(int)Ds_KV->degreesOfFreedom());
                      if(Ds_Prob_tmp<0.01){
                     //    std::cout << "continue --> bad vtx prob from Ds fit" << endl;
	                 continue;
                      }
         

                      DsVertexFitTree->movePointerToTheFirstChild();
                      RefCountedKinematicParticle T1CandMC = DsVertexFitTree->currentParticle();

                      DsVertexFitTree->movePointerToTheNextChild();
                      RefCountedKinematicParticle T2CandMC = DsVertexFitTree->currentParticle();

                      DsVertexFitTree->movePointerToTheNextChild();
                      RefCountedKinematicParticle T3CandMC = DsVertexFitTree->currentParticle(); 
                    
                      DsVertexFitTree->movePointerToTheTop();
                      RefCountedKinematicParticle Ds_KP_withMC = DsVertexFitTree->currentParticle();
          


                    // #####################################
                    // # Bc vertex fit with vtx constraint #
                    // #####################################
 

                      //make Bc particles
                      vector<RefCountedKinematicParticle> BcParticles;
                      BcParticles.push_back(pFactory.particle(muon1TT,MuonMass_,chi,ndf,MuonMassErr_));
                      BcParticles.push_back(pFactory.particle(muon2TT,MuonMass_,chi,ndf,MuonMassErr_));
                      BcParticles.push_back(Ds_KP_withMC);
                     //      BcParticles.push_back(pFactory.particle(pion1TT,KaonMass_,chi,ndf,KaonMassErr_));
                       //    BcParticles.push_back(pFactory.particle(pion2TT,KaonMass_,chi,ndf,KaonMassErr_));
                         //  BcParticles.push_back(pFactory.particle(pion3TT,PionMass_,chi,ndf,PionMassErr_));
 
                      KinematicConstrainedVertexFitter kcvFitter;
                      RefCountedKinematicTree BcVertexFitTree = kcvFitter.fit(BcParticles);
                
                      if(!BcVertexFitTree->isValid()) {
            //             std::cout << "Bc vertex fit tree is not valid" << endl;
                         continue;
                      }
 
                      BcVertexFitTree->movePointerToTheTop();
 		
                      RefCountedKinematicParticle Bc_KP = BcVertexFitTree->currentParticle();	
                      RefCountedKinematicVertex Bc_KV = BcVertexFitTree->currentDecayVertex();
 
                      if(!Bc_KV->vertexIsValid()){
                         std::cout << "Bc MC fit vertex is not valid" << endl;
                         continue;
                      }
		
                      //Bc mass cut
                      if(Bc_KP->currentState().mass()<5.7 || Bc_KP->currentState().mass()>6.8) continue;
		 
                      if(Bc_KV->chiSquared()<0 || Bc_KV->chiSquared()>50 ) 
                       {
	//	        std::cout << " continue from invalid chi2 for Bc_KV = " << Bc_KV->chiSquared() << endl;
		        continue;
                      }
		
                      double Bc_Prob_tmp  = TMath::Prob(Bc_KV->chiSquared(),(int)Bc_KV->degreesOfFreedom());
                      if(Bc_Prob_tmp<0.01)
                      {
		       continue;
                      }

                      std::cout << "Bc mass "<<Bc_KP->currentState().mass()<< std::endl;
 
                      //get childrens from the final Bc fit
                      BcVertexFitTree->movePointerToTheFirstChild();
                      RefCountedKinematicParticle mu1Cand = BcVertexFitTree->currentParticle();
                    
                      BcVertexFitTree->movePointerToTheNextChild();
                      RefCountedKinematicParticle mu2Cand = BcVertexFitTree->currentParticle();

                      BcVertexFitTree->movePointerToTheNextChild();
                      RefCountedKinematicParticle DsCand = BcVertexFitTree->currentParticle();


	              //---------------------------- End of Vertexing -----------------------------------
 
                      //###############################################                     
                      //#      Fill candidate variables now           #
	              //###############################################	
                    
                      if(nBc==0){
		         nMu  = nMu_tmp;
                      }

                      //Bc candidate variables
                      Bc_mass->push_back(Bc_KP->currentState().mass());           
                      Bc_px->push_back(Bc_KP->currentState().globalMomentum().x());
                      Bc_py->push_back(Bc_KP->currentState().globalMomentum().y());
                      Bc_pz->push_back(Bc_KP->currentState().globalMomentum().z());
                      Bc_eta->push_back(Bc_KP->currentState().globalMomentum().eta());
                      Bc_phi->push_back(Bc_KP->currentState().globalMomentum().phi());	
                      Bc_pt->push_back(Bc_KP->currentState().globalMomentum().perp());
         

//		      Bc_Prob->push_back(Bc_Prob_tmp);
                      Bc_chi2->push_back(Bc_KV->chiSquared());
                      Bc_DecayVtxCL->push_back(TMath::Prob(static_cast<double>(Bc_KV->chiSquared()), static_cast<int>(rint(Bc_KV->degreesOfFreedom()))));
		      Bc_DecayVtxX->push_back(Bc_KV->position().x());
                      Bc_DecayVtxY->push_back(Bc_KV->position().y());
                      Bc_DecayVtxZ->push_back(Bc_KV->position().z());
                      Bc_DecayVtxXE->push_back(Bc_KV->error().cxx());
                      Bc_DecayVtxYE->push_back(Bc_KV->error().cyy());
                      Bc_DecayVtxZE->push_back(Bc_KV->error().czz());
		      Bc_DecayVtxXYE->push_back(Bc_KV->error().cyx());
		      Bc_DecayVtxXZE->push_back(Bc_KV->error().czx());
		      Bc_DecayVtxYZE->push_back(Bc_KV->error().czy());

                      
                      //Ds candidate variables
                      Bc_Ds_mass->push_back(Ds_KP->currentState().mass() );
                      Bc_Ds_px->push_back(Ds_KP->currentState().globalMomentum().x());
                      Bc_Ds_py->push_back(Ds_KP->currentState().globalMomentum().x());
                      Bc_Ds_pz->push_back(Ds_KP->currentState().globalMomentum().x());
                      Bc_Ds_eta->push_back(Ds_KP->currentState().globalMomentum().eta());
                      Bc_Ds_phi->push_back(Ds_KP->currentState().globalMomentum().phi());
                      Bc_Ds_pt->push_back(Ds_KP->currentState().globalMomentum().perp());

	     //         Bc_Ds_Prob->push_back(Ds_Prob_tmp);
		      Bc_Ds_chi2->push_back(Ds_KV->chiSquared());
                      Bc_Ds_DecayVtxCL->push_back(TMath::Prob(static_cast<double>(Ds_KV->chiSquared()), static_cast<int>(rint(Ds_KV->degreesOfFreedom()))));
          	      Bc_Ds_DecayVtxX->push_back(Ds_KV->position().x());
                      Bc_Ds_DecayVtxY->push_back(Ds_KV->position().y());
                      Bc_Ds_DecayVtxZ->push_back(Ds_KV->position().z());
                      Bc_Ds_DecayVtxXE->push_back(Ds_KV->error().cxx());
                      Bc_Ds_DecayVtxYE->push_back(Ds_KV->error().cyy());
                      Bc_Ds_DecayVtxZE->push_back(Ds_KV->error().czz());
		      Bc_Ds_DecayVtxXYE->push_back(Ds_KV->error().cyx());
		      Bc_Ds_DecayVtxXZE->push_back(Ds_KV->error().czx());
		      Bc_Ds_DecayVtxYZE->push_back(Ds_KV->error().czy());
                  

                      //Dimuon candidate variables
                      Bc_mumu_mass->push_back(mumu_KP->currentState().mass());
                      Bc_mumu_px->push_back(mumu_KP->currentState().globalMomentum().x());
                      Bc_mumu_py->push_back(mumu_KP->currentState().globalMomentum().y());
                      Bc_mumu_pz->push_back(mumu_KP->currentState().globalMomentum().z());
                      Bc_mumu_eta->push_back(mumu_KP->currentState().globalMomentum().eta());
                      Bc_mumu_phi->push_back(mumu_KP->currentState().globalMomentum().phi());
                      Bc_mumu_pt->push_back(mumu_KP->currentState().globalMomentum().perp());

//		      Bc_mumu_Prob->push_back(mumu_Prob_tmp);
		      Bc_mumu_chi2->push_back(mumu_KV->chiSquared());
                      Bc_mumu_DecayVtxCL->push_back(TMath::Prob(static_cast<double>(mumu_KV->chiSquared()), static_cast<int>(rint(mumu_KV->degreesOfFreedom()))));
         	      Bc_mumu_DecayVtxX->push_back(mumu_KV->position().x());
                      Bc_mumu_DecayVtxY->push_back(mumu_KV->position().y());
                      Bc_mumu_DecayVtxZ->push_back(mumu_KV->position().z());
                      Bc_mumu_DecayVtxXE->push_back(mumu_KV->error().cxx());
                      Bc_mumu_DecayVtxYE->push_back(mumu_KV->error().cyy());
                      Bc_mumu_DecayVtxZE->push_back(mumu_KV->error().czz());
		      Bc_mumu_DecayVtxXYE->push_back(mumu_KV->error().cyx());
		      Bc_mumu_DecayVtxXZE->push_back(mumu_KV->error().czx());
		      Bc_mumu_DecayVtxYZE->push_back(mumu_KV->error().czy());
                      


                      //Ds childrens' variables
                      Bc_Ds_px1_track->push_back(iTrack1->px());
                      Bc_Ds_py1_track->push_back(iTrack1->py());
                      Bc_Ds_pz1_track->push_back(iTrack1->pz());


                      Bc_Ds_px2_track->push_back(iTrack2->px());
                      Bc_Ds_py2_track->push_back(iTrack2->py());
                      Bc_Ds_pz2_track->push_back(iTrack2->pz());


                      Bc_Ds_px3_track->push_back(iTrack3->px());
                      Bc_Ds_py3_track->push_back(iTrack3->py());
                      Bc_Ds_pz3_track->push_back(iTrack3->pz());



                      nBc++;
                      muonParticles.clear();
                      DsParticles.clear();
                      BcParticles.clear();



                      std::cout << "\n" << __LINE__ << " :I am here" << std::endl;                    

                     }//closing track3 loop

//                      std::cout << "\n" << __LINE__ << " :closed track3 loop" << std::endl;                    
                 }//closing track2 loop

//                    std::cout << "\n" << __LINE__ << " :closed track2 loop" << std::endl;                    
             }//closing track1 loop
//             trklist.clear();

                      std::cout << "\n" << __LINE__ << " :closed track loops" << std::endl;                    

        }//closing muon2 loop

                 //     std::cout << "\n" << __LINE__ << " :closed muon2 loop" << std::endl;                    
    }//closing muon1 loop

                      std::cout << "\n" << __LINE__ << " :closed muon loops" << std::endl;                    



  //fill the tree and clear the vectors

  if (nBc > 0 ) 
    {
      std::cout << "filling tree, nBc = " << nBc<<std::endl;
      tree_->Fill();
    }


//  tree_->Fill();

  nBc = 0; nMu = 0;


  Bc_mass->clear();          Bc_px->clear();           Bc_py->clear();           Bc_pz->clear(); 
  Bc_eta->clear();           Bc_phi->clear();          Bc_pt->clear();           Bc_chi2->clear();
  Bc_DecayVtxCL->clear();    Bc_DecayVtxX->clear();    Bc_DecayVtxY->clear();    Bc_DecayVtxZ->clear();
  Bc_DecayVtxXE->clear();    Bc_DecayVtxYE->clear();   Bc_DecayVtxZE->clear();
  Bc_DecayVtxXYE->clear();   Bc_DecayVtxXZE->clear();  Bc_DecayVtxYZE->clear();

  Bc_Ds_mass->clear();        Bc_Ds_px->clear();          Bc_Ds_py->clear();          Bc_Ds_pz->clear();
  Bc_Ds_eta->clear();         Bc_Ds_phi->clear();         Bc_Ds_pt->clear();          Bc_Ds_chi2->clear();
  Bc_Ds_DecayVtxCL->clear();  Bc_Ds_DecayVtxX->clear();   Bc_Ds_DecayVtxY->clear();   Bc_Ds_DecayVtxZ->clear();
  Bc_Ds_DecayVtxXE->clear();  Bc_Ds_DecayVtxYE->clear();  Bc_Ds_DecayVtxZE->clear();
  Bc_Ds_DecayVtxXYE->clear(); Bc_Ds_DecayVtxXZE->clear(); Bc_Ds_DecayVtxYZE->clear();

  Bc_mumu_mass->clear();         Bc_mumu_px->clear();          Bc_mumu_py->clear();          Bc_mumu_pz->clear();  
  Bc_mumu_eta->clear();          Bc_mumu_phi->clear();         Bc_mumu_pt->clear();          Bc_mumu_chi2->clear();
  Bc_mumu_DecayVtxCL->clear();   Bc_mumu_DecayVtxX->clear();   Bc_mumu_DecayVtxY->clear();   Bc_mumu_DecayVtxZ->clear();
  Bc_mumu_DecayVtxXE->clear();   Bc_mumu_DecayVtxYE->clear();  Bc_mumu_DecayVtxZE->clear();
  Bc_mumu_DecayVtxXYE->clear();  Bc_mumu_DecayVtxXZE->clear(); Bc_mumu_DecayVtxYZE->clear();

  Bc_Ds_px1_track->clear();   Bc_Ds_py1_track->clear();  Bc_Ds_pz1_track->clear();
  Bc_Ds_px2_track->clear();   Bc_Ds_py2_track->clear();  Bc_Ds_pz2_track->clear();
  Bc_Ds_px3_track->clear();   Bc_Ds_py2_track->clear();  Bc_Ds_pz3_track->clear();




}	  


bool DsMuMu::IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

/*

bool DsMuMu::IsTheSame(const reco::Track& tk, const pat::Muon& mu)
{
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}
*/

bool DsMuMu::IsTheSamePFtk(const pat::PackedCandidate& tk, const pat::Muon& mu){
//  double DeltaEta = fabs(mu.eta()-tk.eta());
//  double DeltaP   = fabs(mu.p()-tk.p());
  if(deltaR(tk.eta(),tk.phi(),mu.eta(),mu.phi()) < 0.0001 )return true;
  //if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}



// ------------ method called once each job just before starting event loop  ------------

void DsMuMu::beginJob()
{
  event_counter_ = 0;
  
  edm::Service<TFileService> fs;
  
  tree_ = fs->make<TTree>("BcToDsMuMuTree", "Bc->DsMuMu Tree");
 
  tree_->Branch("nBc",&nBc,"nBc/i"); 

  tree_->Branch("Bc_mass", &Bc_mass); 
  tree_->Branch("Bc_px", &Bc_px);
  tree_->Branch("Bc_py", &Bc_py);
  tree_->Branch("Bc_pz", &Bc_pz);
  tree_->Branch("Bc_eta", &Bc_eta); 
  tree_->Branch("Bc_phi", &Bc_phi);
  tree_->Branch("Bc_pt", &Bc_pt);
  tree_->Branch("Bc_chi2", &Bc_chi2);
  tree_->Branch("Bc_DecayVtxCL", &Bc_DecayVtxCL);
  tree_->Branch("Bc_DecayVtxX", &Bc_DecayVtxX);
  tree_->Branch("Bc_DecayVtxY", &Bc_DecayVtxY);
  tree_->Branch("Bc_DecayVtxZ", &Bc_DecayVtxZ);
  tree_->Branch("Bc_DecayVtxXE", &Bc_DecayVtxXE);
  tree_->Branch("Bc_DecayVtxYE", &Bc_DecayVtxYE);
  tree_->Branch("Bc_DecayVtxZE", &Bc_DecayVtxZE);
  tree_->Branch("Bc_DecayVtxXYE", &Bc_DecayVtxXYE);
  tree_->Branch("Bc_DecayVtxXZE", &Bc_DecayVtxXZE);
  tree_->Branch("Bc_DecayVtxYZE", &Bc_DecayVtxYZE);



  tree_->Branch("Bc_Ds_mass", &Bc_Ds_mass);
  tree_->Branch("Bc_Ds_px", &Bc_Ds_px);
  tree_->Branch("Bc_Ds_py", &Bc_Ds_py);
  tree_->Branch("Bc_Ds_pz", &Bc_Ds_pz);
  tree_->Branch("Bc_Ds_eta", &Bc_Ds_eta); 
  tree_->Branch("Bc_Ds_phi", &Bc_Ds_phi);
  tree_->Branch("Bc_Ds_pt", &Bc_Ds_pt);
  tree_->Branch("Bc_Ds_chi2", &Bc_Ds_chi2);
  tree_->Branch("Bc_Ds_DecayVtxCL", &Bc_Ds_DecayVtxCL);
  tree_->Branch("Bc_Ds_DecayVtxX", &Bc_Ds_DecayVtxX);
  tree_->Branch("Bc_Ds_DecayVtxY", &Bc_Ds_DecayVtxY);
  tree_->Branch("Bc_Ds_DecayVtxZ", &Bc_Ds_DecayVtxZ);
  tree_->Branch("Bc_Ds_DecayVtxXE", &Bc_Ds_DecayVtxXE);
  tree_->Branch("Bc_Ds_DecayVtxYE", &Bc_Ds_DecayVtxYE);
  tree_->Branch("Bc_Ds_DecayVtxZE", &Bc_Ds_DecayVtxZE);
  tree_->Branch("Bc_Ds_DecayVtxXYE", &Bc_Ds_DecayVtxXYE);
  tree_->Branch("Bc_Ds_DecayVtxXZE", &Bc_Ds_DecayVtxXZE);
  tree_->Branch("Bc_Ds_DecayVtxYZE", &Bc_Ds_DecayVtxYZE);



  tree_->Branch("Bc_mumu_mass", &Bc_mumu_mass);
  tree_->Branch("Bc_mumu_px", &Bc_mumu_px);
  tree_->Branch("Bc_mumu_py", &Bc_mumu_py);
  tree_->Branch("Bc_mumu_pz", &Bc_mumu_pz);
  tree_->Branch("Bc_mumu_eta", &Bc_mumu_eta); 
  tree_->Branch("Bc_mumu_phi", &Bc_mumu_phi);
  tree_->Branch("Bc_mumu_pt", &Bc_mumu_pt);
  tree_->Branch("Bc_mumu_chi2", &Bc_mumu_chi2);
  tree_->Branch("Bc_mumu_DecayVtxCL", &Bc_mumu_DecayVtxCL);
  tree_->Branch("Bc_mumu_DecayVtxX", &Bc_mumu_DecayVtxX);
  tree_->Branch("Bc_mumu_DecayVtxY", &Bc_mumu_DecayVtxY);
  tree_->Branch("Bc_mumu_DecayVtxZ", &Bc_mumu_DecayVtxZ);
  tree_->Branch("Bc_mumu_DecayVtxXE", &Bc_mumu_DecayVtxXE);
  tree_->Branch("Bc_mumu_DecayVtxYE", &Bc_mumu_DecayVtxYE);
  tree_->Branch("Bc_mumu_DecayVtxZE", &Bc_mumu_DecayVtxZE);
  tree_->Branch("Bc_mumu_DecayVtxXYE", &Bc_mumu_DecayVtxXYE);
  tree_->Branch("Bc_mumu_DecayVtxXZE", &Bc_mumu_DecayVtxXZE);
  tree_->Branch("Bc_mumu_DecayVtxYZE", &Bc_mumu_DecayVtxYZE);
  



  tree_->Branch("Bc_Ds_px1_track", &Bc_Ds_px1_track);
  tree_->Branch("Bc_Ds_py1_track", &Bc_Ds_py1_track);
  tree_->Branch("Bc_Ds_pz1_track", &Bc_Ds_pz1_track);


  tree_->Branch("Bc_Ds_px2_track", &Bc_Ds_px2_track);
  tree_->Branch("Bc_Ds_py2_track", &Bc_Ds_py2_track);
  tree_->Branch("Bc_Ds_pz2_track", &Bc_Ds_pz2_track);


  tree_->Branch("Bc_Ds_px3_track", &Bc_Ds_px3_track);
  tree_->Branch("Bc_Ds_py3_track", &Bc_Ds_py3_track);
  tree_->Branch("Bc_Ds_pz3_track", &Bc_Ds_pz3_track);





  tree_->Branch("run",        &run,       "run/I");
  tree_->Branch("event",      &event,     "event/I");
  tree_->Branch("lumiblock",  &lumiblock, "lumiblock/I");
  

  for(int i=0; i<nHistNameSize; i++) {
    histos[i] = fs->make<TH1D>(hist_args[i].name, hist_args[i].title, hist_args[i].n_bins,hist_args[i].x_min, hist_args[i].x_max);
  }

}


// ------------ method called once each job just after ending the event loop  ------------
void DsMuMu::endJob() 
{
  tree_->GetDirectory()->cd();
  tree_->Write();
  
  printf(" \n ---------- End Job ---------- \n" ) ;
  
  cout << "Total number of Events processed : " << event_counter_ << endl;
  
}


//define this as a plug-in
DEFINE_FWK_MODULE(DsMuMu);

