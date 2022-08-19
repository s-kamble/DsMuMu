#ifndef _DsMuMu_h
#define _DsMuMu_h

// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"


#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"


#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <utility>
#include <string>


#include "DataFormats/Math/interface/Vector3D.h"

//
// struct decleration
//


//class definition
class DsMuMu : public edm::EDAnalyzer {
public:
  explicit DsMuMu(const edm::ParameterSet&);
  ~DsMuMu();
  
  bool IsTheSametk(const pat::GenericParticle& tk, const pat::Muon& mu);
  bool IsTheSamePFtk(const pat::PackedCandidate& tk, const pat::Muon& mu);
//  bool IsTheSame(const reco::Track& tk, const pat::Muon& mu);
  bool IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu);


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void printout(const RefCountedKinematicVertex& myVertex) const;
  void printout(const RefCountedKinematicParticle& myParticle) const;
  void printout(const RefCountedKinematicTree& myTree) const;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<pat::Muon>> dimuon_Label;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trackCollection_label;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<reco::BeamSpot> BSLabel_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  //edm::EDGetTokenT<pat::PackedTriggerPrescales>            triggerPrescales_;
  //edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;

  //std::vector<std::string> trigTable_;
  //edm::EDGetTokenT<std::vector< PileupSummaryInfo>> puToken_;
  //edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> v0PtrCollection_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticles_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGenToken_;


  bool OnlyBest_;
  bool isMC_;
  bool OnlyGen_;
  bool doMC_;


  //particle properties 
  ParticleMass MuonMass_;
  float MuonMassErr_;
  ParticleMass KaonMass_;
  float KaonMassErr_;
  ParticleMass PionMass_;
  float PionMassErr_;
  ParticleMass DsMass_;
  float DsMassErr_;
  //double BcMass_;

  //pre-selection cuts
  double MuonMinPt_;
  double MuonMaxEta_;
/*  double TrkMinPt_;
  double TrkMaxEta_;
  double MuMuMinPt_;
  double MuMuMinInvMass_;
  double MuMuMaxInvMass_;
  double DsMinMass_;
  double DsMaxMass_;
  double BcMinMass_;
  double BcMaxMass_;*/


  const MagneticField   *fMagneticField;

  vector<string> TriggerNames_;

  TTree*   tree_;
  int mupCategory;
  int mumCategory;
  int mupME1Clean;
  int mumME1Clean;


  unsigned int             nBc;  
  unsigned int             nMu;
 

  //dimuon
  std::vector<float>       *mumC2;
  std::vector<int>         *mumNHits, *mumNPHits; 
  std::vector<float>       *mupC2;
  std::vector<int>         *mupNHits, *mupNPHits;
  std::vector<float>       *mumdxy, *mupdxy, *mumdz, *mupdz;
  std::vector<float>       *muon_dca;
  std::vector<double>      *mumumass, *mumumasserr;
  std::vector<bool>        *mumisgoodmuon, *mupisgoodmuon ;
  std::vector<int>         *mumcharge, *mupcharge ;
  std::vector<bool>        *mumloosemuon, *muploosemuon ; 
  
  std::vector<int>         *tri_Dim25, *tri_JpsiTk, *tri_JpsiTkTk;
  std::vector<float>       *dr0, *dr1, *dpt0, *dpt1;
  std::vector<bool>        *mu1soft, *mu2soft, *mu1tight, *mu2tight;  
  std::vector<bool>        *mu1PF, *mu2PF, *mu1loose, *mu2loose;  

  int                      muAcc, muTrig, weight;
 
  //**********************************

  //Bc candidate variables 
  std::vector<float>       *Bc_mass, *Bc_px, *Bc_py, *Bc_pz, *Bc_pt, *Bc_eta, *Bc_phi, *Bc_Prob;
  std::vector<short>       *Bc_charge;
  std::vector<float>       *Bc_chi2, *Bc_DecayVtxCL, *Bc_DecayVtxX, *Bc_DecayVtxY, *Bc_DecayVtxZ;
  std::vector<float>       *Bc_DecayVtxXE, *Bc_DecayVtxYE, *Bc_DecayVtxZE;
  std::vector<float>       *Bc_DecayVtxXYE, *Bc_DecayVtxXZE, *Bc_DecayVtxYZE;

  //Ds candidate variables
  std::vector<float>       *Bc_Ds_mass, *Bc_Ds_px, *Bc_Ds_py, *Bc_Ds_pz, *Bc_Ds_pt, *Bc_Ds_eta, *Bc_Ds_phi, *Bc_Ds_Prob;
  std::vector<int>         *Bc_Ds_charge;
  std::vector<float>       *Bc_Ds_chi2, *Bc_Ds_DecayVtxCL, *Bc_Ds_DecayVtxX, *Bc_Ds_DecayVtxY, *Bc_Ds_DecayVtxZ;
  std::vector<float>       *Bc_Ds_DecayVtxXE, *Bc_Ds_DecayVtxYE, *Bc_Ds_DecayVtxZE;
  std::vector<float>       *Bc_Ds_DecayVtxXYE, *Bc_Ds_DecayVtxXZE, *Bc_Ds_DecayVtxYZE;
 
  //Dimuon candidate variables
  std::vector<float>       *Bc_mumu_mass, *Bc_mumu_px, *Bc_mumu_py, *Bc_mumu_pz, *Bc_mumu_pt, *Bc_mumu_eta, *Bc_mumu_phi, *Bc_mumu_Prob;
  std::vector<int>         *Bc_mumu_charge;
  std::vector<float>       *Bc_mumu_chi2, *Bc_mumu_DecayVtxCL, *Bc_mumu_DecayVtxX, *Bc_mumu_DecayVtxY, *Bc_mumu_DecayVtxZ;
  std::vector<float>       *Bc_mumu_DecayVtxXE, *Bc_mumu_DecayVtxYE, *Bc_mumu_DecayVtxZE;
  std::vector<float>       *Bc_mumu_DecayVtxXYE, *Bc_mumu_DecayVtxXZE, *Bc_mumu_DecayVtxYZE;

  //Ds childrens' variables
  std::vector<float>       *Bc_Ds_pt1, *Bc_Ds_px1, *Bc_Ds_py1, *Bc_Ds_pz1;
  std::vector<int>         *Bc_Ds_charge1;
  std::vector<float>       *Bc_Ds_pt2, *Bc_Ds_px2, *Bc_Ds_py2, *Bc_Ds_pz2;
  std::vector<int>         *Bc_Ds_charge2;
  std::vector<float>       *Bc_Ds_pt3, *Bc_Ds_px3, *Bc_Ds_py3, *Bc_Ds_pz3;
  std::vector<int>         *Bc_Ds_charge3;
  
  std::vector<float>       *Bc_Ds_px1_track, *Bc_Ds_py1_track, *Bc_Ds_pz1_track;
  std::vector<float>       *Bc_Ds_px2_track, *Bc_Ds_py2_track, *Bc_Ds_pz2_track;
  std::vector<float>       *Bc_Ds_px3_track, *Bc_Ds_py3_track, *Bc_Ds_pz3_track;

  //Dimuon childrens' variables
  std::vector<float>       *Bc_mumu_pt1, *Bc_mumu_px1, *Bc_mumu_py1, *Bc_mumu_pz1;
  std::vector<int>         *Bc_mumu_charge1;
  std::vector<float>       *Bc_mumu_pt2, *Bc_mumu_px2, *Bc_mumu_py2, *Bc_mumu_pz2;
  std::vector<int>         *Bc_mumu_charge2;
  

  //vertice primario CON mayor Pt
  unsigned int             nVtx;
  float                    priVtxX, priVtxY, priVtxZ, priVtxXE, priVtxYE, priVtxZE, priVtxCL;
  float                    priVtxXYE, priVtxXZE, priVtxYZE;
 

  std::vector<float>       *pVtxIPX,  *pVtxIPY, *pVtxIPZ, *pVtxIPXE, *pVtxIPYE, *pVtxIPZE, *pVtxIPCL;
  std::vector<float>       *pVtxIPXYE,  *pVtxIPXZE, *pVtxIPYZE;

  
  int  run, event;
  int  lumiblock;


  std::vector<double>      *trk1px, *trk1py, *trk1pz, *trk1chg;
  std::vector<double>      *trk2px, *trk2py, *trk2pz, *trk2chg;

  
  //variables to monitor
  int event_counter_;

  TLorentzVector gen_b_p4,gen_ks_p4,gen_pion1_p4,gen_pion2_p4,gen_jpsi_p4,gen_muon1_p4,gen_muon2_p4;
  TVector3       gen_b_vtx,gen_jpsi_vtx;
  float          gen_b_ct;


};


#endif
