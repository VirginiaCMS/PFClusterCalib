#ifndef PFCLUSTERNTUPLIZER_H
#define PFCLUSTERNTUPLIZER_H

#include <TTree.h>

#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/ParticleFlowReco/interface/PFCluster.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>

using namespace std;

class PFClusterNtuplizer : public edm::EDAnalyzer {

 public:

   PFClusterNtuplizer(const edm::ParameterSet& cfg);
   virtual ~PFClusterNtuplizer() { };

 private:

   virtual void analyze(const edm::Event& e, const edm::EventSetup& es);

   // tokens to collections
   edm::EDGetTokenT<vector<PileupSummaryInfo> >         mTokenGenPileup;
   edm::EDGetTokenT<vector<reco::GenParticle> >         mTokenGenParticles;
   edm::EDGetTokenT<vector<reco::PFCluster> >           mTokenPfClusters;
   edm::EDGetTokenT<reco::PFCluster::EEtoPSAssociation> mTokenPsClusters;
   edm::EDGetTokenT<EcalRecHitCollection>               mTokenRecHitsEB;
   edm::EDGetTokenT<EcalRecHitCollection>               mTokenRecHitsEE;
   edm::EDGetTokenT<vector<reco::Vertex> >              mTokenVertices;

   TTree* mTree;

   // variables associated with tree branches
   // MC pileup info
   vector<Int_t> mMcPUBunchCross;
   vector<Int_t> mMcPUNumIntObs;
   vector<float> mMcPUNumIntTrue;

   // MC truth
   float mMcVtxX;
   float mMcVtxY;
   float mMcVtxZ;
   float mMcPt;
   float mMcEta;
   float mMcPhi;
   float mMcE;

   // deltaR between MC photon and its matched PFCluster
   float mPfPhoDeltaR;

   // number of reconstructed primary vertices
   Int_t mNVtx;

   // PFClusters
   Int_t mPfSize;
   Int_t mPfSize5x5_ZS;
   Int_t mPfSize5x5_noZS;
   float mPfPt;
   float mPfEta;
   float mPfPhi;
   float mPfE;

   // EcalClusterLazyTools
   float mPfE1x3;
   float mPfE2x2;
   float mPfE2x5Max;
   float mPfE3x3;
   float mPfE5x5;

   // noZS::EcalClusterLazyTools
   float mPfE1x3_noZS;
   float mPfE2x2_noZS;
   float mPfE2x5Max_noZS;
   float mPfE3x3_noZS;
   float mPfE5x5_noZS;

   // number of preshower clusters, sums of preshower energies (2 layers)
   Int_t mPs1N;
   Int_t mPs2N;
   float mPs1E;
   float mPs2E;
};

#endif
