#include <FWCore/ServiceRegistry/interface/Service.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h>

#include <TVector3.h>
#include <TLorentzVector.h>

#include "ncuAnalysis/PFClusterCalib/interface/PFClusterNtuplizer.h"

using namespace std;

PFClusterNtuplizer::PFClusterNtuplizer(const edm::ParameterSet& cfg)
{

   // initialize tokens to collections
   mTokenGenPileup    = consumes<vector<PileupSummaryInfo> >        (cfg.getParameter<edm::InputTag>("pileupSummaryLabel"));
   mTokenGenParticles = consumes<vector<reco::GenParticle> >        (cfg.getParameter<edm::InputTag>("genParticlesLabel"));
   mTokenPfClusters   = consumes<vector<reco::PFCluster> >          (cfg.getParameter<edm::InputTag>("pfClustersLabel"));
   mTokenPsClusters   = consumes<reco::PFCluster::EEtoPSAssociation>(cfg.getParameter<edm::InputTag>("pfClustersLabel"));
   mTokenRecHitsEB    = consumes<EcalRecHitCollection>              (cfg.getParameter<edm::InputTag>("recHitsEBLabel"));
   mTokenRecHitsEE    = consumes<EcalRecHitCollection>              (cfg.getParameter<edm::InputTag>("recHitsEELabel"));
   mTokenVertices     = consumes<vector<reco::Vertex> >             (cfg.getParameter<edm::InputTag>("verticesLabel"));

   // initialize output TTree
   edm::Service<TFileService> fs;
   mTree = fs->make<TTree>("PFClusterTree", "PFClusters from photons");

   mTree->Branch("mcPUBunchCross", &mMcPUBunchCross);
   mTree->Branch("mcPUNumIntObs",  &mMcPUNumIntObs);
   mTree->Branch("mcPUNumIntTrue", &mMcPUNumIntTrue);

   mTree->Branch("mcVtxX", &mMcVtxX);
   mTree->Branch("mcVtxY", &mMcVtxY);
   mTree->Branch("mcVtxZ", &mMcVtxZ);
   mTree->Branch("mcPt",   &mMcPt);
   mTree->Branch("mcEta",  &mMcEta);
   mTree->Branch("mcPhi",  &mMcPhi);
   mTree->Branch("mcE",    &mMcE);

   mTree->Branch("pfPhoDeltaR",    &mPfPhoDeltaR);

   mTree->Branch("nVtx", &mNVtx);

   mTree->Branch("pfSize",         &mPfSize);
   mTree->Branch("pfSize5x5_ZS",   &mPfSize5x5_ZS);
   mTree->Branch("pfSize5x5_noZS", &mPfSize5x5_noZS);

   mTree->Branch("pfPt",   &mPfPt);
   mTree->Branch("pfEta",  &mPfEta);
   mTree->Branch("pfPhi",  &mPfPhi);
   mTree->Branch("pfE",    &mPfE);

   mTree->Branch("pfE1x3",    &mPfE1x3);
   mTree->Branch("pfE2x2",    &mPfE2x2);
   mTree->Branch("pfE2x5Max", &mPfE2x5Max);
   mTree->Branch("pfE3x3",    &mPfE3x3);
   mTree->Branch("pfE5x5",    &mPfE5x5);

   mTree->Branch("pfE1x3_noZS",    &mPfE1x3_noZS);
   mTree->Branch("pfE2x2_noZS",    &mPfE2x2_noZS);
   mTree->Branch("pfE2x5Max_noZS", &mPfE2x5Max_noZS);
   mTree->Branch("pfE3x3_noZS",    &mPfE3x3_noZS);
   mTree->Branch("pfE5x5_noZS",    &mPfE5x5_noZS);

   mTree->Branch("ps1N", &mPs1N);
   mTree->Branch("ps2N", &mPs2N);
   mTree->Branch("ps1E", &mPs1E);
   mTree->Branch("ps2E", &mPs2E);
}

bool SortByKey(const reco::PFCluster::EEtoPSAssociation::value_type& a,
               const reco::PFCluster::EEtoPSAssociation::value_type& b)
{
   return a.first < b.first;
}

void PFClusterNtuplizer::analyze(const edm::Event& e, const edm::EventSetup& es)
{

   // get collections of objects
   edm::Handle<vector<PileupSummaryInfo> >         handleGenPileup;
   edm::Handle<vector<reco::GenParticle> >         handleGenParticles;
   edm::Handle<vector<reco::PFCluster> >           handlePfClusters;
   edm::Handle<reco::PFCluster::EEtoPSAssociation> handlePsClusters;
   edm::Handle<vector<reco::Vertex> >              handleVertices;

   e.getByToken(mTokenGenPileup,    handleGenPileup);
   e.getByToken(mTokenGenParticles, handleGenParticles);
   e.getByToken(mTokenPfClusters,   handlePfClusters);
   e.getByToken(mTokenPsClusters,   handlePsClusters);
   e.getByToken(mTokenVertices,     handleVertices);

   // pileup
   mMcPUBunchCross.clear();
   mMcPUNumIntObs.clear();
   mMcPUNumIntTrue.clear();
   for (auto pu = handleGenPileup->begin(); pu != handleGenPileup->end(); pu++) {
      mMcPUBunchCross.push_back(pu->getBunchCrossing());
      mMcPUNumIntObs.push_back(pu->getPU_NumInteractions());
      mMcPUNumIntTrue.push_back(pu->getTrueNumInteractions());
   }

   // number of reconstructed primary vertices
   mNVtx = 0;
   for (auto vtx = handleVertices->begin(); vtx != handleVertices->end(); vtx++)
      if (!vtx->isFake()) mNVtx++;

   EcalClusterLazyTools       lazyTool     (e, es, mTokenRecHitsEB, mTokenRecHitsEE);
   noZS::EcalClusterLazyTools lazyTool_noZS(e, es, mTokenRecHitsEB, mTokenRecHitsEE);

   // final-state photons
   for (auto p = handleGenParticles->begin(); p != handleGenParticles->end(); p++) {
      if (p->status() != 1 || p->mother() || p->pdgId() != 22)
         continue;

      mMcVtxX = p->vx();
      mMcVtxY = p->vy();
      mMcVtxZ = p->vz();
      mMcPt   = p->pt();
      mMcEta  = p->eta();
      mMcPhi  = p->phi();
      mMcE    = p->energy();

      TLorentzVector pho;
      pho.SetPtEtaPhiE(p->pt(), p->eta(), p->phi(), p->energy());

      mPfPhoDeltaR = 0.10001;

      const vector<reco::PFCluster>* pfClusters = handlePfClusters.product();

      // find PFCluster which corresponds to the photon direction
      for (size_t i = 0; i < pfClusters->size(); i++) {
         reco::PFCluster c = (*pfClusters)[i];

         TVector3 clus;
         clus.SetXYZ(c.x() - p->vx(), c.y() - p->vy(), c.z() - p->vz());
         float delta = clus.DeltaR(pho.Vect());

         if (delta < mPfPhoDeltaR) {
            mPfPhoDeltaR = delta;

            mPfSize = c.size();
            mPfSize5x5_ZS = lazyTool.n5x5(c);
            mPfSize5x5_noZS = lazyTool_noZS.n5x5(c);

            mPfPt  = c.pt();
            mPfEta = c.eta();
            mPfPhi = c.phi();
            mPfE   = c.energy();

            mPfE1x3    = lazyTool.e1x3(c);
            mPfE2x2    = lazyTool.e2x2(c);
            mPfE3x3    = lazyTool.e3x3(c);
            mPfE5x5    = lazyTool.e5x5(c);
            mPfE2x5Max = lazyTool.e2x5Max(c);

            mPfE1x3_noZS    = lazyTool_noZS.e1x3(c);
            mPfE2x2_noZS    = lazyTool_noZS.e2x2(c);
            mPfE3x3_noZS    = lazyTool_noZS.e3x3(c);
            mPfE5x5_noZS    = lazyTool_noZS.e5x5(c);
            mPfE2x5Max_noZS = lazyTool_noZS.e2x5Max(c);

            mPs1N = 0;
            mPs2N = 0;
            mPs1E = 0;
            mPs2E = 0;

            // search for preshower clusters associated with given PFCluster
            auto key = std::make_pair(i, edm::Ptr<reco::PFCluster>());
            auto range = std::equal_range(handlePsClusters.product()->begin(),
                                          handlePsClusters.product()->end(),
                                          key, SortByKey);

            // loop over matched preshower clusters
            for (auto pair = range.first; pair != range.second; pair++) {
               edm::Ptr<reco::PFCluster> ps = pair->second;

               if (ps->layer() == PFLayer::PS1) {
                  mPs1N++;
                  mPs1E += ps->energy();
               }
               else if (ps->layer() == PFLayer::PS2) {
                  mPs2N++;
                  mPs2E += ps->energy();
               }
            }

         } // delta < minDeltaR

      } // PFClusters

      if (mPfPhoDeltaR < 0.1)
         mTree->Fill();

   } // gen-level photons

}

DEFINE_FWK_MODULE(PFClusterNtuplizer);
