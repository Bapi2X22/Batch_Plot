//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Aug 16 18:22:34 2018 by ROOT version 6.06/01
// from TTree hits/HGC rechits
// found on file: muon_v10.root
//////////////////////////////////////////////////////////

#ifndef HGCNtupleVariables_h
#define HGCNtupleVariables_h

#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

using namespace std;

class HGCNtupleVariables {
public:
  HGCNtupleVariables(TTree * /*tree*/ = 0) : fChain(0) {}
  ~HGCNtupleVariables() {}
  // void    Init(TTree *tree);
  void Init(TTree *tree, TTree *tree2);
  Bool_t Notify();
  Int_t GetEntry(Long64_t entry, Int_t getall = 0) {
    return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
  }

  TTree *fChain;   //! pointer to the analyzed TTree or TChain
  TTree *fChain2;  //! pointer to the analyzed TTree or TChain
  Int_t fCurrent;  //! current Tree number in a TChain
  Int_t fCurrent2; //! current Tree number in a TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types
  //
vector<float>   *iEtaPho1;
   vector<float>   *iPhiPho1;
   vector<float>   *Hit_ES_Eta_Pho1;
   vector<float>   *Hit_ES_Phi_Pho1;
   vector<float>   *Hit_ES_X_Pho1;
   vector<float>   *Hit_ES_Y_Pho1;
   vector<float>   *Hit_ES_Z_Pho1;
   vector<float>   *ES_RecHitEnPho1;
   vector<float>   *Hit_Eta_Pho1;
   vector<float>   *Hit_Phi_Pho1;
   vector<float>   *Hit_X_Pho1;
   vector<float>   *Hit_Y_Pho1;
   vector<float>   *Hit_Z_Pho1;
   vector<float>   *RecHitEnPho1;
   vector<float>   *RecHitFracPho1;
   vector<int>     *RecHitGain1;
   vector<bool>    *RecHitQuality1;
   vector<float>   *HitNoisePho1;
   vector<float>   *dRHit_Eta_Pho1;
   vector<float>   *dRHit_Phi_Pho1;
   vector<float>   *RawHit_Eta_Pho1;
   vector<float>   *RawHit_Phi_Pho1;
   vector<float>   *RawRecHitEnPho1;
   vector<float>   *dRHit_X_Pho1;
   vector<float>   *dRHit_Y_Pho1;
   vector<float>   *dRHit_Z_Pho1;
   vector<float>   *dRRecHitEnPho1;
   vector<float>   *dRRecHitFracPho1;
   vector<float>   *dRHit_Eta_Pho2;
   vector<float>   *dRHit_Phi_Pho2;
   vector<float>   *dRHit_X_Pho2;
   vector<float>   *dRHit_Y_Pho2;
   vector<float>   *dRHit_Z_Pho2;
   vector<float>   *dRRecHitEnPho2;
   vector<float>   *dRRecHitFracPho2;
   vector<bool>    *RecHitFlag_kGood_pho1;
   vector<bool>    *RecHitFlag_kPoorReco_pho1;
   vector<bool>    *RecHitFlag_kOutOfTime_pho1;
   vector<bool>    *RecHitFlag_kFaultyHardware_pho1;
   vector<bool>    *RecHitFlag_kNoisy_pho1;
   vector<bool>    *RecHitFlag_kPoorCalib_pho1;
   vector<bool>    *RecHitFlag_kSaturated_pho1;
   vector<bool>    *RecHitFlag_kLeadingEdgeRecovered_pho1;
   vector<bool>    *RecHitFlag_kNeighboursRecovered_pho1;
   vector<bool>    *RecHitFlag_kTowerRecovered_pho1;
   vector<bool>    *RecHitFlag_kDead_pho1;
   vector<bool>    *RecHitFlag_kKilled_pho1;
   vector<bool>    *RecHitFlag_kTPSaturated_pho1;
   vector<bool>    *RecHitFlag_kL1SpikeFlag_pho1;
   vector<bool>    *RecHitFlag_kWeird_pho1;
   vector<bool>    *RecHitFlag_kDiWeird_pho1;
   vector<bool>    *RecHitFlag_kHasSwitchToGain6_pho1;
   vector<bool>    *RecHitFlag_kHasSwitchToGain1_pho1;
   vector<bool>    *RecHitFlag_kESGood_pho1;
   vector<bool>    *RecHitFlag_kESDead_pho1;
   vector<bool>    *RecHitFlag_kESHot_pho1;
   vector<bool>    *RecHitFlag_kESPassBX_pho1;
   vector<bool>    *RecHitFlag_kESTwoGoodRatios_pho1;
   vector<bool>    *RecHitFlag_kESBadRatioFor12_pho1;
   vector<bool>    *RecHitFlag_kESBadRatioFor23Upper_pho1;
   vector<bool>    *RecHitFlag_kESBadRatioFor23Lower_pho1;
   vector<bool>    *RecHitFlag_kESTS1Largest_pho1;
   vector<bool>    *RecHitFlag_kESTS3Largest_pho1;
   vector<bool>    *RecHitFlag_kESTS3Negative_pho1;
   vector<bool>    *RecHitFlag_kESSaturated_pho1;
   vector<bool>    *RecHitFlag_kESTS2Saturated_pho1;
   vector<bool>    *RecHitFlag_kESTS3Saturated_pho1;
   vector<bool>    *RecHitFlag_kESTS13Sigmas_pho1;
   vector<bool>    *RecHitFlag_kESTS15Sigmas_pho1;
   vector<float>   *iEtaPho2;
   vector<float>   *iPhiPho2;
   vector<float>   *Hit_ES_Eta_Pho2;
   vector<float>   *Hit_ES_Phi_Pho2;
   vector<float>   *Hit_ES_X_Pho2;
   vector<float>   *Hit_ES_Y_Pho2;
   vector<float>   *Hit_ES_Z_Pho2;
   vector<float>   *ES_RecHitEnPho2;
   vector<float>   *Hit_Eta_Pho2;
   vector<float>   *Hit_Phi_Pho2;
   vector<float>   *Hit_X_Pho2;
   vector<float>   *Hit_Y_Pho2;
   vector<float>   *Hit_Z_Pho2;
   vector<float>   *RecHitEnPho2;
   vector<float>   *RecHitFracPho2;
   vector<int>     *RecHitGain2;
   vector<bool>    *RecHitQuality2;
   vector<float>   *HitNoisePho2;
   vector<bool>    *RecHitFlag_kGood_pho2;
   vector<bool>    *RecHitFlag_kPoorReco_pho2;
   vector<bool>    *RecHitFlag_kOutOfTime_pho2;
   vector<bool>    *RecHitFlag_kFaultyHardware_pho2;
   vector<bool>    *RecHitFlag_kNoisy_pho2;
   vector<bool>    *RecHitFlag_kPoorCalib_pho2;
   vector<bool>    *RecHitFlag_kSaturated_pho2;
   vector<bool>    *RecHitFlag_kLeadingEdgeRecovered_pho2;
   vector<bool>    *RecHitFlag_kNeighboursRecovered_pho2;
   vector<bool>    *RecHitFlag_kTowerRecovered_pho2;
   vector<bool>    *RecHitFlag_kDead_pho2;
   vector<bool>    *RecHitFlag_kKilled_pho2;
   vector<bool>    *RecHitFlag_kTPSaturated_pho2;
   vector<bool>    *RecHitFlag_kL1SpikeFlag_pho2;
   vector<bool>    *RecHitFlag_kWeird_pho2;
   vector<bool>    *RecHitFlag_kDiWeird_pho2;
   vector<bool>    *RecHitFlag_kHasSwitchToGain6_pho2;
   vector<bool>    *RecHitFlag_kHasSwitchToGain1_pho2;
   vector<bool>    *RecHitFlag_kESGood_pho2;
   vector<bool>    *RecHitFlag_kESDead_pho2;
   vector<bool>    *RecHitFlag_kESHot_pho2;
   vector<bool>    *RecHitFlag_kESPassBX_pho2;
   vector<bool>    *RecHitFlag_kESTwoGoodRatios_pho2;
   vector<bool>    *RecHitFlag_kESBadRatioFor12_pho2;
   vector<bool>    *RecHitFlag_kESBadRatioFor23Upper_pho2;
   vector<bool>    *RecHitFlag_kESBadRatioFor23Lower_pho2;
   vector<bool>    *RecHitFlag_kESTS1Largest_pho2;
   vector<bool>    *RecHitFlag_kESTS3Largest_pho2;
   vector<bool>    *RecHitFlag_kESTS3Negative_pho2;
   vector<bool>    *RecHitFlag_kESSaturated_pho2;
   vector<bool>    *RecHitFlag_kESTS2Saturated_pho2;
   vector<bool>    *RecHitFlag_kESTS3Saturated_pho2;
   vector<bool>    *RecHitFlag_kESTS13Sigmas_pho2;
   vector<bool>    *RecHitFlag_kESTS15Sigmas_pho2;
   vector<float>   *iEtaPho3;
   vector<float>   *iPhiPho3;
   vector<float>   *Hit_ES_Eta_Pho3;
   vector<float>   *Hit_ES_Phi_Pho3;
   vector<float>   *Hit_ES_X_Pho3;
   vector<float>   *Hit_ES_Y_Pho3;
   vector<float>   *Hit_ES_Z_Pho3;
   vector<float>   *ES_RecHitEnPho3;
   vector<float>   *Hit_Eta_Pho3;
   vector<float>   *Hit_Phi_Pho3;
   vector<float>   *Hit_X_Pho3;
   vector<float>   *Hit_Y_Pho3;
   vector<float>   *Hit_Z_Pho3;
   vector<float>   *RecHitEnPho3;
   vector<float>   *RecHitFracPho3;
   vector<int>     *RecHitGain3;
   vector<bool>    *RecHitQuality3;
   vector<float>   *HitNoisePho3;
   vector<float>   *iEtaPho4;
   vector<float>   *iPhiPho4;
   vector<float>   *Hit_ES_Eta_Pho4;
   vector<float>   *Hit_ES_Phi_Pho4;
   vector<float>   *Hit_ES_X_Pho4;
   vector<float>   *Hit_ES_Y_Pho4;
   vector<float>   *Hit_ES_Z_Pho4;
   vector<float>   *ES_RecHitEnPho4;
   vector<float>   *Hit_Eta_Pho4;
   vector<float>   *Hit_Phi_Pho4;
   vector<float>   *Hit_X_Pho4;
   vector<float>   *Hit_Y_Pho4;
   vector<float>   *Hit_Z_Pho4;
   vector<float>   *RecHitEnPho4;
   vector<float>   *RecHitFracPho4;
   vector<int>     *RecHitGain4;
   vector<bool>    *RecHitQuality4;
   vector<float>   *HitNoisePho4;
   Int_t           nPhotons;
   vector<int>     *A_flags;
   vector<float>   *pt;
   vector<float>   *eta;
   vector<float>   *phi;
   vector<float>   *Pho_cluster_seed_x;
   vector<float>   *Pho_cluster_seed_y;
   vector<float>   *Pho_cluster_seed_z;
   vector<float>   *Pho_cluster_seed_eta;
   vector<float>   *Pho_cluster_seed_phi;
   vector<float>   *energy;
   vector<float>   *energy_ecal_mustache;
   vector<int>     *passMediumId;
   vector<int>     *passTightId;
   vector<int>     *passMVAMediumId;
   vector<float>   *Pho_R9;
   vector<float>   *Pho_S4;
   vector<float>   *Pho_SigIEIE;
   vector<float>   *Pho_SigIPhiIPhi;
   vector<float>   *Pho_SCEtaW;
   vector<float>   *Pho_SCPhiW;
   vector<float>   *Pho_CovIEtaIEta;
   vector<float>   *Pho_CovIEtaIPhi;
   vector<float>   *Pho_ESSigRR;
   vector<float>   *Pho_SCRawE;
   vector<float>   *Pho_SC_ESEnByRawE;
   vector<float>   *Pho_HadOverEm;
   vector<float>   *Pho_PFChIso;
   vector<float>   *Pho_PFChPVIso;
   vector<float>   *Pho_PFPhoIso;
   vector<float>   *Pho_PFNeuIso;
   vector<float>   *Pho_PFChWorstVetoIso;
   vector<float>   *Pho_PFChWorstIso;
   vector<float>   *Pho_EcalPFClusterIso;
   vector<float>   *Pho_HcalPFClusterIso;
   vector<float>   *Pho_CorrectedEnergy;
   vector<float>   *Pho_CorrectedEnergyError;
   vector<float>   *A_lead_Gen_mass;
   vector<float>   *A_lead_Gen_pt;
   vector<float>   *A_lead_Gen_eta;
   vector<float>   *A_lead_Gen_phi;
   vector<float>   *A_sublead_Gen_mass;
   vector<float>   *A_sublead_Gen_pt;
   vector<float>   *A_sublead_Gen_eta;
   vector<float>   *A_sublead_Gen_phi;
   vector<float>   *H_Gen_mass;
   vector<float>   *H_Gen_pt;
   vector<float>   *H_Gen_eta;
   vector<float>   *H_Gen_phi;
   vector<float>   *A_lead_Pho_Gen_Pt;
   vector<float>   *A_lead_Pho_Gen_Eta;
   vector<float>   *A_lead_Pho_Gen_Phi;
   vector<float>   *A_lead_Pho_Gen_E;
   vector<float>   *A_sublead_Pho_Gen_Pt;
   vector<float>   *A_sublead_Pho_Gen_Eta;
   vector<float>   *A_sublead_Pho_Gen_Phi;
   vector<float>   *A_sublead_Pho_Gen_E;
   Float_t         rho;
   Int_t           run;
   Int_t           event;
   Int_t           lumi;   
// List of branches
vector<float>   *iEtaEle1;
vector<float>   *iPhiEle1;
vector<float>   *Hit_ES_Eta_Ele1;
vector<float>   *Hit_ES_Phi_Ele1;
vector<float>   *Hit_ES_X_Ele1;
vector<float>   *Hit_ES_Y_Ele1;
vector<float>   *Hit_ES_Z_Ele1;
vector<float>   *ES_RecHitEnEle1;
vector<float>   *Hit_Eta_Ele1;
vector<float>   *Hit_Phi_Ele1;
vector<float>   *Hit_X_Ele1;
vector<float>   *Hit_Y_Ele1;
vector<float>   *Hit_Z_Ele1;
vector<float>   *RecHitEnEle1;
vector<float>   *RecHitFracEle1;
vector<int>     *RecHitGain1;
vector<bool>    *RecHitQuality1;
vector<bool>    *RecHitFlag_kGood_ele1;
vector<bool>    *RecHitFlag_kPoorReco_ele1;
vector<bool>    *RecHitFlag_kOutOfTime_ele1;
vector<bool>    *RecHitFlag_kFaultyHardware_ele1;
vector<bool>    *RecHitFlag_kNoisy_ele1;
vector<bool>    *RecHitFlag_kPoorCalib_ele1;
vector<bool>    *RecHitFlag_kSaturated_ele1;
vector<bool>    *RecHitFlag_kLeadingEdgeRecovered_ele1;
vector<bool>    *RecHitFlag_kNeighboursRecovered_ele1;
vector<bool>    *RecHitFlag_kTowerRecovered_ele1;
vector<bool>    *RecHitFlag_kDead_ele1;
vector<bool>    *RecHitFlag_kKilled_ele1;
vector<bool>    *RecHitFlag_kTPSaturated_ele1;
vector<bool>    *RecHitFlag_kL1SpikeFlag_ele1;
vector<bool>    *RecHitFlag_kWeird_ele1;
vector<bool>    *RecHitFlag_kDiWeird_ele1;
vector<bool>    *RecHitFlag_kHasSwitchToGain6_ele1;
vector<bool>    *RecHitFlag_kHasSwitchToGain1_ele1;
vector<bool>    *RecHitFlag_kESGood_ele1;
vector<bool>    *RecHitFlag_kESDead_ele1;
vector<bool>    *RecHitFlag_kESHot_ele1;
vector<bool>    *RecHitFlag_kESPassBX_ele1;
vector<bool>    *RecHitFlag_kESTwoGoodRatios_ele1;
vector<bool>    *RecHitFlag_kESBadRatioFor12_ele1;
vector<bool>    *RecHitFlag_kESBadRatioFor23Upper_ele1;
vector<bool>    *RecHitFlag_kESBadRatioFor23Lower_ele1;
vector<bool>    *RecHitFlag_kESTS1Largest_ele1;
vector<bool>    *RecHitFlag_kESTS3Largest_ele1;
vector<bool>    *RecHitFlag_kESTS3Negative_ele1;
vector<bool>    *RecHitFlag_kESSaturated_ele1;
vector<bool>    *RecHitFlag_kESTS2Saturated_ele1;
vector<bool>    *RecHitFlag_kESTS3Saturated_ele1;
vector<bool>    *RecHitFlag_kESTS13Sigmas_ele1;
vector<bool>    *RecHitFlag_kESTS15Sigmas_ele1;
vector<float>   *iEtaEle2;
vector<float>   *iPhiEle2;
vector<float>   *Hit_ES_Eta_Ele2;
vector<float>   *Hit_ES_Phi_Ele2;
vector<float>   *Hit_ES_X_Ele2;
vector<float>   *Hit_ES_Y_Ele2;
vector<float>   *Hit_ES_Z_Ele2;
vector<float>   *ES_RecHitEnEle2;
vector<float>   *Hit_Eta_Ele2;
vector<float>   *Hit_Phi_Ele2;
vector<float>   *Hit_X_Ele2;
vector<float>   *Hit_Y_Ele2;
vector<float>   *Hit_Z_Ele2;
vector<float>   *RecHitEnEle2;
vector<float>   *RecHitFracEle2;
vector<int>     *RecHitGain2;
vector<bool>    *RecHitQuality2;
vector<bool>    *RecHitFlag_kGood_ele2;
vector<bool>    *RecHitFlag_kPoorReco_ele2;
vector<bool>    *RecHitFlag_kOutOfTime_ele2;
vector<bool>    *RecHitFlag_kFaultyHardware_ele2;
vector<bool>    *RecHitFlag_kNoisy_ele2;
vector<bool>    *RecHitFlag_kPoorCalib_ele2;
vector<bool>    *RecHitFlag_kSaturated_ele2;
vector<bool>    *RecHitFlag_kLeadingEdgeRecovered_ele2;
vector<bool>    *RecHitFlag_kNeighboursRecovered_ele2;
vector<bool>    *RecHitFlag_kTowerRecovered_ele2;
vector<bool>    *RecHitFlag_kDead_ele2;
vector<bool>    *RecHitFlag_kKilled_ele2;
vector<bool>    *RecHitFlag_kTPSaturated_ele2;
vector<bool>    *RecHitFlag_kL1SpikeFlag_ele2;
vector<bool>    *RecHitFlag_kWeird_ele2;
vector<bool>    *RecHitFlag_kDiWeird_ele2;
vector<bool>    *RecHitFlag_kHasSwitchToGain6_ele2;
vector<bool>    *RecHitFlag_kHasSwitchToGain1_ele2;
vector<bool>    *RecHitFlag_kESGood_ele2;
vector<bool>    *RecHitFlag_kESDead_ele2;
vector<bool>    *RecHitFlag_kESHot_ele2;
vector<bool>    *RecHitFlag_kESPassBX_ele2;
vector<bool>    *RecHitFlag_kESTwoGoodRatios_ele2;
vector<bool>    *RecHitFlag_kESBadRatioFor12_ele2;
vector<bool>    *RecHitFlag_kESBadRatioFor23Upper_ele2;
vector<bool>    *RecHitFlag_kESBadRatioFor23Lower_ele2;
vector<bool>    *RecHitFlag_kESTS1Largest_ele2;
vector<bool>    *RecHitFlag_kESTS3Largest_ele2;
vector<bool>    *RecHitFlag_kESTS3Negative_ele2;
vector<bool>    *RecHitFlag_kESSaturated_ele2;
vector<bool>    *RecHitFlag_kESTS2Saturated_ele2;
vector<bool>    *RecHitFlag_kESTS3Saturated_ele2;
vector<bool>    *RecHitFlag_kESTS13Sigmas_ele2;
vector<bool>    *RecHitFlag_kESTS15Sigmas_ele2;
Int_t           nElectrons;
vector<float>   *pt;
vector<float>   *eta;
vector<float>   *phi;
vector<float>   *energy;
vector<float>   *energy_error;
vector<float>   *energy_ecal_mustache;
vector<int>     *passLooseId;
vector<int>     *passMediumId;
vector<int>     *passTightId;
vector<float>   *Ele_R9;
vector<float>   *Ele_S4;
vector<float>   *Ele_SigIEIE;
vector<float>   *Ele_SigIPhiIPhi;
vector<float>   *Ele_SCEtaW;
vector<float>   *Ele_SCPhiW;
vector<float>   *Ele_CovIEtaIEta;
vector<float>   *Ele_CovIEtaIPhi;
vector<float>   *Ele_ESSigRR;
vector<float>   *Ele_SCRawE;
vector<float>   *Ele_SC_ESEnByRawE;
vector<float>   *Ele_HadOverEm;
vector<float>   *sumChargedHadronPt;
vector<float>   *sumChargedParticlePt;
vector<float>   *sumEcalClusterEt;
vector<float>   *sumHcalClusterEt;
vector<float>   *sumNeutralHadronEt;
vector<float>   *sumPhotonEt;
vector<float>   *sumPUPt;
vector<float>   *Ele_EcalPFClusterIso;
vector<float>   *Ele_HcalPFClusterIso;
vector<float>   *matchedGenEta_;
vector<float>   *matchedGenphi_;
vector<float>   *matchedGenpt_;
vector<float>   *matchedGenEnergy_;
vector<float>   *recoDRNEnergy;
Float_t         rho;
Int_t           run;
Int_t           event;
Int_t           lumi;
Bool_t          isRefinedSC;

// List of branches
TBranch        *b_iEtaEle1;   //!
TBranch        *b_iPhiEle1;   //!
TBranch        *b_Hit_ES_Eta_Ele1;   //!
TBranch        *b_Hit_ES_Phi_Ele1;   //!
TBranch        *b_Hit_ES_X_Ele1;   //!
TBranch        *b_Hit_ES_Y_Ele1;   //!
TBranch        *b_Hit_ES_Z_Ele1;   //!
TBranch        *b_ES_RecHitEnEle1;   //!
TBranch        *b_Hit_Eta_Ele1;   //!
TBranch        *b_Hit_Phi_Ele1;   //!
TBranch        *b_Hit_X_Ele1;   //!
TBranch        *b_Hit_Y_Ele1;   //!
TBranch        *b_Hit_Z_Ele1;   //!
TBranch        *b_RecHitEnEle1;   //!
TBranch        *b_RecHitFracEle1;   //!
TBranch        *b_RecHitGain1;   //!
TBranch        *b_RecHitQuality1;   //!
TBranch        *b_RecHitFlag_kGood_ele1;   //!
TBranch        *b_RecHitFlag_kPoorReco_ele1;   //!
TBranch        *b_RecHitFlag_kOutOfTime_ele1;   //!
TBranch        *b_RecHitFlag_kFaultyHardware_ele1;   //!
TBranch        *b_RecHitFlag_kNoisy_ele1;   //!
TBranch        *b_RecHitFlag_kPoorCalib_ele1;   //!
TBranch        *b_RecHitFlag_kSaturated_ele1;   //!
TBranch        *b_RecHitFlag_kLeadingEdgeRecovered_ele1;   //!
TBranch        *b_RecHitFlag_kNeighboursRecovered_ele1;   //!
TBranch        *b_RecHitFlag_kTowerRecovered_ele1;   //!
TBranch        *b_RecHitFlag_kDead_ele1;   //!
TBranch        *b_RecHitFlag_kKilled_ele1;   //!
TBranch        *b_RecHitFlag_kTPSaturated_ele1;   //!
TBranch        *b_RecHitFlag_kL1SpikeFlag_ele1;   //!
TBranch        *b_RecHitFlag_kWeird_ele1;   //!
TBranch        *b_RecHitFlag_kDiWeird_ele1;   //!
TBranch        *b_RecHitFlag_kHasSwitchToGain6_ele1;   //!
TBranch        *b_RecHitFlag_kHasSwitchToGain1_ele1;   //!
TBranch        *b_RecHitFlag_kESGood_ele1;   //!
TBranch        *b_RecHitFlag_kESDead_ele1;   //!
TBranch        *b_RecHitFlag_kESHot_ele1;   //!
TBranch        *b_RecHitFlag_kESPassBX_ele1;   //!
TBranch        *b_RecHitFlag_kESTwoGoodRatios_ele1;   //!
TBranch        *b_RecHitFlag_kESBadRatioFor12_ele1;   //!
TBranch        *b_RecHitFlag_kESBadRatioFor23Upper_ele1;   //!
TBranch        *b_RecHitFlag_kESBadRatioFor23Lower_ele1;   //!
TBranch        *b_RecHitFlag_kESTS1Largest_ele1;   //!
TBranch        *b_RecHitFlag_kESTS3Largest_ele1;   //!
TBranch        *b_RecHitFlag_kESTS3Negative_ele1;   //!
TBranch        *b_RecHitFlag_kESSaturated_ele1;   //!
TBranch        *b_RecHitFlag_kESTS2Saturated_ele1;   //!
TBranch        *b_RecHitFlag_kESTS3Saturated_ele1;   //!
TBranch        *b_RecHitFlag_kESTS13Sigmas_ele1;   //!
TBranch        *b_RecHitFlag_kESTS15Sigmas_ele1;   //!
TBranch        *b_iEtaEle2;   //!
TBranch        *b_iPhiEle2;   //!
TBranch        *b_Hit_ES_Eta_Ele2;   //!
TBranch        *b_Hit_ES_Phi_Ele2;   //!
TBranch        *b_Hit_ES_X_Ele2;   //!
TBranch        *b_Hit_ES_Y_Ele2;   //!
TBranch        *b_Hit_ES_Z_Ele2;   //!
TBranch        *b_ES_RecHitEnEle2;   //!
TBranch        *b_Hit_Eta_Ele2;   //!
TBranch        *b_Hit_Phi_Ele2;   //!
TBranch        *b_Hit_X_Ele2;   //!
TBranch        *b_Hit_Y_Ele2;   //!
TBranch        *b_Hit_Z_Ele2;   //!
TBranch        *b_RecHitEnEle2;   //!
TBranch        *b_RecHitFracEle2;   //!
TBranch        *b_RecHitGain2;   //!
TBranch        *b_RecHitQuality2;   //!
TBranch        *b_RecHitFlag_kGood_ele2;   //!
TBranch        *b_RecHitFlag_kPoorReco_ele2;   //!
TBranch        *b_RecHitFlag_kOutOfTime_ele2;   //!
TBranch        *b_RecHitFlag_kFaultyHardware_ele2;   //!
TBranch        *b_RecHitFlag_kNoisy_ele2;   //!
TBranch        *b_RecHitFlag_kPoorCalib_ele2;   //!
TBranch        *b_RecHitFlag_kSaturated_ele2;   //!
TBranch        *b_RecHitFlag_kLeadingEdgeRecovered_ele2;   //!
TBranch        *b_RecHitFlag_kNeighboursRecovered_ele2;   //!
TBranch        *b_RecHitFlag_kTowerRecovered_ele2;   //!
TBranch        *b_RecHitFlag_kDead_ele2;   //!
TBranch        *b_RecHitFlag_kKilled_ele2;   //!
TBranch        *b_RecHitFlag_kTPSaturated_ele2;   //!
TBranch        *b_RecHitFlag_kL1SpikeFlag_ele2;   //!
TBranch        *b_RecHitFlag_kWeird_ele2;   //!
TBranch        *b_RecHitFlag_kDiWeird_ele2;   //!
TBranch        *b_RecHitFlag_kHasSwitchToGain6_ele2;   //!
TBranch        *b_RecHitFlag_kHasSwitchToGain1_ele2;   //!
TBranch        *b_RecHitFlag_kESGood_ele2;   //!
TBranch        *b_RecHitFlag_kESDead_ele2;   //!
TBranch        *b_RecHitFlag_kESHot_ele2;   //!
TBranch        *b_RecHitFlag_kESPassBX_ele2;   //!
TBranch        *b_RecHitFlag_kESTwoGoodRatios_ele2;   //!
TBranch        *b_RecHitFlag_kESBadRatioFor12_ele2;   //!
TBranch        *b_RecHitFlag_kESBadRatioFor23Upper_ele2;   //!
TBranch        *b_RecHitFlag_kESBadRatioFor23Lower_ele2;   //!
TBranch        *b_RecHitFlag_kESTS1Largest_ele2;   //!
TBranch        *b_RecHitFlag_kESTS3Largest_ele2;   //!
TBranch        *b_RecHitFlag_kESTS3Negative_ele2;   //!
TBranch        *b_RecHitFlag_kESSaturated_ele2;   //!
TBranch        *b_RecHitFlag_kESTS2Saturated_ele2;   //!
TBranch        *b_RecHitFlag_kESTS3Saturated_ele2;   //!
TBranch        *b_RecHitFlag_kESTS13Sigmas_ele2;   //!
TBranch        *b_RecHitFlag_kESTS15Sigmas_ele2;   //!
TBranch        *b_nEle;   //!
TBranch        *b_pt;   //!
TBranch        *b_eta;   //!
TBranch        *b_phi;   //!
TBranch        *b_energy;   //!
TBranch        *b_energy_error;   //!
TBranch        *b_energy_ecal_mustache;   //!
TBranch        *b_passLooseId;   //!
TBranch        *b_passMediumId;   //!
TBranch        *b_passTightId;   //!
TBranch        *b_Ele_R9;   //!
TBranch        *b_Ele_S4;   //!
TBranch        *b_Ele_SigIEIE;   //!
TBranch        *b_Ele_SigIPhiIPhi;   //!
TBranch        *b_Ele_SCEtaW;   //!
TBranch        *b_Ele_SCPhiW;   //!
TBranch        *b_Ele_CovIEtaIEta;   //!
TBranch        *b_Ele_CovIEtaIPhi;   //!
TBranch        *b_Ele_ESSigRR;   //!
TBranch        *b_Ele_SCRawE;   //!
TBranch        *b_Ele_SC_ESEnByRawE;   //!
TBranch        *b_Ele_HadOverEm;   //!
TBranch        *b_sumChargedHadronPt;   //!
TBranch        *b_sumChargedParticlePt;   //!
TBranch        *b_sumEcalClusterEt;   //!
TBranch        *b_sumHcalClusterEt;   //!
TBranch        *b_sumNeutralHadronEt;   //!
TBranch        *b_sumPhotonEt;   //!
TBranch        *b_sumPUPt;   //!
TBranch        *b_Ele_EcalPFClusterIso;   //!
TBranch        *b_Ele_HcalPFClusterIso;   //!
TBranch        *b_matchedGenEta_;   //!
TBranch        *b_matchedGenphi_;   //!
TBranch        *b_matchedGenpt_;   //!
TBranch        *b_matchedGenEnergy_;   //!
TBranch        *b_recoDRNEnergy;   //!
TBranch        *b_rho;   //!
TBranch        *b_run;   //!
TBranch        *b_event;   //!
TBranch        *b_lumi;   //!
TBranch        *b_isRefinedSC;   //!
}; //Modified by Somanko

#endif
