#include "Analysis/VVAnalysis/interface/SameSignLept.h"

#include <TStyle.h>
#include <regex>
#include "TParameter.h"


#define Fill1D(NAME, VALUE_) HistFullFill(histMap1D_ NAME, variation.second, VALUE_, weight);

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>LorentzVector; 

enum Lepton {Muon = 13, Electron = 11}; 



void TTTSelector::SlaveBegin(TTree * /*tree*/)
{
  return;
  pileupSF_ = (ScaleFactor *) GetInputList()->FindObject("pileupSF");
  if (pileupSF_ == nullptr ) 
    Abort("Must pass pileup weights SF");
  eIdSF_ = (ScaleFactor *) GetInputList()->FindObject("electronTightIdSF");
  if (eIdSF_ == nullptr ) 
    Abort("Must pass electron ID SF");
  eGsfSF_ = (ScaleFactor *) GetInputList()->FindObject("electronGsfSF");
  if (eGsfSF_ == nullptr ) 
    Abort("Must pass electron GSF SF");
  mIdSF_ = (ScaleFactor *) GetInputList()->FindObject("muonTightIdSF");
  if (mIdSF_ == nullptr ) 
    Abort("Must pass muon ID SF");
  mIsoSF_ = (ScaleFactor *) GetInputList()->FindObject("muonIsoSF");
  if (mIsoSF_ == nullptr ) 
    Abort("Must pass muon Iso SF");

  prefireEff_ = (TEfficiency*) GetInputList()->FindObject("prefireEfficiencyMap");
  if (prefireEff_ == nullptr ) 
    Abort("Must pass prefiring efficiency map");

}


void SameSignLept::Init(TTree *tree)
{

  /* if (addSumweights_) {
    TFile* file = fChain->GetTree()->GetDirectory()->GetFile(); 
    TTree* metaInfo = dynamic_cast<TTree*>(file->Get("metaInfo/metaInfo"));
    if (metaInfo == nullptr)
      std::cerr << "WARNING: Failed to add sumWeights histogram" << std::endl;
    else {
      metaInfo->Draw("1>>sumweights", "summedWeights");
    }
  }
  */
  b.SetTree(tree);
  allChannels_ = {"ee", "mm", "em", "all"};
  hists1D_ = {"CutFLow", "ZMass", "ptl1", "etal1", "ptl2", "etal1", "etal2", "SR", "bjetpt", "jetpt", "nbjet", "njet"};
  SelectorBase::Init(tree);
}


void SameSignLept::SetBranchesNanoAOD() {
  /*
  fChain->SetBranchAddress("nElectron", &nElectron, &b_nElectron);
  fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
  fChain->SetBranchAddress("Electron_pt", &Electron_pt, &b_Electron_pt);
  fChain->SetBranchAddress("Electron_eta", &Electron_eta, &b_Electron_eta);
  fChain->SetBranchAddress("Electron_phi", &Electron_phi, &b_Electron_phi);
  fChain->SetBranchAddress("Muon_pt", &Muon_pt, &b_Muon_pt);
  fChain->SetBranchAddress("Muon_eta", &Muon_eta, &b_Muon_eta);
  fChain->SetBranchAddress("Muon_phi", &Muon_phi, &b_Muon_phi);
  fChain->SetBranchAddress("Electron_cutBased", &Electron_cutBased, &b_Electron_cutBased);
  fChain->SetBranchAddress("Electron_mvaFall17V2Iso", &Electron_mvaFall17V2Iso, &b_Electron_mvaFall17V2Iso);
  fChain->SetBranchAddress("Muon_tightId", &Muon_tightId, &b_Muon_tightId);
  fChain->SetBranchAddress("Muon_mediumId", &Muon_mediumId, &b_Muon_mediumId);
  fChain->SetBranchAddress("Muon_pfIsoId", &Muon_pfIsoId, &b_Muon_pfIsoId);
  fChain->SetBranchAddress("Muon_pfRelIso04_all", &Muon_pfRelIso04_all, &b_Muon_pfRelIso04_all);
  fChain->SetBranchAddress("MET_pt", &MET, &b_MET);
  fChain->SetBranchAddress("MET_phi", &type1_pfMETPhi, &b_type1_pfMETPhi);
  fChain->SetBranchAddress("Electron_charge", &Electron_charge, &b_Electron_charge);
  fChain->SetBranchAddress("Muon_charge", &Muon_charge, &b_Muon_charge);
  fChain->SetBranchAddress("Electron_mass", &Electron_mass, &b_Electron_mass);
  fChain->SetBranchAddress("Muon_mass", &Muon_mass, &b_Muon_mass);
  fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
  fChain->SetBranchAddress("Jet_pt", &Jet_pt, &b_Jet_pt);
  fChain->SetBranchAddress("Jet_mass", &Jet_mass, &b_Jet_mass);
  fChain->SetBranchAddress("Jet_btagCSVV2", &Jet_btagCSVV2, &b_Jet_btagCSVV2);
  fChain->SetBranchAddress("Jet_btagDeepB", &Jet_btagDeepB, &b_Jet_btagDeepB);
  fChain->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);
  fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);

  if (isMC_) {
    fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
    fChain->SetBranchAddress("Pileup_nPU", &numPU, &b_numPU);
    } */

  b.CleanUp();

  b.SetBranch("nElectron", nElectron);
  b.SetBranch("Electron_pt", Electron_pt);
  b.SetBranch("Electron_eta", Electron_eta);
  b.SetBranch("Electron_phi", Electron_phi);
  b.SetBranch("Electron_charge", Electron_charge);
  b.SetBranch("Electron_mass", Electron_mass);
  b.SetBranch("Electron_miniPFRelIso_all", Electron_miniPFRelIso_all);
  b.SetBranch("Electron_mvaFall17V1noIso", Electron_MVA);
  b.SetBranch("Electron_cutBased_Fall17_V1", Electron_cutBased);
  

  b.SetBranch("nMuon", nMuon);
  b.SetBranch("Muon_pt", Muon_pt);
  b.SetBranch("Muon_eta", Muon_eta);
  b.SetBranch("Muon_phi", Muon_phi);
  b.SetBranch("Muon_tightId", Muon_tightId);
  b.SetBranch("Muon_mediumId", Muon_mediumId);
  b.SetBranch("Muon_pfRelIso04_all", Muon_pfRelIso04_all);
  b.SetBranch("Muon_miniPFRelIso_all", Muon_miniPFRelIso_all);
  b.SetBranch("Muon_charge", Muon_charge);
  b.SetBranch("Muon_mass", Muon_mass);

  b.SetBranch("nJet", nJet);
  b.SetBranch("Jet_btagCSVV2", Jet_btagCSVV2);
  b.SetBranch("Jet_btagDeepB", Jet_btagDeepB);
  b.SetBranch("Jet_eta", Jet_eta);
  b.SetBranch("Jet_phi", Jet_phi);
  b.SetBranch("Jet_pt", Jet_pt);
  b.SetBranch("Jet_mass", Jet_mass);

  b.SetBranch("Jet_neHEF", Jet_neHEF);
  b.SetBranch("Jet_neEmEF", Jet_neEmEF);
  b.SetBranch("Jet_nConstituents", Jet_nConstituents);
  b.SetBranch("Jet_chHEF", Jet_chHEF);
  b.SetBranch("Jet_chEmEF", Jet_chEmEF);

  b.SetBranch("MET_pt", MET);
  b.SetBranch("MET_phi", type1_pfMETPhi);
  
  if (isMC_) {
    b.SetBranch("genWeight", genWeight);
    b.SetBranch("Pileup_nPU", numPU);
  }
}
void SameSignLept::SetBranchesUWVV() {
  return;
}
void SameSignLept::LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation) {
  return;
}

void SameSignLept::clearValues() {
  weight = 1;
  Ht = 0;
  nTightJet = 0;
  nbJet = 0;
  goodParts.clear();
}
    
void SameSignLept::LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) { 
  /*
  weight = 1;
  b_nElectron->GetEntry(entry);
  b_nMuon->GetEntry(entry);
  b_Electron_pt->GetEntry(entry);
  b_Electron_eta->GetEntry(entry);
  b_Electron_phi->GetEntry(entry);
  b_Muon_pt->GetEntry(entry);
  b_Muon_eta->GetEntry(entry);
  b_Muon_phi->GetEntry(entry);
  b_Electron_cutBased->GetEntry(entry);
  b_Electron_mvaFall17V2Iso->GetEntry(entry);
  b_Muon_tightId->GetEntry(entry);
  b_Muon_mediumId->GetEntry(entry);
  b_Muon_pfRelIso04_all->GetEntry(entry);
  b_Electron_charge->GetEntry(entry);
  b_Muon_charge->GetEntry(entry);
  b_Electron_mass->GetEntry(entry);
  b_Muon_mass->GetEntry(entry);
  b_MET->GetEntry(entry);
  b_nJet->GetEntry(entry);
  b_Jet_pt->GetEntry(entry);
  b_Jet_mass->GetEntry(entry);
  b_Jet_btagCSVV2->GetEntry(entry);
  b_Jet_btagDeepB->GetEntry(entry);
  b_Jet_eta->GetEntry(entry); 
  b_genWeight->GetEntry(entry);
  */

  clearValues();
  b.SetEntry(entry);


  if (nElectron > N_KEEP_MU_E_ || nMuon > N_KEEP_MU_E_) {
    std::string message = "Found more electrons or muons than max read number.\n    Found ";
    message += std::to_string(nElectron);
    message += " electrons.\n    Found ";
    message += std::to_string(nMuon);
    message += " Muons\n  --> Max read number was ";
    message += std::to_string(N_KEEP_MU_E_);
    message += "\nExiting because this can cause problems. Increase N_KEEP_MU_E_ to avoid this error.\n";
    throw std::domain_error(message);
  }
  else if (nJet > N_KEEP_JET) {
    std::string message = "Found more jets  than max read number.\n    Found ";
    message += std::to_string(nJet);
    message += " jets.\n    --> Max read number was ";
    message += std::to_string(N_KEEP_JET);
    message += "\nExiting because this can cause problems. Increase N_KEEP_JET to avoid this error.\n";
    throw std::domain_error(message); 
  }
    
  /*
  CombMass = 0;
  l1Pt = 0;
  l2Pt = 0;
  l1Eta = 0;
  l2Eta = 0;
  l1Phi = 0;
  l2Phi = 0;
  l1Mass = 0;
  l2Mass = 0;
  jetMass = 0;
  jetPt = 0;
  jetEta =0;

*/
  // cut-based ID Fall17 V2 (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
  /*   
  nCBVIDTightElec = std::count(Electron_cutBased, Electron_cutBased+nElectron, 4);
  nCBVIDVetoElec = std::count(Electron_cutBased, Electron_cutBased+nElectron, 1);
  nTightIdMuon = std::count(Muon_tightId, Muon_tightId+nMuon, true);
  nMediumIdMuon = std::count(Muon_mediumId, Muon_mediumId+nMuon, true);
  */
  ///////////////////////////// counting mva electrons  /////////////////////////////////////////////////////////
  /*
  nMVAElec =0;
  std::vector<size_t>Electron_mva ={};
  for (size_t i=0; i<nElectron; i++) {
    if (std::abs(Electron_eta[i]<0.8)) {
      if ((10<Electron_pt[i]) && (Electron_pt[i]<25)) {
	if (Electron_mvaFall17V2Iso[i]>(0.2 +0.032*(Electron_pt[i]-10))) {  
	  nMVAElec +=1;
	  Electron_mva.push_back(1);
	}
	else 
	  Electron_mva.push_back(0);
      } 

      else if (Electron_pt[i]>25) {
	if (Electron_mvaFall17V2Iso[i]>0.68) {
	  nMVAElec +=1;
	  Electron_mva.push_back(1);
	}
	else 
	  Electron_mva.push_back(0);
      }
      else 
	Electron_mva.push_back(0);
    }
    else if ((0.8<std::abs(Electron_eta[i])) && (std::abs(Electron_eta[i])<1.479)) {
      if ((10<Electron_pt[i]) && (Electron_pt[i]<25)) {
	if (Electron_mvaFall17V2Iso[i]>(0.1+0.025*(Electron_pt[i]-10))) {
	  nMVAElec +=1;
	  Electron_mva.push_back(1);
	}
	else 
	  Electron_mva.push_back(0);
      }
      else if (Electron_pt[i]>25) {
	if (Electron_mvaFall17V2Iso[i]>0.475) {
	  nMVAElec +=1; 
	  Electron_mva.push_back(1);
	}
	else 
	  Electron_mva.push_back(0);
      }
      else 
	Electron_mva.push_back(0);
    }
    else if ((1.479<std::abs(Electron_eta[i])) && (std::abs(Electron_eta[i])<2.5)) {
      if ((10<Electron_pt[i]) && (Electron_pt[i]<25)) {
	if (Electron_mvaFall17V2Iso[i]>(-0.1+0.028*(Electron_pt[i]-10))) {
	  nMVAElec +=1; 
	  Electron_mva.push_back(1);
	}
	else 
	  Electron_mva.push_back(0);
      }
      else if (Electron_pt[i]>25) {
	if (Electron_mvaFall17V2Iso[i]>0.320) {
	  nMVAElec +=1; 
	  Electron_mva.push_back(1);
	}
	else 
	  Electron_mva.push_back(0);
      }
      else 
	Electron_mva.push_back(0);
    }
    else
      Electron_mva.push_back(0);
  }
  */
  for (size_t i=0; i<nElectron; i++) {
    if(IsGoodMVAElec(i)) {
      goodParts.push_back(GoodPart(Electron_pt[i], Electron_eta[i], Electron_phi[i], Electron_mass[i]));
      goodParts.back().SetPdgId(Electron*Electron_charge[i]);
    }
  }

  for (size_t i=0; i<nMuon; i++) {
    if(IsGoodMuon(i)) {
      goodParts.push_back(GoodPart(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]));
      goodParts.back().SetPdgId(Muon*Muon_charge[i]);
    }
  }
  for(size_t i = 0; i < nJet; i++) {
    if(goodParts.size() != 2) break;

    //common to both jets, may put more here
    if(!isLooseJetId(i)) continue;

    // regular jets
    if(IsGoodJet(i)) {
      nTightJet++;
      Ht += Jet_pt[i];
    }

    // bjet 
    if(IsGoodBJet(i)) nbJet++;
  }

  channel_ = channelMap_[channelName_];
  if(goodParts.size() != 2) {
    channel_ = Unknown;
    channelName_ = "Unknown";
  } else if(goodParts[0].Id() == Muon && goodParts[1].Id() == Muon) {
    channel_ = mm;
    channelName_ = "mm";
  } else if(goodParts[0].Id() == Electron && goodParts[1].Id() == Electron) {
    channel_ = ee;
    channelName_ = "ee";
  } else {
    channel_ = em;
    channelName_ ="em";
    /// fix order of leptons by pt
    if(goodParts[0].Pt() < goodParts[1].Pt()) {
      std::swap(goodParts[0], goodParts[1]);
    }
  }


  //nTightIsoMuon = std::count(Muon_pfIsoId, Muon_pfIsoId+nMuon, 4);
  //nLooseIsoMuon = std::count(Muon_pfIsoId, Muon_pfIsoId+nMuon, 1);
   
  // addition of Ht variable; need to make good Jet array for jets that pass cuts;
  // if it doesn't pass, to keep same indices, make it 0
  /*
  for (size_t i =0; i<N_KEEP_JET; i++) {
    if (Jet_pt[i] > 40) {
      if (std::abs(Jet_eta[i]) < 2.4) {  
	goodJet[i] = Jet_pt[i];
	jetMass = Jet_mass[i];
	jetPt = Jet_pt[i];
	jetEta = Jet_eta[i];

      }
    }  
    else
      goodJet[i]=0;
  }
  // calculating Ht, number of good jets, number of bJets (CSVV2 and deepB)
 
  ngoodJet = 0;
  for (auto& num : goodJet) {
    Ht += num;  
    if (num != 0) 
      ngoodJet += 1; 
  }

  //nbJet =0;
  for (size_t id =0; id <N_KEEP_JET; id++) {
    if (Jet_btagCSVV2[id] > 0.9693) {
      if (Jet_pt[id] > 40) { 
	if (std:: abs(Jet_eta[id]) <2.4) { 
	  nbJet += 1; 
	}
      } 
    }
  }
    
  ndeepbJet =0;
  for (size_t id =0; id <N_KEEP_JET; id++) {
    if (Jet_btagDeepB[id] > 0.7489) {
      if (Jet_pt[id]>40) {
	if (std::abs(Jet_eta[id]) <2.4) {
	  ndeepbJet +=1; 
	}
      }
    }
  }


  /////////////////////////////// creating channels em, me, ee, mm  //////////////////////////////////////////////////////
  channel_ = channelMap_[channelName_];
  std::vector<size_t> goodIndices = {};
   
  if (nMediumIdMuon == 2) {
    channel_=mm;
    if(!(Muon_mediumId[0] && Muon_mediumId[1])) {
      for (size_t i=0; i<nMuon; i++) {
	if (Muon_mediumId[i]) 
	  goodIndices.push_back(i);
      } 
      if (goodIndices.size()<2) {
	return;
      }
    }
    else 
      goodIndices = {0,1};
    if (Muon_pt[goodIndices[1]]>Muon_pt[goodIndices[0]]) {
      std::reverse(goodIndices.begin(), goodIndices.end());
    }
      
    if (Muon_charge[goodIndices[0]] == Muon_charge[goodIndices[1]]) {
      l1Pt = Muon_pt[goodIndices[0]];
      l2Pt = Muon_pt[goodIndices[1]];
      l1Eta = Muon_eta[goodIndices[0]];
      l2Eta = Muon_eta[goodIndices[1]];
      l1Phi = Muon_phi[goodIndices[0]];
      l2Phi = Muon_phi[goodIndices[1]];
      l1Mass = Muon_mass[goodIndices[0]];
      l2Mass = Muon_mass[goodIndices[1]];
      l1IsMed = (Muon_mediumId[goodIndices[0]] && (Muon_pfRelIso04_all[goodIndices[0]] < 0.15));
      l2IsMed = (Muon_mediumId[goodIndices[1]] && (Muon_pfRelIso04_all[goodIndices[1]] < 0.15));
    }
  }
  // currently doing MVA Electrons --not cut based 
  else if (nMVAElec == 2) {
    channel_ = ee;
    if(!(Electron_mva[0]==1 && Electron_mva[1]==1)) {
      for (size_t i=0; i<nElectron; i++) {
	if (Electron_mva[i]==1) 
	  goodIndices.push_back(i);
      }
      if (goodIndices.size() <2)
	return; 
    }
    else 
      goodIndices = {0,1};

    if (Electron_pt[goodIndices[1]]>Electron_pt[goodIndices[0]]) {
      std::reverse(goodIndices.begin(), goodIndices.end());
    }

    if (Electron_charge[goodIndices[0]] == Electron_charge[goodIndices[1]]) {
      l1Pt = Electron_pt[goodIndices[0]];
      l2Pt = Electron_pt[goodIndices[1]];
      l1Eta = Electron_eta[goodIndices[0]];
      l2Eta = Electron_eta[goodIndices[1]];
      l1Phi = Electron_phi[goodIndices[0]];
      l2Phi = Electron_phi[goodIndices[1]];
      l1Mass = Electron_mass[goodIndices[0]];
      l2Mass = Electron_mass[goodIndices[1]];
      l1IsTight = (Electron_mva[goodIndices[0]] == 1);
      l2IsTight = (Electron_mva[goodIndices[1]] == 1);
      
    }
  }
      
  // Cut Based Electrons currently commited out
  /*
    else if (nCBVIDTightElec == 2) {
      channel_ = ee;
      if( Electron_cutBased[0]==4 && Electron_cutBased[1]==4) {
      for (size_t i=0; i<nElectron; i++) {
        if (Electron_cutBased[i]==4) 
	    goodIndices.push_back(i);
	    }
	    if (goodIndices.size() <2)
	      return; 
      }
      else 
      goodIndices = {0,1};
      if (Electron_pt[goodIndices[1]]>Electron_pt[goodIndices[0]]) {
      std::reverse(goodIndices.begin(), goodIndices.end());
      }
      if (Electron_charge[goodIndices[0]] == Electron_charge[goodIndices[1]]) {
      l1Pt = Electron_pt[goodIndices[0]];
      l2Pt = Electron_pt[goodIndices[1]];
      l1Eta = Electron_eta[goodIndices[0]];
      l2Eta = Electron_eta[goodIndices[1]];
      l1Phi = Electron_phi[goodIndices[0]];
      l2Phi = Electron_phi[goodIndices[1]];
      l1Mass = Electron_mass[goodIndices[0]];
      l2Mass = Electron_mass[goodIndices[1]];
      l1IsTight = (Electron_cutBased[goodIndices[0]] == 4);
      l2IsTight = (Electron_cutBased[goodIndices[1]] == 4);
      
      }
    }
  

  // changed to mva electron here as well
  else if ((nMVAElec == 1) && (nMediumIdMuon ==1)) {
    if(!(Electron_mva[0] == 1 && Muon_mediumId[0])) {
      for (size_t i=0; i<nElectron; i++) {
	if (Electron_mva[0]==1) {
	  goodIndices.push_back(i);
	}
      }
      for (size_t i=0; i<nMuon; i++) {
	if (Muon_mediumId[i]) {
	  goodIndices.push_back(i);
	}
      }
      if (goodIndices.size()<2)
	return; 
    }
    else 
      goodIndices = {0,0}; 

    if (Muon_pt[goodIndices[1]] > Electron_pt[goodIndices[0]] ) {
      channel_= me;
      std::reverse(goodIndices.begin(), goodIndices.end());
    }

    else if (Electron_pt[goodIndices[0]] > Muon_pt[goodIndices[1]]) {
      channel_ = em;
    }
    
    if (Electron_charge[goodIndices[0]] == Muon_charge[goodIndices[1]]) {
      l1Pt = Electron_pt[goodIndices[0]];
      l2Pt = Muon_pt[goodIndices[1]];
      l1Eta = Electron_eta[goodIndices[0]];
      l2Eta = Muon_eta[goodIndices[1]];
      l1Phi = Electron_phi[goodIndices[0]];
      l2Phi = Muon_phi[goodIndices[1]];
      l1Mass = Electron_mass[goodIndices[0]];
      l2Mass = Muon_mass[goodIndices[1]];
      l1IsTight = (Electron_mva[goodIndices[0]] == 1);
      l2IsMed = (Muon_mediumId[goodIndices[1]] && (Muon_pfRelIso04_all[goodIndices[1]] < 0.15) );
    }
  }
*/
  if (isMC_) {
    ApplyScaleFactors();

    //  genWeight->GetEntry(entry);
    //TODO: add scale factors
    //b_numPU->GetEntry(entry);
    //b_l1GenPt->GetEntry(entry);
    //b_l2GenPt->GetEntry(entry);
    //b_l3GenPt->GetEntry(entry);
    //ApplyScaleFactors();
  }
  // else {
    //TODO: add MET filters
    //b_Flag_duplicateMuonsPass->GetEntry(entry);          
    //b_Flag_badMuonsPass->GetEntry(entry);          
  
  
  // originally has nMediumIdMuon but using medium id...should we change to tight here? 
  // passesLeptonVeto = (std::min(nTightIdMuon, nLooseIsoMuon) + nCBVIDVetoElec) == 2;
   
}
void SameSignLept::ApplyScaleFactors() {
  weight *= (genWeight>0) ? 1.0: -1.0;
  return;
}

bool SameSignLept::IsGoodMuon(size_t index) {
  return ( (Muon_pt[index] > 20) &&
	   (abs(Muon_eta[index]) <2.4) &&
	   Muon_mediumId[index] && 
	   (Muon_miniPFRelIso_all[index]<0.16) );
} 
// tight electron == 4 ; med=3; loose = 2; 
bool SameSignLept::IsGoodCBElec(size_t index) {
  return ((Electron_pt[index] > 20) &&
	  (abs(Electron_eta[index]) < 2.5) && 
	  (Electron_miniPFRelIso_all[index] < 0.12) && 
	  (Electron_cutBased[index] ==4));
} 
//currently 2017
bool SameSignLept::IsGoodMVAElec(size_t index) {
  bool mvaRec = false;
  if(abs(Electron_eta[index]) < 0.8)
    mvaRec = std::max(0.52, 0.77-0.025*(Electron_pt[index]-15));
  else if(abs(Electron_eta[index]) < 1.479)
    mvaRec = std::max(0.11, 0.56-0.045*(Electron_pt[index]-15));
  else if(abs(Electron_eta[index]) < 2.5)
    mvaRec = std::max(-0.01, 0.48-0.049*(Electron_pt[index]-15));
  return ((Electron_pt[index] > 20) &&
	  mvaRec);
}

bool SameSignLept::IsGoodJet(size_t index) {
  return ((Jet_pt[index] > 40.0) &&
	  (abs(Jet_eta[index]) < 2.4) &&
	  // isOverlap(index)
	  );
}

/// TODO: add toggle for different btag stuff
bool SameSignLept::IsGoodBJet(size_t index) {
  return ((Jet_pt[index] > 25.0) &&
	  (abs(Jet_eta[index]) < 2.4) &&
	  (Jet_btagCSVV2[index] > 0.8484) &&  
	  // (Jet_btagDeepB[index] > 0.6324) &&
	  // isOverlap(index)
	  );
}


bool SameSignLept::isLooseJetId(size_t index) {
  return (Jet_neHEF[index] < 0.99 &&
	    Jet_neEmEF[index] < 0.99 &&
	    Jet_nConstituents[index] > 1 &&
	    Jet_chHEF[index] > 0 &&
	    Jet_chEmEF[index] < 0.99
	  );
}

bool SameSignLept::isTightJetId(size_t index) {
  return (Jet_neHEF[index] < 0.9 &&
	    Jet_neEmEF[index] < 0.9 &&
	    Jet_nConstituents[index] > 1 &&
	    Jet_chHEF[index] > 0
	  );
}


///  Filling Histograms
void SameSignLept::FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) { 
  int step = 0;
  Fill1D("CutFlow", 0);
  
  /// 2 good leptons
  if(goodParts.size() != 2) return;
  Fill1D("CutFlow", ++step);

  // first lep requirement
  if(goodParts[0].Pt() < 25) return;
  Fill1D("CutFlow", ++step);
  
  // same sign requirement
  if(goodParts[0].Charge() * goodParts[1].Charge() <= 0) return;
  Fill1D("CutFlow", ++step);

  // met cut
  if (MET < 50) return;
  Fill1D("CutFlow", ++step);

  // ht cut
  if(Ht < 300 ) return;
  Fill1D("CutFlow", ++step);
  
  // jet cut
  if(nTightJet < 4) return;
  Fill1D("CutFlow", ++step);
  
  // bjet cut
  // if(nBJets < 2) return;
  // Fill1D("CutFlow", ++step);
  
  // // veto cut
  // if(!passesLeptonVeto)
  //   Fill1D("CutFlow", ++step);
  // in SR stuff

  // met cut
  if (MET < 50) return;
  SafeHistFill(histMap1D_, getHistName("CutFlow", variation.second), step++, weight);

  Fill1D("ptl1", goodParts[0].Pt());
  Fill1D("ptl2", goodParts[1].Pt());
  //  Fill1D("SR", getSRBin());
  Fill1D("njet", nTightJet);
  Fill1D("nbjet", nbJet);


  
  for(size_t i = 0; i < nJet; i++) {
    if(IsGoodJet(i)) {
      Fill1D("jetpt", Jet_pt[i]);
    }
    if(IsGoodBJet(i)) {
      Fill1D("bjetpt", Jet_pt[i]);
    }
  }
}


void SameSignLept::SetupNewDirectory() {
  SelectorBase::SetupNewDirectory();

  InitializeHistogramsFromConfig();

  /*
  AddObject<TH1D>(cutflow_ee_, "cutflow_ee", "Tight leptons; Cut flow", 6, 0, 6);
  AddObject<TH1D>(cutflow_mm_, "cutflow_mm", "Tight leptons; Cut flow", 6, 0, 6);
  AddObject<TH1D>(cutflow_em_, "cutflow_em", "Tight e, Med. m; Cut flow", 6, 0, 6);
  AddObject<TH1D>(cutflow_me_, "cutflow_me", "Tight e, Med. m; Cut flow", 6, 0, 6);
  AddObject<TH1D>(CombMass_ee_, "CombMass_ee", "Tight leptons; m_{ee} [GeV]", 80, 102, 22);
  AddObject<TH1D>(CombMass_mm_, "CombMass_mm", "Tight leptons; m_{#mu#mu} [GeV]", 80, 102, 22);
  AddObject<TH1D>(CombMass_em_, "CombMass_em", "Tight leptons; m_{e#mu} [GeV]", 80, 102, 22);
  AddObject<TH1D>(CombMass_me_, "CombMass_me", "Tight leptons; m_{#mu e} [GeV]", 80, 102, 22);
  
  AddObject<TH1D>(elect_mass, "Electron mass", "Tight electrons; m_{e} [GeV]",50, 300, 0);
  AddObject<TH1D>(muon_mass, "Muon mass", "Medium muons; m_{#mu} [GeV]",50, 300, 0);
  AddObject<TH1D>(jet_mass, "Jet mass", "Jets; m_{j} [GeV]", 30, 300, 0);
  
  AddObject<TH1D>(ptl1_ee_, "ptl1_ee", "Tight leptons; p_{T}(e_{1}) [GeV]", 100, 0, 200);
  AddObject<TH1D>(ptl1_mm_, "ptl1_mm", "Tight leptons; p_{T}(#mu_{1}) [GeV]", 100, 0, 100);
  AddObject<TH1D>(ptl2_ee_, "ptl2_ee", "Tight leptons; p_{T}(e_{2}) [GeV]", 100, 0, 200);
  AddObject<TH1D>(ptl2_mm_, "ptl2_mm", "Tight leptons; p_{T}(#mu_{2}) [GeV]", 100, 0, 100);
  
  AddObject<TH1D>(elect_pt, "Electron pt", "Tight electrons; p_{T}(e) [GeV]", 100, 0, 250);
  AddObject<TH1D>(muon_pt, "Muon pt", "Medium muons; p_{T}(#mu) [GeV]", 100, 0 ,250);
  AddObject<TH1D>(jet_pt, "Jet pt", "Jets; p_{T}(j) [Gev]", 100, 0, 250);

  AddObject<TH1D>(l1eta_mm_, "l1eta_mm", "Tight leptons; eta(#mu_{1})", 100, -5, 5);
  AddObject<TH1D>(l2eta_mm_, "l2eta_mm", "Tight leptons; eta(#mu_{2})", 100, -5, 5);
  AddObject<TH1D>(l1eta_ee_, "l1eta_ee", "Tight leptons; eta(e_{1})", 100, -5, 5);
  AddObject<TH1D>(l2eta_ee_, "l2eta_ee", "Tight leptons; eta(e_{2})", 100, -5, 5);
  
  AddObject<TH1D>(elect_eta, "Electron eta", "Tight electrons; eta(e)", 100, -5, 5);
  AddObject<TH1D>(muon_eta, "Muon eta", "Medium muons; eta(#mu)", 100, -5, 5);
  AddObject<TH1D>(jet_eta, "Jet eta", "Jets; eta(j)", 100, -5, 5);
  
  AddObject<TH1D>(l1phi_mm_, "l1phi_mm", "Tight leptons; #phi(#mu_{1})", 100, -5, 5);
  AddObject<TH1D>(l2phi_mm_, "l2phi_mm", "Tight leptons; phi(#mu_{2})", 100, -5, 5);
  AddObject<TH1D>(l1phi_ee_, "l1phi_ee", "Tight leptons; phi(e_{1})", 100, -5, 5);
  AddObject<TH1D>(l2phi_ee_, "l2phi_ee", "Tight leptons; phi(e_{2})", 100, -5, 5);
  
  AddObject<TH1D>(elect_phi, "Electron phi", "Tight electrons; phi(e)", 100, -5, 5); 
  AddObject<TH1D>(muon_phi, "Muon phi", "Medium muons; phi(#mu)", 100, -5,5);
  
  AddObject<TH1D>(SR_uw, "nEvents", "Signal Regions; regions", 8, 0, 8);
  AddObject<TH1D>(n_bjet_, "nbJets", "B Jets; number", 10, 0, 10);
  AddObject<TH1D>(n_deepbjet_, "ndeepbJets", "B Jets; number", 10, 0, 10);
  AddObject<TH1D>(n_cut_elec, "nCutElec", "Cut Based Elec; number", 10, 0, 10);
  AddObject<TH1D>(n_mva_elec, "nMVAElec", "MVA Elec; number", 10, 0, 10);
  */
} 
