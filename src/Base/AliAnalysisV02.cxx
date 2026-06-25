/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

#include "AliAnalysisV02.h"

#include "Riostream.h"

using std::cout;
using std::endl;

#define Fatal(d, s) printf("***  FATAL: %s: %s\n",d,s)

//ClassImp(AliAnalysisV02)

// ------------------------------------------------------------------------- //
AliAnalysisV02::AliAnalysisV02():
  fInputList(NULL),
  fHistList(NULL),
  fDebugLevel(0),
  fCentrality(0.),
  fUseMultiplicityFlowWeights(kTRUE),
  fptSubMax(0.8), fptSubMin(0.4),
  fPtBins({0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.5, 5})
{
// Dummy constructor of the class.
  printf("AliAnalysisV02::AliAnalysisV02()\n");
  this->InitializeArrays();
}

// ------------------------------------------------------------------------- //
AliAnalysisV02::AliAnalysisV02(const char *name):
  fInputList(NULL),
  fHistList(NULL),
  fDebugLevel(0),
  fCentrality(0.),
  fUseMultiplicityFlowWeights(kTRUE),
  fptSubMax(0.8), fptSubMin(0.4),
  fPtBins({0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.5, 5})
{
// Constructor of the class.
  printf("AliAnalysisV02::AliAnalysisV02(const char *name)\n");

  // Define the base list.
  fHistList = new TList();
  fHistList->SetName(name);
  fHistList->SetOwner(kTRUE);

  // Initialize all arrays within the data members.
  InitializeArrays();

  //DefineOutput(1, TList::Class());
}


// ------------------------------------------------------------------------- //
AliAnalysisV02::~AliAnalysisV02()
{
// Destructor of the class.
  if (fHistList) {delete fHistList;}
}

void AliAnalysisV02::InitializeArrays() {
  for (Int_t i = 0; i < nSpecies; ++i) {
    fV02Profile[i] = NULL;
    fPtNchProfile[i] = NULL;
    fV22Profile[i] = NULL;
  }
}

void AliAnalysisV02::BookAndNestAllLists() {
  const char *speciesNames[nSpecies] = {"hadron", "pion", "kaon", "proton"};

  for (Int_t i = 0; i < nSpecies; ++i) {
    TString nameV02 = Form("V02corr_%s", speciesNames[i]);
    TString titleV02 = Form("V02corr %s", speciesNames[i]);
    fV02Profile[i] = new TProfile2D(nameV02.Data(), titleV02.Data(),
                                    nHistCent, 0., float(nHistCent),
                                    nPtBins, fPtBins);

    TString namePtNch = Form("PtNch_%s", speciesNames[i]);
    TString titlePtNch = Form("PtNch %s", speciesNames[i]);
    fPtNchProfile[i] = new TProfile2D(namePtNch.Data(), titlePtNch.Data(),
                                      nHistCent, 0., float(nHistCent),
                                      nPtBins, fPtBins);

    TString nameV22 = Form("V22_%s", speciesNames[i]);
    TString titleV22 = Form("V22 %s", speciesNames[i]);
    fV22Profile[i] = new TProfile(nameV22.Data(), titleV22.Data(),
                                  nHistCent, 0., float(nHistCent));

    fHistList->Add(fV02Profile[i]);
    fHistList->Add(fPtNchProfile[i]);
    fHistList->Add(fV22Profile[i]);
  }

}
// ------------------------------------------------------------------------- //
void AliAnalysisV02::UserCreateOutputObjects() {
  this->BookAndNestAllLists();
}

// ------------------------------------------------------------------------- //
void AliAnalysisV02::WriteLists(TString name) {
// Write the corresponding outputlists.
  fHistList->Write(name.Data(),1);
  delete fHistList;
}

Int_t AliAnalysisV02::getPIDIndex(AliJBaseTrack *aTrack) {
  Int_t pidIndex = 0;
  if (aTrack->GetParticleType() == 0x211 || aTrack->GetParticleType() == -0x211) {pidIndex = 1;}
  else if (aTrack->GetParticleType() == 0x321 || aTrack->GetParticleType() == -0x321) {pidIndex = 2;}
  else if (aTrack->GetParticleType() == 0x2212 || aTrack->GetParticleType() == -0x2212) {pidIndex = 3;}
  else {pidIndex = 0;}
  return pidIndex;
}
// ------------------------------------------------------------------------- //
void AliAnalysisV02::UserExec(Option_t *option) {
  Int_t centralityBin = SelectCentrality(fCentrality);
  TH1F *hPtMid[nSpecies];
  for (Int_t iSpec = 0; iSpec < nSpecies; ++iSpec) {
    hPtMid[iSpec] = new TH1F(Form("hPtMid_%d", iSpec), "hPtMid", nPtBins, &fPtBins[0]);
  }
  Int_t nTracks = fInputList->GetEntriesFast();
  
  TComplex Qpos[nSpecies][2];
  TComplex Qneg[nSpecies][2];
  for (Int_t iSpec = 0; iSpec < nSpecies; ++iSpec) {
    Qpos[iSpec][0] = TComplex(0., 0.);
    Qpos[iSpec][1] = TComplex(0., 0.);
    Qneg[iSpec][0] = TComplex(0., 0.);
    Qneg[iSpec][1] = TComplex(0., 0.);
  }

  if (fDebugLevel>10) printf("Centrality: %d, nTracks: %d\n", centralityBin, nTracks);

  // Track loop
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++)
  {
    Double_t iPhi = 0.;
    Double_t iPt = 0.;
    Double_t iEta = 0.;

    // Parse the selected track and get their information into the arrays.
    AliJBaseTrack *aTrack = (AliJBaseTrack*)fInputList->At(iTrack);
    if (!aTrack) {continue;}

    iPhi = aTrack->Phi();
    iPt = aTrack->Pt();
    iEta = aTrack->Eta();

    Int_t pidIndex = getPIDIndex(aTrack); // 1 for pion, 2 for kaon, 3 for proton
    if(std::abs(iEta)<0.4){
      // Currently fill only hadron species (index 0); others can be used
      // once PID selection is added.
      hPtMid[0]->Fill(iPt);
      if(pidIndex > 0) {
        hPtMid[pidIndex]->Fill(iPt);
      }
    }
    if(iEta < -fptSubMin){
      Qneg[0][0] += 1;
      Qneg[0][1] += TComplex(TMath::Cos(-2*iPhi), TMath::Sin(-2*iPhi));
      if(pidIndex > 0) {
        Qneg[pidIndex][0] += 1;
        Qneg[pidIndex][1] += TComplex(TMath::Cos(-2*iPhi), TMath::Sin(-2*iPhi));
      }
    }
    if(iEta > fptSubMin){
      Qpos[0][0] += 1;
      Qpos[0][1] += TComplex(TMath::Cos(2*iPhi), TMath::Sin(2*iPhi));
      if(pidIndex > 0) {
        Qpos[pidIndex][0] += 1;
        Qpos[pidIndex][1] += TComplex(TMath::Cos(2*iPhi), TMath::Sin(2*iPhi));
      }
    }
  } // Go to the next iTrack.

  // For now, fill only the "hadron" species (index 0).
  const Int_t iSpeciesHadron = 0;

  Double_t corr[nSpecies];
  for (Int_t iSpec = 0; iSpec < nSpecies; ++iSpec) {
    corr[iSpec] = (Qpos[iSpec][1]*Qneg[iSpec][1]).Re()/((Qpos[iSpec][0]*Qneg[iSpec][0]).Re()); //! <v22 gap>

    for (uint i=0; i<nPtBins; i++) {
      double ptFraction = 0.;
      if (hPtMid[iSpec]->Integral() > 0) {
        ptFraction = hPtMid[iSpec]->GetBinContent(i) / hPtMid[0]->Integral();
      }
      double ptVal = (fPtBins[i] + fPtBins[i+1])/2.;
      double v02 = corr[0]*ptFraction; // Use integrated v22 for all species
      if (TMath::Finite(v02)) {
        fV02Profile[iSpec]->Fill((Float_t)(centralityBin)+.5, ptVal, v02, (fUseMultiplicityFlowWeights) ? (Qpos[0][0]*Qneg[0][0]).Re() : 1.0); // v22 integrated and ptfraction from species
        fPtNchProfile[iSpec]->Fill((Float_t)(centralityBin)+.5, ptVal, ptFraction); // Fill for all species
      }
    }
    if (TMath::Finite(corr[iSpec])) fV22Profile[iSpec]->Fill((Float_t)(centralityBin)+.5, corr[iSpec], (fUseMultiplicityFlowWeights) ? (Qpos[iSpec][0]*Qneg[iSpec][0]).Re() : 1.0); // FIll for all species
  }

  // Reset event-by-event objects.
  for (Int_t iSpec = 0; iSpec < nSpecies; ++iSpec) {
    hPtMid[iSpec]->Reset();
    delete hPtMid[iSpec];
  }
  nTracks = 0;
}

// ------------------------------------------------------------------------- //
void AliAnalysisV02::Terminate(Option_t *option)
{
// Execute the additional steps after the loop over the events is over.
}


// ------------------------------------------------------------------------- //
Int_t AliAnalysisV02::SelectCentrality(Double_t centValue)
{
// Get the centrality bin corresponding to the given centrality value.
  Float_t centralityArray[10] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80.};
  // Float_t centralityArray[12] = {0., 2., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50.};
  for (Int_t iCent = 0; iCent < nHistCent; iCent++)
  {
    if (centralityArray[iCent] < 0) {return -1;}
    if (centValue >= centralityArray[iCent]) {continue;}
    else {return iCent-1;}
  }

  // We went through all centrality edges without returning. This means that the measured value
  // is bigger than the maximum centrality we want for our analysis.
  return -1;
}
