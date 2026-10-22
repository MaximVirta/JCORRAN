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

#include "AliAnalysisPtVnWeighted.h"

#include "Riostream.h"
#include "TMath.h"

using std::cout;
using std::endl;

#define Fatal(d, s) printf("***  FATAL: %s: %s\n",d,s)

//ClassImp(AliAnalysisPtVnWeighted)

// ------------------------------------------------------------------------- //
AliAnalysisPtVnWeighted::AliAnalysisPtVnWeighted():
  fInputList(NULL),
  fHistList(NULL),
  fDebugLevel(0),
  fCentrality(0.),
  fUseMultWeight(kTRUE),
  fptSubMax(0.8), fptSubMin(0.4),
  fCentralityBins(16)
{
// Dummy constructor of the class.
  printf("AliAnalysisPtVnWeighted::AliAnalysisPtVnWeighted()\n");
  this->InitializeArrays();
}

// ------------------------------------------------------------------------- //
AliAnalysisPtVnWeighted::AliAnalysisPtVnWeighted(const char *name):
  fInputList(NULL),
  fHistList(NULL),
  fDebugLevel(0),
  fCentrality(0.),
  fUseMultWeight(kTRUE),
  fptSubMax(0.8), fptSubMin(0.4),
  fCentralityBins(9)
{
// Constructor of the class.
  printf("AliAnalysisPtVnWeighted::AliAnalysisPtVnWeighted(const char *name)\n");

  // Define the base list.
  fHistList = new TList();
  fHistList->SetName(name);
  fHistList->SetOwner(kTRUE);

  // Initialize all arrays within the data members.
  InitializeArrays();

  //DefineOutput(1, TList::Class());
}


// ------------------------------------------------------------------------- //
AliAnalysisPtVnWeighted::~AliAnalysisPtVnWeighted()
{
// Destructor of the class.
  if (fHistList) {delete fHistList;}
}

void AliAnalysisPtVnWeighted::InitializeArrays() {
  for (int i = 0; i < nCorr; ++i) {
    corrProfiles[i] = NULL;
  }


}

void AliAnalysisPtVnWeighted::BookAndNestAllLists() {
  for (int i = 0; i < nCorr; ++i) {
    corrProfiles[i] = new TProfile(Form("profile_%s", corrlistnames[i].Data()), corrlistnames[i].Data(),
      nHistCent, 0, float(nHistCent), "s");
    fHistList->Add(corrProfiles[i]);
  }
}

// ------------------------------------------------------------------------- //
void AliAnalysisPtVnWeighted::UserCreateOutputObjects() {
  this->BookAndNestAllLists();
}

// ------------------------------------------------------------------------- //
void AliAnalysisPtVnWeighted::WriteLists(TString name) {
// Write the corresponding outputlists.
  fHistList->Write(name.Data(),1);
  delete fHistList;
}

// ------------------------------------------------------------------------- //
void AliAnalysisPtVnWeighted::UserExec(Option_t *option) {
  Int_t centralityBin = SelectCentrality(fCentrality);

  Int_t nTracks = fInputList->GetEntriesFast();

  double psubA[3] = {0.,0.,0.};
  double psubB[3] = {0.,0.,0.};
  double psubMid[3] = {0.,0.,0.};
  TComplex QA(0., 0.), QpA(0., 0.), Q4A(0., 0.), Q3A(0., 0.), Q6A(0., 0.);
  TComplex QB(0., 0.), QpB(0., 0.), Q4B(0., 0.), Q3B(0., 0.), Q6B(0., 0.);
  TComplex Qpos[2] = {TComplex(0., 0.)};
  TComplex Qneg[2] = {TComplex(0., 0.)};
  TComplex Q3pos[2] = {TComplex(0., 0.)};
  TComplex Q3neg[2] = {TComplex(0., 0.)};

  if (fDebugLevel>10) printf("Centrality: %d, nTracks: %d\n", centralityBin, nTracks);

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

    if(-fptSubMax < iEta && iEta < -fptSubMin){
      psubA[0] += 1;
      psubA[1] += iPt;
      psubA[2] += iPt*iPt;
      const TComplex e2(TMath::Cos(2.*iPhi), TMath::Sin(2.*iPhi));
      QA += e2;
      QpA += e2*iPt;
      Q4A += TComplex(TMath::Cos(4.*iPhi), TMath::Sin(4.*iPhi));
      Q3A += TComplex(TMath::Cos(3.*iPhi), TMath::Sin(3.*iPhi));
      Q6A += TComplex(TMath::Cos(6.*iPhi), TMath::Sin(6.*iPhi));
    }
    if(fptSubMin < iEta && iEta < fptSubMax){
      psubB[0] += 1;
      psubB[1] += iPt;
      psubB[2] += iPt*iPt;
      const TComplex e2(TMath::Cos(2.*iPhi), TMath::Sin(2.*iPhi));
      QB += e2;
      QpB += e2*iPt;
      Q4B += TComplex(TMath::Cos(4.*iPhi), TMath::Sin(4.*iPhi));
      Q3B += TComplex(TMath::Cos(3.*iPhi), TMath::Sin(3.*iPhi));
      Q6B += TComplex(TMath::Cos(6.*iPhi), TMath::Sin(6.*iPhi));
    }
    if(std::abs(iEta)<0.4){
      psubMid[0] += 1;
      psubMid[1] += iPt;
      psubMid[2] += iPt*iPt;
    }
    if(-fptSubMax < iEta && iEta < -fptSubMin){
      Qneg[0] += 1;
      Qneg[1] += {TMath::Cos(-2.*iPhi),TMath::Sin(-2.*iPhi)};
    }
    if(fptSubMin < iEta && iEta < fptSubMax){
      Qpos[0] += 1;
      Qpos[1] += {TMath::Cos(2.*iPhi),TMath::Sin(2.*iPhi)};
    }
    if(-fptSubMax < iEta && iEta < -fptSubMin){
      Q3neg[0] += 1;
      Q3neg[1] += {TMath::Cos(-3.*iPhi),TMath::Sin(-3.*iPhi)};
    }
    if(fptSubMin < iEta && iEta < fptSubMax){
      Q3pos[0] += 1;
      Q3pos[1] += {TMath::Cos(3.*iPhi),TMath::Sin(3.*iPhi)};
    }

  } // Go to the next iTrack.
  const Double_t mA = psubA[0];
  const Double_t mB = psubB[0];
  const Double_t mMid = psubMid[0];
  // 1-, 2-, 3- and 4-particle weights of the A/B subevent pT moments.
  // The two 3-particle topologies are not interchangeable.
  const Double_t w1A = mA;
  const Double_t w1B = mB;
  const Double_t w2AB = mA*mB;
  const Double_t w2AA = mA*(mA - 1.);
  const Double_t w2BB = mB*(mB - 1.);
  const Double_t w2AMid = mA*mMid;
  const Double_t w2MidB = mMid*mB;
  const Double_t w3AAB = mA*(mA - 1.)*mB;
  const Double_t w3ABB = mA*mB*(mB - 1.);
  const Double_t w3AMidB = mA*mB*mMid;
  const Double_t w4AABB = w2AA*w2BB;
  const Double_t w2Mid = mMid*(mMid - 1.);

  const Double_t wGap22 = (Qpos[0]*Qneg[0]).Re();
  const Double_t wGap22pt = wGap22*mMid;
  const Double_t wGap32 = (Q3pos[0]*Q3neg[0]).Re();
  const Double_t wGap32pt = wGap32*mMid;

  const Double_t ptA = psubA[1]/mA;
  const Double_t ptB = psubB[1]/mB;
  const Double_t ptMid = psubMid[1]/mMid;
  const Double_t pt2AB = psubA[1]*psubB[1]/w2AB;
  const Double_t pt2AA = (psubA[1]*psubA[1] - psubA[2])/w2AA;
  const Double_t pt2BB = (psubB[1]*psubB[1] - psubB[2])/w2BB;
  const Double_t pt2AMid = psubA[1]*psubMid[1]/w2AMid;
  const Double_t pt2MidB = psubB[1]*psubMid[1]/w2MidB;
  const Double_t pt3AAB = pt2AA*ptB;
  const Double_t pt3ABB = pt2BB*ptA;
  const Double_t pt3AMidB = ptA*ptB*ptMid;
  const Double_t pt4AABB = pt2AA*pt2BB;

  Double_t corr[nCorr] = {0.};
  Double_t evtWeight[nCorr] = {0.};

  // <pt> with 1-, 2-, 3- and 4-particle weights. 3-particle: AAB and ABB.
  corr[kPtA] = corr[kPtA_w2] = corr[kPtA_w3AAB] = corr[kPtA_w3ABB] = corr[kPtA_w4] = ptA;
  evtWeight[kPtA] = w1A;
  evtWeight[kPtA_w2] = w2AB;
  evtWeight[kPtA_w3AAB] = w3AAB;
  evtWeight[kPtA_w3ABB] = w3ABB;
  evtWeight[kPtA_w4] = w4AABB;

  corr[kPtB] = corr[kPtB_w2] = corr[kPtB_w3AAB] = corr[kPtB_w3ABB] = corr[kPtB_w4] = ptB;
  evtWeight[kPtB] = w1B;
  evtWeight[kPtB_w2] = w2AB;
  evtWeight[kPtB_w3AAB] = w3AAB;
  evtWeight[kPtB_w3ABB] = w3ABB;
  evtWeight[kPtB_w4] = w4AABB;

  corr[kPtMid] = corr[kPtMid_w3] = ptMid;
  evtWeight[kPtMid] = mMid;
  evtWeight[kPtMid_w3] = wGap22pt;

  // <pt^2> with 2-, 3- and 4-particle weights. 3-particle: AAB and ABB.
  corr[kPt2AB] = corr[kPt2AB_w3AAB] = corr[kPt2AB_w3ABB] = corr[kPt2AB_w4] = pt2AB;
  evtWeight[kPt2AB] = w2AB;
  evtWeight[kPt2AB_w3AAB] = w3AAB;
  evtWeight[kPt2AB_w3ABB] = w3ABB;
  evtWeight[kPt2AB_w4] = w4AABB;

  corr[kPt2AA] = corr[kPt2AA_w3AAB] = corr[kPt2AA_w3ABB] = corr[kPt2AA_w4] = pt2AA;
  evtWeight[kPt2AA] = w2AA;
  evtWeight[kPt2AA_w3AAB] = w3AAB;
  evtWeight[kPt2AA_w3ABB] = w3ABB;
  evtWeight[kPt2AA_w4] = w4AABB;

  corr[kPt2BB] = corr[kPt2BB_w3AAB] = corr[kPt2BB_w3ABB] = corr[kPt2BB_w4] = pt2BB;
  evtWeight[kPt2BB] = w2BB;
  evtWeight[kPt2BB_w3AAB] = w3AAB;
  evtWeight[kPt2BB_w3ABB] = w3ABB;
  evtWeight[kPt2BB_w4] = w4AABB;

  corr[kPt2AMid] = pt2AMid;
  evtWeight[kPt2AMid] = w2AMid;
  corr[kPt2MidB] = pt2MidB;
  evtWeight[kPt2MidB] = w2MidB;

  // <pt^3> with 3- and 4-particle weights. 3-particle fills are AAB and ABB.
  corr[kPt3AAB] = corr[kPt3AAB_w4] = pt3AAB;
  evtWeight[kPt3AAB] = w3AAB;
  evtWeight[kPt3AAB_w4] = w4AABB;
  corr[kPt3ABB] = corr[kPt3ABB_w4] = pt3ABB;
  evtWeight[kPt3ABB] = w3ABB;
  evtWeight[kPt3ABB_w4] = w4AABB;

  corr[kPt3AMidB] = pt3AMidB;
  evtWeight[kPt3AMidB] = w3AMidB;

  // <pt^4> with the 4-particle weight.
  corr[kPt4AABB] = pt4AABB;
  evtWeight[kPt4AABB] = w4AABB;

  // <v_2^2>, <v_2^2 ptA>, <v_2^2 ptB>, <v_2^2 ptA ptB> with the same w4.
  // Flow and pT particles are distinct: one flow leg and one pT leg in A,
  // and the same in B.
  const TComplex flowPtA = QA*psubA[1] - QpA;
  const TComplex flowPtB = QB*psubB[1] - QpB;
  corr[kV22_w4] = (QA*TComplex::Conjugate(QB)).Re()/w2AB;
  corr[kV22ptA_w4] = (flowPtA*TComplex::Conjugate(QB)).Re()/(w2AA*mB);
  corr[kV22ptB_w4] = (QA*TComplex::Conjugate(flowPtB)).Re()/(mA*w2BB);
  corr[kV22ptAptB_w4] = (flowPtA*TComplex::Conjugate(flowPtB)).Re()/w4AABB;
  evtWeight[kV22_w4] = evtWeight[kV22ptA_w4] = evtWeight[kV22ptB_w4] = evtWeight[kV22ptAptB_w4] = w4AABB;

  corr[kPt2mid] = (psubMid[1]*psubMid[1] - psubMid[2])/w2Mid;
  evtWeight[kPt2mid] = w2Mid;

  corr[kV22pt_gap] = psubMid[1]*(Qpos[1]*Qneg[1]).Re()/wGap22pt;
  evtWeight[kV22pt_gap] = wGap22pt;
  corr[kV22_gap] = corr[kV22_w3] = (Qpos[1]*Qneg[1]).Re()/wGap22;
  evtWeight[kV22_gap] = wGap22;
  evtWeight[kV22_w3] = wGap22pt;

  corr[kV32pt_gap] = psubMid[1]*(Q3pos[1]*Q3neg[1]).Re()/wGap32pt;
  evtWeight[kV32pt_gap] = wGap32pt;
  corr[kV32_gap] = corr[kV32_w3] = (Q3pos[1]*Q3neg[1]).Re()/wGap32;
  evtWeight[kV32_gap] = wGap32;
  evtWeight[kV32_w3] = wGap32pt;

  // Two particles in A and two in B. Same four-particle weight as <pt^4>_AABB.
  const TComplex c2A = QA*QA - Q4A;
  const TComplex c2B = QB*QB - Q4B;
  corr[kV24_gap] = (c2A*TComplex::Conjugate(c2B)).Re()/w4AABB;
  evtWeight[kV24_gap] = w4AABB;
  const TComplex c3A = Q3A*Q3A - Q6A;
  const TComplex c3B = Q3B*Q3B - Q6B;
  corr[kV34_gap] = (c3A*TComplex::Conjugate(c3B)).Re()/w4AABB;
  evtWeight[kV34_gap] = w4AABB;

  for (int i = 0; i < nCorr; ++i) {
    FillCorr(centralityBin, i, corr[i], evtWeight[i]);
  }



  // Reset event-by-event objects.
  nTracks = 0;
}

// ------------------------------------------------------------------------- //
void AliAnalysisPtVnWeighted::FillCorr(Int_t centralityBin, Int_t index, Double_t value, Double_t weight) {
  if (!fUseMultWeight) weight = 1.;
  if (!TMath::Finite(value) || !TMath::Finite(weight) || weight <= 0.) return;
  corrProfiles[index]->Fill(static_cast<Float_t>(centralityBin) + 0.5, value, weight);
}

// ------------------------------------------------------------------------- //
void AliAnalysisPtVnWeighted::Terminate(Option_t *option)
{
// Execute the additional steps after the loop over the events is over.
}


// ------------------------------------------------------------------------- //
Int_t AliAnalysisPtVnWeighted::SelectCentrality(Double_t centValue)
{
// Get the centrality bin corresponding to the given centrality value.
  // Float_t centralityArray[10] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80.};
  Float_t centralityArray[12] = {0., 2., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50.};
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
