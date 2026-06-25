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

#include "AliAnalysisPtVn.h"

#include "Riostream.h"

using std::cout;
using std::endl;

#define Fatal(d, s) printf("***  FATAL: %s: %s\n",d,s)

//ClassImp(AliAnalysisPtVn)

// ------------------------------------------------------------------------- //
AliAnalysisPtVn::AliAnalysisPtVn():
  fInputList(NULL),
  fHistList(NULL),
  fDebugLevel(0),
  fCentrality(0.),
  fptSubMax(0.8), fptSubMin(0.4),
  fCentralityBins(16)
{
// Dummy constructor of the class.
  printf("AliAnalysisPtVn::AliAnalysisPtVn()\n");
  this->InitializeArrays();
}

// ------------------------------------------------------------------------- //
AliAnalysisPtVn::AliAnalysisPtVn(const char *name):
  fInputList(NULL),
  fHistList(NULL),
  fDebugLevel(0),
  fCentrality(0.),
  fptSubMax(0.8), fptSubMin(0.4),
  fCentralityBins(9)
{
// Constructor of the class.
  printf("AliAnalysisPtVn::AliAnalysisPtVn(const char *name)\n");

  // Define the base list.
  fHistList = new TList();
  fHistList->SetName(name);
  fHistList->SetOwner(kTRUE);

  // Initialize all arrays within the data members.
  InitializeArrays();

  //DefineOutput(1, TList::Class());
}


// ------------------------------------------------------------------------- //
AliAnalysisPtVn::~AliAnalysisPtVn()
{
// Destructor of the class.
  if (fHistList) {delete fHistList;}
}

void AliAnalysisPtVn::InitializeArrays() {
  for (int i = 0; i < nCorr; ++i) {
    corrProfiles[i] = NULL;
  }


}

void AliAnalysisPtVn::BookAndNestAllLists() {
  for (int i = 0; i < nCorr; ++i) {
    corrProfiles[i] = new TProfile(Form("profile_%s", corrlistnames[i].Data()), corrlistnames[i].Data(),
      nHistCent, 0, float(nHistCent), "s");
    fHistList->Add(corrProfiles[i]);
  }
}

// ------------------------------------------------------------------------- //
void AliAnalysisPtVn::UserCreateOutputObjects() {
  this->BookAndNestAllLists();
}

// ------------------------------------------------------------------------- //
void AliAnalysisPtVn::WriteLists(TString name) {
// Write the corresponding outputlists.
  fHistList->Write(name.Data(),1);
  delete fHistList;
}

// ------------------------------------------------------------------------- //
void AliAnalysisPtVn::UserExec(Option_t *option) {
  Int_t centralityBin = SelectCentrality(fCentrality);

  Int_t nTracks = fInputList->GetEntriesFast();
  TComplex abcd[3][3][3][3] = {{{{TComplex(0., 0.)}}}};
  TComplex wabcd[3][3][3][3] = {{{{TComplex(0., 0.)}}}};

  TComplex abcd_v3[3][3][3][3] = {{{{TComplex(0., 0.)}}}};
  TComplex wabcd_v3[3][3][3][3] = {{{{TComplex(0., 0.)}}}};

  double psubA[3] = {0.,0.,0.};
  double psubB[3] = {0.,0.,0.};
  double psubMid[3] = {0.,0.,0.};
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
    }
    if(fptSubMin < iEta && iEta < fptSubMax){
      psubB[0] += 1;
      psubB[1] += iPt;
      psubB[2] += iPt*iPt;
    }
    if(std::abs(iEta)<0.4){
      psubMid[0] += 1;
      psubMid[1] += iPt;
      psubMid[2] += iPt*iPt;
    }
    if(iEta < -fptSubMin){
      Qneg[0] += 1;
      Qneg[1] += {TMath::Cos(-2.*iPhi),TMath::Sin(-2.*iPhi)};
    }
    if(iEta > fptSubMin){
      Qpos[0] += 1;
      Qpos[1] += {TMath::Cos(2.*iPhi),TMath::Sin(2.*iPhi)};
    }
    if(iEta < -fptSubMin){
      Q3neg[0] += 1;
      Q3neg[1] += {TMath::Cos(-3.*iPhi),TMath::Sin(-3.*iPhi)};
    }
    if(iEta > fptSubMin){
      Q3pos[0] += 1;
      Q3pos[1] += {TMath::Cos(3.*iPhi),TMath::Sin(3.*iPhi)};
    }

    FillABCDCounter(abcd, TComplex(TMath::Cos(2*iPhi), TMath::Sin(2*iPhi)), TComplex(TMath::Cos(-2*iPhi), TMath::Sin(-2*iPhi)), TComplex(iPt,0.), TComplex(1.,0.));
    FillABCDCounter(wabcd, TComplex(1, 0.), TComplex(1, 0.), TComplex(1, 0.), TComplex(1.,0.));

    FillABCDCounter(abcd_v3, TComplex(TMath::Cos(3.*iPhi), TMath::Sin(3.*iPhi)), TComplex(TMath::Cos(-3.*iPhi), TMath::Sin(-3.*iPhi)), TComplex(iPt,0.), TComplex(1.,0.));
    FillABCDCounter(wabcd_v3, TComplex(1, 0.), TComplex(1, 0.), TComplex(1, 0.), TComplex(1.,0.));

  } // Go to the next iTrack.
  Double_t corr[nCorr] = {0.};
  corr[0] = getStdAABBCC(abcd)/getStdAABBCC(wabcd); //! v_2^4, p_T^2
  corr[1] = getStdAABBC(abcd)/getStdAABBCC(wabcd); //! v_2^4, p_T 6pc
  corr[2] = getStdAABB(abcd)/getStdAABBCC(wabcd); //! v_2^4 6pc

  corr[3] = getStdABCC(abcd)/getStdABCC(wabcd); //! v_2^2, p_T^2
  corr[4] = getStdABC(abcd)/getStdABCC(wabcd); //! v_2^2, p_T 4pc
  corr[5] = getStdAB(abcd)/getStdABCC(wabcd); //! v_2^2 4pc

  corr[6] = getStdABC(abcd)/getStdABC(wabcd); //! v_2^2 p_T
  corr[7] = getStdAB(abcd)/getStdABC(wabcd); //! v_2^2 3pc
  corr[8] = getStdAB(abcd)/getStdAB(wabcd); //! v_2^2

  corr[9] = getStdAABB(abcd)/getStdAABB(wabcd); //! v_2^4
  corr[10] = abcd[0][0][1][0]/wabcd[0][0][1][0]; //! p_T
  corr[11] = getStdCC(abcd)/getStdCC(wabcd); //! p_T^2

  //extra
  corr[12] = getStdAABBCD(abcd)/getStdAABBCD(wabcd); //! v_2^4 p_T 6pc weighted
  corr[13] = getStdAABBDD(abcd)/getStdAABBDD(wabcd); //! v_2^4 6pc weighted
  corr[14] = getStdABCD(abcd)/getStdABCD(wabcd); //! v_2^2 p_T 4pc weighted
  corr[15] = getStdABDD(abcd)/getStdABDD(wabcd); //! v_2^2 4pc weighted
  corr[16] = getStdABD(abcd)/getStdABD(wabcd); //! v_2^2 3pc weighted

  //pt sub
  corr[17] = psubA[1]/psubA[0]; //! [pt]^A
  corr[18] = psubB[1]/psubB[0]; //! [pt]^B
  corr[19] = psubMid[1]/psubMid[0]; //! [pt]^Mid
  corr[20] = psubA[1]*psubB[1]/(psubA[0]*psubB[0]); //! [pt^2]^AB
  corr[21] = (psubA[1]*psubA[1] - psubA[2])/(psubA[0]*(psubA[0]-1)); //! [pt^2]^AA
  corr[22] = (psubB[1]*psubB[1] - psubB[2])/(psubB[0]*(psubB[0]-1)); //! [pt^2]^BB
  corr[23] = psubA[1]*psubMid[1]/(psubA[0]*psubMid[0]); //! [pt^2]^AMid
  corr[24] = psubB[1]*psubMid[1]/(psubB[0]*psubMid[0]); //! [pt^2]^MidB
  corr[25] = psubB[1]*(psubA[1]*psubA[1] - psubA[2])/(psubA[0]*(psubA[0]-1)*psubB[0]); //! [pt^3]^AAB
  corr[26] = psubA[1]*(psubB[1]*psubB[1] - psubB[2])/(psubB[0]*(psubB[0]-1)*psubA[0]); //! [pt^3]^ABB
  corr[27] = psubA[1]*psubB[1]*psubMid[1]/(psubA[0]*psubB[0]*psubMid[0]); //! [pt^3]^AMidB
  corr[28] = (psubA[1]*psubA[1] - psubA[2])*(psubB[1]*psubB[1] - psubB[2])/(psubA[0]*(psubA[0]-1)*psubB[0]*(psubB[0]-1)); //! [pt^4]^AABB

  //rho(v2^2,[pt])
  corr[29] = psubMid[1]/psubMid[0]; //! [pt]^mid
  corr[30] = (psubMid[1]*psubMid[1] - psubMid[2])/(psubMid[0]*(psubMid[0]-1)); //! [pt^2]^mid
  corr[31] = psubMid[1]*(Qpos[1]*Qneg[1]).Re()/((Qpos[0]*Qneg[0]).Re()*psubMid[0]); //! <v22 gap * [pt]^mid>
  corr[32] = (Qpos[1]*Qneg[1]).Re()/((Qpos[0]*Qneg[0]).Re()); //! <v22 gap>

  //rho(v3^2,[pt])
  corr[33] = psubMid[1]*(Q3pos[1]*Q3neg[1]).Re()/((Q3pos[0]*Q3neg[0]).Re()*psubMid[0]); //! <v32 gap * [pt]^mid>
  corr[34] = (Q3pos[1]*Q3neg[1]).Re()/((Q3pos[0]*Q3neg[0]).Re()); //! <v32 gap>
  corr[35] = getStdAB(abcd_v3)/getStdAB(wabcd_v3); //! v_3^2
  corr[36] = getStdAABB(abcd_v3)/getStdAABB(wabcd_v3); //! v_3^4

  for (int i = 0; i < nCorr; i++) {
    double multWeight = 1.;
    if (i==31) {
      multWeight = (Qpos[0]*Qneg[0]).Re()*psubMid[0];
    } else if (i==32) {
      multWeight = (Qpos[0]*Qneg[0]).Re()*psubMid[0];
    }
    if (TMath::Finite(corr[i])) corrProfiles[i]->Fill((Float_t)(centralityBin)+0.5, corr[i],multWeight); //Event weight for v22-pt gap for consistency with pub.
  }



  // Reset event-by-event objects.
  nTracks = 0;
}

// ------------------------------------------------------------------------- //
void AliAnalysisPtVn::Terminate(Option_t *option)
{
// Execute the additional steps after the loop over the events is over.
}


// ------------------------------------------------------------------------- //
Int_t AliAnalysisPtVn::SelectCentrality(Double_t centValue)
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
