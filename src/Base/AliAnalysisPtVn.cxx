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
  fCentrality(0.), fcent_0(0.), fcent_1(0.), fcent_2(0.), fcent_3(0.), fcent_4(0.),
  fcent_5(0.), fcent_6(0.), fcent_7(0.), fcent_8(0.), fcent_9(0.),
  fcent_10(0.), fcent_11(0.), fcent_12(0.), fcent_13(0.), fcent_14(0.),
  fcent_15(0.), fcent_16(0.),
  fCentralityBins(16), bComputeEtaGap(kFALSE), fEtaGap(0.8)
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
  fCentrality(0.), fcent_0(0.), fcent_1(0.), fcent_2(0.), fcent_3(0.), fcent_4(0.),
  fcent_5(0.), fcent_6(0.), fcent_7(0.), fcent_8(0.), fcent_9(0.),
  fcent_10(0.), fcent_11(0.), fcent_12(0.), fcent_13(0.), fcent_14(0.),
  fcent_15(0.), fcent_16(0.),
  fCentralityBins(9), bComputeEtaGap(kFALSE), fEtaGap(0.8)
{
// Constructor of the class.
//  printf("AliAnalysisPtVn::AliAnalysisPtVn(const char *name)\n");

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
    corrLists[i] = NULL;
    corrProfiles[i] = NULL;
    for (Int_t icent = 0; icent < nHistCent; icent++) {
      // corrHists[i][icent] = NULL;
    }
  }


}

void AliAnalysisPtVn::BookAndNestAllLists() {
  for (int i = 0; i < nCorr; ++i) {
    corrLists[i] = new TList();
    corrLists[i]->SetName(Form("list_%s",corrlistnames[i].Data()));
    corrLists[i]->SetOwner(kTRUE);
    // fHistList->Add(corrLists[i]);
    corrProfiles[i] = new TProfile(Form("profile_%s", corrlistnames[i].Data()), corrlistnames[i].Data(),
      nHistCent, 0, float(nHistCent), "s");
    fHistList->Add(corrProfiles[i]);
    for (int j = 0; j < nHistCent; ++j) {
      TString hname = Form("hist_%s_c%02d", corrlistnames[i].Data(), j);
      // corrHists[i][j] = new TH1D(hname.Data(), hname.Data(), 200, minRange[i], maxRange[i]);
      // corrLists[i]->Add(corrHists[i][j]);
    }
  }
}

// ------------------------------------------------------------------------- //
void AliAnalysisPtVn::UserCreateOutputObjects() {
  this->BookAndNestAllLists();
}

// ------------------------------------------------------------------------- //
void AliAnalysisPtVn::WriteLists(char* name) {
// Write the corresponding outputlists.
  fHistList->Write(name,1);
  delete fHistList;
  fHistList=NULL;
}

// ------------------------------------------------------------------------- //
void AliAnalysisPtVn::UserExec(Option_t *option) {
  Int_t centralityBin = SelectCentrality(fCentrality);

  Int_t nTracks = fInputList->GetEntriesFast();
  TComplex abcd[3][3][3][3] = {{{{TComplex(0., 0.)}}}};
  TComplex wabcd[3][3][3][3] = {{{{TComplex(0., 0.)}}}};

  if (fDebugLevel>10) printf("Centrality: %d, nTracks: %d\n", centralityBin, nTracks);

  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++)
  {
    Double_t iPhi = 0.;
    Double_t iPt = 0.;

    // Parse the selected track and get their information into the arrays.
    AliJBaseTrack *aTrack = (AliJBaseTrack*)fInputList->At(iTrack);
    if (!aTrack) {continue;}

    iPhi = aTrack->Phi();
    iPt = aTrack->Pt();

    FillABCDCounter(abcd, TComplex(TMath::Cos(2*iPhi), TMath::Sin(2*iPhi)), TComplex(TMath::Cos(-2*iPhi), TMath::Sin(-2*iPhi)), TComplex(iPt,0.), TComplex(1.,0.));
    FillABCDCounter(wabcd, TComplex(1, 0.), TComplex(1, 0.), TComplex(1, 0.), TComplex(1.,0.));

  } // Go to the next iTrack.
  Double_t corr[nCorr] = {0.};
  Int_t M = nTracks;
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

  for (int i = 0; i < nCorr; i++) {
    // corrHists[i][centralityBin]->Fill(corr[i]);
    corrProfiles[i]->Fill((Float_t)(centralityBin)+0.5, corr[i]);
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
  Float_t centralityArray[10] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80.};
  for (Int_t iCent = 0; iCent < 10; iCent++)
  {
    if (centralityArray[iCent] < 0) {return -1;}
    if (centValue >= centralityArray[iCent]) {continue;} 
    else {return iCent-1;} 
  } 

  // We went through all centrality edges without returning. This means that the measured value
  // is bigger than the maximum centrality we want for our analysis.
  return -1;
}
