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

#include "AliAnalysisAnaTwoMultiCorrelations.h"
#include "Riostream.h"
#include "TList.h"
#include "TComplex.h"
#include "TFile.h"
#include "TH2I.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TProfile.h"
#include "TMath.h"
#include "TDirectoryFile.h"
#include "TSystem.h"
#include "AliJBaseTrack.h"
#include <vector>
#include <TExMap.h>
#include <TClonesArray.h>
#include <sstream>

using std::cout;
using std::endl;


AliAnalysisAnaTwoMultiCorrelations::AliAnalysisAnaTwoMultiCorrelations():
  fInputList(0), fMainList(NULL), //fMPCList(NULL),
  fCentrality(-1), fCentralityBin(-1),
  fMultiplicity(0), fMultiplicityMin(10),
  fPtMin(0.2), fPtMax(5.0),
  fNCombi(7), fGetSC(kTRUE), fGetLowerHarmos(kTRUE)
{
// Dummy constructor of the class.
  //printf("AliAnalysisTaskHOCFA Dummy constructor\n");
  InitialiseArrayMembers();
}

/* ------------------------------------------------------------------------- */
AliAnalysisAnaTwoMultiCorrelations::AliAnalysisAnaTwoMultiCorrelations(const char*name, Bool_t usePWeights):
  fInputList(0), fMainList(NULL), //fMPCList(NULL),
  fCentrality(-1), fCentralityBin(-1),
  fMultiplicity(0), fMultiplicityMin(10),
  fPtMin(0.2), fPtMax(5.0),
  fNCombi(7), fGetSC(kTRUE), fGetLowerHarmos(kTRUE)
{
// Constructor of the class.
  //printf("AliAnalysisTaskHOCFA Constructor\n");
  InitialiseArrayMembers();

// Define the main list where everything is saved.
  fMainList = new TList();
  fMainList->SetName("MPC");
  fMainList->SetOwner(kTRUE);
}

/* ------------------------------------------------------------------------- */
AliAnalysisAnaTwoMultiCorrelations::~AliAnalysisAnaTwoMultiCorrelations()
{
// Destructor of the class.
  if (fMainList) {delete fMainList;}
}

/* ------------------------------------------------------------------------- */
void AliAnalysisAnaTwoMultiCorrelations::UserCreateOutputObjects()
{
// Declare the outputs of the task at the beginning of the analysis.
  this->BookFinalResults();
}

/* ------------------------------------------------------------------------- */
void AliAnalysisAnaTwoMultiCorrelations::UserExec(Option_t *)
{
// Execute the analysis for each event in the input sample.

// Get the centrality and multiplicity of the trimmed events and reject
// the ones with not enough tracks.
  fMultiplicity = fInputList->GetEntriesFast();
  fCentralityBin = GetCentralityBin(fCentrality);
  if (fMultiplicity < fMultiplicityMin) {return;}

// Get the information on the selected tracks.
  Double_t *iPt  = new Double_t[fMultiplicity]();       // Transverse momentum.
  double *iEta = new double[fMultiplicity]();   // Pseudorapidity.
  Double_t *iPhi = new Double_t[fMultiplicity]();       // Azimuthal angles.
  Int_t iIndex = 0;           // Index of the selected track in the final arrays.

  for (Int_t iTrack = 0; iTrack < fMultiplicity; iTrack++) {
    AliJBaseTrack *aTrack = (AliJBaseTrack*)fInputList->At(iTrack);
    if (!aTrack) {continue;}

    iPt[iIndex] = aTrack->Pt();
    iEta[iTrack] = aTrack->Eta();
    iPhi[iTrack] = aTrack->Phi();

    iIndex++;
  }

// Compute the Q-vectors and multiparticle correlations.
  CalculateQvectors(fMultiplicity, iPhi);
  ComputeAllTerms();

// Reset the variables for the next event.
  fMultiplicity = 0;
  delete [] iPt;
  delete [] iEta;
  delete [] iPhi;
}

/* ------------------------------------------------------------------------- */
void AliAnalysisAnaTwoMultiCorrelations::Terminate(Option_t *)
{
// All the additional steps after the loop over the events.
  //printf("AliAnalysisTaskHOCFA is done! \n");
} // End: void Terminate(Option_t *).

/* ------------------------------------------------------------------------- */
Int_t AliAnalysisAnaTwoMultiCorrelations::GetCentralityBin(Float_t cent) {
// The analysis always uses the same division of centralities.
  Float_t centralityArray[10] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80.};
  for (Int_t i = 0; i < 10; i++) {
    if (cent >= centralityArray[i]) {continue;}
    else {return i-1;}
  }
  return -1;
}

/* ------------------------------------------------------------------------- */
void AliAnalysisAnaTwoMultiCorrelations::InitialiseArrayMembers()
{
// Initialise the arrays in the data members.
  for (Int_t i = 0; i < 9; i++) {
    fCentralityList[i] = NULL;

    fCorrel2p[i] = NULL;
    for (int j = 0; j < 7; j++) {fCorrel2h[j][i] = NULL;}
    fCorrel3h[i] = NULL;
  }

  int harmo2h[13][2] = {{3,2}, {4,2}, {4,3}, {5,2}, {5,3}, {5,4}, {6,2},  // Lower: 7 first.
                        {6,3}, {6,4}, {7,2}, {7,3}, {8,2}, {8,3}};
  int powers[13][2] = { {1,1}, {2,0}, {2,1}, {3,0}, {3,1}, {4,0}, {4,1},  // For the 2-h AC.
                               {0,2}, {1,2}, {0,3}, {1,3}, {0,4}, {1,4}};
  for (int i = 0; i < 13; i++) {
    for (int j = 0; j < 2; j++) {
      fHarmoArray2h[i][j] = harmo2h[i][j];
      fPowers[i][j] = powers[i][j];
    }
  }

  int harmo3h[9][3] = { {2,3,4}, {2,3,5}, {2,4,5}, {2,3,6}, {2,3,7},  // Lower: 3 first.
                        {2,3,8}, {2,4,6}, {3,4,5}, {3,4,6}};
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 3; j++) {fHarmoArray3h[i][j] = harmo3h[i][j];}
  }

  for (Int_t iHarmo = 0; iHarmo < 81; iHarmo++) {
    for (Int_t iPower = 0; iPower < 11; iPower++) {
      fQvectors[iHarmo][iPower] = TComplex(0., 0.);
    }
  }
}

/* ------------------------------------------------------------------------- */
void AliAnalysisAnaTwoMultiCorrelations::BookFinalResults()
{
  if (!fMainList) {Fatal("AliAnalysisTaskHOCFA::BookFinalResults()", "FATAL: fMainList not found.");}

  Float_t centralityArray[10] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80.};

  for (Int_t iCent = 0; iCent < 9; iCent++) {
  // Book the needed centrality lists.
    fCentralityList[iCent] = new TList();
    fCentralityList[iCent]->SetName(Form("Centrality_%.1f-%.1f", centralityArray[iCent], centralityArray[iCent+1]));
    fCentralityList[iCent]->SetOwner(kTRUE);
    fMainList->Add(fCentralityList[iCent]);

    // Book the profiles for the resulting correlators.
    fCorrel2p[iCent] = new TProfile("", "", 8, 0., 8.);
    fCorrel2p[iCent]->SetName(Form("fCorrel2p_Bin%d", iCent));
    fCorrel2p[iCent]->SetTitle(Form("2-p terms, no #Delta #eta, bin%d", iCent));
    fCorrel2p[iCent]->SetStats(kTRUE);
    fCentralityList[iCent]->Add(fCorrel2p[iCent]);

    for (int j = 0; j < fNCombi; j++) {   // fNCombi = 1 (SC), 7 (lower AC) or 6 (higher AC).
      fCorrel2h[j][iCent] = new TProfile("", "", 13, 0., 13.);
      if (fGetSC) {
        fCorrel2h[j][iCent]->SetName(Form("fCorrel2h_allSC_Bin%d", iCent));
        fCorrel2h[j][iCent]->SetTitle(Form("2-h terms for SC, bin%d", iCent));
      } else {
        if (fGetLowerHarmos) {  // We measure only the 7 first 2-h terms.
          fCorrel2h[j][iCent]->SetName(Form("fCorrel2h_AC%d%d_Bin%d",
            fHarmoArray2h[j][0], fHarmoArray2h[j][1], iCent));
           fCorrel2h[j][iCent]->SetTitle(Form("2-h terms for AC(%d,%d), bin%d",
            fHarmoArray2h[j][0], fHarmoArray2h[j][1], iCent));
        } else {  // Else we measure only the 6 last 2-h terms.
          fCorrel2h[j][iCent]->SetName(Form("fCorrel2h_AC%d%d_Bin%d",
            fHarmoArray2h[j+7][0], fHarmoArray2h[j+7][1], iCent));
          fCorrel2h[j][iCent]->SetTitle(Form("2-h terms for AC(%d,%d), bin%d",
            fHarmoArray2h[j+7][0], fHarmoArray2h[j+7][1], iCent));
        }
      }   // The names and titles for the 2-h profile have been set.
      fCorrel2h[j][iCent]->Sumw2();
      fCentralityList[iCent]->Add(fCorrel2h[j][iCent]);
    } // Go to the next combination in the list.

    fCorrel3h[iCent] = new TProfile("", "", 9, 0., 9.);
    fCorrel3h[iCent]->SetName(Form("fCorrel3h_Bin%d", iCent));
    fCorrel3h[iCent]->SetTitle(Form("3-h terms for SC, bin%d", iCent));
    fCorrel3h[iCent]->SetStats(kTRUE);
    if (fGetSC) {fCentralityList[iCent]->Add(fCorrel3h[iCent]);}
  } // Go to the next centrality bin.
}

// ------------------------------------------------------------------------- //
void AliAnalysisAnaTwoMultiCorrelations::CalculateQvectors(Long64_t multiplicity, Double_t angles[])
{
// Calculate the needed Q-vectors.
  Double_t iAngle = 0.; // Azimuthal angle of the current track.

// Ensure all the Q-vectors are initially zero.
  for (Int_t iHarmo = 0; iHarmo < 81; iHarmo++){
    for (Int_t iPower = 0; iPower < 11; iPower++){
      fQvectors[iHarmo][iPower] = TComplex(0., 0.);
    }
  }

// Compute the Q-vectors.
  for (Long64_t iTrack = 0; iTrack < multiplicity; iTrack++){
    iAngle = angles[iTrack];
    for (Int_t iHarmo = 0; iHarmo < 81; iHarmo++){
      for (Int_t iPower = 0; iPower < 11; iPower++){
        fQvectors[iHarmo][iPower] += TComplex(TMath::Cos(iHarmo*iAngle), TMath::Sin(iHarmo*iAngle));
      }
    }
  }
// Reset the variables.
  iAngle = 0.;
}

// ------------------------------------------------------------------------- //
TComplex AliAnalysisAnaTwoMultiCorrelations::Q(Int_t n, Int_t p)
{
// Alias for fQvectors to make them more easy to use.
  if (n >= 0) {return fQvectors[n][p];}
  return TComplex::Conjugate(fQvectors[-n][p]); // Use that Q*(n,p) = Q(-n,p).
}

// ------------------------------------------------------------------------- //
TComplex AliAnalysisAnaTwoMultiCorrelations::CalculateRecursion(Int_t n, Int_t *harmonic, Int_t mult, Int_t skip)
{
/* Calculate the multi-particle correlators by using the recursion method.                  /
/ Improved, faster version originally developed by Kristjan Gulbrandsen (gulbrand@nbi.dk). */
  Int_t nm1 = n-1;
  TComplex c(Q(harmonic[nm1], mult));
  if (nm1 == 0) return c;
  c *= CalculateRecursion(nm1, harmonic);
  if (nm1 == skip) return c;

  Int_t multp1 = mult+1;
  Int_t nm2 = n-2;
  Int_t counter1 = 0;
  Int_t hhold = harmonic[counter1];
  harmonic[counter1] = harmonic[nm2];
  harmonic[nm2] = hhold + harmonic[nm1];
  TComplex c2(CalculateRecursion(nm1, harmonic, multp1, nm2));
  Int_t counter2 = n-3;

  while (counter2 >= skip) {
    harmonic[nm2] = harmonic[counter1];
    harmonic[counter1] = hhold;
    ++counter1;
    hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[nm2];
    harmonic[nm2] = hhold + harmonic[nm1];
    c2 += CalculateRecursion(nm1, harmonic, multp1, counter2);
    --counter2;
  }
  harmonic[nm2] = harmonic[counter1];
  harmonic[counter1] = hhold;

  if (mult == 1) return c-c2;
  return c-Double_t(mult)*c2;
}

// ------------------------------------------------------------------------- //
void AliAnalysisAnaTwoMultiCorrelations::ComputeAllTerms()
{
// Compute all the terms needed for the ACs/SCs for all the observables.
  int nPart = 2;        // We start with the case of 2-particle correlations.
  int lHarmo[7] = {0};  // We work with "symmetric" correlators of max 14 particles > max 7 harmonics.
  int lPower[2] = {0};  // Powers of {vm2,vn2} for the AC case.

  // Fill the profile for the 2-p as it is common to all analysis cases.
  for (int iBin = 0; iBin < 8; iBin++){
    lHarmo[0] = iBin+1; // Only the first element needs to be filled in the 2-p case.
    CalculateCorrelator(nPart, lHarmo, fCorrel2p[fCentralityBin], iBin, lPower);
  }

  // Fill the 2-harmonic profile according to the analysis configuration.
  // fNCombi = 1 (for SC, 13 combis in 1 profile), 7 (for lower AC) or 6 (for higher AC).
  for (int iProf = 0; iProf < fNCombi; iProf++) {
    for (int iBin = 0; iBin < 13; iBin++) {
      // Define the number of particles and the harmonics for AC/SC individually.
      if (fGetSC) {
        nPart = 4;  // 2-h terms in SC are always 4-particle correlators.
        lHarmo[0] = fHarmoArray2h[iBin][0];
        lHarmo[1] = fHarmoArray2h[iBin][1];
      } else {
        nPart = 0;  // Reset the number of particles for this bin.
        int cHarmo = 0; // Current index in the harmonic array.
        for (int iHarmo = 0; iHarmo < 2; iHarmo++) {
          if (fPowers[iBin][iHarmo] == 0) {continue;}   // Skip the unneeded harmonics.
          lPower[iHarmo] = fPowers[iBin][iHarmo];   // We change the power element only if non-zero.

          for (int jPower = 0; jPower < fPowers[iBin][iHarmo]; jPower++) {
            lHarmo[cHarmo] = fHarmoArray2h[iProf][iHarmo];   // Write the harmo as many times as its power is.
            cHarmo++;
          }

          nPart += 2*fPowers[iBin][iHarmo]; // 2-h terms in AC have twice as many particles as their cumulant order.
        } // Go to the next harmonic of the pair.
      } // We have the harmonic array and the number of particles at this point.

      // Calculate the multiparticle correlator itself using the recursion method.
      CalculateCorrelator(nPart, lHarmo, fCorrel2h[iProf][fCentralityBin], iBin, lPower);

      // Reset the variables for safety purposes.
      lPower[0] = 0.; lPower[1] = 0.;
    } // Go to the next bin in the current profile.
  } // Go to the next profile in the array.

  // Fill the 3-harmonic profile only if the analysis is configured for SC.
  if (fGetSC) {
    nPart = 6;  // 3-h terms in SC are always 6-particle correlators.
    for (int iBin = 0; iBin < 9; iBin++) {
      // Define the harmonic array for this bin.
      for (int iH = 0; iH < 3; iH++) {lHarmo[iH] = fHarmoArray3h[iBin][iH];}

      // Calculate the multiparticle correlator for this combination of harmonics.
      CalculateCorrelator(nPart, lHarmo, fCorrel3h[fCentralityBin], iBin, lPower);
    } // Go to the next combination of 3 harmonics.
  }

}

// ------------------------------------------------------------------------- //
void AliAnalysisAnaTwoMultiCorrelations::CalculateCorrelator(int myMulti, int myHarmos[], TProfile *myProfile, int myBin, int myPowers[])
{
// Calculate the multiparticle correlator corresponding to harmonics[].
  TComplex cCorrel = TComplex(0., 0.);
  double eventWeight = 0.;  // Event weight for this correlator.
  double rCorrel = 0.;      // Final value to fill in the profile.

  int numer2h[2] = {myHarmos[0], -1*myHarmos[0]};
  int denom2h[2] = {0};
  int numer4h[4] = {myHarmos[0], myHarmos[1],
                    -1*myHarmos[0], -1*myHarmos[1]};
  int denom4h[4] = {0};
  int numer6h[6] = {myHarmos[0], myHarmos[1], myHarmos[2],
                    -1*myHarmos[0], -1*myHarmos[1], -1*myHarmos[2]};
  int denom6h[6] = {0};
  int numer8h[8] = {myHarmos[0], myHarmos[1], myHarmos[2], myHarmos[3],
                    -1*myHarmos[0], -1*myHarmos[1], -1*myHarmos[2], -1*myHarmos[3]};
  int denom8h[8] = {0};
  int numer10h[10] = {myHarmos[0], myHarmos[1], myHarmos[2], myHarmos[3], myHarmos[4],
                    -1*myHarmos[0], -1*myHarmos[1], -1*myHarmos[2], -1*myHarmos[3], -1*myHarmos[4]};
  int denom10h[10] = {0};

  // Compute the denominator (= event weight) and numerator with the provided harmonics.
  switch(myMulti) {
  case 2 :
    eventWeight = ( CalculateRecursion(2, denom2h) ).Re();
    cCorrel = ( CalculateRecursion(2, numer2h) )/eventWeight;
    rCorrel = cCorrel.Re();
    break;
  case 4 :
    eventWeight = ( CalculateRecursion(4, denom4h) ).Re();
    cCorrel = ( CalculateRecursion(4, numer4h) )/eventWeight;
    rCorrel = cCorrel.Re();
    break;
  case 6 :    
    eventWeight = ( CalculateRecursion(6, denom6h) ).Re();
    cCorrel = ( CalculateRecursion(6, numer6h) )/eventWeight;
    rCorrel = cCorrel.Re();
    break;
  case 8 :
    eventWeight = ( CalculateRecursion(8, denom8h) ).Re();
    cCorrel = ( CalculateRecursion(8, numer8h) )/eventWeight;
    rCorrel = cCorrel.Re();
    break;
  case 10 :
    eventWeight = ( CalculateRecursion(10, denom10h) ).Re();
    cCorrel = ( CalculateRecursion(10, numer10h) )/eventWeight;
    rCorrel = cCorrel.Re();
    break;
  default :
    printf("Error: invalid number of particles.\n");
    break;
  }

  // Fill the corresponding bin in the right TProfile.
  myProfile->Fill( (float)myBin + 0.5, rCorrel, eventWeight );

  if (myMulti == 2) {   // Bins are filled with v1-v8.
    myProfile->GetXaxis()->SetBinLabel(myBin+1, Form("v_{%d}", myBin+1));
  } else if (fGetSC) {
    if (myMulti == 4) {  // Bins are filled with 2-h combinations.
      myProfile->GetXaxis()->SetBinLabel(myBin+1, Form("(%d,%d)", myHarmos[0], myHarmos[1]));
    } else {
      myProfile->GetXaxis()->SetBinLabel(myBin+1, Form("(%d,%d,%d)", myHarmos[0], myHarmos[1], myHarmos[2]));
    }
  } else {  // Bins are filled with the power patterns.
    myProfile->GetXaxis()->SetBinLabel(myBin+1, Form("{%d,%d}", myPowers[0], myPowers[1]));
  }

/*// Reset the local variables for the next call.
  cCorrel = TComplex(0., 0.);
  rCorrel = 0.;
  eventWeight = 0.;*/
}
