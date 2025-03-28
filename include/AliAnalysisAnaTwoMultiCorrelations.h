/* -------------------------------------------------------------------------- /
/ Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.      /
/ See cxx source for full Copyright notice                                    /
/ --------------------------------------------------------------------------- /
/ Analysis task for the computation of the correlators needed for the CFA.    /
/                                                                             /
/ Author: Cindy Mordasini (cindy.mordasini@cern.ch)                           /
/ -------------------------------------------------------------------------- */

#ifndef AliAnalysisAnaTwoMultiCorrelations_H
#define AliAnalysisAnaTwoMultiCorrelations_H

#include "TSystem.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TList.h"
#include "TComplex.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TProfile.h"
#include <TExMap.h>
#include <sstream>
class TClonesArray;

class AliAnalysisAnaTwoMultiCorrelations{
public:
// Methods inherited from AliAnalysisTaskSE.
  AliAnalysisAnaTwoMultiCorrelations();
  AliAnalysisAnaTwoMultiCorrelations(const char *name, Bool_t usePWeights=kFALSE);
  virtual ~AliAnalysisAnaTwoMultiCorrelations();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);

// Methods specific for this class.
  void SetInputList(TClonesArray *inputarray) {fInputList = inputarray;}
  TClonesArray* GetInputList() const {return fInputList;}
  void SetEventCentrality(float cent) {fCentrality = cent;}
  Int_t GetCentralityBin(Float_t cent);   // Centrality array is hard-coded with 9 bins.

  virtual void InitialiseArrayMembers();
  virtual void BookFinalResults();
  TList* GetTList() const {return fMainList;} // Get the list for external task.

  void SetMinMultiplicity(int minMult) {fMultiplicityMin = minMult;}
  void SetPtRange(double minPt, double maxPt) {fPtMin = minPt; fPtMax = maxPt;}
  void SetObservable(bool thisObs, bool thisOrder) {fGetSC = thisObs; fGetLowerHarmos = thisOrder;
    if (thisObs) {fNCombi = 1;}       // All 13 combinations are obtained in one go for the SC.
    else {
      if (thisOrder) {fNCombi = 7;}   // We get the AC only for the 7 lowest 2-h combinations.
      else {fNCombi = 6;}             // We get the AC only for the 6 higher 2-h combinations.
    }
  }

  virtual void CalculateQvectors(Long64_t multiplicity, Double_t angles[]);
  TComplex Q(Int_t n, Int_t p);
  TComplex CalculateRecursion(Int_t n, Int_t *harmonic, Int_t mult=1, Int_t skip=0);
  virtual void ComputeAllTerms(); // TBC: Do I need to pass the angles and weights?
  virtual void CalculateCorrelator(int myMulti, int myHarmos[], TProfile *myProfile, int myBin, int myPowers[]);//(Int_t combi, Int_t bin, Int_t nParticles, Int_t harmonics[], Double_t *errorTerms);
  
/* ------------------------------------------------------------------------- */
private:
  AliAnalysisAnaTwoMultiCorrelations(const AliAnalysisAnaTwoMultiCorrelations& aattmc);
  AliAnalysisAnaTwoMultiCorrelations& operator=(const AliAnalysisAnaTwoMultiCorrelations& aattmc);
  TClonesArray *fInputList;
  TList *fMainList;         // Main list where all results are saved.
  TList *fCentralityList[9];   //! Results per each centrality bin.
  
  Float_t fCentrality;      // Centrality of the current event.
  Int_t fCentralityBin;     //! Corresponding centrality bin.
  Long64_t fMultiplicity;   //! Multiplicity of the event after full selection.
  Int_t fMultiplicityMin;   // Minimum multiplicity to calculate the correlators.
  double fPtMin;                // Minimum transverse momentum.
  double fPtMax;                // Maximum transverse momentum.


  int fHarmoArray2h[13][2];     // Combinations of 2-harmonics for AC/SC.
    // 13: max number of combinations of 2 harmonics.
  int fHarmoArray3h[9][3];      // Combinations of 3-harmonics for SC only.
  int fPowers[13][2];           // List of 2-h terms by their power for AC only.
  TComplex fQvectors[81][11];   // All the needed combinations of Q-vectors.
    // Size: [(v8*10part)+1][10part+1].
  int fNCombi;                  // Number of combinations of harmonics in one wagon (max 7).
  bool fGetSC;                  // kTRUE: Measure 2-h and 3-h SC, else 2-h AC (default: true).
  bool fGetLowerHarmos;         // kTRUE: Measure the terms for the lower harmonics (default: true).

  TProfile *fCorrel2p[9];        //! 2-p terms for v1-v8 without eta gap (8 bins).
  TProfile *fCorrel2h[7][9];     //! 2-harmonic terms for SC/AC (13 bins).
    // 7: max number of combi in one wagon (all in 1 for SC), 16: max number of centrality bins.
  TProfile *fCorrel3h[9];        //! 3-harmonic terms for SC (9 bins, not enabled for AC).
};  // End of the class.

#endif
