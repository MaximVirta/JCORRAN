/* -------------------------------------------------------------------------- /
/ Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.      /
/ See cxx source for full Copyright notice                                    /
/ --------------------------------------------------------------------------- /
/ Weighted copy of AliAnalysisPtVn. Profile fills use the combinatorial       /
/ event weight of each correlator. Mean-pT moments are also filled with the   /
/ higher-order weights needed to build cumulants. v_2^4 p_T^n is omitted.     /
/                                                                             /
/ Authors: Cindy Mordasini (cindy.mordasini@cern.ch)                          /
/          Maxim Virta                                                        /
/ -------------------------------------------------------------------------- */
#ifndef ALIANALYSISPTVNWEIGHTED_H
#define ALIANALYSISPTVNWEIGHTED_H

//#include "AliHeader.h"
#include "TComplex.h"
#include "TList.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TProfile.h"
#include "TSystem.h"
#include "AliJBaseTrack.h"
#include "TClonesArray.h"
#include <vector>
#include <complex>

class TClonesArray;
class AliAnalysisPtVnWeighted {
 public:
  // Methods inherited from AliAnalysisTaskSE.
  AliAnalysisPtVnWeighted();
  AliAnalysisPtVnWeighted(const char *name);
  virtual ~AliAnalysisPtVnWeighted();
  virtual void WriteLists(TString name);
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  // General methods specific to this analysis class.
  virtual void InitializeArrays();
  virtual void BookAndNestAllLists();
  virtual void SetUseMultWeight(Bool_t useWeight) {fUseMultWeight = useWeight;}
  void FillCorr(Int_t centralityBin, Int_t index, Double_t value, Double_t weight);

  /// General setters/getters.
  void SetInputList(TClonesArray *inputarray) {fInputList = inputarray;}
  TClonesArray *GetInputList() const {return fInputList;}
  TList* GetMainList() const{return fHistList;}
  void SetDebugLevel(Int_t debuglevel) {
    fDebugLevel = debuglevel;
    cout << "Setting Debug Level = " << fDebugLevel << endl;
  }

  /// Centrality-related methods.
  void SetEventCentrality(Double_t cent) {fCentrality = cent;}
  Int_t SelectCentrality(Double_t centValue);
  void SetPtSubRange(Float_t ptSubMin, Float_t ptSubMax) {
    fptSubMin = ptSubMin;
    fptSubMax = ptSubMax;
  }

 private:
  AliAnalysisPtVnWeighted(const AliAnalysisPtVnWeighted& aat);
  AliAnalysisPtVnWeighted& operator=(const AliAnalysisPtVnWeighted& aat);

  TClonesArray *fInputList;         // Input tracks selected in the JCatalyst.
  TList *fHistList;                 // Base list to hold all output objects.
  Int_t fDebugLevel;                // Verbosity of the class in the terminal.
  Bool_t fUseMultWeight;

  Double_t fCentrality;             // Centrality of the current event.

  Int_t fCentralityBins;            //! Number of centrality bins (Size(array)-1).
  Float_t fptSubMin;
  Float_t fptSubMax;

  // Natural-weight profiles keep the original names. A suffix is the event
  // weight used for that copy of the same per-event value:
  //   w2     = M_A * M_B
  //   w3AAB  = M_A * (M_A - 1) * M_B
  //   w3ABB  = M_A * M_B * (M_B - 1)
  //   w4     = M_A * (M_A - 1) * M_B * (M_B - 1)
  // ptA/ptB themselves are the 1-particle fills (weights M_A and M_B).
  // pt2* without a suffix are the 2-particle fills, pt3* the 3-particle fills,
  // and pt4AABB the 4-particle fill.
  enum ECorr {
    kPtA,
    kPtA_w2,
    kPtA_w3AAB,
    kPtA_w3ABB,
    kPtA_w4,
    kPtB,
    kPtB_w2,
    kPtB_w3AAB,
    kPtB_w3ABB,
    kPtB_w4,
    kPtMid,
    kPtMid_w3,
    kPt2AB,
    kPt2AB_w3AAB,
    kPt2AB_w3ABB,
    kPt2AB_w4,
    kPt2AA,
    kPt2AA_w3AAB,
    kPt2AA_w3ABB,
    kPt2AA_w4,
    kPt2BB,
    kPt2BB_w3AAB,
    kPt2BB_w3ABB,
    kPt2BB_w4,
    kPt2AMid,
    kPt2MidB,
    kPt3AAB,
    kPt3ABB,
    kPt3AAB_w4,
    kPt3ABB_w4,
    kPt3AMidB,
    kPt4AABB,
    kPt2mid,
    kV22pt_gap,
    kV22_gap,
    kV22_w3,
    kV32pt_gap,
    kV32_gap,
    kV32_w3,
    kV24_gap,
    kV34_gap,
    // A–B subevent v_2^2, with the four-particle weight w4.
    // v22_w4 is one particle in A and one in B. The pT legs are distinct
    // from those flow particles.
    kV22_w4,
    kV22ptA_w4,
    kV22ptB_w4,
    kV22ptAptB_w4,
    kNCorr
  };
  const static Int_t nCorr = kNCorr;
  const static Int_t nHistCent = 12;

  TProfile *corrProfiles[nCorr];

  TString corrlistnames[nCorr] = {
    "ptA", "ptA_w2", "ptA_w3AAB", "ptA_w3ABB", "ptA_w4",
    "ptB", "ptB_w2", "ptB_w3AAB", "ptB_w3ABB", "ptB_w4",
    "ptMid", "ptMid_w3",
    "pt2AB", "pt2AB_w3AAB", "pt2AB_w3ABB", "pt2AB_w4",
    "pt2AA", "pt2AA_w3AAB", "pt2AA_w3ABB", "pt2AA_w4",
    "pt2BB", "pt2BB_w3AAB", "pt2BB_w3ABB", "pt2BB_w4",
    "pt2AMid", "pt2MidB",
    "pt3AAB", "pt3ABB", "pt3AAB_w4", "pt3ABB_w4",
    "pt3AMidB",
    "pt4AABB",
    "pt2mid",
    "v22pt_gap", "v22_gap", "v22_w3",
    "v32pt_gap", "v32_gap", "v32_w3",
    "v24_gap", "v34_gap",
    "v22_w4", "v22ptA_w4", "v22ptB_w4", "v22ptAptB_w4"};

  //ClassDef(AliAnalysisPtVnWeighted, 1);
};

#endif  // ALIANALYSISPTVNWEIGHTED_H
