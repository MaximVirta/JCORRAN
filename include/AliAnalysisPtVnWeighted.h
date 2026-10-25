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

  // The suffix is the event weight. C is the mid-rapidity subevent, |eta| < 0.4.
  //   w1A    = M_A
  //   w1B    = M_B
  //   w1C    = M_C
  //   w2AB   = M_A * M_B
  //   w2AA   = M_A * (M_A - 1)
  //   w2BB   = M_B * (M_B - 1)
  //   w2AC   = M_A * M_C
  //   w2CB   = M_C * M_B
  //   w2CC   = M_C * (M_C - 1)
  //   w3AAB  = M_A * (M_A - 1) * M_B
  //   w3ABB  = M_A * M_B * (M_B - 1)
  //   w3ABC  = M_A * M_B * M_C
  //   w4AABB = M_A * (M_A - 1) * M_B * (M_B - 1)
  enum ECorr {
    kPtA_w1A,
    kPtA_w2AB,
    kPtA_w2AA,
    kPtA_w3AAB,
    kPtA_w3ABB,
    kPtA_w4AABB,
    kPtB_w1B,
    kPtB_w2AB,
    kPtB_w2BB,
    kPtB_w3AAB,
    kPtB_w3ABB,
    kPtB_w4AABB,
    kPtC_w1C,
    kPtC_w3ABC,
    kPt2AB_w2AB,
    kPt2AB_w3AAB,
    kPt2AB_w3ABB,
    kPt2AB_w4AABB,
    kPt2AA_w2AA,
    kPt2AA_w3AAB,
    kPt2AA_w3ABB,
    kPt2AA_w4AABB,
    kPt2BB_w2BB,
    kPt2BB_w3AAB,
    kPt2BB_w3ABB,
    kPt2BB_w4AABB,
    kPt2AC_w2AC,
    kPt2CB_w2CB,
    kPt3AAB_w3AAB,
    kPt3ABB_w3ABB,
    kPt3AAB_w4AABB,
    kPt3ABB_w4AABB,
    kPt3ABC_w3ABC,
    kPt4AABB_w4AABB,
    kPt2CC_w2CC,
    kV22ptC_w3ABC,
    kV22_w2AB,
    kV22_w3ABC,
    kV32ptC_w3ABC,
    kV32_w2AB,
    kV32_w3ABC,
    kV32_w4AABB,
    kV24_w4AABB,
    kV34_w4AABB,
    // A–B subevent v_2^2, filled with w4AABB.
    // v22_w4AABB is one particle in A and one in B. The pT legs are distinct
    // from those flow particles.
    kV22_w4AABB,
    kV22ptA_w4AABB,
    kV22ptB_w4AABB,
    kV22ptAptB_w4AABB,
    kNCorr
  };
  const static Int_t nCorr = kNCorr;
  const static Int_t nHistCent = 12;

  TProfile *corrProfiles[nCorr];

  TString corrlistnames[nCorr] = {
    "ptA_w1A", "ptA_w2AB", "ptA_w2AA", "ptA_w3AAB", "ptA_w3ABB", "ptA_w4AABB",
    "ptB_w1B", "ptB_w2AB", "ptB_w2BB", "ptB_w3AAB", "ptB_w3ABB", "ptB_w4AABB",
    "ptC_w1C", "ptC_w3ABC",
    "pt2AB_w2AB", "pt2AB_w3AAB", "pt2AB_w3ABB", "pt2AB_w4AABB",
    "pt2AA_w2AA", "pt2AA_w3AAB", "pt2AA_w3ABB", "pt2AA_w4AABB",
    "pt2BB_w2BB", "pt2BB_w3AAB", "pt2BB_w3ABB", "pt2BB_w4AABB",
    "pt2AC_w2AC", "pt2CB_w2CB",
    "pt3AAB_w3AAB", "pt3ABB_w3ABB", "pt3AAB_w4AABB", "pt3ABB_w4AABB",
    "pt3ABC_w3ABC",
    "pt4AABB_w4AABB",
    "pt2CC_w2CC",
    "v22ptC_w3ABC", "v22_w2AB", "v22_w3ABC",
    "v32ptC_w3ABC", "v32_w2AB", "v32_w3ABC", "v32_w4AABB",
    "v24_w4AABB", "v34_w4AABB",
    "v22_w4AABB", "v22ptA_w4AABB", "v22ptB_w4AABB", "v22ptAptB_w4AABB"};

  //ClassDef(AliAnalysisPtVnWeighted, 1);
};

#endif  // ALIANALYSISPTVNWEIGHTED_H
