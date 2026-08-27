/* -------------------------------------------------------------------------- /
/ Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.      /
/ See cxx source for full Copyright notice                                    /
/ --------------------------------------------------------------------------- /
/ Analysis task for the computation of the numerator and denominators for the /
/ various SPC measured for Run2 Pb-Pb data                                    /
/                                                                             /
/ Authors: Cindy Mordasini (cindy.mordasini@cern.ch)                          /
/          Maxim Virta                                                        /
/ -------------------------------------------------------------------------- */
#ifndef ALIANALYSISV02_H
#define ALIANALYSISV02_H

//#include "AliHeader.h"
#include "TComplex.h"
#include "TList.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TSystem.h"
#include "AliJBaseTrack.h"
#include "TClonesArray.h"
#include <vector>
#include <complex>

class TClonesArray;
class AliAnalysisV02 {
 public:
  // Methods inherited from AliAnalysisTaskSE.
  AliAnalysisV02();
  AliAnalysisV02(const char *name);
  virtual ~AliAnalysisV02();
  virtual void WriteLists(TString name);
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  // General methods specific to this analysis class.
  virtual void InitializeArrays();
  virtual void BookAndNestAllLists();
  Int_t getPIDIndex(AliJBaseTrack *aTrack);
  /// General setters/getters.
  void SetInputList(TClonesArray *inputarray) {fInputList = inputarray;}
  void SetUseFlowWeights(bool useFlowWeights) {fUseMultiplicityFlowWeights = useFlowWeights;}
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
  AliAnalysisV02(const AliAnalysisV02& aat);
  AliAnalysisV02& operator=(const AliAnalysisV02& aat);

  TClonesArray *fInputList;         // Input tracks selected in the JCatalyst.
  TList *fHistList;                 // Base list to hold all output objects.
  Int_t fDebugLevel;                // Verbosity of the class in the terminal.

  Double_t fCentrality;             // Centrality of the current event.

  bool fUseMultiplicityFlowWeights;
  Float_t fptSubMin;
  Float_t fptSubMax;
  const static Int_t nPtBins = 26;
  const static Int_t nHistCent = 10;
  const static Int_t nSpecies = 4;
  Double_t fPtBins[nPtBins+1];

  TProfile2D *fV02Profile[nSpecies];
  TProfile2D *fPtNchProfile[nSpecies];
  TProfile *fV22Profile[nSpecies];



  //ClassDef(AliAnalysisV02, 1);
};

#endif  // ALIANALYSISV02
