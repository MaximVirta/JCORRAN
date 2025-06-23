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
#ifndef ALIANALYSISPtVn_H
#define ALIANALYSISPtVn_H

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

#include "AliJEfficiency.h"

class TClonesArray;
class AliAnalysisPtVn {
 public:
  // Methods inherited from AliAnalysisTaskSE.
  AliAnalysisPtVn();
  AliAnalysisPtVn(const char *name);
  virtual ~AliAnalysisPtVn();
  virtual void WriteLists(char* name);
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  // General methods specific to this analysis class.
  virtual void InitializeArrays();
  virtual void BookAndNestAllLists();


//// Correlation calculations ////
  double getStdAABBCC(TComplex abcdvec[3][3][3][3]){
    TComplex a = abcdvec[1][0][0][0];
    TComplex b = abcdvec[0][1][0][0];
    TComplex c = abcdvec[0][0][1][0];
    TComplex aa = abcdvec[2][0][0][0];
    TComplex bb = abcdvec[0][2][0][0];
    TComplex cc = abcdvec[0][0][2][0];
    TComplex ab = abcdvec[1][1][0][0];
    TComplex ac = abcdvec[1][0][1][0];
    TComplex bc = abcdvec[0][1][1][0];
    TComplex aab = abcdvec[2][1][0][0];
    TComplex aac = abcdvec[2][0][1][0];
    TComplex abb = abcdvec[1][2][0][0];
    TComplex acc = abcdvec[1][0][2][0];
    TComplex abc = abcdvec[1][1][1][0];
    TComplex bbc = abcdvec[0][2][1][0];
    TComplex bcc = abcdvec[0][1][2][0];
    TComplex aabb = abcdvec[2][2][0][0];
    TComplex aacc = abcdvec[2][0][2][0];
    TComplex aabc = abcdvec[2][1][1][0];
    TComplex abbc = abcdvec[1][2][1][0];
    TComplex abcc = abcdvec[1][1][2][0];
    TComplex bbcc = abcdvec[0][2][2][0];
    TComplex aabbc = abcdvec[2][2][1][0];
    TComplex aabcc = abcdvec[2][1][2][0];
    TComplex abbcc = abcdvec[1][2][2][0];
    TComplex aabbcc = abcdvec[2][2][2][0];
    return (a*a*b*b*c*c - aa*b*b*c*c - a*a*bb*c*c - a*a*b*b*cc - 4.*a*ab*b*c*c -
   4.*a*ac*b*b*c - 4.*a*a*b*bc*c + 4.*aab*b*c*c + 4.*aac*b*b*c +
   4.*a*abb*c*c + 4.*a*acc*b*b + 4.*a*a*bbc*c + 4.*a*a*b*bcc +
   16.*a*abc*b*c + aa*bb*c*c + aa*b*b*cc + a*a*bb*cc + 2.*ab*ab*c*c +
   2.*ac*ac*b*b + 2.*a*a*bc*bc + 4.*aa*b*bc*c + 4.*a*ac*bb*c +
   4.*a*ab*b*cc + 8.*ab*ac*b*c + 8.*a*ab*bc*c + 8.*a*ac*b*bc - 6.*aabb*c*c -
   24.*aabc*b*c - 6.*aacc*b*b - 24.*abbc*a*c - 24.*abcc*a*b - 6.*bbcc*a*a -
   8.*aab*bc*c - 8.*aac*b*bc - 4.*aac*bb*c - 4.*aab*b*cc - 8.*abb*ac*c -
   4.*abb*a*cc - 8.*acc*ab*b - 4.*acc*a*bb - 8.*bbc*a*ac - 4.*bbc*aa*c -
   8.*bcc*a*ab - 4.*bcc*aa*b - 16.*abc*ab*c - 16.*abc*ac*b - 16.*abc*a*bc -
   aa*bb*cc - 2.*ab*ab*cc - 2.*ac*ac*bb - 2.*bc*bc*aa - 8.*ab*ac*bc +
   48.*aabbc*c + 48.*aabcc*b + 48.*abbcc*a + 6.*aabb*cc + 6.*aacc*bb +
   6.*bbcc*aa + 24.*aabc*bc + 24.*abbc*ac + 24.*abcc*ab + 8.*aab*bcc +
   8.*aac*bbc + 8.*abb*acc + 16.*abc*abc - 120.*aabbcc).Re();
  }
  double getStdAABBC(TComplex abcdvec[3][3][3][3]){
    TComplex a = abcdvec[1][0][0][0];
    TComplex b = abcdvec[0][1][0][0];
    TComplex c = abcdvec[0][0][1][0];
    TComplex aa = abcdvec[2][0][0][0];
    TComplex ab = abcdvec[1][1][0][0];
    TComplex ac = abcdvec[1][0][1][0];
    TComplex bb = abcdvec[0][2][0][0];
    TComplex bc = abcdvec[0][1][1][0];
    TComplex aab = abcdvec[2][1][0][0];
    TComplex aac =  abcdvec[2][0][1][0];
    TComplex abb = abcdvec[1][2][0][0];
    TComplex abc = abcdvec[1][1][1][0];
    TComplex bbc = abcdvec[0][2][1][0];
    TComplex aabb = abcdvec[2][2][0][0];
    TComplex aabc = abcdvec[2][1][1][0];
    TComplex abbc = abcdvec[1][2][1][0];
    TComplex aabbc = abcdvec[2][2][1][0];
    return (a*a*b*b*c - aa*b*b*c - a*a*bb*c - 4.*ab*a*b*c - 2.*a*ac*b*b - 2.*a*a*bc*b
    + 2.*ab*ab*c + 4.*ab*ac*b + 4.*ab*bc*a + 8.*abc*a*b + 4.*aab*b*c + 2.*aac*b*b + 4.*abb*a*c + 2.*bbc*a*a + aa*bb*c + 2.*aa*b*bc + 2.*bb*a*ac
    - 12.*aabc*b - 12.*abbc*a - 6.*aabb*c - 8.*abc*ab - 2.*bbc*aa - 2.*aac*bb - 4.*aab*bc - 4.*abb*ac
    + 24.*aabbc).Re(); }
  double getStdAABB(TComplex abcdvec[3][3][3][3]){
  TComplex a = abcdvec[1][0][0][0];
  TComplex b = abcdvec[0][1][0][0];
  TComplex aa = abcdvec[2][0][0][0];
  TComplex bb = abcdvec[0][2][0][0];
  TComplex ab = abcdvec[1][1][0][0];
  TComplex aab = abcdvec[2][1][0][0];
  TComplex abb = abcdvec[1][2][0][0];
  TComplex aabb = abcdvec[2][2][0][0];
    return (a*a*b*b - aa*b*b - a*a*bb - 4*a*ab*b + aa*bb + 2*ab*ab  + 4*aab*b + 4*a*abb - 6*aabb).Re(); }
  double getStdABCC(TComplex abcdvec[3][3][3][3]){
    TComplex a = abcdvec[1][0][0][0];
    TComplex b = abcdvec[0][1][0][0];
    TComplex c = abcdvec[0][0][1][0];
    TComplex ab  = abcdvec[1][1][0][0];
    TComplex ac = abcdvec[1][0][1][0];
    TComplex bc = abcdvec[0][1][1][0];
    TComplex cc = abcdvec[0][0][2][0];
    TComplex abc = abcdvec[1][1][1][0];
    TComplex acc = abcdvec[1][0][2][0];
    TComplex bcc = abcdvec[0][1][2][0];
    TComplex abcc = abcdvec[1][1][2][0];
    return (a*b*c*c - a*b*cc - 2.*a*bc*c - 2.*ac*b*c - ab*c*c
    + 2.*acc*b + 2.*a*bcc + 4.*abc*c + ab*cc + 2.*ac*bc
    - 6.*abcc).Re(); }
  double getStdABC(TComplex abcdvec[3][3][3][3]){
  TComplex a = abcdvec[1][0][0][0];
  TComplex b = abcdvec[0][1][0][0];
  TComplex c = abcdvec[0][0][1][0];
  TComplex ab = abcdvec[1][1][0][0];
  TComplex ac = abcdvec[1][0][1][0];
  TComplex bc = abcdvec[0][1][1][0];
  TComplex abc = abcdvec[1][1][1][0];
    return (a*b*c - ab*c - ac*b - a*bc + 2.*abc).Re(); }
  double getStdAB(TComplex abcdvec[3][3][3][3]){
  TComplex a = abcdvec[1][0][0][0];
  TComplex b = abcdvec[0][1][0][0];
  TComplex ab = abcdvec[1][1][0][0];
    return (a*b-ab).Re(); }
  double getStdCC(TComplex abcdvec[3][3][3][3]){
  TComplex c = abcdvec[0][0][1][0];
  TComplex cc = abcdvec[0][0][2][0];
    return (c*c-cc).Re(); }
  template<typename T>
  void FillABCDCounter(T inarr[3][3][3][3], T a, T b, T c, T d)
  {
    for(int i = 0; i<3; ++i)
      for(int j = 0; j < 3; ++j)
        for(int k = 0; k < 3; ++k)
          for(int l = 0; l < 3; ++l)
            inarr[i][j][k][l] += TComplex::Power(a,i)*TComplex::Power(b,j)*TComplex::Power(c,k)*TComplex::Power(d,l);
    return;
  }
//// Correlation calculations ////




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
  void SetCentrality(Float_t cen0, Float_t cen1, Float_t cen2, Float_t cen3, Float_t cen4,
      Float_t cen5, Float_t cen6, Float_t cen7, Float_t cen8, Float_t cen9, Float_t cen10,
      Float_t cen11, Float_t cen12, Float_t cen13, Float_t cen14, Float_t cen15, Float_t cen16) {
    fcent_0 = cen0; fcent_1 = cen1; fcent_2 = cen2; fcent_3 = cen3; fcent_4 = cen4;
    fcent_5 = cen5; fcent_6 = cen6; fcent_7 = cen7; fcent_8 = cen8; fcent_9 = cen9;
    fcent_10 = cen10; fcent_11 = cen11; fcent_12 = cen12; fcent_13 = cen13; fcent_14 = cen14;
    fcent_15 = cen15; fcent_16 = cen16;
  }

  void SetEtaGaps(Bool_t ComputeEtaGap, Float_t EtaGap) {
    bComputeEtaGap = ComputeEtaGap;
    fEtaGap = EtaGap;
  }

 private:
  AliAnalysisPtVn(const AliAnalysisPtVn& aat);
  AliAnalysisPtVn& operator=(const AliAnalysisPtVn& aat);

  TClonesArray *fInputList;         // Input tracks selected in the JCatalyst.
  TList *fHistList;                 // Base list to hold all output objects.
  Int_t fDebugLevel;                // Verbosity of the class in the terminal.

  Double_t fCentrality;             // Centrality of the current event.
  Float_t fcent_0, fcent_1, fcent_2, fcent_3, fcent_4, fcent_5, fcent_6, fcent_7, fcent_8,
    fcent_9, fcent_10, fcent_11, fcent_12, fcent_13, fcent_14, fcent_15, fcent_16;
      // fcent_i holds the edge of a centrality bin.
  Int_t fCentralityBins;            //! Number of centrality bins (Size(array)-1).

  const static Int_t nCorr = 12;
  const static Int_t nHistCent = 7;

  TList *corrLists[nCorr];
  TProfile *corrProfiles[nCorr];
  Float_t minRange[nCorr] = {-.1, -.1, -.1, -.1, -.1, -.1, -.1, -.1, -.1, -0.1, 0.2, 0.2};
  Float_t maxRange[nCorr] = {.1, .1, .1, .1, .1, .1, .1, .1, .1, 0.1, 1., 1.};

  TString corrlistnames[nCorr] = {
    "v24pt2", "v24pt_6pc", "v24_6pc",
    "v22pt2", "v22pt_4pc", "v22_4pc",
    "v22pt", "v22_3pc", "v22",
    "v24", "pt", "pt2"}; 

  TH1D *corrHists[nCorr][nHistCent];

  Bool_t bComputeEtaGap;            // kTRUE: Calculate 2-particle eta gaps (default: kFALSE).
  Float_t fEtaGap;                  // Value of the eta gap itself.


  //ClassDef(AliAnalysisPtVn, 1); 
};

#endif  // ALIANALYSISPtVn
