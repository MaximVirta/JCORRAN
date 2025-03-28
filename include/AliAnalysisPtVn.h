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
class AliAnalysisPtVn {
 public:
  // Methods inherited from AliAnalysisTaskSE.
  AliAnalysisPtVn();
  AliAnalysisPtVn(const char *name);
  virtual ~AliAnalysisPtVn();
  virtual void WriteLists(TString name);
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
  double getStdAABBCD(TComplex abcdvec[3][3][3][3]){
  TComplex a = abcdvec[1][0][0][0];
  TComplex b = abcdvec[0][1][0][0];
  TComplex c = abcdvec[0][0][1][0];
  TComplex d = abcdvec[0][0][0][1];
  TComplex aa = abcdvec[2][0][0][0];
  TComplex bb = abcdvec[0][2][0][0];
  TComplex ab = abcdvec[1][1][0][0];
  TComplex ac = abcdvec[1][0][1][0];
  TComplex ad = abcdvec[1][0][0][1];
  TComplex bc = abcdvec[0][1][1][0];
  TComplex bd = abcdvec[0][1][0][1];
  TComplex cd = abcdvec[0][0][1][1];
  TComplex aab = abcdvec[2][1][0][0];
  TComplex aac = abcdvec[2][0][1][0];
  TComplex aad = abcdvec[2][0][0][1];
  TComplex abb = abcdvec[1][2][0][0];
  TComplex abc = abcdvec[1][1][1][0];
  TComplex abd = abcdvec[1][1][0][1];
  TComplex acd = abcdvec[1][0][1][1];
  TComplex bbc = abcdvec[0][2][1][0];
  TComplex bbd = abcdvec[0][2][0][1];
  TComplex bcd = abcdvec[0][1][1][1];
  TComplex aabb = abcdvec[2][2][0][0];
  TComplex aabc = abcdvec[2][1][1][0];
  TComplex aabd = abcdvec[2][1][0][1];
  TComplex aacd = abcdvec[2][0][1][1];
  TComplex abbc = abcdvec[1][2][1][0];
  TComplex abbd = abcdvec[1][2][0][1];
  TComplex abcd = abcdvec[1][1][1][1];
  TComplex bbcd = abcdvec[0][2][1][1];
  TComplex aabbc = abcdvec[2][2][1][0];
  TComplex aabbd = abcdvec[2][2][0][1];
  TComplex aabcd = abcdvec[2][1][1][1];
  TComplex abbcd = abcdvec[1][2][1][1];
  TComplex aabbcd = abcdvec[2][2][1][1];
  return (-120.*aabbcd + 48.*a*abbcd + 24.*ab*abcd + 16.*abc*abd + 12.*abbd*ac +
 8.*abb*acd + 12.*abbc*ad + 48.*aabcd*b - 24.*a*abcd*b - 8.*abd*ac*b -
 8.*ab*acd*b - 8.*abc*ad*b - 6.*aacd*b*b + 4.*a*acd*b*b + 2.*ac*ad*b*b +
 6.*aacd*bb - 4.*a*acd*bb - 2.*ac*ad*bb + 4.*aad*bbc - 4.*a*ad*bbc -
 6.*a*a*bbcd + 6.*aa*bbcd + 4.*aac*bbd - 4.*a*ac*bbd + 12.*aabd*bc -
 8.*a*abd*bc - 4.*ab*ad*bc - 4.*aad*b*bc + 4.*a*ad*b*bc + 8.*aab*bcd -
 8.*a*ab*bcd + 4.*a*a*b*bcd - 4.*aa*b*bcd + 12.*aabc*bd - 8.*a*abc*bd -
 4.*ab*ac*bd - 4.*aac*b*bd + 4.*a*ac*b*bd + 2.*a*a*bc*bd - 2.*aa*bc*bd +
 24.*aabbd*c - 12.*a*abbd*c - 8.*ab*abd*c - 4.*abb*ad*c - 12.*aabd*b*c +
 8.*a*abd*b*c + 4.*ab*ad*b*c + 2.*aad*b*b*c - 2.*a*ad*b*b*c -
 2.*aad*bb*c + 2.*a*ad*bb*c + 2.*a*a*bbd*c - 2.*aa*bbd*c - 4.*aab*bd*c +
 4.*a*ab*bd*c - 2.*a*a*b*bd*c + 2.*aa*b*bd*c + 6.*aabb*cd - 2.*ab*ab*cd -
 4.*a*abb*cd - 4.*aab*b*cd + 4.*a*ab*b*cd - a*a*b*b*cd + aa*b*b*cd +
 a*a*bb*cd - aa*bb*cd + 24.*aabbc*d - 12.*a*abbc*d - 8.*ab*abc*d -
 4.*abb*ac*d - 12.*aabc*b*d + 8.*a*abc*b*d + 4.*ab*ac*b*d + 2.*aac*b*b*d -
 2.*a*ac*b*b*d - 2.*aac*bb*d + 2.*a*ac*bb*d + 2.*a*a*bbc*d - 2.*aa*bbc*d -
 4.*aab*bc*d + 4.*a*ab*bc*d - 2.*a*a*b*bc*d + 2.*aa*b*bc*d - 6.*aabb*c*d +
 2.*ab*ab*c*d + 4.*a*abb*c*d + 4.*aab*b*c*d - 4.*a*ab*b*c*d +
 a*a*b*b*c*d - aa*b*b*c*d - a*a*bb*c*d + aa*bb*c*d).Re();
}
  double getStdAABBDD(TComplex abcdvec[3][3][3][3]){
  TComplex a = abcdvec[1][0][0][0];
  TComplex b = abcdvec[0][1][0][0];
  TComplex d = abcdvec[0][0][0][1];
  TComplex aa = abcdvec[2][0][0][0];
  TComplex bb = abcdvec[0][2][0][0];
  TComplex dd = abcdvec[0][0][0][2];
  TComplex ab = abcdvec[1][1][0][0];
  TComplex ad = abcdvec[1][0][0][1];
  TComplex bd = abcdvec[0][1][0][1];
  TComplex aab = abcdvec[2][1][0][0];
  TComplex aad = abcdvec[2][0][0][1];
  TComplex abb = abcdvec[1][2][0][0];
  TComplex add = abcdvec[1][0][0][2];
  TComplex abd = abcdvec[1][1][0][1];
  TComplex bbd = abcdvec[0][2][0][1];
  TComplex bdd = abcdvec[0][1][0][2];
  TComplex aabb = abcdvec[2][2][0][0];
  TComplex aadd = abcdvec[2][0][0][2];
  TComplex aabd = abcdvec[2][1][0][1];
  TComplex abbd = abcdvec[1][2][0][1];
  TComplex abdd = abcdvec[1][1][0][2];
  TComplex bbdd = abcdvec[0][2][0][2];
  TComplex aabbd = abcdvec[2][2][0][1];
  TComplex aabdd = abcdvec[2][1][0][2];
  TComplex abbdd = abcdvec[1][2][0][2];
  TComplex aabbdd  = abcdvec[2][2][0][2];
  return (-120.*aabbdd + 48.*a*abbdd + 16.*abd*abd + 24.*ab*abdd + 24.*abbd*ad +
 8.*abb*add + 48.*aabdd*b - 24.*a*abdd*b - 16.*abd*ad*b - 8.*ab*add*b -
 6.*aadd*b*b + 2.*ad*ad*b*b + 4.*a*add*b*b + 6.*aadd*bb - 2.*ad*ad*bb -
 4.*a*add*bb + 8.*aad*bbd - 8.*a*ad*bbd - 6.*a*a*bbdd + 6.*aa*bbdd +
 24.*aabd*bd - 16.*a*abd*bd - 8.*ab*ad*bd - 8.*aad*b*bd + 8.*a*ad*b*bd +
 2.*a*a*bd*bd - 2.*aa*bd*bd + 8.*aab*bdd - 8.*a*ab*bdd + 4.*a*a*b*bdd -
 4.*aa*b*bdd + 48.*aabbd*d - 24.*a*abbd*d - 16.*ab*abd*d - 8.*abb*ad*d -
 24.*aabd*b*d + 16.*a*abd*b*d + 8.*ab*ad*b*d + 4.*aad*b*b*d -
 4.*a*ad*b*b*d - 4.*aad*bb*d + 4.*a*ad*bb*d + 4.*a*a*bbd*d - 4.*aa*bbd*d -
 8.*aab*bd*d + 8.*a*ab*bd*d - 4.*a*a*b*bd*d + 4.*aa*b*bd*d - 6.*aabb*d*d +
 2.*ab*ab*d*d + 4.*a*abb*d*d + 4.*aab*b*d*d - 4.*a*ab*b*d*d +
 a*a*b*b*d*d - aa*b*b*d*d - a*a*bb*d*d + aa*bb*d*d + 6.*aabb*dd -
 2.*ab*ab*dd - 4.*a*abb*dd - 4.*aab*b*dd + 4.*a*ab*b*dd - a*a*b*b*dd +
 aa*b*b*dd + a*a*bb*dd - aa*bb*dd).Re();
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
    return (a*a*b*b - aa*b*b - a*a*bb - 4.*a*ab*b + aa*bb + 2.*ab*ab  + 4.*aab*b + 4.*a*abb - 6.*aabb).Re(); }
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
  double getStdABCD(TComplex abcdvec[3][3][3][3]){
  TComplex a = abcdvec[1][0][0][0];
  TComplex b = abcdvec[0][1][0][0];
  TComplex c = abcdvec[0][0][1][0];
  TComplex d = abcdvec[0][0][0][1];
  TComplex ab  = abcdvec[1][1][0][0];
  TComplex ac = abcdvec[1][0][1][0];
  TComplex ad = abcdvec[1][0][0][1];
  TComplex bc = abcdvec[0][1][1][0];
  TComplex bd = abcdvec[0][1][0][1];
  TComplex cd = abcdvec[0][0][1][1];
  TComplex abc = abcdvec[1][1][1][0];
  TComplex abd = abcdvec[1][1][0][1];
  TComplex acd = abcdvec[1][0][1][1];
  TComplex bcd = abcdvec[0][1][1][1];
  TComplex abcd = abcdvec[1][1][1][1];
  return (-6.*abcd + 2.*acd*b + ad*bc + 2.*a*bcd + ac*bd + 2.*abd*c - ad*b*c -
 a*bd*c + ab*cd - a*b*cd + 2.*abc*d - ac*b*d - a*bc*d - ab*c*d +
 a*b*c*d).Re();
 }
  double getStdABDD(TComplex abcdvec[3][3][3][3]){
  TComplex a = abcdvec[1][0][0][0];
  TComplex b = abcdvec[0][1][0][0];
  TComplex d = abcdvec[0][0][0][1];
  TComplex ab  = abcdvec[1][1][0][0];
  TComplex ad = abcdvec[1][0][0][1];
  TComplex bd = abcdvec[0][1][0][1];
  TComplex dd = abcdvec[0][0][0][2];
  TComplex abd = abcdvec[1][1][0][1];
  TComplex add = abcdvec[1][0][0][2];
  TComplex bdd = abcdvec[0][1][0][2];
  TComplex abdd = abcdvec[1][1][0][2];
  return (a*b*d*d - a*b*dd - 2.*a*bd*d - 2.*ad*b*d - ab*d*d
  + 2.*add*b + 2.*a*bdd + 4.*abd*d + ab*dd + 2.*ad*bd
  - 6.*abdd).Re(); }
  double getStdABC(TComplex abcdvec[3][3][3][3]){
  TComplex a = abcdvec[1][0][0][0];
  TComplex b = abcdvec[0][1][0][0];
  TComplex c = abcdvec[0][0][1][0];
  TComplex ab = abcdvec[1][1][0][0];
  TComplex ac = abcdvec[1][0][1][0];
  TComplex bc = abcdvec[0][1][1][0];
  TComplex abc = abcdvec[1][1][1][0];
    return (a*b*c - ab*c - ac*b - a*bc + 2.*abc).Re(); }
  double getStdABD(TComplex abcdvec[3][3][3][3]){
  TComplex a = abcdvec[1][0][0][0];
  TComplex b = abcdvec[0][1][0][0];
  TComplex d = abcdvec[0][0][0][1];
  TComplex ab = abcdvec[1][1][0][0];
  TComplex ad = abcdvec[1][0][0][1];
  TComplex bd = abcdvec[0][1][0][1];
  TComplex abd = abcdvec[1][1][0][1];
  return (a*b*d - ab*d - ad*b - a*bd + 2.*abd).Re(); }
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

  const static Int_t nCorr = 29;
  const static Int_t nHistCent = 7;

  TList *corrLists[nCorr];
  TProfile *corrProfiles[nCorr];
  
  TString corrlistnames[nCorr] = {
    "v24pt2", "v24pt_6pc", "v24_6pc",
    "v22pt2", "v22pt_4pc", "v22_4pc",
    "v22pt", "v22_3pc", "v22",
    "v24", "pt", "pt2",
    "v24pt_6pc_w","v24_6pc_w","v22pt_4pc_w","v22_4pc_w","v22_3pc_w",
    "ptA", "ptB", "pt2AB", "pt2AA","pt2BB",
    "pt3AAB", "pt3ABB", "pt4AABB",
    "ptmid", "pt2mid",
    "v22pt_gap", "v22_gap"};

  TH1D *corrHists[nCorr][nHistCent];

  Bool_t bComputeEtaGap;            // kTRUE: Calculate 2-particle eta gaps (default: kFALSE).
  Float_t fEtaGap;                  // Value of the eta gap itself.


  //ClassDef(AliAnalysisPtVn, 1);
};

#endif  // ALIANALYSISPtVn
