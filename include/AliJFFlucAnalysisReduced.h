#ifndef AliJFFlucAnalysisReduced_cxx
#define AliJFFlucAnalysisReduced_cxx

#include "AliJHistManager.h"
#include <TComplex.h>

class TClonesArray;

class AliJFFlucAnalysisReduced{// : public AliAnalysisTaskSE {
public:
	AliJFFlucAnalysisReduced();
	AliJFFlucAnalysisReduced(const char *name);
	AliJFFlucAnalysisReduced(const AliJFFlucAnalysisReduced& a); // not implemented
	AliJFFlucAnalysisReduced& operator=(const AliJFFlucAnalysisReduced& ap); // not implemented

	virtual ~AliJFFlucAnalysisReduced();
	virtual void UserCreateOutputObjects(const std::string managerName = "jfluc");
	virtual void UserExec(Option_t *option);
	virtual void Terminate(Option_t *);

	void SetInputList(TClonesArray *inputarray){fInputList = inputarray;}
	TClonesArray * GetInputList() const{return fInputList;}
	void SetEventCentrality( float cent ){fCent = cent;}
	float GetEventCentrality() const{return fCent;}
	void SetEventImpactParameter( float ip ){ fImpactParameter = ip; }
	void SetEtaGapMins(double etaGapMins[]){
		for(int i=0; i<kEtaGaps; i++){
			EtaGapMin[i] = etaGapMins[i];
		}
	}

	void SetEtaRange( double eta_min, double eta_max){fEta_min = eta_min; fEta_max = eta_max; }
	void Fill_QA_plot(double etamax );

	// new function for QC method //
	void CalculateQvectorsQC(double);
	TComplex Q(int n, int p);
	TComplex Two( int n1, int n2);
	TComplex Four( int n1, int n2, int n3, int n4);

	enum SUBEVENT{
		SUBEVENT_A = 0x1,
		SUBEVENT_B = 0x2
	};
	void SelectSubevents(UInt_t _subeventMask){
		subeventMask = _subeventMask;
	}
	enum BINNING{
		BINNING_CENT_PbPb,
		BINNING_CENT_OO,
		BINNING_MULT_PbPb_1,
		BINNING_MULT_pPb_1
	};
	void SetBinning(BINNING _binning){
		binning = _binning;
	}
	enum{
		FLUC_EBE_WEIGHTING = 0x8
	};
	void AddFlags(UInt_t nflags){
		flags |= nflags;
	}

	static Double_t CentBin_PbPb_default[][2];
	static Double_t CentBin_OO_central[][2];
	static Double_t MultBin_PbPb_1[][2];
	static Double_t MultBin_pPb_1[][2];
	static Double_t (*pBin[4])[2];
	static UInt_t NBin[4];

	static int GetBin(Double_t, BINNING);

	enum{kH0, kH1, kH2, kH3, kH4, kH5, kH6, kH7, kH8, kH9, kH10, kH11, kH12, kH13, kH14, kH15, kH16, kNH}; // Q-vector harmonics; v4 autocorrelations need n <= 12
	enum{kK0, kK1, kK2, kK3, kK4, nKL}; // order
	enum{kMaxVn = kH4}; // filled flow harmonics: v2, v3, v4
	// AliJFFlucAnalysis.h and AliJFFlucAnalysisTProfile.h #define kcNH as kH6.
#pragma push_macro("kcNH")
#undef kcNH
	enum{kcNH = kMaxVn + 1}; // max second-harmonic index + 1
#pragma pop_macro("kcNH")
	enum{kNoGap = 0, kEtaGap00 = 1, kEtaGap02 = 2, kEtaGap10 = 3, kEtaGap14 = 4, kEtaGaps = 5}; // eta gap methods
	enum{kCharged, kPion, kKaon, kProton, kSpecies};

	//private:


	bool fDebug;
	TClonesArray *fInputList;
	Float_t	fCent;
	Float_t	fImpactParameter;
	int fCBin;
	UInt_t subeventMask;
	BINNING binning;
	UInt_t flags;

	double fEta_min;
	double fEta_max;
	Double_t EtaGapMin[kEtaGaps];

	TComplex QvectorQC[kEtaGaps][2][kNH][nKL];

	AliJHistManager * fHMG;//!

	AliJBin fBin_Subset;//!
	AliJBin fBin_k;//!
	AliJBin fBin_hh;//!
	AliJBin fBin_etaGap;//!
	AliJBin fHistCentBin;//!

	AliJTH1D fh_cent;//! // for cent dist
	AliJTH1D fh_ImpactParameter;//! // for impact parameter for mc
	AliJTProfile fh_ecc;//! // for eccentricity
	AliJTProfile fh_pt_eta05;//! // for pt dist at eta = 0.5
	AliJTProfile fh_pt_rap05;//! // for pt dist at rap = 0.5
	AliJTH1D fh_eta;//! // for eta dist of tracks
	AliJTProfile fh_mult_eta05;//! // for mult dist at eta = 0.5
	AliJTProfile fh_mult_rap05;//! // for mult dist at rap = 0.5
	AliJTH1D fh_phi;//! // for phi dist [ic][isub]

	AliJTH1D fh_ntracks;//! // for number of tracks dist
	AliJTProfile fh_vna;//! // single vn^k with autocorrelation removed (up to a limited order)

	// additinal variables for SC with QC
	AliJTProfile fh_SC_with_QC_4corr;//! // for <vn^2 vm^2>
	//ClassDef(AliJFFlucAnalysisReduced, 1); // example of analysis
};

#endif
