#ifndef AliJFFlucAnalysisTProfile_cxx
#define AliJFFlucAnalysisTProfile_cxx

#include "AliJHistManager.h"
#include <TComplex.h>
#include <TFile.h>
#include <TF3.h>

class TClonesArray;

class AliJFFlucAnalysisTProfile{// : public AliAnalysisTaskSE {
public:
	AliJFFlucAnalysisTProfile();
	AliJFFlucAnalysisTProfile(const char *name);
	AliJFFlucAnalysisTProfile(const AliJFFlucAnalysisTProfile& a); // not implemented
	AliJFFlucAnalysisTProfile& operator=(const AliJFFlucAnalysisTProfile& ap); // not implemented

	virtual ~AliJFFlucAnalysisTProfile();
	virtual void UserCreateOutputObjects(const std::string managerName = "jfluc");
	virtual void UserExec(Option_t *option);
	virtual void Terminate(Option_t *);

	void SetInputList(TClonesArray *inputarray){fInputList = inputarray;}
	TClonesArray * GetInputList() const{return fInputList;}
	void SetEventCentrality( float cent ){fCent = cent;}
	float GetEventCentrality() const{return fCent;}
	void SetEventImpactParameter( float ip ){ fImpactParameter = ip; }
	void SetEventVertex( const double *vtx ){ fVertex = vtx; }
	//void SetIsPhiModule( Bool_t isphi ){//{ IsPhiModule = isphi ; }
	//void SetPhiModuleHistos( int cent, int sub, TH1D *hModuledPhi);

	void SetEtaRange( double eta_min, double eta_max){fEta_min = eta_min; fEta_max = eta_max; }
	void SetEffConfig( int Mode, int FilterBit ){ fEffMode = Mode; fEffFilterBit = FilterBit; cout << "fEffMode set = " << fEffMode << endl;}
	//void SetIsSCptdep( Bool_t isSCptdep ){ IsSCptdep = isSCptdep; cout << "doing addtional loop to check SC pt dep = "<< IsSCptdep << endl;}
	//void SetSCwithQC(Bool_t isSCwithQC){ IsSCwithQC = isSCwithQC; cout << "doing additinal loop for SC results with QC method = " << IsSCwithQC << endl;}
	//void SetEbEWeight(Bool_t isEbEWeighted){ IsEbEWeighted = isEbEWeighted; cout << "use event weight = " << IsEbEWeighted << endl;}
	void SetQCEtaCut( Double_t QC_eta_cut_min, Double_t QC_eta_cut_max, Double_t QC_eta_gap_half){
		fQC_eta_cut_min = QC_eta_cut_min;
		fQC_eta_cut_max = QC_eta_cut_max;
		fQC_eta_gap_half = QC_eta_gap_half;
		cout<<"setting eta range for QC" << fQC_eta_cut_min << "~" << fQC_eta_cut_max << endl;
	}
	void SetPhiWeights(TH1 *p){
		pPhiWeights = p;
	}
	void SetPhiWeights(TF3 *p){
		pPhiWeightsAna = p;
	}

	void SetEventTracksQA(unsigned int tpc, unsigned int glb){ fTPCtrks = (float)tpc; fGlbtrks = (float)glb;}
	void SetEventFB32TracksQA(unsigned int fb32, unsigned int fb32tof){ fFB32trks = (float)fb32; fFB32TOFtrks = (float)fb32tof;}
	
	//TComplex CalculateQnSP( double eta1, double eta2, int harmonics);

	TComplex Get_Qn_pt(double eta1, double eta2, int harmonics, int ipt, double pt_min, double pt_max);
	double Get_QC_Vn( double QnA_real, double QnA_img, double QnB_real, double QnB_img);
	void Fill_QA_plot(double eta1, double eta2 );

	double Get_ScaledMoments( int k, int harmonics);

	// new function for QC method //
	void CalculateQvectorsQC(double, double);
	TComplex Q(int n, int p);
	TComplex Two( int n1, int n2);
	TComplex Four( int n1, int n2, int n3, int n4);

	// Getter for single vn
	Double_t Get_vn( int ih, int imethod ){ return fSingleVn[ih][imethod]; } // method 0:SP, 1:QC(with eta gap), 2:QC(without eta gap)

	enum SUBEVENT{
		SUBEVENT_A = 0x1,
		SUBEVENT_B = 0x2
	};
	void SelectSubevents(UInt_t _subeventMask){
		subeventMask = _subeventMask;
	}
	enum BINNING{
		BINNING_CENT_PbPb,
		BINNING_MULT_PbPb_1,
		BINNING_MULT_pPb_1
	};
	void SetBinning(BINNING _binning){
		binning = _binning;
	}
	enum{
		FLUC_PHI_CORRECTION = 0x2,
		FLUC_SCPT = 0x4,
		FLUC_EBE_WEIGHTING = 0x8
	};
	void AddFlags(UInt_t nflags){
		flags |= nflags;
	}

	static Double_t CentBin_PbPb_default[][2];
	static Double_t MultBin_PbPb_1[][2];
	static Double_t MultBin_pPb_1[][2];
	static Double_t (*pBin[3])[2];
	static Double_t pttJacek[74];
	//static UInt_t CentralityTranslationMap[CENTN_NAT];
	//static UInt_t NCentBin;
	static UInt_t NBin[3];
	static UInt_t NpttJacek;

	//static int GetCentralityClass(Double_t);
	//static int GetMultiplicityBin(Double_t, BINNING);
	static int GetBin(Double_t, BINNING);

	enum{kH0, kH1, kH2, kH3, kH4, kH5, kH6, kH7, kH8, kH9, kH10, kH11, kH12, kNH}; //harmonics
	enum{kK0, kK1, kK2, kK3, kK4, nKL}; // order
#define kcNH kH6 //max second dimension + 1
//private:

	bool b_vnpt_graphs;
	bool fDebug;
	TClonesArray *fInputList;
	const double *fVertex;//!
	TH1 *pPhiWeights;//!
	TF3 *pPhiWeightsAna;//!
	Float_t	fCent;
	Float_t	fImpactParameter;
	int fCBin;
	int fEffMode;
	int fEffFilterBit;
	float fTPCtrks;
	float fGlbtrks;
	float fFB32trks;
	float fFB32TOFtrks;
	UInt_t subeventMask;
	BINNING binning;
	UInt_t flags;
	Double_t fSingleVn[kNH][3]; // 3 methods

	double fEta_min;
	double fEta_max;
	double NSubTracks[2];

	Double_t fQC_eta_cut_min;
	Double_t fQC_eta_cut_max;
	Double_t fQC_eta_gap_half;

	TComplex QvectorQC[kNH][nKL];
	TComplex QvectorQCeta10[2][kNH][nKL]; // ksub

	AliJHistManager * fHMG;//!

	AliJBin fBin_Subset;//!
	AliJBin fBin_h;//!
	AliJBin fBin_k;//!
	AliJBin fBin_hh;//!
	AliJBin fBin_kk;//!
	AliJBin fHistCentBin;//!
	AliJBin fVertexBin;//! // x, y, z
	AliJBin fCorrBin;//!

	AliJTH1D fh_cent;//! // for cent dist
	AliJTH1D fh_ImpactParameter;//! // for impact parameter for mc
	AliJTH1D fh_ecc;//!
	AliJTH1D fh_pt;//! // for pt dist of tracks
	AliJTH1D fh_pt_Pion;//! // for pt dist of tracks
	AliJTH1D fh_pt_Kaon;//! // for pt dist of tracks
	AliJTH1D fh_pt_Proton;//! // for pt dist of tracks
	AliJTH1D fh_eta;//! // for eta dist of tracks
	AliJTH1D fh_eta_Pion;//! // for eta dist of tracks
	AliJTH1D fh_eta_Kaon;//! // for eta dist of tracks
	AliJTH1D fh_eta_Proton;//! // for eta dist of tracks
	AliJTH1D fh_mult_eta05;//!
	AliJTH1D fh_mult_eta05_Pion;//!
	AliJTH1D fh_mult_eta05_Kaon;//!
	AliJTH1D fh_mult_eta05_Proton;//!
	AliJTH1D fh_mult_rap05;//!
	AliJTH1D fh_mult_rap05_Pion;//!
	AliJTH1D fh_mult_rap05_Kaon;//!
	AliJTH1D fh_mult_rap05_Proton;//!
	AliJTH1D fh_eta_NoPtCut;//! // for eta dist of tracks
	AliJTH1D fh_phi;//! // for phi dist [ic][isub]

	AliJTH1D fh_ntracks;//! // for number of tracks dist
	AliJTProfile fh_vn;//!  // single vn^k  array [ih][ik][iCent]
	AliJTProfile fh_vna;//! // single vn^k with autocorrelation removed (up to a limited order)
	AliJTProfile fh_vn_vn;//! // combination for <vn*vn> [ih][ik][ihh][ikk][iCent]

	AliJTProfile fh_correlator;//! // some more complex correlators

	// additional variables for ptbins(Standard Candles only)
	enum{kPt0, kPt1, kPt2, kPt3, kPt4, kPt5, kPt6, kPt7, N_ptbins};
	double NSubTracks_pt[2][N_ptbins];

	AliJBin fBin_Nptbins;//!
	AliJTH1D fh_SC_ptdep_4corr;//! // for < vn^2 vm^2 >
	AliJTH1D fh_SC_ptdep_2corr;//!  // for < vn^2 >
	// additinal variables for SC with QC
	AliJTH1D fh_SC_with_QC_4corr;//! // for <vn^2 vm^2>
	AliJTH1D fh_SC_with_QC_2corr;//! // for <vn^2>
	AliJTH1D fh_SC_with_QC_2corr_eta10;//!
	//AliJTH2D fh_QvectorQC;//! // check for Q-vec dist for [ic][ih]
	//AliJTH1D fh_QvectorQCphi;//!
	AliJTH1D fh_evt_SP_QC_ratio_2p;//! // check SP QC evt by evt ratio
	AliJTH1D fh_evt_SP_QC_ratio_4p;//! // check SP QC evt by evt ratio
	//ClassDef(AliJFFlucAnalysisTProfile, 1); // example of analysis

	static const int kNPtBinsQA = 22;  // Number of QA pT bins from the data
	Double_t fPtBinsQA[kNPtBinsQA+1] = {0.0, 0.2, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.125, 1.375, 1.625, 1.875, 2.125, 2.375, 2.75, 3.25, 3.75, 4.5, 5.5, 7.0, 9.0};  // Binning matching the v2 measurement in Pb-Pb at 5.02 TeV
};

#endif
