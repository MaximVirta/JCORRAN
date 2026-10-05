#include <cmath>
#include <TH1D.h>
#include <TMath.h>
#include <TComplex.h>
#include <TClonesArray.h>
#include "AliJBaseTrack.h"
#include "AliJFFlucAnalysisReduced.h"
#pragma GCC diagnostic warning "-Wall"

namespace {
inline bool IsFinite(double x) { return std::isfinite(x); }
}

//ClassImp(AliJFFlucAnalysisReduced)

		//________________________________________________________________________
AliJFFlucAnalysisReduced::AliJFFlucAnalysisReduced() :
	//: AliAnalysisTaskSE(),
	fInputList(0),
	fCent(0),
	fCBin(0),
	fHMG(0),
	fBin_Subset(),
	fBin_k(),
	fBin_hh(),
	fBin_etaGap(),
	fHistCentBin(),
	fh_cent(),
	fh_ImpactParameter(),
	fh_ecc(),
	fh_pt_eta05(),
	fh_pt_rap05(),
	fh_eta(),
	fh_mult_eta05(),
	fh_mult_rap05(),
	fh_phi(),
	fh_ntracks(),
	fh_vna(),
	fh_SC_with_QC_4corr()
{
	subeventMask = SUBEVENT_A|SUBEVENT_B;
	binning = BINNING_CENT_PbPb;
	flags = 0;
	fEta_min = 0;
	fEta_max = 0.8;
	for(int i=0; i<kEtaGaps; i++)
		EtaGapMin[i] = 0;
	fImpactParameter = -1;
	fDebug = false;
}

//________________________________________________________________________
AliJFFlucAnalysisReduced::AliJFFlucAnalysisReduced(const char *name) :
	//: AliAnalysisTaskSE(name),
	fInputList(0),
	fCent(0),
	fCBin(0),
	fHMG(0),
	fBin_Subset(),
	fBin_k(),
	fBin_hh(),
	fBin_etaGap(),
	fHistCentBin(),
	fh_cent(),
	fh_ImpactParameter(),
	fh_ecc(),
	fh_pt_eta05(),
	fh_pt_rap05(),
	fh_eta(),
	fh_mult_eta05(),
	fh_mult_rap05(),
	fh_phi(),
	fh_ntracks(),
	fh_vna(),
	fh_SC_with_QC_4corr()
{
	cout << "analysis task created " << endl;

	subeventMask = SUBEVENT_A|SUBEVENT_B;
	binning = BINNING_CENT_PbPb;
	flags = 0;
	fEta_min = 0;
	fEta_max = 0.8;
	for(int i=0; i<kEtaGaps; i++)
		EtaGapMin[i] = 0;
	fImpactParameter = -1;
	fDebug = false;
}

Double_t AliJFFlucAnalysisReduced::CentBin_PbPb_default[][2] = {{0,1},{1,2},{2,5},{5,10},{10,20},{20,30},{30,40},{40,50},{50,60}};
Double_t AliJFFlucAnalysisReduced::CentBin_OO_central[][2] = {{0,1},{1,2},{2,5},{5,10},{10,15},{15,20},{20,25},{25,30},{30,35},{35,40},{40,45},{45,50}};
Double_t AliJFFlucAnalysisReduced::MultBin_PbPb_1[][2] = {{25.245,25.611},{31.555,31.982},{37.898,38.387},{44.251,44.803},{50.584,51.198},{56.942,57.62},{66.341,67.113},{79.03,79.93},{91.764,92.792},{104.431,105.587},{117.103,118.387},{136.063,137.541},{161.3,163.033},{186.707,188.698},{212.126,214.374},{281.085,284.042},{344.727,348.332},{408.197,412.447},{471.798,476.696},{535.041,540.583},{598.475,604.663},{662.018,668.854},{725.484,732.967},{788.636,796.763},{852.345,861.121},{915.632,925.054},{979.134,989.204},{1042.376,1053.091},{1105.751,1117.114},{1169.205,1181.215},{1232.423,1245.079},{1295.89,1309.194},{1359.578,1373.533},{1422.675,1437.274},{1485.99,1501.236},{1549.608,1565.503},{1612.809,1629.351},{1676.158,1693.347},{1739.647,1757.486},{1802.865,1821.351},{1866.328,1885.462},{1929.777,1949.561},{1993.119,2013.553},{2056.422,2077.504},{2119.716,2141.448},{2182.935,2205.317},{2246.261,2269.296},{2309.523,2333.209},{2373.038,2397.381},{2436.146,2461.145},{2499.398,2525.057},{2563.07,2589.39},{2626.482,2653.462},{2689.435,2717.075},{2753.073,2781.382},{2815.911,2844.889},{2879.821,2909.486},{2943.2,2973.535},{3006.237,3037.292},{3068.718,3100.47},{3132.669,3165.259}};
Double_t AliJFFlucAnalysisReduced::MultBin_pPb_1[][2] = {{30.651,31.125},{42.627,43.263},{54.598,55.399},{66.613,67.581},{78.665,79.802},{90.789,92.098},{102.805,104.286},{114.903,116.569}};
Double_t (*AliJFFlucAnalysisReduced::pBin[4])[2] = {&CentBin_PbPb_default[0],&CentBin_OO_central[0],&MultBin_PbPb_1[0],&MultBin_pPb_1[0]};
UInt_t AliJFFlucAnalysisReduced::NBin[4] = {
	sizeof(AliJFFlucAnalysisReduced::CentBin_PbPb_default)/sizeof(AliJFFlucAnalysisReduced::CentBin_PbPb_default[0]),
	sizeof(AliJFFlucAnalysisReduced::CentBin_OO_central)/sizeof(AliJFFlucAnalysisReduced::CentBin_OO_central[0]),
	sizeof(AliJFFlucAnalysisReduced::MultBin_PbPb_1)/sizeof(AliJFFlucAnalysisReduced::MultBin_PbPb_1[0]),
	sizeof(AliJFFlucAnalysisReduced::MultBin_pPb_1)/sizeof(AliJFFlucAnalysisReduced::MultBin_pPb_1[0])
};

//________________________________________________________________________
AliJFFlucAnalysisReduced::AliJFFlucAnalysisReduced(const AliJFFlucAnalysisReduced& a):
	//AliAnalysisTaskSE(a.GetName()),
	fInputList(a.fInputList),
	fCent(a.fCent),
	fCBin(a.fCBin),
	fHMG(a.fHMG),
	fBin_Subset(a.fBin_Subset),
	fBin_k(a.fBin_k),
	fBin_hh(a.fBin_hh),
	fBin_etaGap(a.fBin_etaGap),
	fHistCentBin(a.fHistCentBin),
	fh_cent(a.fh_cent),
	fh_ImpactParameter(a.fh_ImpactParameter),
	fh_ecc(a.fh_ecc),
	fh_pt_eta05(a.fh_pt_eta05),
	fh_pt_rap05(a.fh_pt_rap05),
	fh_eta(a.fh_eta),
	fh_mult_eta05(a.fh_mult_eta05),
	fh_mult_rap05(a.fh_mult_rap05),
	fh_phi(a.fh_phi),
	fh_ntracks(a.fh_ntracks),
	fh_vna(a.fh_vna),
	fh_SC_with_QC_4corr(a.fh_SC_with_QC_4corr)
{
	//copy constructor
	//	DefineOutput(1, TList::Class() );
}
//________________________________________________________________________
AliJFFlucAnalysisReduced& AliJFFlucAnalysisReduced::operator = (const AliJFFlucAnalysisReduced& ap){
	// assignment operator
	this->~AliJFFlucAnalysisReduced();
	new(this) AliJFFlucAnalysisReduced(ap);
	return *this;
}

//________________________________________________________________________
void AliJFFlucAnalysisReduced::UserCreateOutputObjects(const std::string managerName){
	
	fHMG = new AliJHistManager("AliJFFlucHistManager",managerName.data());
	// set AliJBin here //
	fBin_Subset .Set("Sub","Sub","Sub:%d", AliJBin::kSingle).SetBin(2);
	fBin_k .Set("K","K","K:%d", AliJBin::kSingle).SetBin(nKL);

	fBin_hh .Set("NHH","NHH","NHH:%d", AliJBin::kSingle).SetBin(kcNH);

	fBin_etaGap .Set("EtaGap","EtaGap","EtaGap:%d", AliJBin::kSingle).SetBin(kEtaGaps);

	//TODO: index with binning the array of pointers
	if(binning != BINNING_CENT_PbPb && binning != BINNING_CENT_OO)
		fHistCentBin.Set("MultBin","MultBin","Cent:%d",AliJBin::kSingle).SetBin(NBin[binning]);
	else fHistCentBin.Set("CentBin","CentBin","Cent:%d",AliJBin::kSingle).SetBin(NBin[binning]);

	// set AliJTH1D here //
	fh_cent
		<< TH1D("h_cent","h_cent", 200, 0, 100)
		<< "END" ;

	fh_ImpactParameter
		<< TH1D("h_IP", "h_IP", 400, -2, 20)
		<< "END" ;
	
	fh_ecc
		<< TProfile("h_ecc","h_ecc",kcNH,0.,static_cast<float>(kcNH))
		<< fHistCentBin
		<< "END" ;
	
	fh_pt_eta05
		<< TProfile("hPtJacekEta05", "", kSpecies, -0.5, static_cast<float>(kSpecies)-0.5)
		<< fHistCentBin
		<< "END" ;

	fh_pt_rap05
		<< TProfile("hPtJacekRap05", "", kSpecies, -0.5, static_cast<float>(kSpecies)-0.5)
		<< fHistCentBin
		<< "END" ;

	fh_eta
		<< TH1D("h_eta", "h_eta", 40, -2.0, 2.0 )
		<< fHistCentBin
		<< "END" ;
	
	fh_mult_eta05
		<< TProfile("hMultEta05", "", kSpecies, -0.5, static_cast<float>(kSpecies)-0.5)
		<< fHistCentBin
		<< "END" ;

	fh_mult_rap05
		<< TProfile("hMultRap05", "", kSpecies, -0.5, static_cast<float>(kSpecies)-0.5)
		<< fHistCentBin
		<< "END" ;
	fh_phi
		<< TH1D("h_phi", "h_phi", 50,-TMath::Pi(),TMath::Pi())
		<< fHistCentBin << fBin_Subset
		<< "END" ;

	fh_ntracks
		<< TH1D("h_tracks", "h_tracks", 500, 0, 500)
		<< fHistCentBin
		<< "END" ;
	// Changed to TProfiles
	fh_vna
		<< TProfile("hvna","hvna", kcNH, 0., static_cast<float>(kcNH))
		<< fBin_k << fBin_etaGap
		<< fHistCentBin
		<< "END";   // histogram of vn_h^k values for [ik][iEtaGap][iCent]

	fh_SC_with_QC_4corr
		<< TProfile("hQC_SC4p", "hQC_SC4p", kcNH, 0., static_cast<float>(kcNH))
		<< fBin_hh << fBin_etaGap
		<< fHistCentBin
		<< "END" ;

	//AliJTH1D set done.

	if (fDebug) fHMG->Print();
	//fHMG->WriteConfig();

}

//________________________________________________________________________
AliJFFlucAnalysisReduced::~AliJFFlucAnalysisReduced() {
	delete fHMG;
}

#define A i
#define B (1-i)
#define C(u) TComplex::Conjugate(u)
//TODO: conjugate macro
inline TComplex TwoGap(const TComplex (*pQq)[AliJFFlucAnalysisReduced::kNH][AliJFFlucAnalysisReduced::nKL], uint i, uint a, uint b){
	return pQq[A][a][1]*C(pQq[B][b][1]);
}

inline TComplex FourGap22(const TComplex (*pQq)[AliJFFlucAnalysisReduced::kNH][AliJFFlucAnalysisReduced::nKL], uint i, uint a, uint b, uint c, uint d){
	return pQq[A][a][1]*pQq[A][b][1]*C(pQq[B][c][1]*pQq[B][d][1])-pQq[A][a+b][2]*C(pQq[B][c][1]*pQq[B][d][1])-pQq[A][a][1]*pQq[A][b][1]*C(pQq[B][c+d][2])+pQq[A][a+b][2]*C(pQq[B][c+d][2]);
}

inline TComplex SixGap33(const TComplex (*pQq)[AliJFFlucAnalysisReduced::kNH][AliJFFlucAnalysisReduced::nKL], uint i, uint n1, uint n2, uint n3, uint n4, uint n5, uint n6){
	return pQq[A][n1][1]*pQq[A][n2][1]*pQq[A][n3][1]*C(pQq[B][n4][1]*pQq[B][n5][1]*pQq[B][n6][1])-pQq[A][n1][1]*pQq[A][n2][1]*pQq[A][n3][1]*C(pQq[B][n4+n5][2]*pQq[B][n6][1])-pQq[A][n1][1]*pQq[A][n2][1]*pQq[A][n3][1]*C(pQq[B][n4+n6][2]*pQq[B][n5][1])-pQq[A][n1][1]*pQq[A][n2][1]*pQq[A][n3][1]*C(pQq[B][n5+n6][2]*pQq[B][n4][1])+2.0*pQq[A][n1][1]*pQq[A][n2][1]*pQq[A][n3][1]*C(pQq[B][n4+n5+n6][3])-pQq[A][n1+n2][2]*pQq[A][n3][1]*C(pQq[B][n4][1]*pQq[B][n5][1]*pQq[B][n6][1])+pQq[A][n1+n2][2]*pQq[A][n3][1]*C(pQq[B][n4+n5][2]*pQq[B][n6][1])+pQq[A][n1+n2][2]*pQq[A][n3][1]*C(pQq[B][n4+n6][2]*pQq[B][n5][1])+pQq[A][n1+n2][2]*pQq[A][n3][1]*C(pQq[B][n5+n6][2]*pQq[B][n4][1])-2.0*pQq[A][n1+n2][2]*pQq[A][n3][1]*C(pQq[B][n4+n5+n6][3])-pQq[A][n1+n3][2]*pQq[A][n2][1]*C(pQq[B][n4][1]*pQq[B][n5][1]*pQq[B][n6][1])+pQq[A][n1+n3][2]*pQq[A][n2][1]*C(pQq[B][n4+n5][2]*pQq[B][n6][1])+pQq[A][n1+n3][2]*pQq[A][n2][1]*C(pQq[B][n4+n6][2]*pQq[B][n5][1])+pQq[A][n1+n3][2]*pQq[A][n2][1]*C(pQq[B][n5+n6][2]*pQq[B][n4][1])-2.0*pQq[A][n1+n3][2]*pQq[A][n2][1]*C(pQq[B][n4+n5+n6][3])-pQq[A][n2+n3][2]*pQq[A][n1][1]*C(pQq[B][n4][1]*pQq[B][n5][1]*pQq[B][n6][1])+pQq[A][n2+n3][2]*pQq[A][n1][1]*C(pQq[B][n4+n5][2]*pQq[B][n6][1])+pQq[A][n2+n3][2]*pQq[A][n1][1]*C(pQq[B][n4+n6][2]*pQq[B][n5][1])+pQq[A][n2+n3][2]*pQq[A][n1][1]*C(pQq[B][n5+n6][2]*pQq[B][n4][1])-2.0*pQq[A][n2+n3][2]*pQq[A][n1][1]*C(pQq[B][n4+n5+n6][3])+2.0*pQq[A][n1+n2+n3][3]*C(pQq[B][n4][1]*pQq[B][n5][1]*pQq[B][n6][1])-2.0*pQq[A][n1+n2+n3][3]*C(pQq[B][n4+n5][2]*pQq[B][n6][1])-2.0*pQq[A][n1+n2+n3][3]*C(pQq[B][n4+n6][2]*pQq[B][n5][1])-2.0*pQq[A][n1+n2+n3][3]*C(pQq[B][n5+n6][2]*pQq[B][n4][1])+4.0*pQq[A][n1+n2+n3][3]*C(pQq[B][n4+n5+n6][3]);
}
#undef C

//________________________________________________________________________
void AliJFFlucAnalysisReduced::UserExec(Option_t *) {
	// find Centrality
	int trk_number = fInputList->GetEntriesFast();
	fCBin = (binning != BINNING_CENT_PbPb && binning != BINNING_CENT_OO)?
		GetBin((double)trk_number,binning):
		GetBin(fCent,binning); //--- similarly in Task

	if(fCBin == -1)
		return;
	
	fh_ntracks[fCBin]->Fill( trk_number ) ;
	fh_cent->Fill(fCent) ;
	fh_ImpactParameter->Fill( fImpactParameter);
	Fill_QA_plot( fEta_max );

	CalculateQvectorsQC(fEta_max);

	// v2^2 :  k=1  /// remember QnQn = vn^(2k) not k
	// harmonics filled only through v4

	TComplex ncorr[kcNH][nKL];

	// Calculate results for each eta gap method
	for(int ietaGap=kEtaGap00; ietaGap<kEtaGaps; ietaGap++){
		const TComplex (*pQq)[kNH][nKL] = QvectorQC[ietaGap];

		for(int i = 0; i < 2; ++i){
			if((subeventMask & (1<<i)) == 0)
				continue;

			Double_t ref_2p = TwoGap(pQq,i,0,0).Re();
			Double_t ref_4p = FourGap22(pQq,i,0,0,0,0).Re();
			Double_t ref_6p = SixGap33(pQq,i,0,0,0,0,0,0).Re();

			Double_t ebe_2p_weight = 1.0;
			Double_t ebe_4p_weight = 1.0;
			Double_t ebe_6p_weight = 1.0;

			if(flags & FLUC_EBE_WEIGHTING){
				ebe_2p_weight = ref_2p;
				ebe_4p_weight = ref_4p;
				ebe_6p_weight = ref_6p;
			}
			Double_t ref_2Np[nKL-1] = {
				ref_2p,
				ref_4p,
				ref_6p
			};
			Double_t ebe_2Np_weight[nKL-1] = {
				ebe_2p_weight,
				ebe_4p_weight,
				ebe_6p_weight
			};

			for(int ih=2; ih<=kMaxVn; ih++){
				ncorr[ih][1] = TwoGap(pQq,i,ih,ih); // vn{2}
				ncorr[ih][2] = FourGap22(pQq,i,ih,ih,ih,ih); // vn{4}
				ncorr[ih][3] = SixGap33(pQq,i,ih,ih,ih,ih,ih,ih); // vn{6}
			}
			
			for(int ih=2; ih<=kMaxVn; ih++){
				for(int ik=1; ik<nKL-1; ik++){ // 2k(0) =1, 2k(1) =2, 2k(2)=4....
					Double_t vna_val = ncorr[ih][ik].Re()/ref_2Np[ik-1]; // vn{2k}(eta gap)
					if(IsFinite(vna_val) && IsFinite(ebe_2Np_weight[ik-1]))
						fh_vna[ik][ietaGap][fCBin]->Fill(static_cast<float>(ih) + 0.5, vna_val,ebe_2Np_weight[ik-1]);
				}

				for(int ihh=2; ihh<=kMaxVn; ihh++){ // SC for n,m <= 4
					TComplex scfour = FourGap22(pQq,i,ih,ihh,ih,ihh) / ref_4p;
					if(IsFinite(scfour.Re()) && IsFinite(ebe_4p_weight))
						fh_SC_with_QC_4corr[ihh][ietaGap][fCBin]->Fill(static_cast<float>(ih) + 0.5, scfour.Re(), ebe_4p_weight );
				}
			}


		}
	}

	// Calculate results for no gap

	Double_t ref_2p = Two(0,0).Re();
	Double_t ref_4p = Four(0,0,0,0).Re();

	Double_t ebe_2p_weight = 1.0;
	Double_t ebe_4p_weight = 1.0;

	if(flags & FLUC_EBE_WEIGHTING){
		ebe_2p_weight = ref_2p;
		ebe_4p_weight = ref_4p;
	}
	Double_t ref_2Np[nKL-1] = {
		ref_2p,
		ref_4p
	};
	Double_t ebe_2Np_weight[nKL-1] = {
		ebe_2p_weight,
		ebe_4p_weight
	};


	for(int ih=2; ih<kMaxVn; ih++){
		ncorr[ih][1] = Two(ih, -ih); // vn{2}
		ncorr[ih][2] = Four(ih,ih, -ih, -ih); // vn{4}
	}
	
	for(int ih=2; ih<kMaxVn; ih++){
		for(int ik=1; ik<nKL-2; ik++){ // 2k(0) =1, 2k(1) =2, 2k(2)=4....
			Double_t vna_val = ncorr[ih][ik].Re()/ref_2Np[ik-1]; // vn{2k}(no gap)
			if(IsFinite(vna_val) && IsFinite(ebe_2Np_weight[ik-1]))
				fh_vna[ik][kNoGap][fCBin]->Fill(static_cast<float>(ih) + 0.5, vna_val, ebe_2Np_weight[ik-1]);
		}
	}
}

//________________________________________________________________________
void AliJFFlucAnalysisReduced::Terminate(Option_t *)
{
	cout<<"Sucessfully Finished"<<endl;
}
//________________________________________________________________________
//________________________________________________________________________
void AliJFFlucAnalysisReduced::Fill_QA_plot( Double_t etamax )
{
	Long64_t ntracks = fInputList->GetEntriesFast();

	for( Long64_t it=0; it < ntracks; it++){
		AliJBaseTrack *itrack = (AliJBaseTrack*)fInputList->At(it); // load track
		Double_t eta = itrack->Eta();

		if(TMath::Abs(eta) > etamax)
			continue;

		Double_t phi = itrack->Phi();

		fh_eta[fCBin]->Fill(eta);
		fh_phi[fCBin][(int)(eta > 0.0)]->Fill(phi);
	}
	// end of track loop
}
//________________________________________________________________________
void AliJFFlucAnalysisReduced::CalculateQvectorsQC(double etamax){
	// calcualte Q-vector for QC method ( no subgroup )
	//init
	for(int ietaGap=0; ietaGap<kEtaGaps; ietaGap++){
		for(int isub=0; isub<2; isub++){
			for(int ih=0; ih<kNH; ih++){
				for(int ik=0; ik<nKL; ++ik){
					QvectorQC[ietaGap][isub][ih][ik] = TComplex(0,0);
				}
			}
		}
	}
	//Calculate Q-vector with particle loop
	Long64_t ntracks = fInputList->GetEntriesFast(); // all tracks from Task input
	for( Long64_t it=0; it<ntracks; it++){
		AliJBaseTrack *itrack = (AliJBaseTrack*)fInputList->At(it); // load track
		Double_t eta = itrack->Eta();

		if( TMath::Abs(eta) > etamax )
			continue;
		/////////////////////////////////////////////////

		int isub = (int)(eta > 0.0);
		Double_t phi = itrack->Phi();
		Double_t pt = itrack->Pt();
		if(pt < 0.2 || pt > 5.0)
			continue;

		for(int ih=0; ih<kNH; ih++){
			TComplex q[nKL];
			for(int ik=0; ik<nKL; ik++){
				q[ik] = TComplex(TMath::Cos(ih*phi),TMath::Sin(ih*phi));

				// Fill both subevents for no gap
				QvectorQC[kNoGap][0][ih][ik] += q[ik];
				QvectorQC[kNoGap][1][ih][ik] += q[ik];

				// Fill subevents for eta gap
				for(int ietaGap=kEtaGap00; ietaGap<kEtaGaps; ietaGap++){
					if(TMath::Abs(eta) > EtaGapMin[ietaGap])
						QvectorQC[ietaGap][isub][ih][ik] += q[ik];
				}
			}
		}
	} // track loop done.
}
//________________________________________________________________________
TComplex AliJFFlucAnalysisReduced::Q(int n, int p){
	// Return QvectorQC
	// Q{-n, p} = Q{n, p}*
	if(n >= 0)
		return QvectorQC[kNoGap][0][n][p];
	return TComplex::Conjugate(QvectorQC[kNoGap][0][-n][p]);
}
//________________________________________________________________________
TComplex AliJFFlucAnalysisReduced::Two(int n1, int n2){
	// two-particle correlation <exp[i(n1*phi1 + n2*phi2)]>
	//	cout << "TWO FUNCTION " << Q(n1,1) << "*" << Q(n2,1) << " - " << Q(n1+n2 , 2) << endl;
	TComplex two = Q(n1, 1) * Q(n2, 1) - Q( n1+n2, 2);
	return two;
}
//________________________________________________________________________
TComplex AliJFFlucAnalysisReduced::Four( int n1, int n2, int n3, int n4){
	TComplex four =
		Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
		- Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
		+ Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
		+ 2.*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
		+ 2.*Q(n2,1)*Q(n1+n3+n4,3)+2.*Q(n1,1)*Q(n2+n3+n4,3)-6.*Q(n1+n2+n3+n4,4);
	return four;
}
//__________________________________________________________________________
int AliJFFlucAnalysisReduced::GetBin(Double_t fq, BINNING _binning){
	for(UInt_t iMbin = 0; iMbin < NBin[_binning]; iMbin++){
		if(fq >= pBin[_binning][iMbin][0] && fq < pBin[_binning][iMbin][1])
			return iMbin;
	}
	return -1;
}

