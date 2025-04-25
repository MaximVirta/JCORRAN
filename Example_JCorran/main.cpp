
#include <stdio.h>
#include <vector>
#include <exception>
#include <fstream>

#include <TRandom3.h>
#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TGraphErrors.h>
#include <TComplex.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TROOT.h>
#include <TStyle.h>


#include "args.hxx"

#include <AliJBaseTrack.h>
#include <AliJFFlucAnalysis.h>
#include <AliAnalysisAnaTwoMultiCorrelations.h>
#include <AliAnalysisSPCRun2.h>
#include <AliAnalysisPtVn.h>
#include <AliAnalysisMeanPt.h>


typedef unsigned int uint;

#define NC 8

int main(int argc, char **pargv){
	args::ArgumentParser parser("JFLUC HIC analyzer","");
	args::Group group_io(parser,"Output files",args::Group::Validators::DontCare);
	args::ValueFlag<std::string> outfilen(group_io,"dst","Output root file name",{'o',"outfile"},"AnalysisResults.root");

	args::Group group_ana(parser,"Analysis",args::Group::Validators::DontCare);

	args::ValueFlag<double> minPt(group_ana,"pTmin","Minimum track pT",{"pTmin"},0.2);
	args::ValueFlag<double> maxPt(group_ana,"pTmax","Maximum track pT",{"pTmax"},5.0);
	args::ValueFlag<double> absEtaMax(group_ana,"absEtaMax","Maximum track |eta|",{"absEtaMax"},0.8);

	try{
		parser.ParseCLI(argc,pargv);

	}catch(args::Help &e){
		std::cout<<parser;
		return 0;

	}catch(args::ParseError &e){
		std::cout<<e.what()<<std::endl;
		std::cout<<parser;
		return 1;

	}

	// SUPPRESSES ROOT WARNINGS!
	gROOT->ProcessLine( "gErrorIgnoreLevel = 6001;");

	printf("----------------------------\npT range:\t%.3f - %.3f\neta range:\t%.3f - %.3f\n----------------------------\n",minPt.Get(),maxPt.Get(),-absEtaMax.Get(),absEtaMax.Get());

	AliJFFlucAnalysis::BINNING binning = AliJFFlucAnalysis::BINNING_CENT_PbPb;

	AliJFFlucAnalysis *pfa = new AliJFFlucAnalysis("bayesian-hydro-jfluc");
	pfa->SetBinning(binning);

	bool bSC = true;
	bool bAC = true;
	bool bSPC = true;
	bool bVNPT = true;

	AliAnalysisAnaTwoMultiCorrelations *pca_SC = new AliAnalysisAnaTwoMultiCorrelations("TwoMultiCorrelations_SC",kFALSE);
	AliAnalysisAnaTwoMultiCorrelations *pca_AC = new AliAnalysisAnaTwoMultiCorrelations("TwoMultiCorrelations_AC",kFALSE);

	AliAnalysisPtVn *pca_PtVn = new AliAnalysisPtVn("PtVnCorr");
	AliAnalysisSPCRun2 *pspc = new AliAnalysisSPCRun2("SPC");
	AliAnalysisMeanPt *pmpt = new AliAnalysisSPCRun2("meanPt");

	if (bSC) {
		cout << "Computing SC" << endl;
		// Compute the higher order SC (Anna's analysis)
		pca_SC->SetJWeights(kFALSE,0);
		// SetHarmoArray not used anymore in the latest version of the analysis script.
		//pca->SetHarmoArray("2 3 4 2 3 5 2 4 6 3 4 5"); //<-- BA
		//pca->SetHarmoArray("2 3 0 2 4 0 3 2 0 4 2 0"); //<-- Cindy
		pca_SC->SetMinMultiplicity(10);
		pca_SC->SetPtRange(0.2, 5.0);
		pca_SC->SetObservable(true, false);
	}
	if (bAC) {
		cout << "Computing AC" << endl;
		// Compute the AC for the lower order combinations of harmonics (Cindy's analysis)
		pca_AC->SetJWeights(kFALSE,0);
		pca_AC->SetMinMultiplicity(10);
		pca_AC->SetPtRange(0.2, 5.0);
		pca_AC->SetObservable(false, true);
	}
	
	if (bVNPT) cout << "Computing v_n p_T correlations" << endl;

	if (bSPC) cout << "Computing SPC" << endl;


	TFile *pfo = new TFile(outfilen.Get().c_str(),"recreate");
	TDirectory *phydroDir = pfo->mkdir("hydro");
	pfo->cd("hydro");

	// Centrality
	const double vertex[3] = {0.0,0.0,0.0};
	const double pcent_vs_mult[] = {2348.0,1928.0,1843.7,1632.2,1342.9,917.3,615.3,396.1,238.7,129.9};

	pfa->AddFlags(AliJFFlucAnalysis::FLUC_EBE_WEIGHTING);
	pfa->SetEtaRange(0.0,0.8);
	pfa->SetEventVertex(vertex);
	pfa->UserCreateOutputObjects();

	pfo->cd();
	if (bSC) {
		phydroDir->mkdir("TwoMultiCorr_SC");
		phydroDir->cd("TwoMultiCorr_SC");
		pca_SC->UserCreateOutputObjects();
	}
	if (bAC) {
		phydroDir->mkdir("TwoMultiCorr_AC");
		phydroDir->cd("TwoMultiCorr_AC");
		pca_AC->UserCreateOutputObjects();
	}
	if (bVNPT) {
		phydroDir->mkdir("vnptCorr");
		phydroDir->cd("vnptCorr");
		pca_PtVn->UserCreateOutputObjects();
	}
	if (bSPC) {
		phydroDir->mkdir("SPC");
		phydroDir->cd("SPC");
		pspc->UserCreateOutputObjects();
	}

	TClonesArray *pinputList = new TClonesArray("AliJBaseTrack",2500);

	// RNG for particles
	TRandom3 *rng = new TRandom3(21);



	for(uint iEvt = 0; iEvt<1000; iEvt++){
		int nTracks = static_cast<int>(rng->Gaus(1000, 600));
		if (nTracks<1) continue;

		double cent = 0.0;
		for(uint ci = 0; ci < sizeof(pcent_vs_mult)/sizeof(pcent_vs_mult[0])-1; ++ci){
			if(nTracks < pcent_vs_mult[ci+1])
				continue;
			cent = std::max((pcent_vs_mult[ci]-nTracks)/(pcent_vs_mult[ci]-pcent_vs_mult[ci+1])*(AliJFFlucAnalysis::CentBin_PbPb_default[ci][1]-AliJFFlucAnalysis::CentBin_PbPb_default[ci][0])+AliJFFlucAnalysis::CentBin_PbPb_default[ci][0],0.0);
			break;
		}

		// Get centrality bin
		int fCBin = AliJFFlucAnalysis::GetBin(cent,binning);
		if(fCBin == -1) continue;
		uint nch05 = 0;
		// Loop over tracks
		for(uint iTrk = 0; iTrk < nTracks; ++iTrk){

			pinputList->Clear();
			double eta = rng->Uniform(-1., 1.);
			double phi = rng->Uniform(0, TMath::TwoPi());
			double pt = rng->PoissonD(2.);


			if (TMath::Abs(eta)<0.5) nch05++;

			pfa->fh_eta_NoPtCut[fCBin]->Fill(eta,1.0);

			//XXX pT cut 0.2 - 5.0 applied in the modified JFluc analysis!

			bool rejectTrack = false;
			if(pt < minPt.Get() || pt > maxPt.Get() || eta < -absEtaMax.Get() || eta > absEtaMax.Get() || rejectTrack)
				continue;

			double px = pt*cos(phi);
			double py = pt*sin(phi);
			double pz = pt*sinh(eta);
			new ((*pinputList)[pinputList->GetEntriesFast()]) AliJBaseTrack(px, py, pz, 1, iTrk, 0, -1);

		}
		
		pfa->fh_mult_eta05[fCBin]->Fill((double)nch05,1.0);

		// Analysis inputs
		pfa->Init(0, "");
		pfa->SetInputList(pinputList);
		pfa->SetEventCentrality(cent);
		pfa->SetEventImpactParameter(-1);

		pfa->SetPhiWeights((TH1*)0);
		pfa->SetPhiWeights((TF3*)0);

		pfa->UserExec("");

		// Optional analysis setup (SPC, NSC, AC, VNPT)
		if (bSC) {
			pca_SC->SetInputList(pinputList);
			pca_SC->SetEventCentrality(cent);
	                                pca_SC->UserExec("");
		}
		if (bAC) {
			pca_AC->SetInputList(pinputList);
			pca_AC->SetEventCentrality(cent);
			pca_AC->UserExec("");
		}
		if (bVNPT) {
			phydroDir->cd("vnptCorr");
			pca_PtVn->SetInputList(pinputList);
			pca_PtVn->SetEventCentrality(cent);
			pca_PtVn->UserExec("");
		}
		if (bSPC) {
			phydroDir->cd("SPC");
			pspc->SetInputList(pinputList);
			pspc->SetEventCentrality(cent);
			pspc->SetSPC(0);
			pspc->UserExec("");
		}
	}

	printf("analysis DONE!\n");

	if (bVNPT) {
		phydroDir->cd("vnptCorr");
		pca_PtVn->WriteLists("OutputPtVnCorr");
	}
	if (bSPC) {
		phydroDir->cd(Form("SPC"));
		pspc->WriteLists(Form("OutputListSPC"));
	}
	
	pfo->Write();
	pfo->Close();
	delete pfo;

	delete pinputList;

	delete pfa;

	return 0;
}

