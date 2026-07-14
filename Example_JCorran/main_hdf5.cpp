#include <stdio.h>
#include <vector>
#include <exception>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>
#include <iostream>
#include <string>
#include <cstdlib>

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

#include <H5Cpp.h>
using namespace H5;

#include "args.hxx"

#include "AliJBaseTrack.h"
#include "AliJFFlucAnalysisTProfile.h"
#include "AliAnalysisAnaTwoMultiCorrelations.h"
#include "AliAnalysisSPCRun2.h"
#include "AliAnalysisPtVn.h"
#include "AliAnalysisV02.h"

typedef unsigned int uint;

#define NC 8

namespace PidUtil {
	inline bool isChargedPion(int64_t pid) {
		return std::llabs(pid) == 211;
	}
	inline bool isChargedKaon(int64_t pid) {
		return std::llabs(pid) == 321;
	}
	inline bool isProton(int64_t pid) {
		return std::llabs(pid) == 2212;
	}
}

void InitialiseDataTypes(CompType &particleDataType);
int InitialiseCentrality(std::string systemStr, int paramID);
int InitialiseForwardCentrality(std::string systemStr, int paramID);
void InitialiseAnalyses(bool forwardCent, double etaMin, double etaMax, double etaMaxVnPt);
void FillAndRunAnalyses(int fCBin, int fCBin_forward, float cent, float cent_forward, bool bCent, bool forwardCent, bool bSPC, TClonesArray *pinputList);

struct ParticleType{
	int64_t sample;
	int64_t pid;
	int64_t charge;
	double pT;
	double ET;
	double mT;
	double phi;
	double y;
	double eta;
};

const double cent_vs_mult2760[] = {1835.7,1528.4,1461.2,1292.3,1064.2,726.2,486.2,312.3,187.0,100.6};
const double cent_vs_mult5020[] = {2348.0,1928.0,1843.7,1632.2,1342.9,917.3,615.3,396.1,238.7,129.9};
const double cent_vs_mult_XeXe5440[] = {1700.5,1300.8,1238.3,1096.4,904.2,616.3,412.7,266.3,161.7,90.0};
const double cent_vs_mult_pPb5020[] = {77.7,44.5,41.5,36.9,32.8,27.7,23.7,19.9,16.1,12.3};
const double cent_vs_mult_AuAu200[] = {1019.7,683.3,610.1,481.5,362.7,229.4,147.7,93.4,55.1,29.0};
const double *pcent_vs_mult;
const double *pcent_forward;
static double centBins[13];
static double centBinsForward[13];

bool bV02;
bool bSPC;
bool bCent;
bool bForward;

uint dNch_deta_PID[4];
uint dNch_dy_PID[4];
AliJFFlucAnalysisTProfile::BINNING binning;
AliJFFlucAnalysisTProfile *pfa;
AliJFFlucAnalysisTProfile *pfa_forward;
AliAnalysisSPCRun2 *pspc;
AliAnalysisSPCRun2 *pspc_forward;
AliAnalysisPtVn *pca_PtVn;
AliAnalysisPtVn *pca_PtVn_forward;
AliAnalysisV02 *pca_V02;
TDirectory *phydroDir;
TFile *pfo;

int main(int argc, char **pargv){
	args::ArgumentParser parser("JFLUC HIC analyzer (HDF5 input)","");
	args::Group group_io(parser,"Output files",args::Group::Validators::DontCare);
	args::ValueFlag<std::string> outfilen(group_io,"dst","Output root file name",{'o',"outfile"},"AnalysisResults.root");
	args::PositionalList<std::string> inputFiles(group_io,"input.hdf","HDF5 particle input file list");

	args::Group group_ana(parser,"Analysis",args::Group::Validators::DontCare);
	args::ValueFlag<std::string> binningMode(group_ana,"mode","Set JFFlucAnalysis binning mode. Accepted modes: PbPb_cent (default), OO_cent, PbPb_mult_1, pPb_mult_1",{'b',"binning"},"PbPb_cent");
	args::ValueFlag<std::string> system(group_ana,"system","Collision system or centrality CSV: PbPb2760, PbPb5020 (default), pPb5020, XeXe5440, AuAu200, or path to CSV",{'s',"system"},"PbPb5020");

	args::ValueFlag<double> minPt(group_ana,"pTmin","Minimum track pT",{"pTmin"}, 0.2);
	args::ValueFlag<double> maxPt(group_ana,"pTmax","Maximum track pT",{"pTmax"}, 5.0);
	args::ValueFlag<double> absEtaMin(group_ana,"absEtaMin","Minimum |eta| for eta gap",{"absEtaMin"}, 0.4);
	args::ValueFlag<double> absEtaMax(group_ana,"absEtaMax","Maximum track |eta|",{"absEtaMax"}, 0.8);
	args::ValueFlag<double> absEtaMaxVnPt(group_ana,"absEtaMaxVnPt","Maximum track |eta| for vn-pt calulations",{"absEtaMaxVnPt"}, 0.8);

	args::ValueFlag<bool> forwardCent(group_ana,"forwardCent","Use additional forward centrality.",{"forwardCent"}, false);
	args::ValueFlag<std::string> forwardCentBinning(group_ana,"forwardCentBinning","Binning file for the additional forward centrality.",{"forwardCentBinning"}, "");

	args::ValueFlag<int> param(group_ana,"param","Parametrization ID (line in centrality CSV)",{"param"},-1);
	args::ValueFlag<std::string> addObs(group_ana, "addObs", "Additional observables: spc, v02", {"addObs"}, "none");

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

	if(inputFiles.end()-inputFiles.begin() <= 0){
		printf("No input files: exiting.\n");
		return 1;
	}

	gROOT->ProcessLine( "gErrorIgnoreLevel = 6001;");

	printf("----------------------------\npT range:\t%.3f - %.3f\neta range:\t%.3f - %.3f\n----------------------------\n",minPt.Get(),maxPt.Get(),-absEtaMax.Get(),absEtaMax.Get());

	Bool_t useOO = kFALSE;

	if(binningMode.Get().compare("PbPb_cent") == 0)
		binning = AliJFFlucAnalysisTProfile::BINNING_CENT_PbPb;
	else if(binningMode.Get().compare("OO_cent") == 0) {
		useOO = kTRUE;
		binning = AliJFFlucAnalysisTProfile::BINNING_CENT_OO;
	} else if(binningMode.Get().compare("PbPb_mult_1") == 0)
		binning = AliJFFlucAnalysisTProfile::BINNING_MULT_PbPb_1;
	else if(binningMode.Get().compare("pPb_mult_1") == 0)
		binning = AliJFFlucAnalysisTProfile::BINNING_MULT_pPb_1;
	else binning = AliJFFlucAnalysisTProfile::BINNING_CENT_PbPb;

	bSPC = addObs.Get().find("spc") != std::string::npos;
	bV02 = addObs.Get().find("v02") != std::string::npos;

	pfo = new TFile(outfilen.Get().c_str(),"recreate");
	phydroDir = pfo->mkdir("hydro");
	pfo->cd("hydro");

	if (InitialiseCentrality(system.Get(), param.Get())) return 1;
	if (forwardCent.Get()) {
		if (InitialiseForwardCentrality(forwardCentBinning.Get(), param.Get())) return 1;
	}

	InitialiseAnalyses(forwardCent.Get(), absEtaMin.Get(), absEtaMax.Get(), absEtaMaxVnPt.Get());

	TClonesArray *pinputList = new TClonesArray("AliJBaseTrack",2500);

	CompType particleDataType(sizeof(ParticleType));
	InitialiseDataTypes(particleDataType);

	for(auto &mf : args::get(inputFiles)){
		H5File *pfile;
		try{
			pfile = new H5File(mf,H5F_ACC_RDONLY);
		}catch(H5::FileIException &e){
			printf("Unable to open file %s\n",mf.c_str());
			continue;
		}
		Group rootGroup(pfile->openGroup("/"));

		printf("Looping over %lld events...\n", rootGroup.getNumObjs());
		for(uint gi = 0; gi < rootGroup.getNumObjs(); ++gi){
			bCent = true;
			bForward = true;
			H5std_string objName = rootGroup.getObjnameByIdx(gi);
			DataSet dataSet = pfile->openDataSet(objName);

			DataSpace dataSpace = dataSet.getSpace();
			size_t dsSize = dataSpace.getSimpleExtentNpoints();
			std::vector<ParticleType> hadrons(dsSize);
			dataSet.read(hadrons.data(),particleDataType);
			std::sort(hadrons.begin(),hadrons.end(),[](ParticleType &a, ParticleType &b)->bool{
				return a.sample < b.sample;
			});

			Attribute attr;
			double dNch_deta;
			double dNch_deta_forward;
			try{
				attr = dataSet.openAttribute("dNch_deta");
				attr.read(attr.getDataType(),&dNch_deta);
				if (forwardCent.Get()) {
					attr = dataSet.openAttribute("dNch_deta_forward");
					attr.read(attr.getDataType(),&dNch_deta_forward);
				}
			}catch(H5::AttributeIException &e){
				printf("No dNch_deta attribute: empty event\n");
				continue;
			}

			double cent = -1.0;
			for(uint ci = 0; ci < sizeof(centBins)/sizeof(centBins[0])-1; ++ci){
				if(dNch_deta < pcent_vs_mult[ci+1])
					continue;
				if (useOO) {
					if (ci>11) continue;
					cent = std::max((pcent_vs_mult[ci]-dNch_deta)/(pcent_vs_mult[ci]-pcent_vs_mult[ci+1])*(AliJFFlucAnalysisTProfile::CentBin_OO_central[ci][1]-AliJFFlucAnalysisTProfile::CentBin_OO_central[ci][0])+AliJFFlucAnalysisTProfile::CentBin_OO_central[ci][0],0.0);
				} else {
					if (ci>8) continue;
					cent = std::max((pcent_vs_mult[ci]-dNch_deta)/(pcent_vs_mult[ci]-pcent_vs_mult[ci+1])*(AliJFFlucAnalysisTProfile::CentBin_PbPb_default[ci][1]-AliJFFlucAnalysisTProfile::CentBin_PbPb_default[ci][0])+AliJFFlucAnalysisTProfile::CentBin_PbPb_default[ci][0],0.0);
				}
				break;
			}

			int fCBin = AliJFFlucAnalysisTProfile::GetBin(cent,binning);
			printf("dNch/deta: %.4f, Centrality: %.2f, CentBin: %d\n", dNch_deta, cent, fCBin);
			if(fCBin == -1) bCent = false;

			int fCBin_forward = -1;
			double cent_forward = -1.0;
			if (forwardCent.Get()) {
				for(uint ci = 0; ci < sizeof(centBins)/sizeof(centBins[0])-1; ++ci){
					if(dNch_deta_forward < pcent_forward[ci+1])
						continue;
					if (useOO) {
						if (ci>11) continue;
						cent_forward = std::max((pcent_forward[ci]-dNch_deta_forward)/(pcent_forward[ci]-pcent_forward[ci+1])*(AliJFFlucAnalysisTProfile::CentBin_OO_central[ci][1]-AliJFFlucAnalysisTProfile::CentBin_OO_central[ci][0])+AliJFFlucAnalysisTProfile::CentBin_OO_central[ci][0],0.0);
					} else {
						if(ci>8) continue;
						cent_forward = std::max((pcent_forward[ci]-dNch_deta_forward)/(pcent_forward[ci]-pcent_forward[ci+1])*(AliJFFlucAnalysisTProfile::CentBin_PbPb_default[ci][1]-AliJFFlucAnalysisTProfile::CentBin_PbPb_default[ci][0])+AliJFFlucAnalysisTProfile::CentBin_PbPb_default[ci][0],0.0);
					}
					break;
				}

				fCBin_forward = AliJFFlucAnalysisTProfile::GetBin(cent_forward,binning);
				if(fCBin_forward == -1) bForward = false;
			} else {
				bForward = false;
			}

			if (!bCent && !bForward) {
				continue;
			}

			double e2, e3;
			attr = dataSet.openAttribute("e2");
			attr.read(attr.getDataType(),&e2);
			attr = dataSet.openAttribute("e3");
			attr.read(attr.getDataType(),&e3);
			if (bCent) {
				pfa->fh_ecc[2][fCBin]->Fill(e2,1.0);
				pfa->fh_ecc[3][fCBin]->Fill(e3,1.0);
			}
			if (bForward) {
				pfa_forward->fh_ecc[2][fCBin_forward]->Fill(e2,1.0);
				pfa_forward->fh_ecc[3][fCBin_forward]->Fill(e3,1.0);
			}

			for(uint i = 0, cs = hadrons[i].sample; i < hadrons.size(); cs = hadrons[i].sample){
				pinputList->Clear();
				for (uint ij=0; ij<4; ij++) {
					dNch_deta_PID[ij] = 0;
					dNch_dy_PID[ij] = 0;
				}
				for(uint trackIndex = 0; i < hadrons.size(); ++i){
					if(hadrons[i].sample != cs)
						break;

					if(hadrons[i].charge == 0)
						continue;

					int64_t pid = hadrons[i].pid;
					double phi = hadrons[i].phi;
					double pt = hadrons[i].pT;
					double eta = hadrons[i].eta;
					double yrap = hadrons[i].y;

					bool eta05 = TMath::Abs(eta) <= 0.5;
					bool rap05 = TMath::Abs(yrap) <= 0.5;
					if (eta05) ++dNch_deta_PID[0];
					if (rap05) ++dNch_dy_PID[0];

					if(PidUtil::isChargedKaon(pid)){
						if (bCent) {
							if (rap05) pfa->fh_pt_Kaon[fCBin]->Fill(pt,1.0);
							pfa->fh_eta_Kaon[fCBin]->Fill(eta,1.0);
						}
						if (bForward) {
							if (rap05) pfa_forward->fh_pt_Kaon[fCBin_forward]->Fill(pt,1.0);
							pfa_forward->fh_eta_Kaon[fCBin_forward]->Fill(eta,1.0);
						}
						if (eta05) ++dNch_deta_PID[1];
						if (rap05) ++dNch_dy_PID[1];
					}else if(PidUtil::isChargedPion(pid)){
						if (bCent) {
							if (rap05) pfa->fh_pt_Pion[fCBin]->Fill(pt,1.0);
							pfa->fh_eta_Pion[fCBin]->Fill(eta,1.0);
						}
						if (bForward) {
							if (rap05) pfa_forward->fh_pt_Pion[fCBin_forward]->Fill(pt,1.0);
							pfa_forward->fh_eta_Pion[fCBin_forward]->Fill(eta,1.0);
						}
						if (eta05) ++dNch_deta_PID[2];
						if (rap05) ++dNch_dy_PID[2];
					}else if(PidUtil::isProton(pid)){
						if (bCent) {
							if (rap05) pfa->fh_pt_Proton[fCBin]->Fill(pt,1.0);
							pfa->fh_eta_Proton[fCBin]->Fill(eta,1.0);
						}
						if (bForward) {
							if (rap05) pfa_forward->fh_pt_Proton[fCBin_forward]->Fill(pt,1.0);
							pfa_forward->fh_eta_Proton[fCBin_forward]->Fill(eta,1.0);
						}
						if (eta05) ++dNch_deta_PID[3];
						if (rap05) ++dNch_dy_PID[3];
					}

					if (bCent) pfa->fh_eta_NoPtCut[fCBin]->Fill(eta,1.0);
					if (bForward) pfa_forward->fh_eta_NoPtCut[fCBin_forward]->Fill(eta,1.0);

					if(pt < minPt.Get() || pt > maxPt.Get() || eta < -absEtaMax.Get() || eta > absEtaMax.Get())
						continue;

					double px = pt*cos(phi);
					double py = pt*sin(phi);
					double pz = pt*sinh(eta);
					new ((*pinputList)[pinputList->GetEntriesFast()]) AliJBaseTrack(px,py,pz,1,trackIndex++,pid,hadrons[i].charge);
				}

				FillAndRunAnalyses(fCBin, fCBin_forward, cent, cent_forward, bCent, bForward, bSPC, pinputList);
			}

			if(!hadrons.empty())
				printf("scanned %lu particles in %lld samples (%u/%llu).\n",hadrons.size(),(long long)hadrons.back().sample,gi+1,rootGroup.getNumObjs());
		}

		printf("analysis DONE: %s\n",mf.c_str());
		delete pfile;
	}

	phydroDir->cd("vnptCorr");
	pca_PtVn->WriteLists(Form("OutputPtVnCorr"));
	if (forwardCent.Get()) {
		phydroDir->cd("vnptCorr");
		pca_PtVn_forward->WriteLists(Form("OutputPtVnCorr_forward"));
	}
	if (bSPC) {
		phydroDir->cd(Form("SPC"));
		pspc->WriteLists(Form("OutputListSPC"));
		if (forwardCent.Get()) pspc_forward->WriteLists(Form("OutputListSPC_forward"));
	}
	if (bV02) {
		phydroDir->cd(Form("V02"));
		pca_V02->WriteLists(Form("OutputV02"));
	}
	pfo->Write();
	pfo->Close();
	delete pfo;
	if (bSPC) delete pspc;
	delete pinputList;

	return 0;
}

void InitialiseDataTypes(CompType &particleDataType) {
	particleDataType.insertMember("sample",HOFFSET(ParticleType,sample),PredType::NATIVE_INT64);
	particleDataType.insertMember("ID",HOFFSET(ParticleType,pid),PredType::NATIVE_INT64);
	particleDataType.insertMember("charge",HOFFSET(ParticleType,charge),PredType::NATIVE_INT64);
	particleDataType.insertMember("pT",HOFFSET(ParticleType,pT),PredType::NATIVE_DOUBLE);
	particleDataType.insertMember("ET",HOFFSET(ParticleType,ET),PredType::NATIVE_DOUBLE);
	particleDataType.insertMember("mT",HOFFSET(ParticleType,mT),PredType::NATIVE_DOUBLE);
	particleDataType.insertMember("phi",HOFFSET(ParticleType,phi),PredType::NATIVE_DOUBLE);
	particleDataType.insertMember("y",HOFFSET(ParticleType,y),PredType::NATIVE_DOUBLE);
	particleDataType.insertMember("eta",HOFFSET(ParticleType,eta),PredType::NATIVE_DOUBLE);
}

int InitialiseCentrality(std::string systemStr, int paramID) {
	if(systemStr.compare("PbPb2760") == 0){
		printf("Using 2.76 TeV centrality binning.\n");
		pcent_vs_mult = cent_vs_mult2760;
	}else if(systemStr.compare("PbPb5020") == 0){
		printf("Using 5.02 TeV centrality binning.\n");
		pcent_vs_mult = cent_vs_mult5020;
	}else if(systemStr.compare("XeXe5440") == 0){
		printf("Using XeXe 5.44 TeV centrality binning.\n");
		pcent_vs_mult = cent_vs_mult_XeXe5440;
	}else if(systemStr.compare("AuAu200") == 0){
		printf("Using AuAu 200 GeV centrality binning.\n");
		pcent_vs_mult = cent_vs_mult_AuAu200;
	}else if(systemStr.compare("pPb5020") == 0){
		printf("Using pPb 5.02 TeV centrality binning.\n");
		pcent_vs_mult = cent_vs_mult_pPb5020;
	}else{
		std::ifstream file(systemStr.c_str());
		if(!file.is_open()){
			printf("Unable to read param %d from %s.\n",paramID,systemStr.data());
			return 1;
		}
		if(paramID < 0){
			printf("Invalid parameter id %d\n",paramID);
			return 1;
		}
		printf("Using centralities from file: %s\n", systemStr.data());
		file.seekg(std::ios::beg);
		for(int i = 0; i < paramID; ++i)
			file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
		std::string s;
		file >> s;
		file.close();
		printf("%s\n", s.data());
		pcent_vs_mult = centBins;

		std::stringstream ss(s);
		std::string s1;
		for(uint i = 0; std::getline(ss,s1,','); i++)
			centBins[i] = std::stod(s1);
	}
	return 0;
}

int InitialiseForwardCentrality(std::string systemStr, int paramID) {
	std::ifstream file(systemStr.c_str());
	if(!file.is_open()){
		printf("Unable to read param %d from %s.\n",paramID,systemStr.data());
		return 1;
	}
	if(paramID < 0){
		printf("Invalid parameter id %d\n",paramID);
		return 1;
	}
	printf("Using forward centralities from file: %s\n", systemStr.data());
	file.seekg(std::ios::beg);
	for(int i = 0; i < paramID; ++i)
		file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	std::string s;
	file >> s;
	file.close();
	printf("%s\n", s.data());
	pcent_forward = centBinsForward;

	std::stringstream ss(s);
	std::string s1;
	for(uint i = 0; std::getline(ss,s1,','); i++)
		centBinsForward[i] = std::stod(s1);
	return 0;
}

void InitialiseAnalyses(bool forwardCent, double etaMin, double etaMax, double etaMaxVnPt) {
	pfa = new AliJFFlucAnalysisTProfile("bayesian-hydro-jfluc");
	pfa->SetBinning(binning);
	pfa->AddFlags(AliJFFlucAnalysisTProfile::FLUC_EBE_WEIGHTING);
	pfa->SetEtaRange(etaMin, etaMax);
	pfa->UserCreateOutputObjects();

	pfo->cd();
	phydroDir->mkdir("vnptCorr");

	pca_PtVn = new AliAnalysisPtVn("PtVnCorr");
	pca_PtVn->UserCreateOutputObjects();
	pca_PtVn->SetPtSubRange(etaMin, etaMaxVnPt);

	if (forwardCent) {
		pfa_forward = new AliJFFlucAnalysisTProfile("jfluc-forward");
		pfa_forward->SetBinning(binning);
		pfa_forward->AddFlags(AliJFFlucAnalysisTProfile::FLUC_EBE_WEIGHTING);
		pfa_forward->SetEtaRange(etaMin, etaMax);
		pfo->cd("hydro");
		pfa_forward->UserCreateOutputObjects("jfluc-forward");

		phydroDir->cd("vnptCorr");
		pca_PtVn_forward = new AliAnalysisPtVn("PtVnCorr_forward");
		pca_PtVn_forward->UserCreateOutputObjects();
		pca_PtVn_forward->SetPtSubRange(etaMin, etaMaxVnPt);
	}

	if (bSPC) {
		cout << "Computing SPC" << endl;
		pspc = new AliAnalysisSPCRun2("SPC_0");
		phydroDir->mkdir("SPC");
		pspc->UserCreateOutputObjects();

		if (forwardCent) {
			pspc_forward = new AliAnalysisSPCRun2("SPC_forward_0");
			pspc_forward->UserCreateOutputObjects();
		}
	}
	if (bV02) {
		cout << "Computing V02" << endl;
		pca_V02 = new AliAnalysisV02("V02_analysis");
		phydroDir->mkdir("V02");
		pca_V02->SetUseFlowWeights(kFALSE);
		pca_V02->UserCreateOutputObjects();
	}
}

void FillAndRunAnalyses(int fCBin, int fCBin_forward, float cent, float cent_forward, bool bCent, bool forwardCent, bool bSPC, TClonesArray *pinputList) {
	if (bCent) {
		pfa->fh_mult_eta05[fCBin]->Fill((double)dNch_deta_PID[0],1.0);
		pfa->fh_mult_eta05_Kaon[fCBin]->Fill((double)dNch_deta_PID[1],1.0);
		pfa->fh_mult_eta05_Pion[fCBin]->Fill((double)dNch_deta_PID[2],1.0);
		pfa->fh_mult_eta05_Proton[fCBin]->Fill((double)dNch_deta_PID[3],1.0);

		pfa->fh_mult_rap05[fCBin]->Fill((double)dNch_dy_PID[0],1.0);
		pfa->fh_mult_rap05_Kaon[fCBin]->Fill((double)dNch_dy_PID[1],1.0);
		pfa->fh_mult_rap05_Pion[fCBin]->Fill((double)dNch_dy_PID[2],1.0);
		pfa->fh_mult_rap05_Proton[fCBin]->Fill((double)dNch_dy_PID[3],1.0);

		pfa->SetInputList(pinputList);
		pfa->SetEventCentrality(cent);
		pfa->SetEventImpactParameter(-1);
		pfa->SetPhiWeights((TH1*)0);
		pfa->SetPhiWeights((TF3*)0);
		pfa->UserExec("");

		phydroDir->cd("vnptCorr");
		pca_PtVn->SetInputList(pinputList);
		pca_PtVn->SetEventCentrality(cent);
		pca_PtVn->UserExec("");
	}

	if (forwardCent) {
		pfa_forward->fh_mult_eta05[fCBin_forward]->Fill((double)dNch_deta_PID[0],1.0);
		pfa_forward->fh_mult_eta05_Kaon[fCBin_forward]->Fill((double)dNch_deta_PID[1],1.0);
		pfa_forward->fh_mult_eta05_Pion[fCBin_forward]->Fill((double)dNch_deta_PID[2],1.0);
		pfa_forward->fh_mult_eta05_Proton[fCBin_forward]->Fill((double)dNch_deta_PID[3],1.0);

		pfa_forward->fh_mult_rap05[fCBin_forward]->Fill((double)dNch_dy_PID[0],1.0);
		pfa_forward->fh_mult_rap05_Kaon[fCBin_forward]->Fill((double)dNch_dy_PID[1],1.0);
		pfa_forward->fh_mult_rap05_Pion[fCBin_forward]->Fill((double)dNch_dy_PID[2],1.0);
		pfa_forward->fh_mult_rap05_Proton[fCBin_forward]->Fill((double)dNch_dy_PID[3],1.0);

		pfa_forward->SetInputList(pinputList);
		pfa_forward->SetEventCentrality(cent);
		pfa_forward->SetEventImpactParameter(-1);
		pfa_forward->SetPhiWeights((TH1*)0);
		pfa_forward->SetPhiWeights((TF3*)0);
		pfa_forward->UserExec("");

		phydroDir->cd("vnptCorr");
		pca_PtVn_forward->SetInputList(pinputList);
		pca_PtVn_forward->SetEventCentrality(cent);
		pca_PtVn_forward->UserExec("");
	}

	if (bSPC && bCent) {
		phydroDir->cd("SPC");
		pspc->SetInputList(pinputList);
		pspc->SetEventCentrality(cent);
		pspc->SetSPC(0);
		pspc->UserExec("");
	}
	if (bV02 && bCent) {
		phydroDir->cd("V02");
		pca_V02->SetInputList(pinputList);
		pca_V02->SetEventCentrality(cent);
		pca_V02->UserExec("");
	}
	if (bSPC && forwardCent) {
		phydroDir->cd("SPC");
		pspc_forward->SetInputList(pinputList);
		pspc_forward->SetEventCentrality(cent);
		pspc_forward->SetSPC(0);
		pspc_forward->UserExec("");
	}
}
