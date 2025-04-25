#include "AliJAnalysisMultiplicityMeanPt.h"
#include <TFile.h>
#include <TDirectory.h>
#include <TMath.h>
#include <iostream>

ClassImp(AliJAnalysisMultiplicityMeanPt);

AliJAnalysisMultiplicityMeanPt::AliJAnalysisMultiplicityMeanPt(const char* name) : 
    fName(name)
{
    // Default centrality bins: 0-5-10-20-30-40-50-60-70-80
    fCentBins = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};
    
    // Initialize PID name mapping
    fPidNames[kCharged] = "Charged";
    fPidNames[kPion] = "Pion";
    fPidNames[kKaon] = "Kaon";
    fPidNames[kProton] = "Proton";
    
    InitializeMaps();
}

AliJAnalysisMultiplicityMeanPt::~AliJAnalysisMultiplicityMeanPt() {
    // Delete all created graphs
    for (auto& pair : fMultGraphs) {
        if (pair.second) delete pair.second;
    }
    for (auto& pair : fMultRapGraphs) {
        if (pair.second) delete pair.second;
    }
    for (auto& pair : fMeanPtGraphs) {
        if (pair.second) delete pair.second;
    }
}

void AliJAnalysisMultiplicityMeanPt::InitializeMaps() {
    // Initialize vectors for each centrality bin
    int nCentBins = fCentBins.size() - 1;
    fEventCounts.resize(nCentBins, 0);
    
    // Initialize for each particle type
    for (const auto& pidPair : fPidNames) {
        int pid = pidPair.first;
        fParticleCounts[pid] = std::vector<int>(nCentBins, 0);
        fMeanPtSums[pid] = std::vector<double>(nCentBins, 0.0);
        fParticleCountsRap[pid] = std::vector<int>(nCentBins, 0);
    }
}

void AliJAnalysisMultiplicityMeanPt::SetSystemParameters(const char* collisionSystem) {
    // Set centrality bins based on collision system
    if (strcmp(collisionSystem, "PbPb2760") == 0) {
        // PbPb 2.76 TeV centrality bins
        fCentBins = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90};
    } else if (strcmp(collisionSystem, "PbPb5020") == 0) {
        // PbPb 5.02 TeV centrality bins
        fCentBins = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90};
    } else if (strcmp(collisionSystem, "XeXe5440") == 0) {
        // XeXe 5.44 TeV centrality bins
        fCentBins = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90};
    } else if (strcmp(collisionSystem, "pPb5020") == 0) {
        // pPb 5.02 TeV centrality bins
        fCentBins = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90};
    } else if (strcmp(collisionSystem, "AuAu200") == 0) {
        // AuAu 200 GeV centrality bins
        fCentBins = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90};
    } else {
        std::cout << "Unknown collision system: " << collisionSystem << ". Using default centrality bins." << std::endl;
    }
    
    // Reinitialize maps with new binning
    InitializeMaps();
}

void AliJAnalysisMultiplicityMeanPt::SetCentralityBins(const std::vector<double>& centBins) {
    fCentBins = centBins;
    InitializeMaps();
}

int AliJAnalysisMultiplicityMeanPt::GetCentBin(double centrality) {
    // Find the centrality bin index
    for (size_t i = 0; i < fCentBins.size() - 1; i++) {
        if (centrality >= fCentBins[i] && centrality < fCentBins[i+1]) {
            return i;
        }
    }
    return -1; // Out of range
}

void AliJAnalysisMultiplicityMeanPt::ProcessEvent(double centrality, const std::vector<double>& pt_values, 
                                                  const std::vector<int>& pid_codes) {
    int centBin = GetCentBin(centrality);
    if (centBin < 0) return; // Skip if centrality is out of range
    
    // Count this event for the centrality bin
    fEventCounts[centBin]++;
    
    // Track counts for each particle type
    std::map<int, int> trackCounts;          // Full eta range
    std::map<int, int> trackCountsRap;       // Mid-rapidity
    std::map<int, double> ptSums;            // For mean pT calculation
    
    // Initialize counters
    for (const auto& pidPair : fPidNames) {
        trackCounts[pidPair.first] = 0;
        trackCountsRap[pidPair.first] = 0;
        ptSums[pidPair.first] = 0.0;
    }
    
    // Count all charged particles
    trackCounts[kCharged] = pt_values.size();
    
    // Process all tracks
    for (size_t i = 0; i < pt_values.size(); i++) {
        double pt = pt_values[i];
        int pid = pid_codes[i];
        
        // Sum pT for all charged particles
        ptSums[kCharged] += pt;
        
        // Count identified particles and sum their pT
        if (abs(pid) == 211) {  // Pions
            trackCounts[kPion]++;
            ptSums[kPion] += pt;
            trackCountsRap[kPion]++; // Simplified - in real code would check rapidity
        } 
        else if (abs(pid) == 321) {  // Kaons
            trackCounts[kKaon]++;
            ptSums[kKaon] += pt;
            trackCountsRap[kKaon]++; // Simplified
        } 
        else if (abs(pid) == 2212) {  // Protons
            trackCounts[kProton]++;
            ptSums[kProton] += pt;
            trackCountsRap[kProton]++; // Simplified
        }
    }
    
    // Charged particles in mid-rapidity - simplified, would check eta in real code
    for (const auto& pidPair : fPidNames) {
        int pid = pidPair.first;
        fParticleCounts[pid][centBin] += trackCounts[pid];
        fParticleCountsRap[pid][centBin] += trackCountsRap[pid];
        fMeanPtSums[pid][centBin] += ptSums[pid];
    }
}

void AliJAnalysisMultiplicityMeanPt::FinalizeMeasurements() {
    int nCentBins = fCentBins.size() - 1;
    
    // Calculate centrality bin centers and widths
    std::vector<double> centBinCenters(nCentBins);
    std::vector<double> centBinWidths(nCentBins);
    
    for (int i = 0; i < nCentBins; i++) {
        centBinCenters[i] = (fCentBins[i] + fCentBins[i+1]) / 2.0;
        centBinWidths[i] = (fCentBins[i+1] - fCentBins[i]) / 2.0;
    }
    
    // Create graphs for each particle type
    for (const auto& pidPair : fPidNames) {
        int pid = pidPair.first;
        std::string name = pidPair.second;
        
        // Create multiplicity graph
        TGraphErrors* gr_mult = new TGraphErrors(nCentBins);
        gr_mult->SetName(TString::Format("gr_mult_%s", name.c_str()));
        gr_mult->SetTitle(TString::Format("%s Multiplicity;Centrality (%%);dN/d#eta", name.c_str()));
        
        // Create multiplicity in rapidity window graph
        TGraphErrors* gr_mult_rap = new TGraphErrors(nCentBins);
        gr_mult_rap->SetName(TString::Format("gr_mult_rap_%s", name.c_str()));
        gr_mult_rap->SetTitle(TString::Format("%s Multiplicity;Centrality (%%);dN/dy", name.c_str()));
        
        // Create mean pT graph
        TGraphErrors* gr_pTmean = new TGraphErrors(nCentBins);
        gr_pTmean->SetName(TString::Format("gr_pTmean_%s", name.c_str()));
        gr_pTmean->SetTitle(TString::Format("%s Mean p_{T};Centrality (%%);#LTp_{T}#GT (GeV/c)", name.c_str()));
        
        // Fill the graphs
        for (int i = 0; i < nCentBins; i++) {
            if (fEventCounts[i] > 0) {
                // Calculate multiplicity per event
                double mult = static_cast<double>(fParticleCounts[pid][i]) / fEventCounts[i];
                double mult_err = TMath::Sqrt(mult / fEventCounts[i]); // Statistical error
                
                // Calculate multiplicity in rapidity window per event
                double mult_rap = static_cast<double>(fParticleCountsRap[pid][i]) / fEventCounts[i];
                double mult_rap_err = TMath::Sqrt(mult_rap / fEventCounts[i]); // Statistical error
                
                // Calculate mean pT
                double mean_pt = 0.0;
                double mean_pt_err = 0.0;
                
                if (fParticleCounts[pid][i] > 0) {
                    mean_pt = fMeanPtSums[pid][i] / fParticleCounts[pid][i];
                    // Simplified error estimation
                    mean_pt_err = mean_pt / TMath::Sqrt(fParticleCounts[pid][i]);
                }
                
                // Set points in graphs
                gr_mult->SetPoint(i, centBinCenters[i], mult);
                gr_mult->SetPointError(i, centBinWidths[i], mult_err);
                
                gr_mult_rap->SetPoint(i, centBinCenters[i], mult_rap);
                gr_mult_rap->SetPointError(i, centBinWidths[i], mult_rap_err);
                
                gr_pTmean->SetPoint(i, centBinCenters[i], mean_pt);
                gr_pTmean->SetPointError(i, centBinWidths[i], mean_pt_err);
            } else {
                // No events in this centrality bin
                gr_mult->SetPoint(i, centBinCenters[i], 0);
                gr_mult->SetPointError(i, centBinWidths[i], 0);
                
                gr_mult_rap->SetPoint(i, centBinCenters[i], 0);
                gr_mult_rap->SetPointError(i, centBinWidths[i], 0);
                
                gr_pTmean->SetPoint(i, centBinCenters[i], 0);
                gr_pTmean->SetPointError(i, centBinWidths[i], 0);
            }
        }
        
        // Store the graphs
        fMultGraphs[name] = gr_mult;
        fMultRapGraphs[name] = gr_mult_rap;
        fMeanPtGraphs[name] = gr_pTmean;
    }
}

TGraphErrors* AliJAnalysisMultiplicityMeanPt::GetMultGraph(const char* particle) {
    std::string key(particle);
    if (fMultGraphs.find(key) != fMultGraphs.end()) {
        return fMultGraphs[key];
    }
    return nullptr;
}

TGraphErrors* AliJAnalysisMultiplicityMeanPt::GetMultRapGraph(const char* particle) {
    std::string key(particle);
    if (fMultRapGraphs.find(key) != fMultRapGraphs.end()) {
        return fMultRapGraphs[key];
    }
    return nullptr;
}

TGraphErrors* AliJAnalysisMultiplicityMeanPt::GetMeanPtGraph(const char* particle) {
    std::string key(particle);
    if (fMeanPtGraphs.find(key) != fMeanPtGraphs.end()) {
        return fMeanPtGraphs[key];
    }
    return nullptr;
}

void AliJAnalysisMultiplicityMeanPt::WriteGraphs(TObject* dir) {
    // Cast to TDirectory if it's a directory
    TDirectory* directory = dynamic_cast<TDirectory*>(dir);
    if (!directory) {
        std::cerr << "Error: Provided object is not a TDirectory." << std::endl;
        return;
    }

    // Change to the directory
    directory->cd();
    
    // Write all multiplicity graphs
    for (auto& pair : fMultGraphs) {
        if (pair.second) pair.second->Write();
    }
    
    // Write all multiplicity in rapidity window graphs
    for (auto& pair : fMultRapGraphs) {
        if (pair.second) pair.second->Write();
    }
    
    // Write all mean pT graphs
    for (auto& pair : fMeanPtGraphs) {
        if (pair.second) pair.second->Write();
    }
} 