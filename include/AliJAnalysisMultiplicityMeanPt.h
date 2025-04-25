#ifndef ALIJANALYSISMULTIPLICITYMEANPT_H
#define ALIJANALYSISMULTIPLICITYMEANPT_H

#include <vector>
#include <string>
#include <map>
#include <TObject.h>
#include <TGraphErrors.h>
#include <TString.h>

class AliJAnalysisMultiplicityMeanPt : public TObject {
public:
    AliJAnalysisMultiplicityMeanPt(const char* name = "MultiplcityMeanPtAnalysis");
    virtual ~AliJAnalysisMultiplicityMeanPt();
    
    // Initialize with system-specific centrality binning
    void SetSystemParameters(const char* collisionSystem);
    
    // Set custom centrality bins
    void SetCentralityBins(const std::vector<double>& centBins);
    
    // Process events for multiplicity and mean pT calculations
    void ProcessEvent(double centrality, const std::vector<double>& pt_values, const std::vector<int>& pid_codes);
    
    // Finalize and create output graphs
    void FinalizeMeasurements();
    
    // Get multiplicity graphs
    TGraphErrors* GetMultGraph(const char* particle);
    TGraphErrors* GetMultRapGraph(const char* particle);
    
    // Get mean pT graphs
    TGraphErrors* GetMeanPtGraph(const char* particle);
    
    // Write all graphs to a ROOT file
    void WriteGraphs(TObject* dir);
    
    // Particle types
    enum ParticleType {
        kCharged = 0,
        kPion = 1,
        kKaon = 2,
        kProton = 3
    };
    
private:
    TString fName;                          // Name of the analysis
    std::vector<double> fCentBins;          // Centrality bins
    std::map<int, std::string> fPidNames;   // Maps PID codes to names
    
    // Maps to store sum of particles per centrality bin
    std::map<int, std::vector<int>> fParticleCounts;    // Particle counts per centrality bin
    std::map<int, std::vector<double>> fMeanPtSums;     // Sum of pT per centrality bin
    std::map<int, std::vector<int>> fParticleCountsRap; // Particle counts in rapidity window
    
    // Track statistics per centrality bin
    std::vector<int> fEventCounts;          // Number of events per centrality bin
    
    // Output graphs
    std::map<std::string, TGraphErrors*> fMultGraphs;       // Multiplicity graphs
    std::map<std::string, TGraphErrors*> fMultRapGraphs;    // Multiplicity in rapidity window graphs
    std::map<std::string, TGraphErrors*> fMeanPtGraphs;     // Mean pT graphs
    
    // Utility functions
    int GetCentBin(double centrality);
    void InitializeMaps();
    
    ClassDef(AliJAnalysisMultiplicityMeanPt, 1);
};

#endif 