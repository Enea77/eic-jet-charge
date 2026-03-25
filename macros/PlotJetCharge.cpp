#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLatex.h" // Added for printing text on the plot

void PlotJetCharge(
    const TString& inputFile = "ep_18x275_minQ2_100_Allfiles.root", //ep_18x275_minQ2_100_100files
    TString jetType = "reco",     // accepts "reco", "gen", or "matched"
    TString trkType = "reco",     // accepts "reco", "gen", or "matched"
    float KAPPA = 0.5,
    float MIN_JET_PT = 20.0,
    float MIN_TRK_PT = 0.5,
    int MIN_TRKS = 2,             // Minimum number of tracks required in the jet cone
    float R_CONE = 0.4,           // Added: Jet cone radius
    float MAX_ETA = 3.          // Added: Maximum absolute jet pseudorapidity
) {
    // -----------------------------------------------------------------------------
    // 1. Input Validation & Setup
    // -----------------------------------------------------------------------------
    jetType.ToLower();
    trkType.ToLower();

    if (jetType != "reco" && jetType != "gen" && jetType != "matched") {
        std::cerr << "Error: jetType must be 'reco', 'gen', or 'matched'." << std::endl;
        return;
    }
    if (trkType != "reco" && trkType != "gen" && trkType != "matched") {
        std::cerr << "Error: trkType must be 'reco', 'gen', or 'matched'." << std::endl;
        return;
    }

    TFile* inFile = TFile::Open(inputFile, "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Could not open input file " << inputFile << "!" << std::endl;
        return;
    }
    TTree* tree = (TTree*)inFile->Get("matched_tree");

    // -----------------------------------------------------------------------------
    // 2. Set Up Input Branches Dynamically
    // -----------------------------------------------------------------------------
    std::vector<float> *jet_pt = nullptr, *jet_eta = nullptr, *jet_phi = nullptr;
    std::vector<float> *trk_pt = nullptr, *trk_eta = nullptr, *trk_phi = nullptr, *trk_charge = nullptr;

    // Map the string argument to the exact branch prefix in the TTree
    TString jetPrefix;
    if (jetType == "reco") jetPrefix = "recoJet";
    else if (jetType == "gen") jetPrefix = "genJet";
    else if (jetType == "matched") jetPrefix = "matchGenJet";

    TString trkPrefix;
    if (trkType == "reco") trkPrefix = "recoTrk";
    else if (trkType == "gen") trkPrefix = "genTrk";
    else if (trkType == "matched") trkPrefix = "matchGenTrk";

    tree->SetBranchAddress(jetPrefix + "_pt",  &jet_pt);
    tree->SetBranchAddress(jetPrefix + "_eta", &jet_eta);
    tree->SetBranchAddress(jetPrefix + "_phi", &jet_phi);

    tree->SetBranchAddress(trkPrefix + "_pt",  &trk_pt);
    tree->SetBranchAddress(trkPrefix + "_eta", &trk_eta);
    tree->SetBranchAddress(trkPrefix + "_phi", &trk_phi);
    tree->SetBranchAddress(trkPrefix + "_charge", &trk_charge);

    // -----------------------------------------------------------------------------
    // 3. Set Up Histogram & Output Naming
    // -----------------------------------------------------------------------------
    TString baseName = inputFile;
    baseName.ReplaceAll(".root", "");
    if (baseName.Contains("/")) {
        baseName = baseName(baseName.Last('/') + 1, baseName.Length());
    }

    // Incorporate all arguments into the output filenames
    TString outSuffix = Form("_%sJet_%sTrk_k%.2f_jpt%.1f_tpt%.1f_minTrk%d_R%.2f_eta%.1f", 
                             jetType.Data(), trkType.Data(), KAPPA, MIN_JET_PT, MIN_TRK_PT, MIN_TRKS, R_CONE, MAX_ETA);
    TString outRootName = "Plots/jet_charge_" + baseName + outSuffix + ".root";
    TString outPngName  = "Plots/JetCharge_" + baseName + outSuffix + ".png";

    TFile* outFile = TFile::Open(outRootName, "RECREATE");

    // Cleaned up title
    TString hTitle = Form("%s Jet + %s Trk;Jet Charge;Counts", jetType.Data(), trkType.Data());
    TH1F* h_Q = new TH1F("h_JetCharge", hTitle, 200, -2, 2);

    // -----------------------------------------------------------------------------
    // 4. Event Loop
    // -----------------------------------------------------------------------------
    Long64_t nEntries = tree->GetEntries();
    std::cout << "Calculating jet charge for " << nEntries << " events..." << std::endl;
    std::cout << "Using: " << jetPrefix << " and " << trkPrefix << std::endl;

    for (Long64_t iEv = 0; iEv < nEntries; iEv++) {
        tree->GetEntry(iEv);

        for (size_t iJ = 0; iJ < jet_pt->size(); iJ++) {
            // Skip invalid/dummy jets, pT cut, and eta cut
            if (jet_pt->at(iJ) == -999.0 || jet_pt->at(iJ) <= MIN_JET_PT) continue;
            if (std::abs(jet_eta->at(iJ)) > MAX_ETA) continue; // New Eta Cut

            TVector3 vJ;
            vJ.SetPtEtaPhi(jet_pt->at(iJ), jet_eta->at(iJ), jet_phi->at(iJ));

            float rawPtSum = 0;
            float sumQWPt = 0;
            int nTrk = 0;

            for (size_t iT = 0; iT < trk_pt->size(); iT++) {
                // Skip invalid/dummy tracks or tracks below pT threshold
                if (trk_pt->at(iT) == -999.0 || trk_pt->at(iT) <= MIN_TRK_PT) continue;
                
                float q = trk_charge->at(iT);
                if (q == 0) continue; // Skip neutral particles 

                TVector3 vT;
                vT.SetPtEtaPhi(trk_pt->at(iT), trk_eta->at(iT), trk_phi->at(iT));

                // Cone matching using the new R_CONE argument
                if (vJ.DeltaR(vT) < R_CONE) {
                    nTrk++;
                    rawPtSum += vT.Pt();
                    sumQWPt  += q * std::pow(vT.Pt(), KAPPA);
                }
            }

            // Fill histogram ONLY if the number of tracks inside the cone passes the cut
            if (nTrk >= MIN_TRKS) {
                h_Q->Fill(sumQWPt / std::pow(rawPtSum, KAPPA));
            }
        }
    }

    // -----------------------------------------------------------------------------
    // 5. Save Results & Plot
    // -----------------------------------------------------------------------------
    outFile->cd();
    h_Q->Write();

    TCanvas* c1 = new TCanvas("c1", "Jet Charge", 800, 600);
    h_Q->SetLineColor(kBlue + 2);
    h_Q->SetLineWidth(2);
    
    // Scale the Y-axis max by 30% so the text box doesn't overlap the peak of the histogram
    h_Q->SetMaximum(h_Q->GetMaximum() * 1.3);
    h_Q->Draw("HIST");

    // Print cuts on the plot using TLatex
    TLatex latex;
    latex.SetNDC();             // Use normalized device coordinates (0 to 1)
    latex.SetTextFont(42);      // Standard root font
    latex.SetTextSize(0.035);   // Adjust size as needed

    // Draw text lines (X = 0.15 is slightly inset from left axis, Y goes top to bottom)
    latex.DrawLatex(0.15, 0.85, Form("Jet p_{T} > %.1f GeV, |#eta| < %.1f", MIN_JET_PT, MAX_ETA));
    latex.DrawLatex(0.15, 0.80, Form("Trk p_{T} > %.1f GeV", MIN_TRK_PT));
    latex.DrawLatex(0.15, 0.75, Form("nTrk #geq %d, R = %.2f", MIN_TRKS, R_CONE));
    latex.DrawLatex(0.15, 0.70, Form("#kappa = %.2f", KAPPA));
    
    c1->SaveAs(outPngName);

    std::cout << "Done! Plot saved as: " << outPngName << std::endl;
    std::cout << "ROOT file saved as: "  << outRootName << std::endl;
}