#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TMath.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TCanvas.h"

void process_eic_trees() {
    // -----------------------------------------------------------------------------
    // 1. Open Input File
    // -----------------------------------------------------------------------------
    TString inFileName = "D:/EIC/EPIC/RECO/23.12.0/epic_craterlake/DIS/NC/18x275/minQ2=100/eic.root";
    TFile* inFile = TFile::Open(inFileName, "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Could not open input file!" << std::endl;
        return;
    }

    // -----------------------------------------------------------------------------
    // 2. Set Up TTreeReader (Bypasses missing dictionaries)
    // -----------------------------------------------------------------------------
    TTreeReader reader("events", inFile);

    // --- JETS ---
    TTreeReaderArray<float> recoJetPx(reader, "ReconstructedChargedJets.momentum.x");
    TTreeReaderArray<float> recoJetPy(reader, "ReconstructedChargedJets.momentum.y");
    TTreeReaderArray<float> recoJetPz(reader, "ReconstructedChargedJets.momentum.z");

    TTreeReaderArray<float> genJetPx(reader, "GeneratedChargedJets.momentum.x");
    TTreeReaderArray<float> genJetPy(reader, "GeneratedChargedJets.momentum.y");
    TTreeReaderArray<float> genJetPz(reader, "GeneratedChargedJets.momentum.z");

    // --- TRACKS / PARTICLES ---
    TTreeReaderArray<float> recoPartPx(reader, "ReconstructedChargedParticles.momentum.x");
    TTreeReaderArray<float> recoPartPy(reader, "ReconstructedChargedParticles.momentum.y");
    TTreeReaderArray<float> recoPartPz(reader, "ReconstructedChargedParticles.momentum.z");
    TTreeReaderArray<float> recoPartCharge(reader, "ReconstructedChargedParticles.charge");

    TTreeReaderArray<float> genPartPx(reader, "GeneratedParticles.momentum.x");
    TTreeReaderArray<float> genPartPy(reader, "GeneratedParticles.momentum.y");
    TTreeReaderArray<float> genPartPz(reader, "GeneratedParticles.momentum.z");
    TTreeReaderArray<float> genPartCharge(reader, "GeneratedParticles.charge");

    // -----------------------------------------------------------------------------
    // 3. Set Up Output File, Tree, and Vectors
    // -----------------------------------------------------------------------------
    TFile* outFile = TFile::Open("small_matched_output.root", "RECREATE");
    TTree* outTree = new TTree("matched_tree", "Tree with matched Jets and Tracks");

    // Output variables (Matched)
    std::vector<float> v_recoJet_pt, v_recoJet_eta, v_recoJet_phi;
    std::vector<float> v_matchGenJet_pt, v_matchGenJet_eta, v_matchGenJet_phi;
    
    std::vector<float> v_recoTrk_pt, v_recoTrk_eta, v_recoTrk_phi, v_recoTrk_charge;
    std::vector<float> v_matchGenTrk_pt, v_matchGenTrk_eta, v_matchGenTrk_phi, v_matchGenTrk_charge;

    // Output variables (All Gen - independent of matching)
    std::vector<float> v_genJet_pt, v_genJet_eta, v_genJet_phi;
    std::vector<float> v_genTrk_pt, v_genTrk_eta, v_genTrk_phi, v_genTrk_charge;

    // Branches (Matched)
    outTree->Branch("recoJet_pt", &v_recoJet_pt);
    outTree->Branch("recoJet_eta", &v_recoJet_eta);
    outTree->Branch("recoJet_phi", &v_recoJet_phi);
    outTree->Branch("matchGenJet_pt", &v_matchGenJet_pt);
    outTree->Branch("matchGenJet_eta", &v_matchGenJet_eta);
    outTree->Branch("matchGenJet_phi", &v_matchGenJet_phi);

    outTree->Branch("recoTrk_pt", &v_recoTrk_pt);
    outTree->Branch("recoTrk_eta", &v_recoTrk_eta);
    outTree->Branch("recoTrk_phi", &v_recoTrk_phi);
    outTree->Branch("recoTrk_charge", &v_recoTrk_charge);
    outTree->Branch("matchGenTrk_pt", &v_matchGenTrk_pt);
    outTree->Branch("matchGenTrk_eta", &v_matchGenTrk_eta);
    outTree->Branch("matchGenTrk_phi", &v_matchGenTrk_phi);
    outTree->Branch("matchGenTrk_charge", &v_matchGenTrk_charge);

    // Branches (All Gen)
    outTree->Branch("genJet_pt", &v_genJet_pt);
    outTree->Branch("genJet_eta", &v_genJet_eta);
    outTree->Branch("genJet_phi", &v_genJet_phi);
    
    outTree->Branch("genTrk_pt", &v_genTrk_pt);
    outTree->Branch("genTrk_eta", &v_genTrk_eta);
    outTree->Branch("genTrk_phi", &v_genTrk_phi);
    outTree->Branch("genTrk_charge", &v_genTrk_charge);

    // Histograms
    TH1F* h_dR_jets = new TH1F("h_dR_jets", "#Delta R (Reco Charged Jet, closest Gen Jet);#Delta R;Counts", 100, 0, 0.5);
    TH1F* h_dR_tracks = new TH1F("h_dR_tracks", "#Delta R (Reco Track, closest Gen Particle);#Delta R;Counts", 100, 0, 0.5);

    const float MAX_DR_JETS = 0.1;
    const float MAX_DR_TRACKS = 0.1;

    // -----------------------------------------------------------------------------
    // 4. Event Loop
    // -----------------------------------------------------------------------------
    int evCount = 0;
    std::cout << "Starting event loop..." << std::endl;

    while (reader.Next()) {
        if (evCount % 1000 == 0) std::cout << "Processing Event " << evCount << std::endl;
        evCount++;

        // Clear vectors for the new event
        v_recoJet_pt.clear(); v_recoJet_eta.clear(); v_recoJet_phi.clear();
        v_matchGenJet_pt.clear(); v_matchGenJet_eta.clear(); v_matchGenJet_phi.clear();
        v_recoTrk_pt.clear(); v_recoTrk_eta.clear(); v_recoTrk_phi.clear(); v_recoTrk_charge.clear();
        v_matchGenTrk_pt.clear(); v_matchGenTrk_eta.clear(); v_matchGenTrk_phi.clear(); v_matchGenTrk_charge.clear();
        
        v_genJet_pt.clear(); v_genJet_eta.clear(); v_genJet_phi.clear();
        v_genTrk_pt.clear(); v_genTrk_eta.clear(); v_genTrk_phi.clear(); v_genTrk_charge.clear();

        // --- SAVE ALL GEN JETS ---
        for (int ig = 0; ig < genJetPx.GetSize(); ig++) {
            TVector3 genP(genJetPx[ig], genJetPy[ig], genJetPz[ig]);
            if (genP.Pt() == 0) continue;
            v_genJet_pt.push_back(genP.Pt());
            v_genJet_eta.push_back(genP.Eta());
            v_genJet_phi.push_back(genP.Phi());
        }

        // --- SAVE ALL GEN TRACKS ---
        for (int ig = 0; ig < genPartPx.GetSize(); ig++) {
            TVector3 genP(genPartPx[ig], genPartPy[ig], genPartPz[ig]);
            if (genP.Pt() == 0) continue;
            v_genTrk_pt.push_back(genP.Pt());
            v_genTrk_eta.push_back(genP.Eta());
            v_genTrk_phi.push_back(genP.Phi());
            v_genTrk_charge.push_back(genPartCharge[ig]);
        }

        // --- MATCH CHARGED JETS TO GEN JETS ---
        for (int ir = 0; ir < recoJetPx.GetSize(); ir++) {
            TVector3 recoP(recoJetPx[ir], recoJetPy[ir], recoJetPz[ir]);
            if (recoP.Pt() == 0) continue; 

            v_recoJet_pt.push_back(recoP.Pt());
            v_recoJet_eta.push_back(recoP.Eta());
            v_recoJet_phi.push_back(recoP.Phi());

            float minDR = 999.0;
            int matchedIdx = -1;

            for (int ig = 0; ig < genJetPx.GetSize(); ig++) {
                TVector3 genP(genJetPx[ig], genJetPy[ig], genJetPz[ig]);
                if (genP.Pt() == 0) continue;
                
                float dR = recoP.DeltaR(genP);
                if (dR < minDR) {
                    minDR = dR;
                    matchedIdx = ig;
                }
            }

            h_dR_jets->Fill(minDR);

            if (minDR < MAX_DR_JETS && matchedIdx != -1) {
                TVector3 genP(genJetPx[matchedIdx], genJetPy[matchedIdx], genJetPz[matchedIdx]);
                v_matchGenJet_pt.push_back(genP.Pt());
                v_matchGenJet_eta.push_back(genP.Eta());
                v_matchGenJet_phi.push_back(genP.Phi());
            } else {
                v_matchGenJet_pt.push_back(-999.0);
                v_matchGenJet_eta.push_back(-999.0);
                v_matchGenJet_phi.push_back(-999.0);
            }
        }

        // --- MATCH RECO TRACKS TO ALL GEN PARTICLES ---
        for (int ir = 0; ir < recoPartPx.GetSize(); ir++) {
            TVector3 recoP(recoPartPx[ir], recoPartPy[ir], recoPartPz[ir]);
            if (recoP.Pt() == 0) continue;

            v_recoTrk_pt.push_back(recoP.Pt());
            v_recoTrk_eta.push_back(recoP.Eta());
            v_recoTrk_phi.push_back(recoP.Phi());
            v_recoTrk_charge.push_back(recoPartCharge[ir]);

            float minDR = 999.0;
            int matchedIdx = -1;

            for (int ig = 0; ig < genPartPx.GetSize(); ig++) {
                TVector3 genP(genPartPx[ig], genPartPy[ig], genPartPz[ig]);
                if (genP.Pt() == 0) continue;

                float dR = recoP.DeltaR(genP);
                if (dR < minDR) {
                    minDR = dR;
                    matchedIdx = ig;
                }
            }

            h_dR_tracks->Fill(minDR);

            if (minDR < MAX_DR_TRACKS && matchedIdx != -1) {
                TVector3 genP(genPartPx[matchedIdx], genPartPy[matchedIdx], genPartPz[matchedIdx]);
                v_matchGenTrk_pt.push_back(genP.Pt());
                v_matchGenTrk_eta.push_back(genP.Eta());
                v_matchGenTrk_phi.push_back(genP.Phi());
                v_matchGenTrk_charge.push_back(genPartCharge[matchedIdx]);
            } else {
                v_matchGenTrk_pt.push_back(-999.0);
                v_matchGenTrk_eta.push_back(-999.0);
                v_matchGenTrk_phi.push_back(-999.0);
                v_matchGenTrk_charge.push_back(-999.0);
            }
        }

        outTree->Fill();
    }

    // -----------------------------------------------------------------------------
    // 5. Save Outputs and Draw Plots
    // -----------------------------------------------------------------------------
    outFile->cd();
    outTree->Write();
    h_dR_jets->Write();
    h_dR_tracks->Write();

    TCanvas* c1 = new TCanvas("c1", "Delta R Plots", 800, 400);
    c1->Divide(2, 1);
    
    c1->cd(1);
    gPad->SetLogy();
    h_dR_jets->Draw();
    
    c1->cd(2);
    gPad->SetLogy();
    h_dR_tracks->Draw();
    
    c1->SaveAs("DeltaR_matching_plots.png");

    outFile->Close();
    inFile->Close();
    std::cout << "Done! Processed " << evCount << " events." << std::endl;
    std::cout << "Saved output to small_matched_output.root and DeltaR_matching_plots.png" << std::endl;
}