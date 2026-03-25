#include <iostream>
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TCanvas.h"

void calc_gen_jet_charge() {
    // -----------------------------------------------------------------------------
    // 1. Open Input File
    // -----------------------------------------------------------------------------
    TFile* inFile = TFile::Open("small_matched_output.root", "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Could not open input file small_matched_output.root!" << std::endl;
        return;
    }
    TTree* tree = (TTree*)inFile->Get("matched_tree");

    // -----------------------------------------------------------------------------
    // 2. Set Up Input Branches
    // -----------------------------------------------------------------------------
    std::vector<float> *genJet_pt = nullptr, *genJet_eta = nullptr, *genJet_phi = nullptr;
    std::vector<float> *genTrk_pt = nullptr, *genTrk_eta = nullptr, *genTrk_phi = nullptr, *genTrk_charge = nullptr;

    tree->SetBranchAddress("genJet_pt", &genJet_pt);
    tree->SetBranchAddress("genJet_eta", &genJet_eta);
    tree->SetBranchAddress("genJet_phi", &genJet_phi);
    
    tree->SetBranchAddress("genTrk_pt", &genTrk_pt);
    tree->SetBranchAddress("genTrk_eta", &genTrk_eta);
    tree->SetBranchAddress("genTrk_phi", &genTrk_phi);
    tree->SetBranchAddress("genTrk_charge", &genTrk_charge);

    // -----------------------------------------------------------------------------
    // 3. Set Up Output File & Histograms
    // -----------------------------------------------------------------------------
    TFile* outFile = TFile::Open("gen_jet_charge_output.root", "RECREATE");

    TH1F* h_Q_GG = new TH1F("h_Q_GG", "Gen Jet Charge (nTrk > 1, p_{T} > 7 GeV, |#eta| < 1.0, #kappa = 0.5);Jet Charge;Counts", 200, -2, 2);
    TH1F* h_topo_GG = new TH1F("h_topo_GG", "Gen Jet Topologies;;Fraction of Jets", 3, 0, 3);

    const char* binLabels[3] = {"Single Trk", "All Same Charge", "Leading Pt >= 90%"};
    for (int i = 1; i <= 3; i++) {
        h_topo_GG->GetXaxis()->SetBinLabel(i, binLabels[i-1]);
    }

    const float R_CONE = 1.0; 
    const float MIN_JET_PT = 7.0;
    const float MAX_JET_ETA = 3.0; // Eta cut parameter
    const float KAPPA = 0.5; // Momentum weight power

    int totalJets_GG = 0;

    // -----------------------------------------------------------------------------
    // 4. Event Loop
    // -----------------------------------------------------------------------------
    Long64_t nEntries = tree->GetEntries();
    std::cout << "Processing " << nEntries << " events for pure Gen Jet Charge..." << std::endl;

    for (Long64_t iEv = 0; iEv < nEntries; iEv++) {
        tree->GetEntry(iEv);

        for (size_t iJ = 0; iJ < genJet_pt->size(); iJ++) {
            
            // --- Jet Kinematic Cuts ---
            if (genJet_pt->at(iJ) <= MIN_JET_PT) continue;
            if (std::abs(genJet_eta->at(iJ)) > MAX_JET_ETA) continue; // Added Eta Cut

            TVector3 vJ_gen;
            vJ_gen.SetPtEtaPhi(genJet_pt->at(iJ), genJet_eta->at(iJ), genJet_phi->at(iJ));

            int nTrk_GG = 0, nPos_GG = 0, nNeg_GG = 0; 
            float sumPt_GG = 0, sumQPt_GG = 0, maxPt_GG = 0;
            float rawPtSum_GG = 0; 
            
            // Track totals for the entire cone (Charged + Neutral)
            TVector3 coneVectorSum(0, 0, 0);
            float coneScalarSum = 0;

            std::vector<int> tracksInCone; 
            std::vector<int> neutralTracksInCone;

            for (size_t iT = 0; iT < genTrk_pt->size(); iT++) {
                TVector3 vT_gen;
                vT_gen.SetPtEtaPhi(genTrk_pt->at(iT), genTrk_eta->at(iT), genTrk_phi->at(iT));
                if (vT_gen.Pt() < 0.5) continue;

                // If particle is inside the jet cone
                if (vJ_gen.DeltaR(vT_gen) < R_CONE) {
                    
                    // Add to total cone momentum trackers
                    coneVectorSum += vT_gen;
                    coneScalarSum += vT_gen.Pt();
                    
                    float q_gen = genTrk_charge->at(iT);
                    
                    // Separate logic based on charge
                    if (q_gen == 0) {
                        neutralTracksInCone.push_back(iT);
                    } else {
                        nTrk_GG++;
                        
                        // Apply power weighting kappa to track pT for Jet Charge
                        float wPt = std::pow(vT_gen.Pt(), KAPPA); 
                        
                        sumPt_GG += wPt;
                        sumQPt_GG += q_gen * wPt;
                        rawPtSum_GG += vT_gen.Pt(); // Unweighted sum for topology logic
                        
                        if (vT_gen.Pt() > maxPt_GG) maxPt_GG = vT_gen.Pt();
                        if (q_gen > 0) nPos_GG++; else if (q_gen < 0) nNeg_GG++;
                        
                        tracksInCone.push_back(iT);
                    }
                }
            }

            if (nTrk_GG > 0) {
                totalJets_GG++;
                
                if (nTrk_GG < 3) {
                    h_topo_GG->Fill(0.5); 
                } else {
                    if (nPos_GG == nTrk_GG || nNeg_GG == nTrk_GG) h_topo_GG->Fill(1.5); 
                    if (maxPt_GG / rawPtSum_GG >= 0.9) h_topo_GG->Fill(2.5); 
                    
                    float jetCharge = sumQPt_GG / pow(rawPtSum_GG, KAPPA); //pow(genJet_pt->at(iJ), KAPPA);        //pow(rawPtSum_GG, KAPPA);
                    h_Q_GG->Fill(jetCharge);

                    // Debugging Printout for magnitude exactly == 1 (using small float tolerance)
                    if (std::abs(std::abs(jetCharge) - 1.0) < 1e-4) {
                        std::cout << "\n======================================================\n";
                        std::cout << "Event " << iEv << " | Jet " << iJ << " | Jet pT = " << genJet_pt->at(iJ) 
                                  << " | Jet Charge = " << jetCharge << " | Charged nTrks = " << nTrk_GG << "\n";
                        
                        std::cout << "--> Vector Sum pT (All Particles) = " << coneVectorSum.Pt() << "\n";
                        std::cout << "--> Scalar Sum pT (All Particles) = " << coneScalarSum << "\n\n";
                        
                        std::cout << "Charged Tracks in Cone:\n";
                        for (int idx : tracksInCone) {
                            std::cout << "  -> Trk pT = " << genTrk_pt->at(idx) 
                                      << "\t| Charge = " << genTrk_charge->at(idx) << "\n";
                        }
                        
                        std::cout << "Neutral Particles in Cone:\n";
                        if (neutralTracksInCone.empty()) {
                            std::cout << "  (None)\n";
                        } else {
                            for (int idx : neutralTracksInCone) {
                                std::cout << "  -> Neu pT = " << genTrk_pt->at(idx) 
                                          << "\t| Charge = 0\n";
                            }
                        }
                    }
                }
            }
        }
    }

    // -----------------------------------------------------------------------------
    // 5. Scale Topologies & Save
    // -----------------------------------------------------------------------------
    if (totalJets_GG > 0) h_topo_GG->Scale(1.0 / totalJets_GG);

    outFile->cd();
    h_Q_GG->SetLineColor(kBlack);
    h_Q_GG->Write();
    h_topo_GG->Write();

    TCanvas* c1 = new TCanvas("c1", "Gen Jet Charge", 800, 600);
    h_Q_GG->Draw("HIST"); 
    c1->SaveAs("GenJetCharge.png");

    TCanvas* c2 = new TCanvas("c2", "Gen Jet Topologies", 800, 600);
    h_topo_GG->SetFillColor(kGray);
    h_topo_GG->Draw("BAR"); 
    h_topo_GG->SetMinimum(0);
    c2->SaveAs("GenJetTopologies.png");
    
    std::cout << "\nDone! Analyzed " << totalJets_GG << " valid jets." << std::endl;
}