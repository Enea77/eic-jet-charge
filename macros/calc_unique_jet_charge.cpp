#include <iostream>
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TMath.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TCanvas.h"

void calc_unique_jet_charge() {
    // -----------------------------------------------------------------------------
    // 1. Open Input File
    // -----------------------------------------------------------------------------
    TFile* inFile = TFile::Open("small_matched_output.root", "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Could not open input file!" << std::endl;
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
    TFile* outFile = TFile::Open("unique_jet_charge_output.root", "RECREATE");

    TH1F* h_Q_GG = new TH1F("h_Q_GG", "Gen Jet Charge (nTrk > 1, p_{T} > 7 GeV);Jet Charge;Counts", 100, -1.5, 1.5);
    TH1F* h_Closure = new TH1F("h_Closure", "Closure Test: Vector Sum p_{T} / Jet p_{T};Ratio;Counts", 100, 0.0, 2.0);
    TH1F* h_topo_GG = new TH1F("h_topo_GG", "Gen Jet Topologies;;Fraction of Jets", 3, 0, 3);

    const char* binLabels[3] = {"Single Trk", "All Same Charge", "Leading Pt >= 90%"};
    for (int i = 1; i <= 3; i++) h_topo_GG->GetXaxis()->SetBinLabel(i, binLabels[i-1]);

    // Anti-kT jet parameter used in EICrecon (Often 1.0 for DIS, change to 0.5 if needed)
    const float R_CONE = 1.0; 
    const float MIN_JET_PT = 7.0;
    const float KAPPA = 0.5; 

    int totalJets_GG = 0;

    // -----------------------------------------------------------------------------
    // 4. Event Loop
    // -----------------------------------------------------------------------------
    Long64_t nEntries = tree->GetEntries();
    std::cout << "Processing " << nEntries << " events using Unique Track Assignment..." << std::endl;

    for (Long64_t iEv = 0; iEv < nEntries; iEv++) {
        tree->GetEntry(iEv);

        int nJets = genJet_pt->size();
        if (nJets == 0) continue;

        // Create vectors to hold the properties assigned to each jet
        std::vector<TVector3> jets(nJets);
        std::vector<TVector3> jet_vectorSum(nJets, TVector3(0,0,0));
        
        std::vector<int> jet_nTrk(nJets, 0);
        std::vector<int> jet_nPos(nJets, 0);
        std::vector<int> jet_nNeg(nJets, 0);
        
        std::vector<float> jet_sumPt(nJets, 0.0);
        std::vector<float> jet_sumQPt(nJets, 0.0);
        std::vector<float> jet_rawPtSum(nJets, 0.0);
        std::vector<float> jet_maxPt(nJets, 0.0);

        // Load jet vectors
        for (int iJ = 0; iJ < nJets; iJ++) {
            jets[iJ].SetPtEtaPhi(genJet_pt->at(iJ), genJet_eta->at(iJ), genJet_phi->at(iJ));
        }

        // --- TRACK ASSIGNMENT LOOP ---
        for (size_t iT = 0; iT < genTrk_pt->size(); iT++) {
            TVector3 vT;
            vT.SetPtEtaPhi(genTrk_pt->at(iT), genTrk_eta->at(iT), genTrk_phi->at(iT));
            float q = genTrk_charge->at(iT);

            // Find the closest jet to this track
            float minDR = R_CONE;
            int bestJetIdx = -1;

            for (int iJ = 0; iJ < nJets; iJ++) {
                float dR = vT.DeltaR(jets[iJ]);
                if (dR < minDR) {
                    minDR = dR;
                    bestJetIdx = iJ;
                }
            }

            // If the track belongs to a jet within R_CONE, assign its properties ONLY to that jet
            if (bestJetIdx != -1) {
                jet_vectorSum[bestJetIdx] += vT;

                if (q != 0) { // Charged track logic
                    jet_nTrk[bestJetIdx]++;
                    
                    float wPt = std::pow(vT.Pt(), KAPPA);
                    jet_sumPt[bestJetIdx] += wPt;
                    jet_sumQPt[bestJetIdx] += q * wPt;
                    jet_rawPtSum[bestJetIdx] += vT.Pt();
                    
                    if (vT.Pt() > jet_maxPt[bestJetIdx]) jet_maxPt[bestJetIdx] = vT.Pt();
                    if (q > 0) jet_nPos[bestJetIdx]++; 
                    else if (q < 0) jet_nNeg[bestJetIdx]++;
                }
            }
        }

        // --- JET EVALUATION LOOP ---
        for (int iJ = 0; iJ < nJets; iJ++) {
            if (jets[iJ].Pt() <= MIN_JET_PT) continue;

            // Fill Closure Test
            float closureRatio = jet_vectorSum[iJ].Pt() / jets[iJ].Pt();
            h_Closure->Fill(closureRatio);

            if (jet_nTrk[iJ] > 0) {
                totalJets_GG++;
                
                if (jet_nTrk[iJ] == 1) {
                    h_topo_GG->Fill(0.5); 
                } else {
                    if (jet_nPos[iJ] == jet_nTrk[iJ] || jet_nNeg[iJ] == jet_nTrk[iJ]) h_topo_GG->Fill(1.5); 
                    if (jet_maxPt[iJ] / jet_rawPtSum[iJ] >= 0.9) h_topo_GG->Fill(2.5); 
                    
                    float jetCharge = jet_sumQPt[iJ] / jet_sumPt[iJ];
                    h_Q_GG->Fill(jetCharge);
                }
            }
        }
    }

    // -----------------------------------------------------------------------------
    // 5. Scale Topologies & Save
    // -----------------------------------------------------------------------------
    if (totalJets_GG > 0) h_topo_GG->Scale(1.0 / totalJets_GG);

    outFile->cd();
    h_Closure->Write();
    h_Q_GG->SetLineColor(kBlack);
    h_Q_GG->Write();
    h_topo_GG->Write();

    TCanvas* c1 = new TCanvas("c1", "Jet Physics", 1000, 400);
    c1->Divide(2, 1);
    
    c1->cd(1);
    h_Closure->SetLineColor(kRed);
    h_Closure->Draw("HIST");
    
    c1->cd(2);
    h_Q_GG->Draw("HIST"); 
    c1->SaveAs("UniqueJetCharge.png");

    //outFile->Close();
    //inFile->Close();
    
    std::cout << "\nDone! Analyzed " << totalJets_GG << " valid jets." << std::endl;
}