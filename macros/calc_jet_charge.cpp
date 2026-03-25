#include <iostream>
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

void calc_jet_charge() {
    // -----------------------------------------------------------------------------
    // 1. Open Input File
    // -----------------------------------------------------------------------------
    TFile* inFile = TFile::Open("ep_18x275_minQ2_100_100files.root", "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Could not open input file small_matched_output.root!" << std::endl;
        return;
    }
    TTree* tree = (TTree*)inFile->Get("matched_tree");

    // -----------------------------------------------------------------------------
    // 2. Set Up Input Branches
    // -----------------------------------------------------------------------------
    std::vector<float> *recoJet_pt = nullptr, *recoJet_eta = nullptr, *recoJet_phi = nullptr;
    std::vector<float> *matchGenJet_pt = nullptr, *matchGenJet_eta = nullptr, *matchGenJet_phi = nullptr;
    std::vector<float> *recoTrk_pt = nullptr, *recoTrk_eta = nullptr, *recoTrk_phi = nullptr, *recoTrk_charge = nullptr;
    std::vector<float> *matchGenTrk_pt = nullptr, *matchGenTrk_eta = nullptr, *matchGenTrk_phi = nullptr, *matchGenTrk_charge = nullptr;

    tree->SetBranchAddress("recoJet_pt", &recoJet_pt);
    tree->SetBranchAddress("recoJet_eta", &recoJet_eta);
    tree->SetBranchAddress("recoJet_phi", &recoJet_phi);
    tree->SetBranchAddress("matchGenJet_pt", &matchGenJet_pt);
    tree->SetBranchAddress("matchGenJet_eta", &matchGenJet_eta);
    tree->SetBranchAddress("matchGenJet_phi", &matchGenJet_phi);
    tree->SetBranchAddress("recoTrk_pt", &recoTrk_pt);
    tree->SetBranchAddress("recoTrk_eta", &recoTrk_eta);
    tree->SetBranchAddress("recoTrk_phi", &recoTrk_phi);
    tree->SetBranchAddress("recoTrk_charge", &recoTrk_charge);
    tree->SetBranchAddress("matchGenTrk_pt", &matchGenTrk_pt);
    tree->SetBranchAddress("matchGenTrk_eta", &matchGenTrk_eta);
    tree->SetBranchAddress("matchGenTrk_phi", &matchGenTrk_phi);
    tree->SetBranchAddress("matchGenTrk_charge", &matchGenTrk_charge);

    // -----------------------------------------------------------------------------
    // 3. Set Up Output File & Histograms
    // -----------------------------------------------------------------------------
    TFile* outFile = TFile::Open("jet_charge_output.root", "RECREATE");

    const float R_CONE    = 0.5;
    const float MIN_JET_PT = 7.0;
    const float MIN_TRK_PT = 0.5;
    const float KAPPA     = 0.5;

    TH1F* h_Q_RR = new TH1F("h_Q_RR", "Reco Jet + Reco Trk | Jet p_{T} > 7 GeV, Trk p_{T} > 0.5 GeV, nTrk > 1, #kappa = 0.5;Jet Charge;Counts", 200, -2, 2);
    TH1F* h_Q_RG = new TH1F("h_Q_RG", "Reco Jet + Gen Trk  | Jet p_{T} > 7 GeV, Trk p_{T} > 0.5 GeV, nTrk > 1, #kappa = 0.5;Jet Charge;Counts", 200, -2, 2);
    TH1F* h_Q_GR = new TH1F("h_Q_GR", "Gen Jet + Reco Trk  | Jet p_{T} > 7 GeV, Trk p_{T} > 0.5 GeV, nTrk > 1, #kappa = 0.5;Jet Charge;Counts", 200, -2, 2);
    TH1F* h_Q_GG = new TH1F("h_Q_GG", "Gen Jet + Gen Trk   | Jet p_{T} > 7 GeV, Trk p_{T} > 0.5 GeV, nTrk > 1, #kappa = 0.5;Jet Charge;Counts", 200, -2, 2);

    TH1F* h_topo_RR = new TH1F("h_topo_RR", "Jet Topologies (Reco Jet, Reco Trk);;Fraction of Jets", 3, 0, 3);
    TH1F* h_topo_RG = new TH1F("h_topo_RG", "Jet Topologies (Reco Jet, Gen Trk);;Fraction of Jets",  3, 0, 3);
    TH1F* h_topo_GR = new TH1F("h_topo_GR", "Jet Topologies (Gen Jet, Reco Trk);;Fraction of Jets",  3, 0, 3);
    TH1F* h_topo_GG = new TH1F("h_topo_GG", "Jet Topologies (Gen Jet, Gen Trk);;Fraction of Jets",   3, 0, 3);

    const char* binLabels[3] = {"Single Trk", "All Same Charge", "Leading Pt >= 90%"};
    for (int i = 1; i <= 3; i++) {
        h_topo_RR->GetXaxis()->SetBinLabel(i, binLabels[i-1]);
        h_topo_RG->GetXaxis()->SetBinLabel(i, binLabels[i-1]);
        h_topo_GR->GetXaxis()->SetBinLabel(i, binLabels[i-1]);
        h_topo_GG->GetXaxis()->SetBinLabel(i, binLabels[i-1]);
    }

    int totalJets_RR = 0, totalJets_RG = 0, totalJets_GR = 0, totalJets_GG = 0;

    // -----------------------------------------------------------------------------
    // 4. Event Loop
    // -----------------------------------------------------------------------------
    Long64_t nEntries = tree->GetEntries();
    std::cout << "Calculating jet charge for " << nEntries << " events..." << std::endl;

    for (Long64_t iEv = 0; iEv < nEntries; iEv++) {
        tree->GetEntry(iEv);

        for (size_t iJ = 0; iJ < recoJet_pt->size(); iJ++) {

            TVector3 vJ_reco, vJ_gen;
            bool hasGenJ = (matchGenJet_pt->at(iJ) != -999.0);

            vJ_reco.SetPtEtaPhi(recoJet_pt->at(iJ), recoJet_eta->at(iJ), recoJet_phi->at(iJ));
            if (hasGenJ) {
                vJ_gen.SetPtEtaPhi(matchGenJet_pt->at(iJ), matchGenJet_eta->at(iJ), matchGenJet_phi->at(iJ));
            }

            bool useRecoJet = (vJ_reco.Pt() > MIN_JET_PT);
            bool useGenJet  = (hasGenJ && vJ_gen.Pt() > MIN_JET_PT);

            if (!useRecoJet && !useGenJet) continue;

            // Per-jet accumulators
            int nTrk_RR = 0, nPos_RR = 0, nNeg_RR = 0;
            float rawPtSum_RR = 0, sumQWPt_RR = 0, maxPt_RR = 0;

            int nTrk_RG = 0, nPos_RG = 0, nNeg_RG = 0;
            float rawPtSum_RG = 0, sumQWPt_RG = 0, maxPt_RG = 0;

            int nTrk_GR = 0, nPos_GR = 0, nNeg_GR = 0;
            float rawPtSum_GR = 0, sumQWPt_GR = 0, maxPt_GR = 0;

            int nTrk_GG = 0, nPos_GG = 0, nNeg_GG = 0;
            float rawPtSum_GG = 0, sumQWPt_GG = 0, maxPt_GG = 0;

            // Track Loop
            for (size_t iT = 0; iT < recoTrk_pt->size(); iT++) {
                bool hasGenT = (matchGenTrk_pt->at(iT) != -999.0);

                TVector3 vT_reco, vT_gen;
                vT_reco.SetPtEtaPhi(recoTrk_pt->at(iT), recoTrk_eta->at(iT), recoTrk_phi->at(iT));
                float q_reco = recoTrk_charge->at(iT);

                if (hasGenT) {
                    vT_gen.SetPtEtaPhi(matchGenTrk_pt->at(iT), matchGenTrk_eta->at(iT), matchGenTrk_phi->at(iT));
                }

                bool validGenT = hasGenT && (matchGenTrk_charge->at(iT) != 0);
                float q_gen = validGenT ? matchGenTrk_charge->at(iT) : 0;

                // Combo 1: RR
                if (useRecoJet && vT_reco.Pt() > MIN_TRK_PT && vJ_reco.DeltaR(vT_reco) < R_CONE) {
                    nTrk_RR++;
                    rawPtSum_RR += vT_reco.Pt();
                    sumQWPt_RR  += q_reco * std::pow(vT_reco.Pt(), KAPPA);
                    if (vT_reco.Pt() > maxPt_RR) maxPt_RR = vT_reco.Pt();
                    if (q_reco > 0) nPos_RR++; else if (q_reco < 0) nNeg_RR++;
                }

                // Combo 2: RG
                if (useRecoJet && validGenT && vT_gen.Pt() > MIN_TRK_PT && vJ_reco.DeltaR(vT_gen) < R_CONE) {
                    nTrk_RG++;
                    rawPtSum_RG += vT_gen.Pt();
                    sumQWPt_RG  += q_gen * std::pow(vT_gen.Pt(), KAPPA);
                    if (vT_gen.Pt() > maxPt_RG) maxPt_RG = vT_gen.Pt();
                    if (q_gen > 0) nPos_RG++; else if (q_gen < 0) nNeg_RG++;
                }

                // Combo 3: GR
                if (useGenJet && vT_reco.Pt() > MIN_TRK_PT && vJ_gen.DeltaR(vT_reco) < R_CONE) {
                    nTrk_GR++;
                    rawPtSum_GR += vT_reco.Pt();
                    sumQWPt_GR  += q_reco * std::pow(vT_reco.Pt(), KAPPA);
                    if (vT_reco.Pt() > maxPt_GR) maxPt_GR = vT_reco.Pt();
                    if (q_reco > 0) nPos_GR++; else if (q_reco < 0) nNeg_GR++;
                }

                // Combo 4: GG
                if (useGenJet && validGenT && vT_gen.Pt() > MIN_TRK_PT && vJ_gen.DeltaR(vT_gen) < R_CONE) {
                    nTrk_GG++;
                    rawPtSum_GG += vT_gen.Pt();
                    sumQWPt_GG  += q_gen * std::pow(vT_gen.Pt(), KAPPA);
                    if (vT_gen.Pt() > maxPt_GG) maxPt_GG = vT_gen.Pt();
                    if (q_gen > 0) nPos_GG++; else if (q_gen < 0) nNeg_GG++;
                }
            }

            // Fill RR
            if (nTrk_RR > 0) {
                totalJets_RR++;
                if (nTrk_RR == 1) h_topo_RR->Fill(0.5);
                if (nTrk_RR > 2) {
                    if (nPos_RR == nTrk_RR || nNeg_RR == nTrk_RR) h_topo_RR->Fill(1.5);
                    if (maxPt_RR / rawPtSum_RR >= 0.9) h_topo_RR->Fill(2.5);
                    h_Q_RR->Fill(sumQWPt_RR / std::pow(rawPtSum_RR, KAPPA));
                }
            }

            // Fill RG
            if (nTrk_RG > 0) {
                totalJets_RG++;
                if (nTrk_RG == 1) h_topo_RG->Fill(0.5);
                if (nTrk_RG > 2) {
                    if (nPos_RG == nTrk_RG || nNeg_RG == nTrk_RG) h_topo_RG->Fill(1.5);
                    if (maxPt_RG / rawPtSum_RG >= 0.9) h_topo_RG->Fill(2.5);
                    h_Q_RG->Fill(sumQWPt_RG / std::pow(rawPtSum_RG, KAPPA));
                }
            }

            // Fill GR
            if (nTrk_GR > 0) {
                totalJets_GR++;
                if (nTrk_GR == 1) h_topo_GR->Fill(0.5);
                if (nTrk_GR > 2) {
                    if (nPos_GR == nTrk_GR || nNeg_GR == nTrk_GR) h_topo_GR->Fill(1.5);
                    if (maxPt_GR / rawPtSum_GR >= 0.9) h_topo_GR->Fill(2.5);
                    h_Q_GR->Fill(sumQWPt_GR / std::pow(rawPtSum_GR, KAPPA));
                }
            }

            // Fill GG
            if (nTrk_GG > 0) {
                totalJets_GG++;
                if (nTrk_GG == 1) h_topo_GG->Fill(0.5);
                if (nTrk_GG > 2) {
                    if (nPos_GG == nTrk_GG || nNeg_GG == nTrk_GG) h_topo_GG->Fill(1.5);
                    if (maxPt_GG / rawPtSum_GG >= 0.9) h_topo_GG->Fill(2.5);
                    h_Q_GG->Fill(sumQWPt_GG / std::pow(rawPtSum_GG, KAPPA));
                }
            }
        }
    }

    // -----------------------------------------------------------------------------
    // 5. Scale Topologies & Save
    // -----------------------------------------------------------------------------
    if (totalJets_RR > 0) h_topo_RR->Scale(1.0 / totalJets_RR);
    if (totalJets_RG > 0) h_topo_RG->Scale(1.0 / totalJets_RG);
    if (totalJets_GR > 0) h_topo_GR->Scale(1.0 / totalJets_GR);
    if (totalJets_GG > 0) h_topo_GG->Scale(1.0 / totalJets_GG);

    outFile->cd();

    h_Q_RR->SetLineColor(kBlue);
    h_Q_RG->SetLineColor(kRed);
    h_Q_GR->SetLineColor(kGreen+2);
    h_Q_GG->SetLineColor(kBlack);

    h_Q_RR->Write(); h_Q_RG->Write(); h_Q_GR->Write(); h_Q_GG->Write();
    h_topo_RR->Write(); h_topo_RG->Write(); h_topo_GR->Write(); h_topo_GG->Write();

    // --- Draw Jet Charge Canvas ---
    TCanvas* c1 = new TCanvas("c1", "Jet Charge Overlay", 800, 600);
    h_Q_GG->Draw("HIST");
    h_Q_RR->Draw("HIST SAME");
    h_Q_RG->Draw("HIST SAME");
    h_Q_GR->Draw("HIST SAME");

    TLegend* leg1 = new TLegend(0.65, 0.7, 0.88, 0.88);
    leg1->AddEntry(h_Q_RR, "Reco Jet, Reco Trk", "l");
    leg1->AddEntry(h_Q_RG, "Reco Jet, Gen Trk",  "l");
    leg1->AddEntry(h_Q_GR, "Gen Jet, Reco Trk",  "l");
    leg1->AddEntry(h_Q_GG, "Gen Jet, Gen Trk",   "l");
    leg1->Draw();
    c1->SaveAs("JetCharge_Overlay.png");

    // --- Draw Topology Canvas ---
    TCanvas* c2 = new TCanvas("c2", "Jet Topologies", 1000, 800);
    c2->Divide(2, 2);

    h_topo_RR->SetFillColor(kBlue-9);
    h_topo_RG->SetFillColor(kRed-9);
    h_topo_GR->SetFillColor(kGreen-9);
    h_topo_GG->SetFillColor(kGray);

    c2->cd(1); h_topo_RR->Draw("BAR"); h_topo_RR->SetMinimum(0);
    c2->cd(2); h_topo_RG->Draw("BAR"); h_topo_RG->SetMinimum(0);
    c2->cd(3); h_topo_GR->Draw("BAR"); h_topo_GR->SetMinimum(0);
    c2->cd(4); h_topo_GG->Draw("BAR"); h_topo_GG->SetMinimum(0);
    c2->SaveAs("JetTopologies_Fractions.png");

    std::cout << "Done! Results saved in jet_charge_output.root" << std::endl;
}