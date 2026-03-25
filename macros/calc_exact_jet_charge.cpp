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

void calc_exact_jet_charge() {
    // -----------------------------------------------------------------------------
    // 1. Open Input File
    // -----------------------------------------------------------------------------
    TFile* inFile = TFile::Open("D:/EIC/EPIC/RECO/23.12.0/epic_craterlake/DIS/NC/18x275/minQ2=100/eic.root", "READ"); // Change to your actual file
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Could not open input file!" << std::endl;
        return;
    }

    TTreeReader reader("events", inFile);

    // -----------------------------------------------------------------------------
    // 2. Set Up Input Branches
    // -----------------------------------------------------------------------------
    // Jets
    TTreeReaderArray<float> jetPx(reader, "GeneratedChargedJets.momentum.x");
    TTreeReaderArray<float> jetPy(reader, "GeneratedChargedJets.momentum.y");
    TTreeReaderArray<float> jetPz(reader, "GeneratedChargedJets.momentum.z");

    // The Relational Branches (The "Magic" Link)
    TTreeReaderArray<unsigned int> partBegin(reader, "GeneratedChargedJets.particles_begin");
    TTreeReaderArray<unsigned int> partEnd(reader,   "GeneratedChargedJets.particles_end");
    TTreeReaderArray<int>          partIdx(reader,   "_GeneratedChargedJets_particles.index");

    // Particles
    TTreeReaderArray<float> trkPx(reader, "GeneratedParticles.momentum.x");
    TTreeReaderArray<float> trkPy(reader, "GeneratedParticles.momentum.y");
    TTreeReaderArray<float> trkPz(reader, "GeneratedParticles.momentum.z");
    TTreeReaderArray<float> trkQ(reader,  "GeneratedParticles.charge");

    // -----------------------------------------------------------------------------
    // 3. Set Up Output File & Histograms
    // -----------------------------------------------------------------------------
    TFile* outFile = TFile::Open("exact_jet_charge_output.root", "RECREATE");

    // Jet Histograms
    TH1F* h_JetPt  = new TH1F("h_JetPt", "Gen Jet p_{T};Jet p_{T} [GeV];Counts", 100, 0, 50);
    TH1F* h_JetEta = new TH1F("h_JetEta", "Gen Jet #eta;Jet #eta;Counts", 100, -4.0, 4.0);
    TH1F* h_Q      = new TH1F("h_Q", "Exact Jet Charge (nTrk > 1, p_{T} > 7 GeV);Jet Charge;Counts", 100, -1.5, 1.5);
    
    // Constituent Histograms
    TH1F* h_TrkPt  = new TH1F("h_TrkPt", "Constituent Particle p_{T};Particle p_{T} [GeV];Counts", 100, 0, 20);
    TH1F* h_TrkEta = new TH1F("h_TrkEta", "Constituent Particle #eta;Particle #eta;Counts", 100, -4.0, 4.0);
    TH1F* h_TrkQ   = new TH1F("h_TrkQ", "Constituent Particle Charge;Charge;Counts", 5, -2.5, 2.5);
    TH1F* h_dR     = new TH1F("h_dR", "#DeltaR (Jet, Particle);#DeltaR;Counts", 100, 0, 1.0);
    
    // Advanced Physics Histograms
    TH1F* h_Zfrac  = new TH1F("h_Zfrac", "Momentum Fraction (z = p_{T,trk} / p_{T,jet});z;Counts", 100, 0, 1.1);
    TH1F* h_jT     = new TH1F("h_jT", "Transverse Momentum w.r.t Jet Axis (j_{T});j_{T} [GeV];Counts", 100, 0, 5.0);
    TH1F* h_Closure= new TH1F("h_Closure", "Closure Test: Vector Sum p_{T} / Jet p_{T};Ratio;Counts", 100, 0.9, 1.1);

    // Topology Histogram
    TH1F* h_topo = new TH1F("h_topo", "Gen Jet Topologies;;Fraction of Jets", 3, 0, 3);
    const char* binLabels[3] = {"Single Trk", "All Same Charge", "Leading Pt >= 90%"};
    for (int i = 1; i <= 3; i++) h_topo->GetXaxis()->SetBinLabel(i, binLabels[i-1]);

    const float MIN_JET_PT = 7.0;
    const float KAPPA = 0.5; 
    int totalJets = 0;

    // -----------------------------------------------------------------------------
    // 4. Event Loop
    // -----------------------------------------------------------------------------
    int evCount = 0;
    std::cout << "Starting Exact Constituent Loop..." << std::endl;

    while (reader.Next()) {
        if (evCount % 1000 == 0) std::cout << "Processing Event " << evCount << std::endl;
        evCount++;

        for (int iJ = 0; iJ < jetPx.GetSize(); iJ++) {
            TVector3 vJ(jetPx[iJ], jetPy[iJ], jetPz[iJ]);
            
            if (vJ.Pt() <= MIN_JET_PT) continue;

            h_JetPt->Fill(vJ.Pt());
            h_JetEta->Fill(vJ.Eta());

            int nTrk = 0, nPos = 0, nNeg = 0; 
            float sumPt_w = 0, sumQPt_w = 0, maxPt = 0, rawPtSum = 0;
            TVector3 vectorSum(0, 0, 0);

            // Get the start and end indices for this jet's particles
            unsigned int startIdx = partBegin[iJ];
            unsigned int endIdx   = partEnd[iJ];

            // Loop strictly over the exact particles assigned to this jet
            for (unsigned int p = startIdx; p < endIdx; p++) {
                int trkId = partIdx[p]; // The actual index in GeneratedParticles

                // Safety check
                if (trkId < 0 || trkId >= trkPx.GetSize()) continue; 

                TVector3 vT(trkPx[trkId], trkPy[trkId], trkPz[trkId]);
                float q = trkQ[trkId];

                // Fill general constituent histograms
                h_TrkPt->Fill(vT.Pt());
                h_TrkEta->Fill(vT.Eta());
                h_TrkQ->Fill(q);
                h_dR->Fill(vJ.DeltaR(vT));
                
                // Advanced Physics
                h_Zfrac->Fill(vT.Pt() / vJ.Pt());
                h_jT->Fill(vT.Perp(vJ)); // Perpendicular momentum relative to the jet axis

                vectorSum += vT;

                // Process for Jet Charge (Charged particles only)
                if (q != 0) {
                    nTrk++;
                    float wPt = std::pow(vT.Pt(), KAPPA); 
                    
                    sumPt_w += wPt;
                    sumQPt_w += q * wPt;
                    rawPtSum += vT.Pt();
                    
                    if (vT.Pt() > maxPt) maxPt = vT.Pt();
                    if (q > 0) nPos++; else if (q < 0) nNeg++;
                }
            }

            // Fill Closure Test (proving we captured the right particles)
            if (vJ.Pt() > 0) h_Closure->Fill(vectorSum.Pt() / vJ.Pt());

            if (nTrk > 0) {
                totalJets++;
                if (nTrk == 1) {
                    h_topo->Fill(0.5); 
                } else {
                    if (nPos == nTrk || nNeg == nTrk) h_topo->Fill(1.5); 
                    if (maxPt / rawPtSum >= 0.9) h_topo->Fill(2.5); 
                    
                    h_Q->Fill(sumQPt_w / sumPt_w);
                }
            }
        }
    }

    // -----------------------------------------------------------------------------
    // 5. Scale Topologies & Save
    // -----------------------------------------------------------------------------
    if (totalJets > 0) h_topo->Scale(1.0 / totalJets);

    outFile->cd();
    h_JetPt->Write(); h_JetEta->Write();
    h_TrkPt->Write(); h_TrkEta->Write(); h_TrkQ->Write(); h_dR->Write();
    h_Zfrac->Write(); h_jT->Write(); h_Closure->Write();
    h_Q->Write(); h_topo->Write();

    // Draw Sanity Check Canvas
    TCanvas* c1 = new TCanvas("c1", "Jet Closure and Charge", 1000, 400);
    c1->Divide(2, 1);
    
    c1->cd(1);
    h_Closure->SetLineColor(kRed);
    h_Closure->Draw("HIST");
    
    c1->cd(2);
    h_Q->SetLineColor(kBlack);
    h_Q->Draw("HIST");
    c1->SaveAs("ExactJetResults.pdf");

    std::cout << "\nDone! Analyzed " << totalJets << " valid jets." << std::endl;
}