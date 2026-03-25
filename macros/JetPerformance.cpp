#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TVector2.h>
#include <TStyle.h>

#include <cmath>
#include <iostream>

using namespace std;

// ---------------- Histograms ----------------
TH1D *hQgen, *hQreco;
TH1D *hDRgen, *hDRreco;
TH1D *hMultiplicity_reco, *hMultiplicity_gen;
TH1D *hJES, *hJER;
TH2D *hMultVsPt_reco, *hMultVsPt_gen;
TH2D *hChargeVsNtrk_gen, *hChargeVsNtrk_reco;

// ---------------- Kinematics ----------------
double calc_pt(double px, double py) { return sqrt(px*px + py*py); }

double calc_p(double px, double py, double pz)
{ return sqrt(px*px + py*py + pz*pz); }

double calc_eta(double px, double py, double pz)
{
    double p = calc_p(px, py, pz);
    if (p == fabs(pz)) return 0.0;
    return 0.5 * log((p + pz)/(p - pz));
}

double calc_phi(double px, double py)
{ return atan2(py, px); }

// ---------------- Jet Charge ----------------
double computeJetCharge(
    double jet_px, double jet_py, double jet_pz,
    const TTreeReaderArray<float>& trk_px,
    const TTreeReaderArray<float>& trk_py,
    const TTreeReaderArray<float>& trk_pz,
    const TTreeReaderArray<float>& trk_charge,
    TH1D* hDR,
    double pt_cut = 5.0,
    double jet_eta_cut = 1.0
) {
    const double k    = 0.5;   
    const double Rmax = 0.5;

    double jet_pt  = calc_pt(jet_px, jet_py);
    double jet_eta = calc_eta(jet_px, jet_py, jet_pz);

    // ---------------- Cuts: Jet ----------------
    if(jet_pt < pt_cut) return 99.0;                 // Jet pT cut > 5 GeV
    if(fabs(jet_eta) < jet_eta_cut) return 99.0;     // Jet eta cut +-1

    double jet_phi = calc_phi(jet_px, jet_py);

    double numerator = 0.0;
    double pt_sum    = 0.0;
    int nTracks = 0;  // Count charged tracks inside jet

    for (int i = 0; i < trk_px.GetSize(); ++i) {
        if (trk_charge[i] == 0.0f) continue;        // ---------------- Cut: neutral particles
        double trk_pt = calc_pt(trk_px[i], trk_py[i]);
        if (trk_pt < 0.2) continue;                // Track pT cut

        double eta = calc_eta(trk_px[i], trk_py[i], trk_pz[i]);
        double phi = calc_phi(trk_px[i], trk_py[i]);

        double dEta = eta - jet_eta;
        double dPhi = TVector2::Phi_mpi_pi(phi - jet_phi);
        double dR   = sqrt(dEta*dEta + dPhi*dPhi);

        hDR->Fill(dR);

        if (dR < Rmax) {
            numerator += trk_charge[i] * pow(trk_pt, k);
            pt_sum    += trk_pt;
            nTracks++;
        }
    }

    if(nTracks < 2) return 99.0;   // ---------------- Cut: remove one-particle jets

    return numerator / pow(pt_sum, k);
}

// ---------------- Main ----------------
void JetPerformance()
{
    gStyle->SetOptStat(0);

    TFile *file = TFile::Open("D:/EIC/EPIC/RECO/23.12.0/epic_craterlake/DIS/NC/18x275/minQ2=100/eic.root");
    if (!file || file->IsZombie()) {
        cout << "Error opening file\n";
        return;
    }

    TTreeReader reader("events", file);

    // GEN
    TTreeReaderArray<float> GenJet_px(reader,"GeneratedJets.momentum.x");
    TTreeReaderArray<float> GenJet_py(reader,"GeneratedJets.momentum.y");
    TTreeReaderArray<float> GenJet_pz(reader,"GeneratedJets.momentum.z");

    TTreeReaderArray<float> GenTrk_px(reader,"GeneratedParticles.momentum.x");
    TTreeReaderArray<float> GenTrk_py(reader,"GeneratedParticles.momentum.y");
    TTreeReaderArray<float> GenTrk_pz(reader,"GeneratedParticles.momentum.z");
    TTreeReaderArray<float> GenTrk_q(reader,"GeneratedParticles.charge");

    // RECO
    TTreeReaderArray<float> RecoJet_px(reader,"ReconstructedJets.momentum.x");
    TTreeReaderArray<float> RecoJet_py(reader,"ReconstructedJets.momentum.y");
    TTreeReaderArray<float> RecoJet_pz(reader,"ReconstructedJets.momentum.z");

    TTreeReaderArray<float> RecoTrk_px(reader,"ReconstructedParticles.momentum.x");
    TTreeReaderArray<float> RecoTrk_py(reader,"ReconstructedParticles.momentum.y");
    TTreeReaderArray<float> RecoTrk_pz(reader,"ReconstructedParticles.momentum.z");
    TTreeReaderArray<float> RecoTrk_q(reader,"ReconstructedParticles.charge");

    // -------- Histograms --------
    // Triple the bins for jet charge histograms
    hQgen  = new TH1D("hQgen","Gen Jet Charge;Q;Jets",180,-3,3);
    hQreco = new TH1D("hQreco","Reco Jet Charge;Q;Jets",180,-3,3);

    hDRgen  = new TH1D("hDRgen","Gen ΔR;ΔR;Tracks",200,0,0.5);
    hDRreco = new TH1D("hDRreco","Reco ΔR;ΔR;Tracks",200,0,0.5);

    hMultiplicity_gen = new TH1D("hMultiplicity_gen","Gen Particles per Jet",20,0,20);
    hMultiplicity_reco     = new TH1D("hMultiplicity_reco","Reco Particles per Jet",20,0,20);

    hJES = new TH1D("hJES","Jet Energy Scale (Reco/Gen pT)",100,0,2);
    hJER = new TH1D("hJER","Jet Energy Resolution ((Reco-Gen)/Gen)",100,-1,1);

    // ---------------- 2D Histograms ----------------
    hMultVsPt_reco = new TH2D("hMultVsPt_reco","Reco Mult vs pT;Jet pT;N_{Tracks}",100,0,50,20,0,20);
    hMultVsPt_gen  = new TH2D("hMultVsPt_gen","Gen Mult vs pT;Jet pT;N_{Tracks}",100,0,50,20,0,20);

    hChargeVsNtrk_gen  = new TH2D("hChargeVsNtrk_gen","Gen Jet Charge vs Ntrk;N_{Tracks};Jet Charge",20,0,20,180,-3,3);
    hChargeVsNtrk_reco = new TH2D("hChargeVsNtrk_reco","Reco Jet Charge vs Ntrk;N_{Tracks};Jet Charge",20,0,20,180,-3,3);

    // -------- Event Loop --------
    while(reader.Next()){

        // ===== GEN =====
        for(int j=0;j<GenJet_px.GetSize();j++){
            double jet_px = GenJet_px[j];
            double jet_py = GenJet_py[j];
            double jet_pz = GenJet_pz[j];

            double jet_pt  = calc_pt(jet_px,jet_py);
            double jet_eta = calc_eta(jet_px,jet_py,jet_pz);
            double jet_phi = calc_phi(jet_px,jet_py);

            // ---------------- Count number of tracks inside jet (with cuts) ----------------
            int nParticles = 0;
            for(int i=0;i<GenTrk_px.GetSize();i++){
                if(GenTrk_q[i]==0) continue;                   // Cut neutral particles
                double pt = calc_pt(GenTrk_px[i],GenTrk_py[i]);
                if(pt < 0.5) continue;
                double eta = calc_eta(GenTrk_px[i],GenTrk_py[i],GenTrk_pz[i]);
                double phi = calc_phi(GenTrk_px[i],GenTrk_py[i]);
                double dR = sqrt(pow(eta-jet_eta,2)+pow(TVector2::Phi_mpi_pi(phi-jet_phi),2));
                if(dR < 0.5) nParticles++;
            }
            if(nParticles < 3) continue;                       // Cut one-particle jets
            if(fabs(jet_eta)<1) continue;                     // Cut jet eta less than +-1
            if(jet_pt < 7) continue;                           // Jet pt cut > 5 GeV

            hMultiplicity_gen->Fill(nParticles);
            hMultVsPt_gen->Fill(jet_pt,nParticles);

            double Qgen = computeJetCharge(jet_px,jet_py,jet_pz,GenTrk_px,GenTrk_py,GenTrk_pz,GenTrk_q,hDRgen);
            if(Qgen != 99.0){
                hQgen->Fill(Qgen);
                hChargeVsNtrk_gen->Fill(nParticles,Qgen);
            }
        }

        // ===== RECO =====
        for(int j=0;j<RecoJet_px.GetSize();j++){
            double jet_px = RecoJet_px[j];
            double jet_py = RecoJet_py[j];
            double jet_pz = RecoJet_pz[j];

            double jet_pt  = calc_pt(jet_px,jet_py);
            double jet_eta = calc_eta(jet_px,jet_py,jet_pz);
            double jet_phi = calc_phi(jet_px,jet_py);

            int nParticles = 0;
            for(int i=0;i<RecoTrk_px.GetSize();i++){
                if(RecoTrk_q[i]==0) continue;                   // Cut neutral particles
                double pt = calc_pt(RecoTrk_px[i],RecoTrk_py[i]);
                if(pt < 0.2) continue;
                double eta = calc_eta(RecoTrk_px[i],RecoTrk_py[i],RecoTrk_pz[i]);
                double phi = calc_phi(RecoTrk_px[i],RecoTrk_py[i]);
                double dR = sqrt(pow(eta-jet_eta,2)+pow(TVector2::Phi_mpi_pi(phi-jet_phi),2));
                if(dR < 0.5) nParticles++;
            }
            if(nParticles < 2) continue;                       // Cut one-particle jets
            if(fabs(jet_eta)<1) continue;                     // Cut jet eta less than +-1
            if(jet_pt < 5) continue;                           // Jet pt cut > 5 GeV

            hMultiplicity_reco->Fill(nParticles);
            hMultVsPt_reco->Fill(jet_pt,nParticles);

            double Qreco = computeJetCharge(jet_px,jet_py,jet_pz,RecoTrk_px,RecoTrk_py,RecoTrk_pz,RecoTrk_q,hDRreco);
            if(Qreco != 99.0){
                hQreco->Fill(Qreco);
                hChargeVsNtrk_reco->Fill(nParticles,Qreco);
            }

            // ===== Jet scale and resolution (match closest GEN jet) =====
            double best_dR = 0.5;
            double matched_Gen_pt = -1.0;
            for(int gj=0;gj<GenJet_px.GetSize();gj++){
                double gen_px = GenJet_px[gj];
                double gen_py = GenJet_py[gj];
                double gen_pz = GenJet_pz[gj];
                double gen_pt = calc_pt(gen_px,gen_py);
                double gen_eta = calc_eta(gen_px,gen_py,gen_pz);
                double gen_phi = calc_phi(gen_px,gen_py);

                double dR = sqrt(pow(jet_eta-gen_eta,2)+pow(TVector2::Phi_mpi_pi(jet_phi-gen_phi),2));
                if(dR < best_dR){
                    best_dR = dR;
                    matched_Gen_pt = gen_pt;
                }
            }

            if(matched_Gen_pt > 0){
                double JES = jet_pt / matched_Gen_pt;
                double JER = (jet_pt - matched_Gen_pt) / matched_Gen_pt;
                hJES->Fill(JES);
                hJER->Fill(JER);
            }
        }
    }

    // ---------------- Save all histograms ----------------
    TFile outFile("JetPerformance.root","RECREATE");
    hQgen->Write();
    hQreco->Write();
    hDRgen->Write();
    hDRreco->Write();
    hMultiplicity_gen->Write();
    hMultiplicity_reco->Write();
    hMultVsPt_gen->Write();
    hMultVsPt_reco->Write();
    hChargeVsNtrk_gen->Write();
    hChargeVsNtrk_reco->Write();
    hJES->Write();
    hJER->Write();
    outFile.Close();

    cout<<"Done. Histograms saved to JetPerformance.root"<<endl;
}