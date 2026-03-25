#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
#include <TString.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <TStopwatch.h>

float DeltaPhi(float phi1, float phi2) {
    float dPhi = phi1 - phi2;
    // TMath::Pi() is pi (3.14159...)
    // TMath::TwoPi() is 2*pi (6.28318...)
    while (dPhi > TMath::Pi()) dPhi -= TMath::TwoPi();
    while (dPhi <= -TMath::Pi()) dPhi += TMath::TwoPi();
    return dPhi;
}

float DeltaR(float eta1, float phi1, float eta2, float phi2) {
    float dPhi = DeltaPhi(phi1, phi2); // Now uses the dedicated wrapping function
    float dEta = eta1 - eta2;
    return std::sqrt(dEta * dEta + dPhi * dPhi);
}

void JetPerformance(const char* fileName, TString outputFile, bool doElecremoval, bool doSingleTrkJetremoval) {

    gInterpreter->GenerateDictionary("vector<vector<float>>", "vector");
    gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
    gInterpreter->GenerateDictionary("vector<vector<bool>>", "vector");

    // 1. Setup File and Tree
    TFile* f = TFile::Open(fileName);
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open file " << fileName << std::endl;
        return;
    }
    TTree* tree = (TTree*)f->Get("JetTree_R1p0");
    if (!tree) {
        std::cerr << "Error: Could not find TTree 'JetTree' in file." << std::endl;
        f->Close();
        return;
    }


    float EventQ2 = 0.0;
    float EventX = 0.0;
 
	// 2a. Pointers for RecoJet branches (all vector<float> except hasElectron)
    std::vector<float>* RecoJet_pt = nullptr;
    std::vector<float>* RecoJet_eta = nullptr;
    std::vector<float>* RecoJet_phi = nullptr;
    std::vector<float>* RecoJet_E = nullptr;
    std::vector<float>* RecoJet_M = nullptr;
    std::vector<bool>* RecoJet_hasElectron = nullptr;        // Note the type: vector<bool>
    std::vector<float>* RecoJet_maxPtPart_pt = nullptr;

	// RecoJet constituents
    std::vector<std::vector<float>> *RecoJet_const_pt = nullptr;
    std::vector<std::vector<float>> *RecoJet_const_eta = nullptr;
    std::vector<std::vector<float>> *RecoJet_const_phi = nullptr;
    std::vector<std::vector<int>>   *RecoJet_const_nhits = nullptr;
    std::vector<std::vector<int>> 	*RecoJet_const_pdgid = nullptr;
    std::vector<std::vector<int>> 	*RecoJet_const_pdgidTruth = nullptr;
    
    // 2b. Pointers for GenJet branches (all vector<float> except hasElectron/hasNeutral)
    std::vector<float>* GenJet_pt = nullptr;
    std::vector<float>* GenJet_eta = nullptr;
    std::vector<float>* GenJet_phi = nullptr;
    std::vector<float>* GenJet_E = nullptr;
    std::vector<float>* GenJet_M = nullptr;
    std::vector<bool>* GenJet_hasElectron = nullptr;         // Note the type: vector<bool>
    std::vector<bool>* GenJet_hasNeutral = nullptr;          // Note the type: vector<bool>
    std::vector<float>* GenJet_maxPtPart_pt = nullptr; 

    // GenJet constituents
    std::vector<std::vector<float>> *GenJet_const_pt = nullptr;
    std::vector<std::vector<float>> *GenJet_const_eta = nullptr;
    std::vector<std::vector<float>> *GenJet_const_phi = nullptr;
    std::vector<std::vector<int>> *GenJet_const_pdgid = nullptr;

    tree->SetBranchAddress("EventQ2", &EventQ2);
    tree->SetBranchAddress("Eventx", &EventX);
	// 2d. Link RecoJet branches
    tree->SetBranchAddress("RecoJet_pt", &RecoJet_pt);
    tree->SetBranchAddress("RecoJet_eta", &RecoJet_eta);
    tree->SetBranchAddress("RecoJet_phi", &RecoJet_phi);
    tree->SetBranchAddress("RecoJet_E", &RecoJet_E);
    tree->SetBranchAddress("RecoJet_M", &RecoJet_M);
    tree->SetBranchAddress("RecoJet_hasElectron", &RecoJet_hasElectron);
    tree->SetBranchAddress("RecoJet_maxPtPart_pt", &RecoJet_maxPtPart_pt);

    // Reco constituents
    tree->SetBranchAddress("RecoJet_constituent_pt",   &RecoJet_const_pt);
    tree->SetBranchAddress("RecoJet_constituent_eta",  &RecoJet_const_eta);
    tree->SetBranchAddress("RecoJet_constituent_phi",  &RecoJet_const_phi);
    tree->SetBranchAddress("RecoJet_constituent_nhits",&RecoJet_const_nhits);
    tree->SetBranchAddress("RecoJet_constituent_pdgid",&RecoJet_const_pdgid);
    tree->SetBranchAddress("RecoJet_constituent_pdgidTruth",&RecoJet_const_pdgidTruth);

    // 2e. Link GenJet branches
    tree->SetBranchAddress("GenJet_pt", &GenJet_pt);
    tree->SetBranchAddress("GenJet_eta", &GenJet_eta);
    tree->SetBranchAddress("GenJet_phi", &GenJet_phi);
    tree->SetBranchAddress("GenJet_E", &GenJet_E);
    tree->SetBranchAddress("GenJet_M", &GenJet_M);
    tree->SetBranchAddress("GenJet_hasElectron", &GenJet_hasElectron);
    tree->SetBranchAddress("GenJet_hasNeutral", &GenJet_hasNeutral);
    tree->SetBranchAddress("GenJet_maxPtPart_pt", &GenJet_maxPtPart_pt);

    // Gen constituents
    tree->SetBranchAddress("GenJet_constituent_pt",  &GenJet_const_pt);
    tree->SetBranchAddress("GenJet_constituent_eta", &GenJet_const_eta);
    tree->SetBranchAddress("GenJet_constituent_phi", &GenJet_const_phi);
    tree->SetBranchAddress("GenJet_constituent_pdgid", &GenJet_const_pdgid);
 
 
/*   // Define logarithmic bins
    const int nBins = 70;             // number of bins
    double xMin = 1e-5;               // minimum x
    double xMax = 1.0;                // maximum x
    double logMin = log10(xMin);
    double logMax = log10(xMax);

    double binEdges[nBins+1];
    for (int i = 0; i <= nBins; ++i) {
        binEdges[i] = pow(10, logMin + i * (logMax - logMin)/nBins);
    }
*/

    double xMin = 1e-5;
    double xMax = 1.0;
    double binWidth = 0.001;
    int nBins = int((xMax - xMin)/binWidth);

    TH1D* h_x = new TH1D("h_x_EPIC", "Bjorken-x (EPIC); x; Events", nBins, xMin, xMax);

    // Create histogram
    //TH1D* h_x = new TH1D("h_x", "Bjorken-x distribution; x; Events", nBins, binEdges);

    
    // Define Histograms
	TH1D *NEvents = new TH1D("NEvents","",2,0.,2.); NEvents->Sumw2();
	TH1D *numRecoJetsEventHist = new TH1D("numRecoJetsEvent","",100,0.,100.); numRecoJetsEventHist->Sumw2();
	TH1D *numRecoJetsEventHistAftCut = new TH1D("numRecoJetsEventHistAftCut","",100,0.,100.); numRecoJetsEventHistAftCut->Sumw2();
	TH1D *numGenJetsEventHist = new TH1D("numGenJetsEventHist","",100,0.,100.); numGenJetsEventHist->Sumw2();
	TH1D *numGenJetsEventHistAftCut = new TH1D("numGenJetsEventHistAftCut","",100,0.,100.); numGenJetsEventHistAftCut->Sumw2();

	// matching
	TH1D *JetdR_all = new TH1D("JetdR_all","",100,0.,10.); JetdR_all->Sumw2();
	TH1D *JetdR_closest = new TH1D("JetdR_closest","",100,0.,10.); JetdR_closest->Sumw2();
	TH1D *JetMindR_closest = new TH1D("JetMindR_closest","",100,0.,10.); JetMindR_closest->Sumw2();

	TH3D *JetMultiplicityReco = new TH3D("JetMultiplicityReco","",100,0.,100.,400,0.0,200.0,100,-5.0,5.0); JetMultiplicityReco->Sumw2();
	TH3D *JetMultiplicityRecoAftCut = new TH3D("JetMultiplicityRecoAftCut","",100,0.,100.,400,0.0,200.0,100,-5.0,5.0); JetMultiplicityRecoAftCut->Sumw2();
	TH3D *JetMultiplicityGen = new TH3D("JetMultiplicityGen","",100,0.,100.,400,0.0,200.0,100,-5.0,5.0); JetMultiplicityGen->Sumw2();
	TH3D *JetMultiplicityGenAftCut = new TH3D("JetMultiplicityGenAftCut","",100,0.,100.,400,0.0,200.0,100,-5.0,5.0); JetMultiplicityGenAftCut->Sumw2();

	// Define multidimensional histograms
	// Jets -> {pT, eta, phi, mass, energy}
	const int NJetAxis = 5;
	int	JetBins[NJetAxis]      =   { 400  ,  100 ,   64		  , 1000  , 2000  };
	double JetXmin[NJetAxis]   =   { 0.0  , -5.0 ,   -3.2 , 0.0 , 0.0  };
	double JetXmax[NJetAxis]   =   { 200.0 ,  5.0 ,  3.2 , 100.0, 1000.0};
	// all
	THnSparseD *mHistJetReco = new THnSparseD("mHistJetReco", "mHistJetReco", NJetAxis, JetBins, JetXmin, JetXmax); mHistJetReco->Sumw2();
	THnSparseD *mHistJetMatch = new THnSparseD("mHistJetMatch", "mHistJetMatch", NJetAxis, JetBins, JetXmin, JetXmax); mHistJetMatch->Sumw2();
	THnSparseD *mHistJetUnMatch = new THnSparseD("mHistJetUnMatch", "mHistJetUnMatch", NJetAxis, JetBins, JetXmin, JetXmax); mHistJetUnMatch->Sumw2();
	THnSparseD *mHistJetGen = new THnSparseD("mHistJetGen", "mHistJetGen", NJetAxis, JetBins, JetXmin, JetXmax); mHistJetGen->Sumw2();

    TH1D* h_GenJet_maxPtRatio = new TH1D("h_GenJet_maxPtRatio", "Leading Trk p_{T} / p_{T,jet};", 1000, 0.0, 10.0); h_GenJet_maxPtRatio->Sumw2();
    TH1D* h_RecoJet_maxPtRatio = new TH1D("h_RecoJet_maxPtRatio", "Leading Trk p_{T} / p_{T,jet};", 1000, 0.0, 10.0); h_RecoJet_maxPtRatio->Sumw2();
    TH1D* h_GenJet_maxPtRatio_AftCut = new TH1D("h_GenJet_maxPtRatio_AftCut", "Leading Trk p_{T} / p_{T,jet};", 1000, 0.0, 10.0); h_GenJet_maxPtRatio_AftCut->Sumw2();
    TH1D* h_RecoJet_maxPtRatio_AftCut = new TH1D("h_RecoJet_maxPtRatio_AftCut", "Leading Trk p_{T} / p_{T,jet};", 1000, 0.0, 10.0); h_RecoJet_maxPtRatio_AftCut->Sumw2();


	// JER/JES -> {..., E, eta }
	const int NJESJERAxis = 3;
	int	JESJERBins[NJESJERAxis]      =   { 2000  , 1000 , 100 };
	double JESJERXmin[NJESJERAxis]   =   { -1.0  , 0.0 , -5.0 };
	double JESJERXmax[NJESJERAxis]   =   { 1.0  , 500.0 , 5.0};
	THnSparseD *mHistJESJERvsE = new THnSparseD("mHistJESJERvsE", "mHistJESJERvsE", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsE->Sumw2();
	THnSparseD *mHistJESJERvsE_DR = new THnSparseD("mHistJESJERvsE_DR", "mHistJESJERvsE_DR", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsE_DR->Sumw2();
//	THnSparseD *mHistJESJERvsE_ratio = new THnSparseD("mHistJESJERvsE_ratio", "mHistJESJERvsE_ratio", NJESJERAxis, JESJERBins, JESJERXmin, JESJERXmax); mHistJESJERvsE_ratio->Sumw2();
	int	JESJERBinsAng[NJESJERAxis]      =   { 20000  , 1000 , 100 };
	double JESJERXminAng[NJESJERAxis]   =   { -0.5  , 0.0 , -5.0 };
	double JESJERXmaxAng[NJESJERAxis]   =   { 0.5  , 500.0 , 5.0};
	THnSparseD *mHistJESJERvsE_DEta = new THnSparseD("mHistJESJERvsE_DEta", "mHistJESJERvsE_DEta", NJESJERAxis, JESJERBinsAng, JESJERXminAng, JESJERXmaxAng); mHistJESJERvsE_DEta->Sumw2();
	THnSparseD *mHistJESJERvsE_DPhi = new THnSparseD("mHistJESJERvsE_DPhi", "mHistJESJERvsE_DPhi", NJESJERAxis, JESJERBinsAng, JESJERXminAng, JESJERXmaxAng); mHistJESJERvsE_DPhi->Sumw2();


	TH2D* hFrag2D = new TH2D("hFrag2D", "Jet Fragmentation; z; j_{T} [GeV]", 50, 0.0, 1.0, 50, 0.0, 5.0);
	TH2D* hFrag2DPi = new TH2D("hFrag2DPi", "Jet Fragmentation; z; j_{T} [GeV]", 50, 0.0, 1.0, 50, 0.0, 5.0);
	TH2D* hFrag2DK = new TH2D("hFrag2DK", "Jet Fragmentation; z; j_{T} [GeV]", 50, 0.0, 1.0, 50, 0.0, 5.0);
	TH2D* hFrag2DP = new TH2D("hFrag2DP", "Jet Fragmentation; z; j_{T} [GeV]", 50, 0.0, 1.0, 50, 0.0, 5.0);
	TH2D* hFrag2DPiTruth = new TH2D("hFrag2DPiTruth", "Jet Fragmentation; z; j_{T} [GeV]", 50, 0.0, 1.0, 50, 0.0, 5.0);
	TH2D* hFrag2DKTruth = new TH2D("hFrag2DKTruth", "Jet Fragmentation; z; j_{T} [GeV]", 50, 0.0, 1.0, 50, 0.0, 5.0);
	TH2D* hFrag2DPTruth = new TH2D("hFrag2DPTruth", "Jet Fragmentation; z; j_{T} [GeV]", 50, 0.0, 1.0, 50, 0.0, 5.0);

	// 4. Loop over Events
    Long64_t nEntries = tree->GetEntries();
    std::cout << "Processing " << nEntries << " events..." << std::endl;
	TStopwatch sw;
	sw.Start();
	
    for (Long64_t i = 0; i < nEntries; i++) {
    
        tree->GetEntry(i);
        NEvents->Fill(1);
		std::vector<bool> reco_matched(RecoJet_pt->size(), false);

 	   if (i % (nEntries/20) == 0) {
	       double elapsed = sw.RealTime();
	       sw.Continue();
 	       double pct = double(i) / nEntries;
 	       double eta = elapsed / pct - elapsed;
  	       std::cout << Form("Progress: %.1f%%  |  ETA: %.1f s", pct*100, eta) << std::endl;
  	       
 	   }

		
		// 4a. GenJet Loop (Efficiency, Response, and GenJet Ratio)
		int ngjets = 0;
		int ngjetscut = 0;
        for (size_t j = 0; j < GenJet_pt->size(); j++) {
            float gen_pt = GenJet_pt->at(j);
            float gen_eta = GenJet_eta->at(j);
            float gen_phi = GenJet_phi->at(j);
            float gen_e = GenJet_E->at(j);
            if (gen_pt <= 1.0) continue; 
            if (doElecremoval) { if(GenJet_hasElectron->at(j)) continue; }
            // Fill GenJet Fragmentation Ratio
            float gen_max_pt = GenJet_maxPtPart_pt->at(j);
            float gen_ratio = gen_max_pt / gen_pt;
            // Loop over constituents
            const auto &ptv  = GenJet_const_pt->at(j);
            const auto &etav = GenJet_const_eta->at(j);
            const auto &phiv = GenJet_const_phi->at(j);
           
            for (size_t k = 0; k < ptv.size(); ++k) {
                float cpt  = ptv[k];
                float ceta = etav[k];
                float cphi = phiv[k];
            }            

            h_GenJet_maxPtRatio->Fill(gen_ratio);
            JetMultiplicityGen->Fill(ptv.size(), gen_pt, gen_eta);
   			if ( ptv.size() > 1 ) h_GenJet_maxPtRatio_AftCut->Fill(gen_ratio);           
            if ( doSingleTrkJetremoval && ptv.size() <= 1 ) continue;
            JetMultiplicityGenAftCut->Fill(ptv.size(), gen_pt, gen_eta);
            ngjets = ngjets + 1;
            if ( gen_pt > 3.0 && gen_eta > -2.0 && gen_eta < 2.5 ) ngjetscut = ngjetscut + 1;

		    double GenJet[5] = {GenJet_pt->at(j), GenJet_eta->at(j), GenJet_phi->at(j), GenJet_M->at(j), GenJet_E->at(j)};
	    	mHistJetGen->Fill(GenJet);

            // Jet Matching logic remains the same)
            float min_dR = 1e6;
            int best_reco_idx = -1;
            for (size_t k = 0; k < RecoJet_pt->size(); k++) {
                float reco_pt = RecoJet_pt->at(k);
                if (reco_pt <= 1.0) continue; 
	            float reco_max_pt = RecoJet_maxPtPart_pt->at(k);
    	        float reco_ratio = reco_max_pt / reco_pt;
                const auto &ptvk  = RecoJet_const_pt->at(k);
	            if ( doSingleTrkJetremoval && ptvk.size() <= 1 ) continue;
                float dR = DeltaR(GenJet_eta->at(j), GenJet_phi->at(j), RecoJet_eta->at(k), RecoJet_phi->at(k));
                JetdR_all->Fill(dR);
                if (dR < min_dR) { min_dR = dR; best_reco_idx = k; }
            }
                            
            if (best_reco_idx >= 0) {
                mHistJetMatch->Fill(GenJet);
                float reco_E_match = RecoJet_E->at(best_reco_idx);
                float reco_Eta_match = RecoJet_eta->at(best_reco_idx);
                float reco_Phi_match = RecoJet_phi->at(best_reco_idx);
                float deltaphi = DeltaPhi(reco_Phi_match, gen_phi);
                float deltaeta = (reco_Eta_match - gen_eta);
                float deltaR = DeltaR(gen_eta, gen_phi, reco_Eta_match, reco_Phi_match);
				JetdR_closest->Fill(deltaR);
				JetMindR_closest->Fill(min_dR);
				if( min_dR < 1.0 ){ // to remove dR > jet radii      
					reco_matched[best_reco_idx] = true;   	
					double JESJER_E[3] = {(reco_E_match - gen_e)/gen_e, gen_e, reco_Eta_match};
					if( gen_pt > 3.0 ) mHistJESJERvsE->Fill(JESJER_E);
					double JESJER_E_ratio[3] = {reco_E_match/gen_e, gen_e, reco_Eta_match};
//					mHistJESJERvsE_ratio->Fill(JESJER_E_ratio);
					double JESJER_DPhi[3] = {deltaphi, gen_e, reco_Eta_match};
					if( gen_pt > 3.0 ) mHistJESJERvsE_DPhi->Fill(JESJER_DPhi);
					double JESJER_DEta[3] = {deltaeta, gen_e, reco_Eta_match};
					if( gen_pt > 3.0 ) mHistJESJERvsE_DEta->Fill(JESJER_DEta);
					double JESJER_DR[3] = {deltaR, gen_e, reco_Eta_match};
					if( gen_pt > 3.0 ) mHistJESJERvsE_DR->Fill(JESJER_DR);
				}
            }
         } // gen loop
         numGenJetsEventHist->Fill(ngjets);       
         numGenJetsEventHistAftCut->Fill(ngjetscut); 
        // 4b. RecoJet Loop (Fake Rate and RecoJet Ratio)
		int nrjets = 0;
		int nrjetscut = 0;		
        for (size_t k = 0; k < RecoJet_pt->size(); k++) {
            float reco_pt = RecoJet_pt->at(k);
            if (reco_pt <= 1.0) continue; 
            if (doElecremoval) { if(RecoJet_hasElectron->at(k)) continue; }
            // Fill RecoJet Fragmentation Ratio
            float reco_max_pt = RecoJet_maxPtPart_pt->at(k);
            float reco_ratio = reco_max_pt / reco_pt;
            // Loop over constituents         
            const auto &ptv  = RecoJet_const_pt->at(k);
            const auto &etav = RecoJet_const_eta->at(k);
            const auto &phiv = RecoJet_const_phi->at(k);
            const auto &nhitsv = RecoJet_const_nhits->at(k);            
            const auto &pidv = RecoJet_const_pdgid->at(k);            
            const auto &pidTv = RecoJet_const_pdgidTruth->at(k);            

            h_RecoJet_maxPtRatio->Fill(reco_ratio);
            JetMultiplicityReco->Fill(ptv.size(), reco_pt, RecoJet_eta->at(k));
			if ( ptv.size() > 1 ) h_RecoJet_maxPtRatio_AftCut->Fill(reco_ratio);       
            if ( doSingleTrkJetremoval && ptv.size() <= 1 ) continue;
            JetMultiplicityRecoAftCut->Fill(ptv.size(), reco_pt, RecoJet_eta->at(k));
           	nrjets = nrjets + 1;
           	if(reco_pt > 3.0 && RecoJet_eta->at(k) > -2.0 && RecoJet_eta->at(k) < 2.5) nrjetscut = nrjetscut + 1;

	    	double RecoJet[5] = {RecoJet_pt->at(k), RecoJet_eta->at(k), RecoJet_phi->at(k), RecoJet_M->at(k), RecoJet_E->at(k)};
	    	mHistJetReco->Fill(RecoJet); // all jets

            if (!reco_matched[k]) mHistJetUnMatch->Fill(RecoJet);

    		// =======================
   			// FRAGMENTATION OBSERVABLES
    		// =======================
   			TVector3 jetVec;
    		jetVec.SetPtEtaPhi(reco_pt, RecoJet_eta->at(k), RecoJet_phi->at(k));
    		double jetMag2 = jetVec.Mag2();
    		TVector3 jetUnit = jetVec.Unit();
			if(reco_pt > 3.0 && RecoJet_eta->at(k) > -2.0 && RecoJet_eta->at(k) < 2.5) {
			
	    		for (size_t c = 0; c < ptv.size(); c++) {
    	    		TVector3 constVec;
        			constVec.SetPtEtaPhi(ptv[c], etav[c], phiv[c]);
	
    	    		// longitudinal fraction
	        		double z = constVec.Dot(jetUnit) / jetVec.Mag();
	   	     		if (z <= 0.0 || z > 1.0) continue;	

    	    		// transverse momentum wrt jet axis
	       	 		double jT = (constVec - (constVec.Dot(jetUnit) * jetUnit)).Mag();

    	    		// fill histograms (define these before loop)
					hFrag2D->Fill(z, jT);
					if(fabs(pidv[c]) == 211) hFrag2DPi->Fill(z, jT);
					if(fabs(pidv[c]) == 321) hFrag2DK->Fill(z, jT);
					if(fabs(pidv[c]) == 2212) hFrag2DP->Fill(z, jT);
					if(fabs(pidTv[c]) == 211) hFrag2DPiTruth->Fill(z, jT);
					if(fabs(pidTv[c]) == 321) hFrag2DKTruth->Fill(z, jT);
					if(fabs(pidTv[c]) == 2212) hFrag2DPTruth->Fill(z, jT);

    			}
    			
    			
    			if(EventX >= xMin && EventX <= xMax) h_x->Fill(EventX);	
    		}


        }
        numRecoJetsEventHist->Fill(nrjets);
        numRecoJetsEventHistAftCut->Fill(nrjetscut);        
    } // event loop
    
    // 5. Calculate Performance Metrics (Efficiency, Fake Rate, JES/JER) - (Same as before)
    TFile* outFile = TFile::Open(Form("results/%s",outputFile.Data()), "RECREATE");
	NEvents->Write();
	numRecoJetsEventHist->Write();
	numRecoJetsEventHistAftCut->Write();
	numGenJetsEventHist->Write();
	numGenJetsEventHistAftCut->Write();
	JetdR_all->Write();
	JetdR_closest->Write();
	JetMindR_closest->Write();
	JetMultiplicityReco->Write();
	JetMultiplicityRecoAftCut->Write();
	JetMultiplicityGen->Write();
	JetMultiplicityGenAftCut->Write();
	h_GenJet_maxPtRatio->Write();
    h_RecoJet_maxPtRatio->Write();
    h_GenJet_maxPtRatio_AftCut->Write();
    h_RecoJet_maxPtRatio_AftCut->Write();
 	mHistJetReco->Write();
	mHistJetMatch->Write();
	mHistJetUnMatch->Write();
	mHistJetGen->Write();
 	mHistJESJERvsE->Write();
// 	mHistJESJERvsE_ratio->Write();
 	mHistJESJERvsE_DEta->Write();
 	mHistJESJERvsE_DPhi->Write();
 	mHistJESJERvsE_DR->Write();
	hFrag2D->Write();
	hFrag2DPi->Write();
	hFrag2DK->Write();
	hFrag2DP->Write();
	hFrag2DPiTruth->Write();
	hFrag2DKTruth->Write();
	hFrag2DPTruth->Write();
	h_x->Write();

    outFile->Close();    

    sw.Stop();

    // --- OUTPUT MESSAGE ---
    std::cout << "===============================" << std::endl;
    std::cout << "Finished processing all events!" << std::endl;
    std::cout << "Total elapsed time: " << sw.RealTime() << " s" << std::endl;
    std::cout << "===============================" << std::endl;
    

}