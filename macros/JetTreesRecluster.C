#include "TreeReader.h"  // call libraries from ROOT and C++

using namespace fastjet;

void JetTreesRecluster(TString InputFileList, TString OutputFile, std::vector<float> R_values, int removeelectrons, int nhitcut){

	typedef ROOT::Math::PxPyPzEVector LorentzVector;
	
    double minCstPt            = 0.2 ;				 // minimum pT of objects
    double maxCstPt            = 100.;  			 // maximum pT of objects
    double minJetPt            = 1.0 ;  			 // minimum jet pT
    double ghostMaxRap         = 3.5;   			 // maximum rapidity of ghosts
    double ghostArea           = 0.01; 				 // area per ghost
    int numGhostRepeat         = 1;                  // reuse count of ghosts


	// Read the list of input file(s)
	fstream FileList;
	FileList.open(Form("%s",InputFileList.Data()), ios::in);
	if(!FileList.is_open()){cout << "List of input files not founded!" << endl; return;}{cout << "List of input files founded! --> " << InputFileList.Data() << endl;}

	// Make a chain and a vector of file names
	std::vector<TString> FileListVector;
	string FileChain;
	while(getline(FileList, FileChain)){FileListVector.push_back(FileChain.c_str());}
	FileList.close();	
	TChain *mychain = new TChain("events");
	for (std::vector<TString>::iterator listIterator = FileListVector.begin(); listIterator != FileListVector.end(); listIterator++){
		TFile *testfile = TFile::Open(*listIterator,"READ");
		if(testfile && !testfile->IsZombie() && !testfile->TestBit(TFile::kRecovered)){ // safety against corrupted files
			cout << "Adding file " << *listIterator << " to the chains" << endl; // adding files to the chains for each step
			mychain->Add(*listIterator);
		}else{cout << "File: " << *listIterator << " failed!" << endl;}
	}

	// Reading trees
	std::unique_ptr<TTreeReader> tree_reader;
	TreeReader(mychain, tree_reader); 
	
	// Make Output
	TFile *OutFile = TFile::Open(Form("%s",OutputFile.Data()),"RECREATE");	

	int NEVENTS = 0;
	int EVETMULTRECO = 0; 
	int EVETMULTGEN = 0;
	int ScatteredERecId = 0;
	int ScatteredEGenId = 0;
	float EventQ2 = -999;
	float Eventx = -999;
	float EventQ2Gen = -999;
	float EventxGen = -999;

	// Vertex
	std::vector<float> Vertex_x;
	std::vector<float> Vertex_y;
	std::vector<float> Vertex_z;
	std::vector<int> Vertex_ndf;
	std::vector<float> Vertex_chi2;
	std::vector<float> VertexErr_xx;
	std::vector<float> VertexErr_yy;
	std::vector<float> VertexErr_zz;

	// Reco Jets (Variable-length vectors for multiple jets per event)
	std::vector<float> RecoJet_pt;
	std::vector<float> RecoJet_eta;
	std::vector<float> RecoJet_phi;
	std::vector<float> RecoJet_E;
	std::vector<float> RecoJet_M;
	std::vector<bool> RecoJet_hasElectron;
	std::vector<float> RecoJet_maxPtPart_pt; 
	// Reco jet constituents
	std::vector<std::vector<float>> RecoJet_constituent_pt;
	std::vector<std::vector<float>> RecoJet_constituent_eta;
	std::vector<std::vector<float>> RecoJet_constituent_phi;
	std::vector<std::vector<int>> RecoJet_constituent_nhits;
	std::vector<std::vector<int>> RecoJet_constituent_pdgid;
	std::vector<std::vector<int>> RecoJet_constituent_pdgidTruth;

	// Gen Jets (Variable-length vectors for multiple jets per event)
	std::vector<float> GenJet_pt;
	std::vector<float> GenJet_eta;
	std::vector<float> GenJet_phi;
	std::vector<float> GenJet_E;
	std::vector<float> GenJet_M;
	std::vector<bool> GenJet_hasElectron;
	std::vector<bool> GenJet_hasNeutral;
	std::vector<float> GenJet_maxPtPart_pt;
	// Gen jet constituents
	std::vector<std::vector<float>> GenJet_constituent_pt;
	std::vector<std::vector<float>> GenJet_constituent_eta;
	std::vector<std::vector<float>> GenJet_constituent_phi;
	std::vector<std::vector<int>> GenJet_constituent_pdgid;
	
	
    // === Create trees ===
    std::vector<TTree*> trees;
    for (auto R : R_values) {

		TString rStr;
		if (R < 1.0) { rStr = Form("JetTree_R0p%d", int(R * 10));   // R=0.1 -> R0p1
		} else { rStr = Form("JetTree_R%dp%d", int(R), int((R - int(R)) * 10));        // R=1.2 -> R1p2
		}
		
        TTree *tree = new TTree(rStr, rStr);
		// Event information
        tree->Branch("NEVENTS", &NEVENTS, "NEVENTS/I");
		tree->Branch("EVETMULTRECO", &EVETMULTRECO, "EVETMULTRECO/I");
		tree->Branch("EVETMULTGEN", &EVETMULTGEN, "EVETMULTGEN/I");
		tree->Branch("EventQ2", &EventQ2, "EventQ2/F");
		tree->Branch("Eventx", &Eventx, "Eventx/F");
		tree->Branch("EventQ2Gen", &EventQ2Gen, "EventQ2Gen/F");
		tree->Branch("EventxGen", &EventxGen, "EventxGen/F");
        // Vertex
        tree->Branch("Vertex_x", &Vertex_x);
        tree->Branch("Vertex_y", &Vertex_y);
        tree->Branch("Vertex_z", &Vertex_z);
        tree->Branch("Vertex_ndf", &Vertex_ndf);
        tree->Branch("Vertex_chi2", &Vertex_chi2);
        tree->Branch("VertexErr_xx", &VertexErr_xx);
        tree->Branch("VertexErr_yy", &VertexErr_yy);
        tree->Branch("VertexErr_zz", &VertexErr_zz);
		// Scattered electron
        tree->Branch("ScatteredERecId", &ScatteredERecId);
        tree->Branch("ScatteredEGenId", &ScatteredEGenId);
        // Reco
        tree->Branch("RecoJet_pt", &RecoJet_pt);
        tree->Branch("RecoJet_eta", &RecoJet_eta);
        tree->Branch("RecoJet_phi", &RecoJet_phi);
        tree->Branch("RecoJet_E", &RecoJet_E);
        tree->Branch("RecoJet_M", &RecoJet_M);
        tree->Branch("RecoJet_hasElectron", &RecoJet_hasElectron);
        tree->Branch("RecoJet_maxPtPart_pt", &RecoJet_maxPtPart_pt);
        tree->Branch("RecoJet_constituent_pt", &RecoJet_constituent_pt);
        tree->Branch("RecoJet_constituent_eta", &RecoJet_constituent_eta);
        tree->Branch("RecoJet_constituent_phi", &RecoJet_constituent_phi);
        tree->Branch("RecoJet_constituent_pdgid", &RecoJet_constituent_pdgid);
        tree->Branch("RecoJet_constituent_pdgidTruth", &RecoJet_constituent_pdgidTruth);
		tree->Branch("RecoJet_constituent_nhits", &RecoJet_constituent_nhits);
        // Gen
        tree->Branch("GenJet_pt", &GenJet_pt);
        tree->Branch("GenJet_eta", &GenJet_eta);
        tree->Branch("GenJet_phi", &GenJet_phi);
        tree->Branch("GenJet_E", &GenJet_E);
        tree->Branch("GenJet_M", &GenJet_M);
        tree->Branch("GenJet_hasElectron", &GenJet_hasElectron);
        tree->Branch("GenJet_hasNeutral", &GenJet_hasNeutral);
        tree->Branch("GenJet_maxPtPart_pt", &GenJet_maxPtPart_pt);
        tree->Branch("GenJet_constituent_pt", &GenJet_constituent_pt);
        tree->Branch("GenJet_constituent_eta", &GenJet_constituent_eta);
        tree->Branch("GenJet_constituent_phi", &GenJet_constituent_phi);
        tree->Branch("GenJet_constituent_pdgid", &GenJet_constituent_pdgid);

        trees.push_back(tree);
    }
	
	// Loop over events	
	int globalEvent = 0;
	while(tree_reader->Next()) {	
	
	    if(globalEvent%100 == 0) cout << "Events Processed: " << globalEvent << endl;
	    globalEvent++;
	    NEVENTS = globalEvent;
	    // For Vertex
	    Vertex_x.clear();
		Vertex_y.clear();
		Vertex_z.clear();
		Vertex_ndf.clear();
		Vertex_chi2.clear();
		VertexErr_xx.clear();
		VertexErr_yy.clear();
		VertexErr_zz.clear();
        for (unsigned int ivtx = 0; ivtx < CTVx->GetSize(); ++ivtx) {
   	    	Vertex_x.push_back((*CTVx)[ivtx]);
   	    	Vertex_y.push_back((*CTVy)[ivtx]);
   	    	Vertex_z.push_back((*CTVz)[ivtx]);
   	    	Vertex_ndf.push_back((*CTVndf)[ivtx]);
   	    	Vertex_chi2.push_back((*CTVchi2)[ivtx]);
   	    	VertexErr_xx.push_back((*CTVerr_xx)[ivtx]);
   	    	VertexErr_yy.push_back((*CTVerr_yy)[ivtx]);
   	    	VertexErr_zz.push_back((*CTVerr_zz)[ivtx]);        
        }
        
		// For Scattered electron
		
		if(ScatElecRecoId->GetSize() > 0){
			int iscatR = (*ScatElecRecoId)[0];
			ScatteredERecId = iscatR;    
		}else{ScatteredERecId = -999;}


		if(ScatElecGenId->GetSize() > 0){
			int iscatS = (*ScatElecGenId)[0];
			ScatteredEGenId = iscatS;    
		}else{ScatteredEGenId = -999;}
		
		// other quantities
		if(EvtQ2->GetSize() > 0){
			float q2valueR = (*EvtQ2)[0];
			EventQ2 = q2valueR;    
		}else{EventQ2 = -999;}
		
		if(Evtx->GetSize() > 0){
			float xvalueR = (*Evtx)[0];
			Eventx = xvalueR;    
		}else{Eventx = -999;}
		
		if(EvtQ2Gen->GetSize() > 0){
			float q2valueG = (*EvtQ2Gen)[0];
			EventQ2Gen = q2valueG;    
		}else{EventQ2Gen = -999;}
		
		if(EvtxGen->GetSize() > 0){
			float xvalueG = (*EvtxGen)[0];
			EventxGen = xvalueG;    
		}else{EventxGen = -999;}


        // Build pseudojets
        std::vector<PseudoJet> particles_reco;
        for (unsigned int i = 0; i < TrkRecoPx->GetSize(); ++i) {
            TVector3 mom((*TrkRecoPx)[i], (*TrkRecoPy)[i], (*TrkRecoPz)[i]);
            if ( mom.Pt() < minCstPt || mom.Pt() > maxCstPt ) continue;
            if ( (*TrkRecoNhits)[i] < nhitcut ) continue;
            if ( removeelectrons == 1 ){
				// Find electron
			    int chargePartIndex = i; 
			    int elecIndex = -1;
			    float elecIndexWeight = -1.0;
		    	for(unsigned int itrkass = 0; itrkass < TrkPartAssocWeight->GetSize(); itrkass++){ // Loop Over All ReconstructedChargedParticleAssociations
					if( TrkPartAssocWeight->GetSize() > 0 ){
						if((*TrkPartAssocRec)[itrkass] == chargePartIndex){ // Select Entry Matching the ReconstructedChargedParticle Index
						    if((*TrkPartAssocWeight)[itrkass] > elecIndexWeight){ // Find Particle with Greatest Weight = Contributed Most Hits to Track
								elecIndex = (*TrkPartAssocSim)[itrkass]; // Get Index of MCParticle Associated with ReconstructedChargedParticle
								elecIndexWeight = (*TrkPartAssocWeight)[itrkass];
			      			}
			  			}
			  		}
		      	}
				if((*TrkMCGenPDG)[elecIndex] == 11) continue;
            }
            if ( removeelectrons == 2 && i == ScatteredERecId ) continue; 
            PseudoJet p((*TrkRecoPx)[i], (*TrkRecoPy)[i], (*TrkRecoPz)[i], (*TrkRecoE)[i]);
            p.set_user_index(i);
            particles_reco.push_back(p);
        }

        std::vector<PseudoJet> particles_gen;
        for (unsigned int i = 0; i < TrkGenPx->GetSize(); ++i) {
            TVector3 mom((*TrkGenPx)[i], (*TrkGenPy)[i], (*TrkGenPz)[i]);
            if ( mom.Pt() < minCstPt || mom.Pt() > maxCstPt ) continue;
            if ( (*TrkGenCharge)[i] == 0 ) continue;
            if ( removeelectrons == 1 && (*TrkGenPDG)[i] == 11 ) continue;
            if ( removeelectrons == 2 && i == ScatteredEGenId ) continue; 
            PseudoJet p((*TrkGenPx)[i], (*TrkGenPy)[i], (*TrkGenPz)[i], (*TrkGenE)[i]);
            p.set_user_index(i);
            particles_gen.push_back(p);
        }

        // --- Loop over R values --- making jets
        for (size_t iR = 0; iR < R_values.size(); ++iR) { 
        
            double R = R_values[iR];
            TTree *tree = trees[iR];
            // Define algorithm
            JetAlgorithm algo = antikt_algorithm;
            RecombinationScheme scheme = E_scheme;
			// Jet definition
            JetDefinition jet_def(algo, R, scheme);
            GhostedAreaSpec ghost_spec(ghostMaxRap, numGhostRepeat, ghostArea);
            AreaType atype = active_area;
            AreaDefinition area_def(atype, ghost_spec);
                        
			// Clear all vectors for the new event
			RecoJet_pt.clear();
			RecoJet_eta.clear();
			RecoJet_phi.clear();
			RecoJet_E.clear();
			RecoJet_M.clear();
			RecoJet_hasElectron.clear();
			RecoJet_maxPtPart_pt.clear();
	        RecoJet_constituent_pt.clear(); 
	        RecoJet_constituent_eta.clear();
	        RecoJet_constituent_phi.clear(); 
	        RecoJet_constituent_nhits.clear();
	        RecoJet_constituent_pdgid.clear();
	        RecoJet_constituent_pdgidTruth.clear();
	        
			GenJet_pt.clear();
			GenJet_eta.clear();
			GenJet_phi.clear();
			GenJet_E.clear();
			GenJet_M.clear();
			GenJet_hasElectron.clear();
			GenJet_hasNeutral.clear();
			GenJet_maxPtPart_pt.clear();
			GenJet_constituent_pt.clear(); 
	        GenJet_constituent_eta.clear();
	        GenJet_constituent_phi.clear(); 
	        GenJet_constituent_pdgid.clear();

			EVETMULTRECO = (int)TrkRecoPx->GetSize();
			EVETMULTGEN = (int)TrkGenPx->GetSize();

            // --- Reco clustering ---
            ClusterSequenceArea cs_reco(particles_reco, jet_def, area_def);
            std::vector<PseudoJet> jets_reco = sorted_by_pt(cs_reco.inclusive_jets(minJetPt));
            for (auto &jet : jets_reco) {
                RecoJet_pt.push_back(jet.pt());
                RecoJet_eta.push_back(jet.eta());
                RecoJet_phi.push_back(jet.phi_std());
                RecoJet_E.push_back(jet.e());
                RecoJet_M.push_back(jet.m());

                bool hasElectron = false;
                float maxPtReco = -1.0;
                std::vector<float> cpt, ceta, cphi;
                std::vector<int> chits, cpdgid, cpdgidtruth;
                cpt.clear(); ceta.clear(); cphi.clear(); chits.clear(); cpdgid.clear(); cpdgidtruth.clear();

                for (auto &c : jet.constituents()) {
                    int idx = c.user_index();
                    TVector3 v3((*TrkRecoPx)[idx], (*TrkRecoPy)[idx], (*TrkRecoPz)[idx]);
                    cpt.push_back(v3.Pt());
                    ceta.push_back(v3.Eta());
                    cphi.push_back(v3.Phi());
                    chits.push_back((*TrkRecoNhits)[idx]);

                    if (v3.Pt() > maxPtReco) maxPtReco = v3.Pt();
					// Find electron
		    		int chargePartIndex = idx; 
		    		int elecIndex = -1;
		    		float elecIndexWeight = -1.0;
	    			for(unsigned int itrkass = 0; itrkass < TrkPartAssocWeight->GetSize(); itrkass++){ // Loop Over All ReconstructedChargedParticleAssociations
						if( TrkPartAssocWeight->GetSize() > 0 ){
							if((*TrkPartAssocRec)[itrkass] == chargePartIndex){ // Select Entry Matching the ReconstructedChargedParticle Index
						    	if((*TrkPartAssocWeight)[itrkass] > elecIndexWeight){ // Find Particle with Greatest Weight = Contributed Most Hits to Track
									elecIndex = (*TrkPartAssocSim)[itrkass]; // Get Index of MCParticle Associated with ReconstructedChargedParticle
									elecIndexWeight = (*TrkPartAssocWeight)[itrkass];
		     	 				}
		  					}
		  				}
					}
					if((*TrkMCGenPDG)[elecIndex] == 11) hasElectron = true;
					cpdgid.push_back((*TrkRecoPDG)[idx]);
					cpdgidtruth.push_back((*TrkMCGenPDG)[elecIndex]);
                }

                RecoJet_constituent_pt.push_back(cpt);
                RecoJet_constituent_eta.push_back(ceta);
                RecoJet_constituent_phi.push_back(cphi);
                RecoJet_constituent_nhits.push_back(chits);
                RecoJet_constituent_pdgid.push_back(cpdgid);
                RecoJet_constituent_pdgidTruth.push_back(cpdgidtruth);
                RecoJet_hasElectron.push_back(hasElectron);
                RecoJet_maxPtPart_pt.push_back(maxPtReco);
            }
            
        	// --- Gen clustering ---
            ClusterSequenceArea cs_gen(particles_gen, jet_def, area_def);
            std::vector<PseudoJet> jets_gen = sorted_by_pt(cs_gen.inclusive_jets(minJetPt));
            for (auto &jet : jets_gen) {
                GenJet_pt.push_back(jet.pt());
                GenJet_eta.push_back(jet.eta());
                GenJet_phi.push_back(jet.phi_std());
                GenJet_E.push_back(jet.e());
                GenJet_M.push_back(jet.m());

                bool hasGenElectron = false;
                bool hasGenNeutral = false;
                float maxPtGen = -1.0;
                std::vector<float> gpt, geta, gphi;
                std::vector<int> gpdgid;
				gpt.clear(); geta.clear(); gphi.clear(); gpdgid.clear();

                for (auto &c : jet.constituents()) {
                    int idx = c.user_index();
                    TVector3 gv((*TrkGenPx)[idx], (*TrkGenPy)[idx], (*TrkGenPz)[idx]);
                    gpt.push_back(gv.Pt());
                    geta.push_back(gv.Eta());
                    gphi.push_back(gv.Phi());
                    gpdgid.push_back((*TrkGenPDG)[idx]);
                    if (gv.Pt() > maxPtGen) maxPtGen = gv.Pt();
                    if ((*TrkGenPDG)[idx] == 11) hasGenElectron = true;
                    if ((*TrkGenCharge)[idx] == 0) hasGenNeutral = true;
                }
                GenJet_constituent_pt.push_back(gpt);
                GenJet_constituent_eta.push_back(geta);
                GenJet_constituent_phi.push_back(gphi);
                GenJet_hasElectron.push_back(hasGenElectron);
                GenJet_hasNeutral.push_back(hasGenNeutral);
                GenJet_maxPtPart_pt.push_back(maxPtGen);
                GenJet_constituent_pdgid.push_back(gpdgid);
            }
			 
			tree->Fill();      
            
		}
		
	}
	
    // Write and close
    for (auto *tree : trees) tree->Write();
    OutFile->Close();

	cout << "Total number of events: " << globalEvent << endl;    
    

}

