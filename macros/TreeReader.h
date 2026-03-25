#include "CallLibraries.h"  // call libraries from ROOT and C++

// Event Kinematics
std::unique_ptr<TTreeReaderArray<float>> EvtQ2;
std::unique_ptr<TTreeReaderArray<float>> Evtx;
std::unique_ptr<TTreeReaderArray<float>> EvtQ2Gen;
std::unique_ptr<TTreeReaderArray<float>> EvtxGen;

// vertex
std::unique_ptr<TTreeReaderArray<float>> CTVx;
std::unique_ptr<TTreeReaderArray<float>> CTVy;
std::unique_ptr<TTreeReaderArray<float>> CTVz;
std::unique_ptr<TTreeReaderArray<int>> CTVndf;
std::unique_ptr<TTreeReaderArray<float>> CTVchi2;
std::unique_ptr<TTreeReaderArray<float>> CTVerr_xx;
std::unique_ptr<TTreeReaderArray<float>> CTVerr_yy;
std::unique_ptr<TTreeReaderArray<float>> CTVerr_zz;

// Scattered electrons
// reconstructed
std::unique_ptr<TTreeReaderArray<int>> ScatElecRecoId;
// truth
std::unique_ptr<TTreeReaderArray<int>> ScatElecGenId;


// Reconstructed charged particles (Tracks)
std::unique_ptr<TTreeReaderArray<float>> TrkRecoE;
std::unique_ptr<TTreeReaderArray<float>> TrkRecoPx;
std::unique_ptr<TTreeReaderArray<float>> TrkRecoPy;
std::unique_ptr<TTreeReaderArray<float>> TrkRecoPz;
std::unique_ptr<TTreeReaderArray<float>> TrkRecoM;
std::unique_ptr<TTreeReaderArray<int>> TrkRecoPDG;
std::unique_ptr<TTreeReaderArray<unsigned int>> TrkRecoNhits;

// Reco to Gen track association
std::unique_ptr<TTreeReaderArray<int>> TrkPartAssocRec;
std::unique_ptr<TTreeReaderArray<int>> TrkPartAssocSim;
std::unique_ptr<TTreeReaderArray<float>> TrkPartAssocWeight;

// Generated particles -> only stable particles from generator level -> Not available ?
std::unique_ptr<TTreeReaderArray<float>> TrkGenPx;
std::unique_ptr<TTreeReaderArray<float>> TrkGenPy;
std::unique_ptr<TTreeReaderArray<float>> TrkGenPz;
std::unique_ptr<TTreeReaderArray<float>> TrkGenE;
std::unique_ptr<TTreeReaderArray<float>> TrkGenM;
std::unique_ptr<TTreeReaderArray<float>> TrkGenCharge;
std::unique_ptr<TTreeReaderArray<int>> TrkGenPDG;

// MCParticles particles -> full generator information + sec decay from GEANT
std::unique_ptr<TTreeReaderArray<double>> TrkMCGenPx;
std::unique_ptr<TTreeReaderArray<double>> TrkMCGenPy;
std::unique_ptr<TTreeReaderArray<double>> TrkMCGenPz;
std::unique_ptr<TTreeReaderArray<double>> TrkMCGenM;
std::unique_ptr<TTreeReaderArray<float>> TrkMCGenCharge;
std::unique_ptr<TTreeReaderArray<int>> TrkMCGenPDG;
std::unique_ptr<TTreeReaderArray<int>> TrkMCGenStatus;


// Reads jet tree variables
/*
Arguments:
chain: TChain of input files
tree_reader: the object for TTreeReader
*/
void TreeReader(TChain* chain, std::unique_ptr<TTreeReader>& tree_reader) {

    tree_reader = std::make_unique<TTreeReader>(chain);
	// Event quantities
    EvtQ2 = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "InclusiveKinematicsElectron.Q2");
    Evtx = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "InclusiveKinematicsElectron.x");
    EvtQ2Gen = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "InclusiveKinematicsTruth.Q2");
    EvtxGen = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "InclusiveKinematicsTruth.x");
	// Vertex
	CTVx = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "CentralTrackVertices.position.x");
	CTVy = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "CentralTrackVertices.position.y");
	CTVz = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "CentralTrackVertices.position.z");
	CTVndf = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "CentralTrackVertices.ndf");
	CTVchi2 = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "CentralTrackVertices.chi2");
	CTVerr_xx = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "CentralTrackVertices.positionError.xx");
	CTVerr_yy = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "CentralTrackVertices.positionError.yy");
	CTVerr_zz = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "CentralTrackVertices.positionError.zz");
	// Scattered electron
    ScatElecRecoId = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "_InclusiveKinematicsElectron_scat.index");
    ScatElecGenId = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "MCScatteredElectrons_objIdx.index");
//    ScatElecGenId = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "_InclusiveKinematicsTruth_scat.index");

	// Reconstructed charged particles (Tracks)
    TrkRecoE = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedChargedParticles.energy");
    TrkRecoPx = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedChargedParticles.momentum.x");
    TrkRecoPy = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedChargedParticles.momentum.y");
    TrkRecoPz = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedChargedParticles.momentum.z");
    TrkRecoM = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedChargedParticles.mass");
    TrkRecoPDG = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "ReconstructedChargedParticles.PDG");
	TrkPartAssocRec = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "_ReconstructedChargedParticleAssociations_rec.index");
	TrkPartAssocSim = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "_ReconstructedChargedParticleAssociations_sim.index");
	TrkPartAssocWeight = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedChargedParticleAssociations.weight");
	TrkRecoNhits = std::make_unique<TTreeReaderArray<unsigned int>>(*tree_reader, "CentralCKFTrajectories.nMeasurements");
	
	// MC Gen -> add GEANT Stuff
	TrkMCGenStatus = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "MCParticles.generatorStatus");
	TrkMCGenPx = std::make_unique<TTreeReaderArray<double>>(*tree_reader, "MCParticles.momentum.x");
	TrkMCGenPy = std::make_unique<TTreeReaderArray<double>>(*tree_reader, "MCParticles.momentum.y");
	TrkMCGenPz = std::make_unique<TTreeReaderArray<double>>(*tree_reader, "MCParticles.momentum.z");
	TrkMCGenM = std::make_unique<TTreeReaderArray<double>>(*tree_reader, "MCParticles.mass");
	TrkMCGenPDG = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "MCParticles.PDG");
    TrkMCGenCharge = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "MCParticles.charge");

	// Full gen level
	TrkGenE = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "GeneratedParticles.energy");
	TrkGenPx = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "GeneratedParticles.momentum.x");
	TrkGenPy = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "GeneratedParticles.momentum.y");
	TrkGenPz = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "GeneratedParticles.momentum.z");
	TrkGenM = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "GeneratedParticles.mass");
	TrkGenPDG = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "GeneratedParticles.PDG");
	TrkGenCharge = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "GeneratedParticles.charge");

}
