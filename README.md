# EIC Jet Charge Analysis

This repository contains the code to analyze jet charge from Electron-Ion Collider (EIC) simulated data. The workflow consists of skimming raw trees to match reconstructed and generated tracks/jets, and then calculating the momentum-weighted jet charge.

## Repository Structure
* **`macros/`**: Contains the core ROOT C++ macros.
  * `process_eic_trees_list.C`: Skims input files, matches Reco/Gen objects using spatial separation (Delta R), and outputs a smaller `matched_tree`.
  * `PlotJetCharge.C`: Reads the `matched_tree` to calculate and plot the jet charge based on track properties within the jet cone.
* **`condor/`**: HTCondor submission scripts for running the skimming step in parallel on batch systems (e.g., lxplus).
* **`inputfiles/`**: Text lists of ROOT files to be processed.

## Prerequisites
* [ROOT](https://root.cern/) (compiled with C++17 or higher recommended)
* HTCondor (if running the batch submission)

## Workflow

### 1. Skimming and Matching
To run the skimming script locally over a list of files:
```bash
root -l -b -q 'macros/process_eic_trees_list.C("inputfiles/your_list.txt", "output_skim.root", 100)'
```

**Running on HTCondor:**
Before submitting, ensure your paths in `condor/read_eic_files_chunks.sh` point to your correct working directory.
```bash
condor_submit condor/read_eic_files_chunks.sub
```

### 2. Plotting Jet Charge
Once you have the skimmed `matched_tree` ROOT file, run the jet charge calculation. The macro takes several arguments so you can easily modify your kinematic cuts and weighting parameters without recompiling.

**Basic execution (uses default cuts):**
```bash
root -l -b -q 'macros/PlotJetCharge.C("output_skim.root")'
```

**Advanced execution (specifying all arguments):**
```bash
root -l -b -q 'macros/PlotJetCharge.C("output_skim.root", "reco", "reco", 0.5, 20.0, 0.5, 2, 0.4, 3.0)'
```

**Argument Breakdown (in order):**
1. **`inputFile`** *(string)*: Path to your skimmed root file containing the `matched_tree`.
2. **`jetType`** *(string)*: The type of jet to analyze. Accepts `"reco"`, `"gen"`, or `"matched"`.
3. **`trkType`** *(string)*: The type of track to use for the charge calculation. Accepts `"reco"`, `"gen"`, or `"matched"`.
4. **`KAPPA`** *(float)*: The momentum-weighting parameter (kappa) for the jet charge formula.
5. **`MIN_JET_PT`** *(float)*: Minimum transverse momentum for jets (GeV).
6. **`MIN_TRK_PT`** *(float)*: Minimum transverse momentum for tracks (GeV).
7. **`MIN_TRKS`** *(int)*: Minimum number of valid tracks required inside the jet cone.
8. **`R_CONE`** *(float)*: The jet cone radius used to associate tracks to the jet.
9. **`MAX_ETA`** *(float)*: Maximum absolute pseudorapidity (|eta|) for the jets.

### Outputs
The `PlotJetCharge` script automatically creates output files named dynamically based on your chosen parameters:
* A `.root` file containing the raw histogram.
* A formatted `.png` plot of the jet charge distribution with your cuts printed directly on the canvas.
Both files are saved in the `Plots/` directory.
