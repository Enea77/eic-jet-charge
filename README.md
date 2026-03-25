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