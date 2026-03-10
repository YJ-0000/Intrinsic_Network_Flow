# Intrinsic Network Flows (INF)

Code repository for the paper on the **Intrinsic Network Flow (INF)** framework 

---
## 📋 Prerequisites

### Datasets

This project uses three fMRI datasets:

- **HCP S1200** — surface-registered (here, we used MSMSulc) ICA-FIX denoised resting-state data (REST1, REST2)
- **MDTB** — Multi-Domain Task Battery ([OpenNeuro ds002105](https://openneuro.org/datasets/ds002105))
- **Audiovisual Speech** — Audiovisual speech perception dataset ([OpenNeuro ds003717](https://openneuro.org/datasets/ds003717))

HCP dataset was spatially smoothed using **Connectome Workbench** with 6mm FWHM.

MDTB and Audiovisual datasets require preprocessing with surface registration and spatial smoothing prior to running these scripts. In this study, we used **fMRIPrep** for preprocessing and surface registration and **Connectome Workbench** for spatial smoothing (8mm FWHM). The scripts assume final outputs in BIDS format with filenames matching `sub-*_desc-8mmSmoothed_bold.dtseries.nii` (and `*8mmSmoothedDenoised*` for ICA-denoised versions).

### HCP Resting-State File List

HCP data directory structures vary across sites and users. To ensure compatibility, all HCP scripts load file paths from `secure_info/hcp_rest_file_list.mat`, which must contain:

- `hcp_rest1_file_list` — `(N_subjects × 2)` cell array of full file paths for REST1 (column 1 = LR, column 2 = RL)
- `hcp_rest2_file_list` — `(N_subjects × 2)` cell array of full file paths for REST2
- `sub_ids_rest1` — `(N_subjects × 1)` numeric array of subject IDs corresponding to REST1
- `sub_ids_rest2` — `(N_subjects × 1)` numeric array of subject IDs corresponding to REST2

⚠️ **Users must prepare this `hcp_rest_file_list.mat` file themselves** to match their local HCP data organization before running any HCP-related scripts.

### HCP Demographic & Behavioral Data

The following `.csv` files from the HCP are required:

- **Behavioral data** (unrestricted dataset)
- **FreeSurfer data** (unrestricted dataset)
- **Restricted data** (includes genetic/zygosity information)

These file paths are stored in `secure_info/path_info.mat` via `Code00_Setup.m`.

### Key Dependencies

- **MATLAB** (tested on R2025a)
- [SPM25](https://github.com/spm/spm)
- [cifti-matlab](https://github.com/Washington-University/cifti-matlab) — CIFTI read/write
- [GIFT](https://trendscenter.org/software/gift/) — `icatb_icaAlgorithm` (Infomax ICA)
- [APACE](https://github.com/NISOx-BDI/APACE) — ACE heritability modeling
- [BrainSpace](https://github.com/MICA-MNI/BrainSpace) — `GradientMaps`, `spin_permutations`, `compute_mem`, `moran_randomization`
- [DataViz](https://github.com/povilaskarvelis/DataViz) - violin plot
---       



## 📁 Repository Structure

### 📊 Prediction Benchmarks

| Script | Description |
|--------|-------------|
| `Code01_INF_G_lev_MIGP_backDR_fbDMD_loop.m` | Predictive validation of INF framework |
| `Code02_INF_S_lev_fbDMD_loop.m` | Predictive validaition Subject level INF |
| `Code03_Compare_Prediction.m` | Compare prediction accuracy across methods (Group-level INF vs. subject-wise INF/with or without fingerpintings) |

### 🧠 Group-level INF Estimation

| Script | Description |
|--------|-------------|
| `Code04_INF_G_lev_MIGP_backDR_fbDMD_ALL.m` | INF mode estimation on **all** subjects (REST1/REST2), used for downstream analyses |

### 🔬 Fingerprinting & Individual Differences

| Script | Description |
|--------|-------------|
| `Code05_Fingerprints_ALL.m` | Extract subject-level spatial modes and temporal fingerprints (amplitude, persistence, progression) |
| `Code06_Temporal_Fingerprints_cognitive_relavance.m` | Partial correlation of temporal fingerprints with behavioral factors (FDR-corrected); heatmap visualization |
| `Code07_Heritability_analysis_Temporal.m` | ACE heritability model (APACE) on amplitude fingerprints; twin/sibling cosine distance analysis |

### 🌊 Mode Characterization

| Script | Description |
|--------|-------------|
| `Code11_INF_mode_consistency.m` | Assess consistency / test-retest of INF modes |
| `Code12_rsfMRI_features.m` | Compute seed-based FC (PCC, SMG, aIns, FEF, TPJ) and Laterality Index |
| `Code13_Gradients_Phasemap.m` | Correlate INF mode phase maps with connectivity gradients (cortex, cerebellum, thalamus, hippocampus, striatum, amygdala, brainstem); spin/Moran permutation tests |
| `Code14_Gradients_rsfMRIfeatures.m` | Correlate INF modes with rsfMRI features (DMN/CEN/SN seed FC, MLI); spin-test corrected; temporal evolution over oscillation cycle |

### 🖼️ Visualization

| Script | Description |
|--------|-------------|
| `Code21_Display_whole.m` | Generating videos for whole-brain temporal evolution of each INF mode |

### 🎯 Task fMRI — MDTB

| Script | Description |
|--------|-------------|
| `Code31_MDTB_denoised.m` | Denoise fMRI data |
| `Code32_MDTB_INF_amp_phase_taskFC.m` | Extract INF and activation betas from MDTB dataset |
| `Code33_MDTB_Task_Classification.m` | Task decoding via LOSO-CV using INF features (phase, amplitude) vs. ICA activations |
| `Code34_MDTB_Task_Activation.m` | Estimate voxel-wise task activation betas (GLM + ReML with motion confounds) |
| `Code35_MDTB_Task_Act_FC_Recon.m` | Display and summarize task activation maps |
| `Code36_MDTB_INF_prediction.m` | BOLD signal prediction on MDTB task data using INF, null model |
| `Code37_MDTB_INF_prediction_results.m` | Aggregate and summarize MDTB prediction results across subjects |

### 🗣️ Task fMRI — Audiovisual Speech

| Script | Description |
|--------|-------------|
| `Code41_Audiovisual_FMRI_denoise.m` | Denoise fMRI data |
| `Code42_Audiovisual_INF_amp_phase.m` | Extract INF and activation betas from audiovisual speech task (AV/A/V/Null conditions). Condition-level analysis and LOSO classification on audiovisual task |

---


