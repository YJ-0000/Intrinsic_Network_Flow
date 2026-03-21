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

MDTB and Audiovisual datasets require preprocessing with surface registration and spatial smoothing prior to running these scripts. In this study, we used **fMRIPrep** for preprocessing and surface registration and **Connectome Workbench** for spatial smoothing (8mm FWHM). The scripts assume final outputs in BIDS format with filenames matching `sub-*_desc-8mmSmoothed_bold.dtseries.nii`.

### HCP Resting-State Image File List

HCP data directory structures can vary across users. To ensure compatibility, all HCP scripts load file paths from `secure_info/hcp_rest_file_list.mat`, which must contain:

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

### Atlas Files

The following files are expected in the `atlas/` directory:

- **Glasser HCP-MMP1.0 parcellation** (360 areas): `Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii`
- **Cole-Anticevic parcellation** (cortical + subcortical): `CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii`
- **Inflated surfaces**: `Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii`, `Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii`
- **Spherical surfaces** (for spin tests): `S1200.L.sphere.32k_fs_LR.surf.gii`, `S1200.R.sphere.32k_fs_LR.surf.gii`

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

Scripts are numbered by analysis domain: `Code0*` for HCP resting-state, `Code1*` for mode characterization, `Code2*` for visualization, `Code3*` for MDTB task, and `Code4*` for audiovisual task.

### INF Estimation (Code01-03: Predictive validation, Code04: Group-level INF mode using all avaliable subjects, Code05-07: Fingerprinting)

| Script | Description |
|--------|-------------|
| `Code01_INF_G_lev_MIGP_backDR_fbDMD_loop.m` | Predictive validation of INF framework |
| `Code02_INF_S_lev_fbDMD_loop.m` | Predictive validaition Subject level INF |
| `Code03_Compare_Prediction.m` | Compare prediction accuracy across methods (Group-level INF vs. subject-wise INF / with or without fingerpintings) |
| `Code04_INF_G_lev_MIGP_backDR_fbDMD_ALL.m` | INF mode estimation on **all** subjects (REST1/REST2), used for downstream analyses |
| `Code05_Fingerprints_ALL.m` | Extract subject-level spatiotemporal fingerprints (amplitude, persistence, progression) |
| `Code06_Temporal_Fingerprints_cognitive_relavance.m` | Correlation of temporal fingerprints with behavioral factors |
| `Code07_Heritability_analysis_Temporal.m` | ACE heritability model (APACE) on amplitude fingerprints|

> **Note:** Although only 1-second-ahead prediction results are reported in the paper, the script attempts predictions up to 8 seconds ahead to examine long-term behavior. As a result, 7 fewer scans from the test-fitting segment (480 scans) are used compared to when only 1-second-ahead prediction is attempted.

> **Note:** `Code06` requires precomputed behavior latent scores from the paper "Exploring the Latent Structure of Behavior Using the Human Connectome Project’s Data" – see [GitHub Repository](https://github.com/connectomicslab/hcp-behavioral-domains) for the analysis code. *For convenience, we provide the precomputed scores as `scores_04_REST1_2_intersect.csv` in the data folder.*

### Mode Characterization (Code11-14)

| Script | Description |
|--------|-------------|
| `Code11_INF_mode_consistency.m` | Assess consistency of INF modes (Discovery vs. Replicaiton sets / REST1 vs. REST2) |
| `Code12_rsfMRI_features.m` | Seed-based FC (PCC, SMG, aIns, FEF, TPJ) and Laterality Index |
| `Code13_Compare_Gradients_Phasemap.m` | Correlate INF mode phase maps with connectivity gradients (cortex, cerebellum, thalamus, hippocampus, striatum, amygdala, brainstem) |
| `Code14_Compare_rsfMRIfeatures_Phasemap.m` | Correlate INF modes with rsfMRI features (DMN/CEN/SN, MLI) |

> **Note:** `Code12` requires precomputed connectivity gradients (via [BrainSpace](https://github.com/MICA-MNI/BrainSpace)) stored in `results/` as `.mat` files. The following files are expected:
> `gradient_rsfmri_cortex`, `gradient_rsfmri_cerebellum`, `gradient_rsfmri_thalamus`, `gradient_rsfmri_hippocampus`, `gradient_rsfmri_striatum`, `gradient_rsfmri_amygdala`, `gradient_rsfmri_brainstem`.
> Each file should contain a `gm` struct with `gm.gradients{1}` as the gradient matrix. **The gradient files used in this study are provided in the repository.**

**The rsfMRI features used in this study are also provided in the repository.** See `results/seed_based_conn.mat` and `MLI.mat`

### Visualization (Code21)

| Script | Description |
|--------|-------------|
| `Code21_Display_whole.m` | Generating videos for whole-brain temporal evolution of each INF mode |

### Task fMRI 1 — MDTB

| Script | Description |
|--------|-------------|
| `Code31_MDTB_denoised.m` | Denoise fMRI data |
| `Code32_MDTB_INF_amp_phase_taskFC.m` | Extract INF and activation betas from MDTB dataset |
| `Code33_MDTB_Task_Classification.m` | Task decoding via LOSO-CV using INF features (phase, amplitude) vs. IN activations |
| `Code34_MDTB_Task_Activation.m` | Estimate voxel-wise task activation betas (GLM) |
| `Code35_MDTB_Task_Act_FC_Recon.m` | Reonstruction of task-activation maps and task FC using INF modes |
| `Code36_MDTB_INF_prediction.m` | BOLD signal prediction on MDTB task data using INF, null model |
| `Code37_MDTB_INF_prediction_results.m` | Summarize MDTB prediction results across subjects |

### Task fMRI 2 — Audiovisual Speech

| Script | Description |
|--------|-------------|
| `Code41_Audiovisual_FMRI_denoise.m` | Denoise fMRI data |
| `Code42_Audiovisual_INF_amp_phase.m` | LOSO classification for four task conditions (AV/A/V/Null) on audiovisual task |

---


