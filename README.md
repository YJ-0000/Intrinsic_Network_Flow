# Intrinsic Network Flows (INF)

Code repository for the paper on the **Intrinsic Network Flow (INF)** framework 

---

## 📁 Repository Structure

### 🧠 Group-level INF Estimation

| Script | Description |
|--------|-------------|
| `Code01_INF_G_lev_MIGP_backDR_fbDMD_loop.m` | Predictive validation of INF framework |
| `Code02_INF_S_lev_fbDMD_loop.m` | Predictive validaition Subject level INF |
| `Code03_INF_G_lev_MIGP_backDR_fbDMD_ALL.m` | INF mode estimation on **all** subjects (REST1/REST2), used for downstream analyses |

### 📊 Prediction Benchmarks

| Script | Description |
|--------|-------------|
| `Code04_Compare_Prediction.m` | Compare prediction accuracy across methods (Group-level INF vs. subject-wise INF/with or without fingerpintings) |

### 🔬 Fingerprinting & Individual Differences

| Script | Description |
|--------|-------------|
| `Code05_Fingerprints_ALL.m` | Extract subject-level spatial modes and temporal fingerprints (amplitude, persistence, progression) |
| `Code06_Temporal_Fingerprints_cognitive_relavance.m` | Partial correlation of temporal fingerprints with behavioral factors (FDR-corrected); heatmap visualization |
| `Code07_Heritability_analysis_Temporal.m` | ACE heritability model (APACE) on amplitude fingerprints; twin/sibling cosine distance analysis |

### 🌊 Mode Characterization

| Script | Description |
|--------|-------------|
| `Code10_INF_mode_consistency.m` | Assess consistency / test-retest of INF modes |
| `Code11_rsfMRI_features.m` | Compute seed-based FC (PCC, SMG, aIns, FEF, TPJ) and Mean Laterality Index (MLI) |
| `Code12_Gradients_Phasemap.m` | Correlate INF mode phase maps with connectivity gradients (cortex, cerebellum, thalamus, hippocampus, striatum, amygdala, brainstem); spin/Moran permutation tests |
| `Code13_Gradients_rsfMRIfeatures.m` | Correlate INF modes with rsfMRI features (DMN/CEN/SN seed FC, MLI); spin-test corrected; temporal evolution over oscillation cycle |

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

## ⚙️ Key Dependencies

- **MATLAB** (tested on R2025a)
- [SPM25](https://github.com/spm/spm)
- [cifti-matlab](https://github.com/Washington-University/cifti-matlab) — CIFTI read/write
- [ICATB / GIFT](https://trendscenter.org/software/gift/) — `icatb_icaAlgorithm` (Infomax ICA)
- [APACE](https://github.com/NISOx-BDI/APACE) — ACE heritability modeling
- [BrainSpace](https://github.com/MICA-MNI/BrainSpace) — `GradientMaps`, `spin_permutations`, `compute_mem`, `moran_randomization`
- [DataViz](https://github.com/povilaskarvelis/DataViz) - violin plot
---


