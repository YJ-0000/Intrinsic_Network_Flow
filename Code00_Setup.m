%% Setup secure_info
% This script should be executed in the project root directory.
% It creates a path_info.mat file in the secure_info folder.
%
% % % ----- MODIFICATION REQUIRED ----- % % %
% Please specify all paths below before running this script.
% % % % % % % % % % % % % % % % % % % % % %

%% ===== HCP Resting-State Data =====
% ICA-FIX denoised data, smoothed to 6 mm via Connectome Workbench.
% NOTE: File paths for individual subjects should be stored separately
%       in secure_info/hcp_rest_file_list.mat (see README for format).
HCP_denoised_path = '';  % Path to the folder containing ICA-FIX denoised resting-state fMRI data

%% ===== HCP Demographic & Behavioral Data =====
behav_data_path      = '';  % Path to the CSV file of behavioral data (unrestricted, HCP S1200)
gene_data_path       = '';  % Path to the CSV file of restricted data (includes genetic info, HCP S1200)
freesurfer_data_path = '';  % Path to the CSV file of FreeSurfer data (unrestricted, HCP S1200)

%% ===== MDTB Task fMRI (OpenNeuro ds002105) =====
% Preprocessed with fMRIPrep, surface-projected and smoothed with Connectome Workbench.
% Expected filename format: sub-*_desc-8mmSmoothed_bold.dtseries.nii
mdtb_preproc_data_path = '';  % Path to fMRIPrep derivatives (e.g., .../ds002105/derivatives/fmriprep)
mdtb_bids_data_path    = '';  % Path to raw BIDS dataset (e.g., .../ds002105)

%% ===== Audiovisual Speech Task fMRI (OpenNeuro ds003717) =====
% Preprocessed with fMRIPrep, surface-projected and smoothed with Connectome Workbench.
% Expected filename format: sub-*_desc-8mmSmoothed_bold.dtseries.nii
audiovisual_preproc_data_path = '';  % Path to fMRIPrep derivatives (e.g., .../ds003717/derivatives/fmriprep)
audiovisual_bids_data_path    = '';  % Path to raw BIDS dataset (e.g., .../ds003717)

%% ===== Save =====
secure_folder = 'secure_info';
if ~exist(secure_folder, 'dir')
    mkdir(secure_folder);
end

file_name = 'path_info.mat';
file_fullpath = fullfile(secure_folder, file_name);

% Check if file already exists
if exist(file_fullpath, 'file')
    choice = questdlg(['The file "' file_name '" already exists in "' secure_folder '". ' ...
                       'Do you want to overwrite it?'], ...
                      'Confirm Overwrite', 'Yes', 'No', 'No');
    if strcmp(choice, 'No')
        disp('File was not saved. Exiting...');
        return;
    end
end

save(file_fullpath, ...
    'HCP_denoised_path', ...
    'behav_data_path', 'gene_data_path', 'freesurfer_data_path', ...
    'mdtb_preproc_data_path', 'mdtb_bids_data_path', ...
    'audiovisual_preproc_data_path', 'audiovisual_bids_data_path');

disp(['Path information saved to: ' file_fullpath]);