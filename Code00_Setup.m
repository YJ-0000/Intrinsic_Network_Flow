%% Setup secure_data
% This script should be executed in the current directory.
% It creates a path_info.mat file in the secure_data folder with paths to HCP datasets.

% % % ----- MODIFICATION REQURIED ----- % % %
% Please specify the path to your downloaded HCP dataset 
% in order to run this script and the other scripts. 
% % % % % % % % % % % % % % % % % % % % % 

%%% Variables for paths to resting-state HCP dataset %%%
% We require both ICA-FIX denoised data and minimally preprocessed data 
% to properly exclude problematic subjects.
HCP_denoised_path = '/data/users1/ysong30/HCP_denoised_MSMAll_2017'; % Path to the folder containing ICA-FIX denoised resting-state fMRI data (downloaded from HCP)
HCP_preprocessed_rest_path = ''; % Path to the folder containing minimally preprocessed resting-state fMRI data (downloaded from HCP)

%%% Variables for paths to task-state HCP dataset (7 tasks) %%%
% The task data are assumed to be stored in separate folders.
% Each folder starts with the same prefix (e.g., 'HCP_3T_') and ends with the task name:
% 'WM', 'EMOTION', 'MOTOR', 'LANGUAGE', 'GAMBLING', 'SOCIAL', 'RELATIONAL'.
% For example, the working memory task data is stored in "HCP_3T_WM".
HCP_preprocessed_task_path = ''; 
% Prefix path for downloaded minimally preprocessed task-state fMRI data from HCP.
% Example: D:\HCP_3T_
% The scripts will automatically append the task name suffix (e.g., WM, EMOTION).

%%% Paths for demographic data files %%%
behav_data_path = 'data/Behavior_Data_8_18_2024_5_18_26.csv'; % Path to the CSV file of behavioral data in HCP S1200
gene_data_path = 'data/RESTRICTED_tbk5452_10_15_2024_18_49_21.csv'; % Path to the CSV file of restricted data (includes genetic info) in HCP S1200
freesurfer_data_path = 'data/HCP_freesurfer_ys1j13_10_16_2024_0_33_16.csv'; % Path to the CSV file of FreeSurfer data in HCP S1200


%% Save path information to secure_data
secure_data_folder = 'secure_data';
if ~exist(secure_data_folder,'dir')
    mkdir(secure_data_folder);
end
file_name = 'path_info.mat';
file_fullpath = fullfile(secure_data_folder,file_name);

% Check if file already exists and ask user before overwriting
if exist(file_fullpath,'file')
    choice = questdlg(['The file "' file_name '" already exists in "' secure_data_folder '". ' ...
                       'Do you want to overwrite it?'], ...
                      'Confirm Overwrite', ...
                      'Yes','No','No');
    if strcmp(choice,'No')
        disp('File was not saved. Exiting...');
        return;
    end
end

save(file_fullpath, ...
    'HCP_denoised_path','HCP_preprocessed_rest_path','HCP_preprocessed_task_path', ...
    'behav_data_path','gene_data_path','freesurfer_data_path');

disp(['Path information saved to: ' file_fullpath]);
