%% setup
clear; clc;

startup;

is_overwrite = true;

current_path = pwd;

path_info = load('secure_info/path_info.mat');

mdtb_data_path = path_info.mdtb_preproc_data_path;
HCP_data_path = path_info.HCP_denoised_path;

% Define the sampling frequency
TR = 1;
Fs = 1/TR;  % Hz

% Define the frequency range for the band-pass filter
lowCutoff = 0.008;  % Hz
highCutoff = 0.15;  % Hz

%% Get env variables
slurm_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));

%%
sub_MDTB_dirs = dir(fullfile(mdtb_data_path,'sub-*'));
sub_MDTB_dirs = sub_MDTB_dirs([sub_MDTB_dirs.isdir]);
num_subjects = length(sub_MDTB_dirs);

HCP_files = dir(fullfile(HCP_data_path,'**','*.dtseries.nii'));
cifti_template = cifti_read(...
    fullfile(HCP_files(1).folder,HCP_files(1).name));

if isnan(slurm_id)
    sub_idx_list = 1:num_subjects;
else
    sub_idx_list = slurm_id;
end

for nsub = sub_idx_list
    cifti_file = dir(...
        fullfile(sub_MDTB_dirs(nsub).folder,sub_MDTB_dirs(nsub).name,'**','func','*space-fsLR_den-91k_desc-8mmSmoothed_bold.dtseries.nii'));
    confounds_file = dir(...
            fullfile(sub_MDTB_dirs(nsub).folder,sub_MDTB_dirs(nsub).name,'**','func','*confounds_timeseries.tsv'));
        
    for nrun = 1:length(cifti_file)
        cifti_file_path = fullfile(cifti_file(nrun).folder,cifti_file(nrun).name);
    
        denoised_cifti_file_name = replace(cifti_file(nrun).name,...
            'desc-8mmSmoothed_bold.dtseries.nii','desc-8mmSmoothedDenoised_bold.dtseries.nii');
        denoised_cifti_file_path = fullfile(cifti_file(nrun).folder,...
            denoised_cifti_file_name);
    
        if ~is_overwrite
            if exist(denoised_cifti_file_path,"file")
                fprintf('%s is already denoised. Since overwrite = false, passing this subject ... \n',sub_MDTB_dirs(nsub).name);
                continue;
            end
        end
    
        data = niftiread(cifti_file_path);
        data = squeeze(data)';
        
        confounds_file_path = fullfile(confounds_file(nrun).folder,confounds_file(nrun).name);
        confounds_table = readtable(confounds_file_path,"FileType","text",...
            "Delimiter","\t");
    
        confounds_names = confounds_table.Properties.VariableNames;
    
        idx = contains(confounds_names, {'trans','rot'}, 'IgnoreCase', true);
        matched_confounds = confounds_names(idx);
    
        motion_confounds = confounds_table{:,matched_confounds};
    
        motion_confounds(isnan(motion_confounds)) = 0;
        
        idx_wm = contains(confounds_names, {'white_matter'}, 'IgnoreCase', true);
        idx_wm = find(idx_wm);
        idx_csf = contains(confounds_names, {'csf'}, 'IgnoreCase', true);
        idx_csf = find(idx_csf);
        matched_confounds = confounds_names(:,[idx_wm(1),idx_csf(1)]);
    
        physio_confounds = confounds_table{:,matched_confounds};
        physio_confounds(isnan(physio_confounds)) = 0;
    
        X = [ones(size(motion_confounds,1),1),...
            (1:size(motion_confounds,1))',((1:size(motion_confounds,1))').^2, ...
            motion_confounds,physio_confounds];
        
        data_denoised = data - (X * (pinv(X) * data'))';
    
        data_denoised_filtered = bandpass(data_denoised',[lowCutoff, highCutoff],Fs)';
    
        data_denoised_cifti = ...
            cifti_struct_create_from_template(cifti_template,data_denoised_filtered, 'dtseries');
    
        cifti_write(data_denoised_cifti,denoised_cifti_file_path);
    
    end
    fprintf('%03d/%03d: %s is denoised! \n',nsub,num_subjects,sub_MDTB_dirs(nsub).name);
    
end
