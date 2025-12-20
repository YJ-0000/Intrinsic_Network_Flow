clear; clc;

startup;

current_path = pwd;

%% Get env variables
slurm_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));
% slurm_id = 2;
%% Run options

% max_dim = 100;
target_dim = 27;
if isnan(slurm_id) || slurm_id == 1 
    % target_dim_list = 10:10:100; 
    REST_num = 'REST1';
elseif slurm_id == 2 
    % target_dim_list = 21:1:29; 
    REST_num = 'REST2';
end

%%

data_load = load('secure_info/path_info.mat');
HCP_data_path = data_load.HCP_denoised_path;
data_load = load('results/split_subjects.mat');
sub_idx_include = data_load.sub_idx_include;
is_cognitive_impaired = data_load.is_cognitive_impaired;
is_excluded_due_movement = data_load.is_excluded_due_movement;
is_RL_processing_errors = data_load.is_RL_processing_errors;

data_folders = dir(HCP_data_path);
isFolder = [data_folders.isdir];
names = {data_folders.name};
mask = isFolder & ~ismember(names, {'.','..'});
data_folders = data_folders(mask);


fprintf('Checking all data.. \n\n')
if strcmp(REST_num,'REST1')
    data_folders_train = data_folders(sub_idx_include);
elseif strcmp(REST_num,'REST2')
    does_not_data_exist = true(length(data_folders),1);
    for nsub = 1:length(data_folders)
        fprintf('.. checking %03d/%03d \n',nsub,length(data_folders));
        ciftifile_LR = fullfile(data_folders(nsub).folder,...
                data_folders(nsub).name,...
                '*REST2_LR*.nii');
        ciftifile_LR = dir(ciftifile_LR);
        ciftifile_RL = fullfile(data_folders(nsub).folder,...
                data_folders(nsub).name,...
                '*REST2_RL*.nii');
        ciftifile_RL = dir(ciftifile_RL);
        if length(ciftifile_LR) == 1 && length(ciftifile_RL) == 1
            data_LR = cifti_read(fullfile(ciftifile_LR(1).folder,ciftifile_LR(1).name));
            if size(data_LR.cdata,2) == 1200
                data_RL = cifti_read(fullfile(ciftifile_RL(1).folder,ciftifile_RL(1).name));
                if size(data_RL.cdata,2) == 1200
                    does_not_data_exist(nsub) = false;
                end
            end
            
        end
    end

    sub_idx_exclude = does_not_data_exist | is_cognitive_impaired | is_excluded_due_movement(:,2) | is_RL_processing_errors;
    sub_idx_include = ~sub_idx_exclude;
    data_folders_train = data_folders(sub_idx_include);
else
    error('Unrecognized session name.')
end


%%
g_ica_result_files = dir(['results',filesep,'INF_G_lev_',REST_num,'_ALL_027_MIGP_results_*.mat']);

load_results = load(fullfile(g_ica_result_files(end).folder,g_ica_result_files(end).name));
lambda = load_results.lambda;
inv_source = load_results.inv_source;
Phi_all = load_results.Phi_all;

%% Loading one file

num_data_used = 1;

comp_est_accum = zeros(num_data_used*2,1);

nsub = 1;

subname = split(data_folders_train(nsub).name,'_');
% REST_num = 'REST1';
subname = subname{1};

datafile_path_LR = fullfile(data_folders_train(nsub).folder, subname,...
                    ['s6_rfMRI_',REST_num,'_LR_Atlas_hp2000_clean.dtseries.nii']);
fprintf(['Loading data ... sub-',subname,' rfMRI LR \n']);
data = cifti_read(datafile_path_LR);
data = data.cdata;
voxel_num = size(data,1);
temporal_dim = size(data,2);

%% Group PCA

num_subjects = length(data_folders_train);
sub_ids = zeros(num_subjects,1);

engagement_level_list = nan(num_subjects,target_dim);
persistent_rate_list = nan(num_subjects,target_dim);
progression_rate_list = nan(num_subjects,target_dim);

for nsub = 1:num_subjects
    fprintf('>> %d / %d \n',nsub,num_subjects);

    subname = split(data_folders_train(nsub).name,'_');
    subname = subname{1};
    sub_ids(nsub) = str2double(subname);

    datafile_path_LR = fullfile(data_folders_train(nsub).folder, ...
                        subname, ...
                        ['s6_rfMRI_',REST_num,'_LR_Atlas_hp2000_clean.dtseries.nii']);
    fprintf(['Loading data ... sub-',subname,' rfMRI LR \n']);
    data = cifti_read(datafile_path_LR);
    data = data.cdata;

    fprintf(['Normalizing data ... sub-',subname,' rfMRI LR \n']);
    data_normalized1 = normalize(data')';
    
    time_course_sub = (inv_source * data_normalized1);
    temp_v = pinv(time_course_sub) * Phi_all;
    Phi_orig_sub1 = data_normalized1 * temp_v;

    datafile_path_RL = fullfile(data_folders_train(nsub).folder, ...
                        subname, ...
                        ['s6_rfMRI_',REST_num,'_RL_Atlas_hp2000_clean.dtseries.nii']);
    fprintf(['Loading data ... sub-',subname,' rfMRI RL \n']);
    data = cifti_read(datafile_path_RL);
    data = data.cdata;

    fprintf(['Normalizing data ... sub-',subname,' rfMRI LR \n']);
    data_normalized2 = normalize(data')';
    
    time_course_sub = (inv_source * data_normalized2);
    temp_v = pinv(time_course_sub) * Phi_all;
    Phi_orig_sub2 = data_normalized2 * temp_v;

    Phi_orig_sub = 0.5 * (Phi_orig_sub1 + Phi_orig_sub2);

    fprintf(['Time course extracted ... sub-',subname,' rfMRI RL \n']);

    X_seg = [data_normalized1(:,2:end),data_normalized2(:,2:end)];
    Y_seg = [data_normalized1(:,1:end-1),data_normalized2(:,1:end-1)];

    [D, B_mean] = computeDMcoefficients([], Phi_orig_sub, [], X_seg, Y_seg);

    engagement_level_list(nsub,:) = B_mean;
    persistent_rate_list(nsub,:) = abs(D(2:end));
    progression_rate_list(nsub,:) = angle(D(2:end));
end



%% 
% Get the current date and time formatted as YYYYMMDD_HHMMSS
dtStr = datestr(now, 'yyyymmdd_HHMMSS');

% Create a filename using the datetime string
filename = ['results/Fingerprints_',REST_num,'_ALL_', num2str(target_dim,'%03d'), '_results_' dtStr '.mat'];

save(filename, 'Phi_all','inv_source','*_list',...
    'sub_ids','data_folders_train','num_subjects');

