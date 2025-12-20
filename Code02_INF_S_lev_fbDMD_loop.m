clear; clc;

startup;

current_path = pwd;

%% Get env variables
slurm_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));

%% Run options

max_dim = 100;
if isnan(slurm_id) || slurm_id == 1 || slurm_id == 3 
    target_dim_list = 10:10:100; 
elseif slurm_id == 2 || slurm_id == 4
    target_dim_list = 21:1:29; 
end

if isnan(slurm_id) || slurm_id == 1 || slurm_id == 2
    sample_set = 'discovery';
elseif slurm_id == 3 || slurm_id == 4
    sample_set = 'replication';
end

fitting_time_window_list =  [1];
predict_time_window_list = [1,2,4,8];

%%

data_load = load('secure_info/path_info.mat');
HCP_data_path = data_load.HCP_denoised_path;
data_load = load('results/split_subjects_105.mat');

data_folders = dir(HCP_data_path);
isFolder = [data_folders.isdir];
names = {data_folders.name};
mask = isFolder & ~ismember(names, {'.','..'});
data_folders = data_folders(mask);

sub_ids = data_load.sub_ids;
if strcmp(sample_set,'discovery')
    sub_ids_set_explore = data_load.sub_ids_set1_explore;
elseif strcmp(sample_set, 'replication')
    sub_ids_set_explore = data_load.sub_ids_set2_explore;
else
    error('Undefined sample set!!');
end
[~,IA,~] = intersect(sub_ids,sub_ids_set_explore);
data_folders_train = data_folders(IA);

sub_ids = data_load.sub_ids;
if strcmp(sample_set,'discovery')
    sub_ids_set_test = data_load.sub_ids_set1_test;
elseif strcmp(sample_set, 'replication')
    sub_ids_set_test = data_load.sub_ids_set2_test;
else
    error('Undefined sample set!!');
end
[~,IA,~] = intersect(sub_ids,sub_ids_set_test);
data_folders_test = data_folders(IA);

%%

num_subjects = length(data_folders_test);

REST_num = 'REST1';

for target_dim = target_dim_list

    %% Test prediction
    TRtarget = 0.72;

    time_len = 1200;

    R2_DM_array_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));
    R2_lin_array_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));
    R2_null_array_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));
    R2_DM_array_cortex_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));
    R2_null_array_cortex_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));
    R2_DM_array_subcortical_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));
    R2_null_array_subcortical_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));

    for nfolder = 1:num_subjects
        fprintf('>> %d / %d \n',nfolder,num_subjects);
        tic

        subname = split(data_folders_test(nfolder).name,'_');
        % REST_num = subname{4};
        subname = subname{1};

        datafile_path_LR = fullfile(data_folders_test(nfolder).folder, ...
                            subname, ...
                            ['s6_rfMRI_',REST_num,'_LR_Atlas_hp2000_clean.dtseries.nii']);
        fprintf(['Loading data ... sub-',subname,' rfMRI LR \n']);
        data = cifti_read(datafile_path_LR);
        data = data.cdata;

        fprintf(['Normalizing data ... sub-',subname,' rfMRI LR \n']);
        data_normalized = normalize(data')';
        data_normalized_centered = data_normalized - mean(data_normalized);

        fprintf(['Performing single subj ICA ... sub-',subname,' rfMRI LR \n']);

        [~,score,lambda_k] = pca(data_normalized_centered(:,1:720),'NumComponents',target_dim,'Centered',true);
        score_whitened = score * diag(1 ./ sqrt(lambda_k(1:target_dim) + eps));
        [~, ~, ~, source_maps] = icatb_icaAlgorithm('infomax', score_whitened', {'verbose', 'off'});

        inv_source = pinv(source_maps');

        time_course_sub = (inv_source * data_normalized(:,1:720));

        [Phi_all,~] = performDMD(time_course_sub(:,2:end),time_course_sub(:,1:end-1),time_len,TRtarget);

    %     [R2_DM_array, R2_null_array, R2_lin_array,~, ~] ...
    %     = benchmark_using_DMs(Phi_all,data_normalized,source_maps',inv_source);

        temp_v = pinv(time_course_sub(:,1:end-1)) * Phi_all;
        Phi_orig_sub = data_normalized(:,2:720) * temp_v;

        [R2_DM_array,R2_DM_array_cortex,R2_DM_array_subcortical, R2_null_array,R2_null_array_cortex,R2_null_array_subcortical] ...
        =benchmark_using_DMs(Phi_orig_sub,data_normalized,1,1,0.6,fitting_time_window_list,predict_time_window_list,false);


        R2_DM_array_list(nfolder,1,:,:) = R2_DM_array;
        R2_null_array_list(nfolder,1,:,:) = R2_null_array;
%         R2_lin_array_list(nfolder,1,:,:) = R2_lin_array;
        R2_DM_array_cortex_list(nfolder,1,:,:) = R2_DM_array_cortex;
        R2_null_array_cortex_list(nfolder,1,:,:) = R2_null_array_cortex;
        R2_DM_array_subcortical_list(nfolder,1,:,:) = R2_DM_array_subcortical;
        R2_null_array_subcortical_list(nfolder,1,:,:) = R2_null_array_subcortical;

        fprintf(['Prediction complete ... sub-',subname,' rfMRI LR \n']);

        datafile_path_RL = fullfile(data_folders_test(nfolder).folder, ...
                            subname, ...
                            ['s6_rfMRI_',REST_num,'_RL_Atlas_hp2000_clean.dtseries.nii']);
        fprintf(['Loading data ... sub-',subname,' rfMRI RL \n']);
        data = cifti_read(datafile_path_RL);
        data = data.cdata;

        fprintf(['Normalizing data ... sub-',subname,' rfMRI LR \n']);
        data_normalized = normalize(data')';
        data_normalized_centered = data_normalized - mean(data_normalized);

        [~,score,lambda_k] = pca(data_normalized_centered(:,1:720),'NumComponents',target_dim,'Centered',true);
        score_whitened = score * diag(1 ./ sqrt(lambda_k(1:target_dim) + eps));
        [~, ~, ~, source_maps] = icatb_icaAlgorithm('infomax', score_whitened', {'verbose', 'off'});

        inv_source = pinv(source_maps');

        time_course_sub = (inv_source * data_normalized(:,1:720));

        [Phi_all,~] = performDMD(time_course_sub(:,2:end),time_course_sub(:,1:end-1),time_len,TRtarget);

    %     [R2_DM_array, R2_null_array, R2_lin_array,~, ~] ...
    %     = benchmark_using_DMs(Phi_all,data_normalized,source_maps',inv_source);

        temp_v = pinv(time_course_sub(:,1:end-1)) * Phi_all;
        Phi_orig_sub = data_normalized(:,2:720) * temp_v;

        [R2_DM_array,R2_DM_array_cortex,R2_DM_array_subcortical, R2_null_array,R2_null_array_cortex,R2_null_array_subcortical] ...
        =benchmark_using_DMs(Phi_orig_sub,data_normalized,1,1,0.6,fitting_time_window_list,predict_time_window_list,false);

        R2_DM_array_list(nfolder,2,:,:) = R2_DM_array;
        R2_null_array_list(nfolder,2,:,:) = R2_null_array;
%         R2_lin_array_list(nfolder,2,:,:) = R2_lin_array;
        R2_DM_array_cortex_list(nfolder,2,:,:) = R2_DM_array_cortex;
        R2_null_array_cortex_list(nfolder,2,:,:) = R2_null_array_cortex;
        R2_DM_array_subcortical_list(nfolder,2,:,:) = R2_DM_array_subcortical;
        R2_null_array_subcortical_list(nfolder,2,:,:) = R2_null_array_subcortical;

        fprintf(['Prediction complete ... sub-',subname,' rfMRI RL \n']);

        toc

    end


    %% display
    figure;
    bar( [...
    squeeze(mean(R2_DM_array_list(:,:,1,:),[1,2])), ...
    squeeze(mean(R2_lin_array_list(:,:,1,:),[1,2])),...
    squeeze(mean(R2_null_array_list(:,:,1,:),[1,2])) ...
    ]');
    legend({'predict 1s ahead','predict 2s ahead','predict 4s ahead','predict 8s ahead'});

    %% 
    % Get the current date and time formatted as YYYYMMDD_HHMMSS
    dtStr = datestr(now, 'yyyymmdd_HHMMSS');

    % Create a filename using the datetime string
    if strcmp(sample_set,'discovery')
        filename = ['results/disc_subject_wise_ica',num2str(target_dim,'%03d'),'_dmd_results_' dtStr '.mat'];
    elseif strcmp(sample_set, 'replication')
        filename = ['results/repl_subject_wise_ica', num2str(target_dim,'%03d'), '_dmd_results_' dtStr '.mat'];
    else
        error('Undefined sample set!!');
    end
    
    save(filename, 'R2*array_list','*time_window_list');
end

%%

%%
function [Phi_all, lambda] = performDMD(X, Y, time_len, TRtarget)
    % Precompute forward and backward operators
    A1 = X * Y';
    A2 = Y * Y';
    A_f = A1 * pinv(A2);

    B1 = Y * X';
    B2 = X * X';
    A_b = B1 * pinv(B2);

    % Initialize warning record
    lastwarn('');

    % Try to compute A = sqrt(A_f / A_b)
    try
        A = (A_f / A_b) ^ 0.5;
    catch ME
        warning('performDMD:MatrixPowerFailed', ...
            'Error computing matrix square root: %s. Falling back to A_f.', ME.message);
        A = A_f;
    end

    % If any warning occurred during matrix power, fall back
    [warnMsg, ~] = lastwarn;
    if ~isempty(warnMsg)
        warning('performDMD:MatrixPowerWarning', ...
            'Warning during matrix power: "%s". Falling back to A_f.', warnMsg);
        A = A_f;
    end

    % Retain only the real part
    A = real(A);

    % Clear warnings before eigen decomposition
    lastwarn('');

    % Perform eigen decomposition on A
    try
        [Phi_all, D] = eig(A);
    catch ME
        warning('performDMD:EigFailed', ...
            'Error during eig(A): %s. Falling back to eig(A_f).', ME.message);
        [Phi_all, D] = eig(A_f);
    end

    % If any warning occurred during eig(A), fall back
    [warnMsg, ~] = lastwarn;
    if ~isempty(warnMsg)
        warning('performDMD:EigWarning', ...
            'Warning during eig(A): "%s". Falling back to eig(A_f).', warnMsg);
        [Phi_all, D] = eig(A_f);
    end

    % Extract eigenvalues
    lambda = diag(D);

    % (Optional) Mode filtering based on magnitude or period
    % idx_exclude = abs(lambda) < 1e-4 | ...
    %               2*pi*TRtarget ./ abs(angle(lambda)) > time_len * TRtarget;
    % lambda(idx_exclude) = [];
    % Phi_all(:,idx_exclude) = [];
end

