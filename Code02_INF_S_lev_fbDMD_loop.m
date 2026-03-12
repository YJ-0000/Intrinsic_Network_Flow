clear; clc;

startup;

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

fitting_time_window_list  = [1];
predict_time_window_list  = [1,2,4,8];

%% Load file list and subject split
data_load = load('results/split_subjects_105.mat');
file_list_data = load('secure_info/hcp_rest_file_list.mat', 'hcp_rest1_file_list', 'sub_ids_rest1');
hcp_rest1_file_list = file_list_data.hcp_rest1_file_list;

sub_ids = file_list_data.sub_ids_rest1;

if strcmp(sample_set, 'discovery')
    sub_ids_test = data_load.sub_ids_set1_test;
elseif strcmp(sample_set, 'replication')
    sub_ids_test = data_load.sub_ids_set2_test;
else
    error('Undefined sample set!!');
end

[~, idx_test] = intersect(sub_ids, sub_ids_test);
file_list_test = hcp_rest1_file_list(idx_test, :);  % (num_test x 2)
num_subjects = size(file_list_test, 1);

%%

for target_dim = target_dim_list

    %% Test prediction
    TRtarget = 0.72;
    time_len = 1200;

    R2_DM_array_list              = zeros(num_subjects, 2, length(fitting_time_window_list), length(predict_time_window_list));
    R2_lin_array_list             = zeros(num_subjects, 2, length(fitting_time_window_list), length(predict_time_window_list));
    R2_null_array_list            = zeros(num_subjects, 2, length(fitting_time_window_list), length(predict_time_window_list));
    R2_DM_array_cortex_list       = zeros(num_subjects, 2, length(fitting_time_window_list), length(predict_time_window_list));
    R2_null_array_cortex_list     = zeros(num_subjects, 2, length(fitting_time_window_list), length(predict_time_window_list));
    R2_DM_array_subcortical_list  = zeros(num_subjects, 2, length(fitting_time_window_list), length(predict_time_window_list));
    R2_null_array_subcortical_list= zeros(num_subjects, 2, length(fitting_time_window_list), length(predict_time_window_list));

    for nfolder = 1:num_subjects
        fprintf('>> %d / %d\n', nfolder, num_subjects);
        tic

        for ndir = 1:2
            fprintf('  Loading %s\n', file_list_test{nfolder, ndir});
            data = cifti_read(file_list_test{nfolder, ndir});
            data_normalized = normalize(data.cdata')';
            data_normalized_centered = data_normalized - mean(data_normalized);

            fprintf('  Performing single subj ICA...\n');
            [~, score, lambda_k] = pca(data_normalized_centered(:,1:720), 'NumComponents', target_dim, 'Centered', true);
            score_whitened = score * diag(1 ./ sqrt(lambda_k(1:target_dim) + eps));
            [~, ~, ~, source_maps] = icatb_icaAlgorithm('infomax', score_whitened', {'verbose', 'off'});

            inv_source = pinv(source_maps');
            time_course_sub = inv_source * data_normalized(:, 1:720);

            [Phi_all, ~] = performDMD(time_course_sub(:, 2:end), time_course_sub(:, 1:end-1), time_len, TRtarget);
            
            %%%% exact DMD
            % temp_v = pinv(time_course_sub(:, 1:end-1)) * Phi_all;
            % Phi_orig_sub = data_normalized(:, 2:720) * temp_v;
            %%%% DR-like
            temp_v = pinv(time_course_sub) * Phi_all;
            Phi_orig_sub = data_normalized(:, 1:720) * temp_v;

            [R2_DM_array, R2_DM_array_cortex, R2_DM_array_subcortical, ...
             R2_null_array, R2_null_array_cortex, R2_null_array_subcortical] ...
                = benchmark_using_DMs(Phi_orig_sub, data_normalized, 1, 1, 0.6, ...
                    fitting_time_window_list, predict_time_window_list, false);

            R2_DM_array_list(nfolder, ndir, :, :)              = R2_DM_array;
            R2_null_array_list(nfolder, ndir, :, :)            = R2_null_array;
            R2_DM_array_cortex_list(nfolder, ndir, :, :)       = R2_DM_array_cortex;
            R2_null_array_cortex_list(nfolder, ndir, :, :)     = R2_null_array_cortex;
            R2_DM_array_subcortical_list(nfolder, ndir, :, :)  = R2_DM_array_subcortical;
            R2_null_array_subcortical_list(nfolder, ndir, :, :)= R2_null_array_subcortical;

            fprintf('  Prediction complete (run %d)\n', ndir);
        end

        toc
    end
    
    %% Save
    dtStr = datestr(now, 'yyyymmdd_HHMMSS');

    if strcmp(sample_set, 'discovery')
        filename = sprintf('results/disc_subject_wise_ica%03d_dmd_results_%s.mat', target_dim, dtStr);
    elseif strcmp(sample_set, 'replication')
        filename = sprintf('results/repl_subject_wise_ica%03d_dmd_results_%s.mat', target_dim, dtStr);
    else
        error('Undefined sample set!!');
    end

    save(filename, 'R2*array_list', '*time_window_list');
end

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
end