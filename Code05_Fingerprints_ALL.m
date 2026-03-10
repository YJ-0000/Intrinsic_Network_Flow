clear; clc;

startup;

%% Get env variables
slurm_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));

%% Run options
target_dim = 27;
if isnan(slurm_id) || slurm_id == 1
    REST_num = 'REST1';
elseif slurm_id == 2
    REST_num = 'REST2';
end

fprintf('Processing ses-%s data in HCP...\n', REST_num);

%% Load file list
file_list_data = load('secure_info/hcp_rest_file_list.mat');

if strcmp(REST_num, 'REST1')
    file_list_all = file_list_data.hcp_rest1_file_list;
    sub_ids       = file_list_data.sub_ids_rest1;
elseif strcmp(REST_num, 'REST2')
    file_list_all = file_list_data.hcp_rest2_file_list;
    sub_ids       = file_list_data.sub_ids_rest2;
else
    error('Unrecognized session name: %s', REST_num);
end

num_subjects = size(file_list_all, 1);
fprintf('Total subjects for %s: %d\n', REST_num, num_subjects);

%% Load group-level results
g_ica_result_files = dir(sprintf('results/INF_G_lev_%s_ALL_027_MIGP_results_*.mat', REST_num));
load_results = load(fullfile(g_ica_result_files(end).folder, g_ica_result_files(end).name));
lambda     = load_results.lambda;
inv_source = load_results.inv_source;
Phi_all    = load_results.Phi_all;

%% Get dimensions from first file
data = cifti_read(file_list_all{1, 1});
data = data.cdata;
voxel_num    = size(data, 1);
temporal_dim = size(data, 2);

%% Subject-level fingerprinting
INF_sub = nan(num_subjects, voxel_num, target_dim, 'single');

engagement_level_list = nan(num_subjects, target_dim, 'single');
persistent_rate_list  = nan(num_subjects, target_dim, 'single');
progression_rate_list = nan(num_subjects, target_dim, 'single');

for nsub = 1:num_subjects
    fprintf('>> %d / %d\n', nsub, num_subjects);

    % Load and process both runs
    Phi_orig_runs = cell(1, 2);
    data_norm_runs = cell(1, 2);

    for ndir = 1:2
        fprintf('  Loading %s\n', file_list_all{nsub, ndir});
        data = cifti_read(file_list_all{nsub, ndir});
        data_normalized = normalize(data.cdata')';

        time_course_sub = inv_source * data_normalized;
        temp_v = pinv(time_course_sub) * Phi_all;
        Phi_orig_runs{ndir} = data_normalized * temp_v;
        data_norm_runs{ndir} = data_normalized;
    end

    % Average spatial modes across runs
    Phi_orig_sub = 0.5 * (Phi_orig_runs{1} + Phi_orig_runs{2});

    % Concatenate time series for temporal fingerprinting
    X_seg = [data_norm_runs{1}(:, 2:end), data_norm_runs{2}(:, 2:end)];
    Y_seg = [data_norm_runs{1}(:, 1:end-1), data_norm_runs{2}(:, 1:end-1)];

    fprintf('  Fingerprinting...\n');
    [D, B_mean] = computeDMcoefficients([], Phi_orig_sub, [], X_seg, Y_seg);

    INF_sub(nsub, :, :) = Phi_orig_sub;

    engagement_level_list(nsub, :) = B_mean;
    persistent_rate_list(nsub, :)  = abs(D(2:end));
    progression_rate_list(nsub, :) = angle(D(2:end));
end

%% Save
dtStr = datestr(now, 'yyyymmdd_HHMMSS');

filename = sprintf('results/Temporal_Fingerprints_%s_ALL_%03d_results_%s.mat', REST_num, target_dim, dtStr);
save(filename, 'Phi_all', 'inv_source', '*_list', 'sub_ids', 'num_subjects');

filename = sprintf('results/Spatial_Fingerprints_%s_ALL_%03d_results_%s.mat', REST_num, target_dim, dtStr);
save(filename, 'INF_sub', '-v7.3');