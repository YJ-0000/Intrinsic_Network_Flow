clear; clc;

startup;

%% Get env variables
slurm_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));

if isnan(slurm_id) || slurm_id == 1
    REST_num = 'REST1';
elseif slurm_id == 2
    REST_num = 'REST2';
end

%% Load file list
file_list_data = load('secure_info/hcp_rest_file_list.mat');

if strcmp(REST_num, 'REST1')
    file_list_all = file_list_data.hcp_rest1_file_list;
elseif strcmp(REST_num, 'REST2')
    file_list_all = file_list_data.hcp_rest2_file_list;
else
    error('Unrecognized session name: %s', REST_num);
end

num_subjects = size(file_list_all, 1);
fprintf('Processing %s: %d subjects\n', REST_num, num_subjects);

%% Get dimensions from first file
data = cifti_read(file_list_all{1, 1});
left_cortex_dim = data.diminfo{1,1}.models{1,2}.start - 1;
cortex_dim      = data.diminfo{1,1}.models{1,3}.start - 1;
spatial_dim     = size(data.cdata, 1);
temporal_dim    = size(data.cdata, 2);

%% Define ROI seeds from atlas
labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

seed_labels = struct( ...
    'PCC',    [35, 35+180], ...
    'aIns_R', 111, ...
    'SMG',    [148, 148+180], ...
    'FEF',    [10, 10+180], ...
    'TPJ',    [139, 139+180]);

seed_names = fieldnames(seed_labels);
num_seeds  = length(seed_names);

seed_idx = false(cortex_dim, num_seeds);
for i = 1:num_seeds
    lbl = seed_labels.(seed_names{i});
    seed_idx(:, i) = ismember(labels.cdata, lbl);
end

%% Seed-based connectivity
seed_conn = zeros(cortex_dim, num_seeds, 'single');

for nsub = 1:num_subjects
    fprintf('Seed conn >> %d / %d\n', nsub, num_subjects); tic;

    for ndir = 1:2
        fprintf('  Loading %s\n', file_list_all{nsub, ndir});
        data = cifti_read(file_list_all{nsub, ndir});
        data_cortex = data.cdata(1:cortex_dim, :);

        data_filtered = bandpass(data_cortex', [0.01, 0.1], 1/0.72)';

        % Global signal regression
        GS = mean(data_filtered);
        GS_mat = [ones(temporal_dim, 1), GS'];
        data_gsr = ((eye(temporal_dim) - GS_mat * pinv(GS_mat)) * data_filtered')';

        % Accumulate seed correlations
        for i = 1:num_seeds
            seed_ts = mean(data_gsr(seed_idx(:,i), :))';
            seed_conn(:, i) = seed_conn(:, i) + corr(seed_ts, data_gsr')';
        end
    end
    toc
end

seed_conn = seed_conn / (2 * num_subjects);

% Unpack into named variables for saving
seed_conn_PCC    = seed_conn(:, 1);
seed_conn_aIns_R = seed_conn(:, 2);
seed_conn_SMG    = seed_conn(:, 3);
seed_conn_FEF    = seed_conn(:, 4);
seed_conn_TPJ    = seed_conn(:, 5);

save('results/seed_based_conn.mat', 'seed_conn_*');

%% MLI (Mean Laterality Index)
frame_dt = 0.72;
time_window_length   = 30;   % seconds
time_window_interval = 10;   % seconds
window_frames = round(time_window_length / frame_dt);
step_frames   = round(time_window_interval / frame_dt);

MLI_all = zeros(num_subjects, 2, cortex_dim);

for nsub = 1:num_subjects
    fprintf('MLI >> %d / %d\n', nsub, num_subjects); tic;

    for ndir = 1:2
        fprintf('  Loading %s\n', file_list_all{nsub, ndir});
        data = cifti_read(file_list_all{nsub, ndir});
        data_cortex = data.cdata(1:cortex_dim, :);

        data_filtered = bandpass(data_cortex', [0.01, 0.1], 1/0.72)';

        GS_L = mean(data_filtered(1:left_cortex_dim, :));
        GS_R = mean(data_filtered(left_cortex_dim+1:cortex_dim, :));

        MLI = zeros(cortex_dim, 1);
        accum_num = 0;
        current_frame = 1;

        while (current_frame + window_frames) <= size(data_filtered, 2)
            accum_num = accum_num + 1;
            frame_range = current_frame:current_frame + window_frames;

            r_L = corr(GS_L(frame_range)', data_filtered(:, frame_range)');
            r_R = corr(GS_R(frame_range)', data_filtered(:, frame_range)');
            MLI = MLI + (r_L - r_R)';

            current_frame = current_frame + step_frames;
        end

        MLI_all(nsub, ndir, :) = MLI / accum_num;
    end
    toc
end

save('results/MLI.mat', 'MLI_all');