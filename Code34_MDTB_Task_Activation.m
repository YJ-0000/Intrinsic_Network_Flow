clear; clc;

startup;

%% Settings
TR = 1;

slurm_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if ~isnan(slurm_id)
    rng('default');
    rng(77 * slurm_id + 123);
    fprintf('SLURM Mode: Seed set to %d\n', 77 * slurm_id + 123);
end

%% Load paths
path_info = load('secure_info/path_info.mat');
mdtb_data_path        = path_info.mdtb_preproc_data_path;
mdtb_data_unprocessed = path_info.mdtb_bids_data_path;

%% Find subjects
sub_dirs       = dir(fullfile(mdtb_data_path, 'sub-*'));
sub_dirs(~[sub_dirs.isdir]) = [];
sub_unproc_dirs = dir(fullfile(mdtb_data_unprocessed, 'sub-*'));

if ~isnan(slurm_id)
    sub_dirs        = sub_dirs(slurm_id);
    sub_unproc_dirs = sub_unproc_dirs(slurm_id);
end

num_subs = length(sub_dirs);

%% Preallocate
beta_act_vals = cell(num_subs, 32);
task_list     = cell(num_subs, 32);

%% Main loop
for nsub = 1:num_subs
    fprintf('\n>> Subject %d / %d\n', nsub, num_subs);

    cifti_files = dir(fullfile(sub_dirs(nsub).folder, sub_dirs(nsub).name, ...
        '**', 'func', '*desc-8mmSmoothed_bold.dtseries.nii'));
    event_files = dir(fullfile(sub_unproc_dirs(nsub).folder, sub_unproc_dirs(nsub).name, ...
        '**', 'func', '*events.tsv'));
    confounds_files = dir(fullfile(sub_dirs(nsub).folder, sub_dirs(nsub).name, ...
        '**', 'func', '*confounds_timeseries.tsv'));

    if isempty(cifti_files), continue; end

    for nrun = 1:length(event_files)
        tic
        fprintf('  Run %d / %d\n', nrun, length(event_files));

        % Load and percent-signal-change normalize
        data = cifti_read(fullfile(cifti_files(nrun).folder, cifti_files(nrun).name));
        data_smooth = data.cdata;
        data_smooth = 100 * data_smooth / mean(data_smooth, 'all');
        num_frames = size(data_smooth, 2);

        % Load events
        event_table = readtable(fullfile(event_files(nrun).folder, event_files(nrun).name), ...
            'FileType', 'text', 'Delimiter', '\t');
        trial_types = event_table.taskName;

        type_list_run = sort(unique(trial_types));
        type_list_run(strcmp(type_list_run, 'instruct')) = [];
        num_type  = length(type_list_run);
        total_time = TR * num_frames;

        % Build HRF-convolved design matrix
        hrf_weights = zeros(num_frames, num_type);
        for nt = 1:num_type
            mask_task = strcmp(trial_types, type_list_run{nt});
            onsets    = event_table.onset(mask_task);
            durations = event_table.duration(mask_task);
            conv_reg  = createTaskRegressor(TR, total_time, onsets, durations, 16);
            hrf_weights(:, nt) = conv_reg(1:num_frames);
        end

        % Load motion confounds
        confounds_table = readtable(fullfile(confounds_files(nrun).folder, confounds_files(nrun).name), ...
            'FileType', 'text', 'Delimiter', '\t');
        motion_cols = contains(confounds_table.Properties.VariableNames, {'trans','rot'}, 'IgnoreCase', true) & ...
                     ~contains(confounds_table.Properties.VariableNames, 'power2', 'IgnoreCase', true);
        motion_confounds = confounds_table{:, motion_cols};
        motion_confounds(isnan(motion_confounds)) = 0;

        % GLM with ReML
        X = [ones(num_frames, 1), hrf_weights, motion_confounds - mean(motion_confounds)];
        beta = REML_1stLV_spm(X, data_smooth', TR, 128);

        beta_act_vals{nsub, nrun} = beta(2:num_type+1, :)';
        task_list{nsub, nrun} = type_list_run;
        toc
    end
end

%% Save
if isnan(slurm_id)
    filename = 'results/MDTB_task_act';
else
    filename = sprintf('results/MDTB_task_act_sub%02d', slurm_id);
end

try
    save(filename, 'beta_act_vals', 'task_list');
catch
    save(filename, 'beta_act_vals', 'task_list', '-v7.3');
end