clear; clc;

startup;

%% Settings
target_dim = 27;
REST_num = 'REST1';
TR = 1;
fine_or_not = false;

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

%% Load group-level ICA-DMD results
g_ica_result_files = dir(sprintf('results/INF_G_lev_%s_ALL_%03d_MIGP_results_*.mat', REST_num, target_dim));
load_results = load(fullfile(g_ica_result_files(end).folder, g_ica_result_files(end).name));

lambda      = load_results.lambda;
source_maps = load_results.source_maps';
inv_source  = load_results.inv_source;
Phi_all     = load_results.Phi_all;

%% Build ROI parcellation
labels_cortical    = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
labels_subcortical = cifti_read('atlas/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii');

label_data = round(labels_cortical.cdata);
try
    subcort_data = round(labels_subcortical.cdata);
    label_data = [label_data; subcort_data(length(label_data)+1:end)];
catch
end

label_idx_list = sort(unique(label_data));
label_idx_list(label_idx_list == 0) = [];
N_ROIs = numel(label_idx_list);

%% Find subject directories
sub_dirs       = dir(fullfile(mdtb_data_path, 'sub-*'));
sub_dirs(~[sub_dirs.isdir]) = [];
sub_unproc_dirs = dir(fullfile(mdtb_data_unprocessed, 'sub-*'));

if ~isnan(slurm_id)
    sub_dirs        = sub_dirs(slurm_id);
    sub_unproc_dirs = sub_unproc_dirs(slurm_id);
end

num_subs = length(sub_dirs);

%% Preallocate output cells
beta_ampl_vals    = cell(num_subs, 32);
beta_sin_vals     = cell(num_subs, 32);
beta_cos_vals     = cell(num_subs, 32);
beta_ic_vals      = cell(num_subs, 32);
FC_ICN_vals       = cell(num_subs, 32);
FC_ROIs_vals      = cell(num_subs, 32);
FC_rest_ROIs_vals = cell(num_subs, 2);
task_list         = cell(num_subs, 32);

%% Main loop
for nsub = 1:num_subs
    fprintf('\n>> Subject %d / %d\n', nsub, num_subs);

    cifti_files = dir(fullfile(sub_dirs(nsub).folder, sub_dirs(nsub).name, ...
        '**', 'func', '*desc-8mmSmoothedDenoised_bold.dtseries.nii'));
    event_files = dir(fullfile(sub_unproc_dirs(nsub).folder, sub_unproc_dirs(nsub).name, ...
        '**', 'func', '*events.tsv'));

    if isempty(cifti_files), continue; end

    %% Task runs
    for nrun = 1:length(event_files)
        tic
        fprintf('  Run %d / %d\n', nrun, length(event_files));

        % Load and normalize data
        data = cifti_read(fullfile(cifti_files(nrun).folder, cifti_files(nrun).name));
        data_raw  = data.cdata;
        data_norm = normalize(data_raw')';
        data_norm(isnan(data_norm)) = 0;

        % Load events
        event_table = readtable(fullfile(event_files(nrun).folder, event_files(nrun).name), ...
            'FileType', 'text', 'Delimiter', '\t');

        if fine_or_not
            trial_types = event_table.trial_type;
        else
            trial_types = event_table.taskName;
        end

        type_list_run = sort(unique(trial_types));
        type_list_run(strcmp(type_list_run, 'instruct')) = [];
        num_type = length(type_list_run);
        num_frames = size(data_norm, 2);
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
        X = [ones(num_frames, 1), hrf_weights];

        % ICA time courses
        IC_tc     = inv_source * data_norm;
        IC_sub_tc = pinv(source_maps) * data_norm;

        % Subject-level spatial modes
        Phi_sub = source_maps * Phi_all;

        % Compute temporal coefficients via sliding DMD
        D = computeDMcoefficients(data_norm, Phi_sub);
        f = 3;
        dd = zeros(target_dim, f);
        for k = 1:f
            dd(:, k) = D(2:end) .^ (k - ceil(f/2));
        end

        B_tc = nan(target_dim, num_frames);
        % Edge frames: direct pseudoinverse
        for t = 1:floor(f/2)
            B_tc(:, t)       = pinv(Phi_sub) * data_norm(:, t);
            B_tc(:, end-t+1) = pinv(Phi_sub) * data_norm(:, end-t+1);
        end
        % Interior frames: sliding window
        for t = 1+floor(f/2):num_frames-floor(f/2)
            t_start = t - floor(f/2);
            B_tc(:, t) = compute_B_from_Y_PHI_D(data_norm(:, t_start:t_start+f-1), Phi_sub, dd);
        end

        % --- INF betas (amplitude, cos, sin) via ReML ---
        fprintf('    ReML: INF...\n');
        b = REML_spm(X, abs(B_tc)');
        beta_ampl_vals{nsub, nrun} = b(2:end, :);

        b = REML_spm(X, cos(angle(B_tc))');
        beta_cos_vals{nsub, nrun} = b(2:end, :);

        b = REML_spm(X, sin(angle(B_tc))');
        beta_sin_vals{nsub, nrun} = b(2:end, :);

        % --- ICA betas via ReML ---
        fprintf('    ReML: ICs...\n');
        try
            b = REML_spm(X, IC_sub_tc');
            beta_ic_vals{nsub, nrun} = b(2:end, :);
        catch e
            warning(sprintf('Error in IC ReML: %s', e.message)); %#ok<SPWRN>
        end

        % --- FC among ICA components (task-weighted) ---
        fprintf('    FC (ICN)...\n');
        lower_idx = find(tril(true(target_dim), -1));
        FCs_ICN = zeros(length(lower_idx), num_type);
        for nt = 1:num_type
            w = hrf_weights(:, nt) .* (hrf_weights(:, nt) > 1);
            r = corr(IC_tc', 'Weights', w);
            FCs_ICN(:, nt) = r(lower_idx);
        end
        FC_ICN_vals{nsub, nrun} = FCs_ICN;

        % --- FC among ROIs (task-weighted) ---
        fprintf('    FC (ROIs)...\n');
        roi_ts = zeros(N_ROIs, num_frames);
        for nr = 1:N_ROIs
            roi_ts(nr, :) = mean(data_raw(label_data == label_idx_list(nr), :));
        end

        FCs_ROIs = zeros(N_ROIs, N_ROIs, num_type);
        for nt = 1:num_type
            w = hrf_weights(:, nt) .* (hrf_weights(:, nt) > 1);
            FCs_ROIs(:, :, nt) = corr(roi_ts', 'Weights', w);
        end
        FC_ROIs_vals{nsub, nrun} = FCs_ROIs;

        task_list{nsub, nrun} = type_list_run;
        toc
    end

    %% Rest runs (if present: runs 33-34)
    if length(cifti_files) > 32
        sm_maps_rest = cell(1, 2);
        for nrest = 1:2
            nrun = 32 + nrest;
            fprintf('  Rest run %d\n', nrest);

            data = cifti_read(fullfile(cifti_files(nrun).folder, cifti_files(nrun).name));
            data_raw  = data.cdata;
            data_norm = normalize(data_raw')';
            data_norm(isnan(data_norm)) = 0;

            % ROI FC
            roi_ts = zeros(N_ROIs, size(data_raw, 2));
            for nr = 1:N_ROIs
                roi_ts(nr, :) = mean(data_raw(label_data == label_idx_list(nr), :));
            end
            FC_rest_ROIs_vals{nsub, nrest} = corr(roi_ts');

            % Subject spatial maps from rest
            IC_tc = inv_source * data_norm;
            sm_maps_rest{nrest} = data_norm * pinv(IC_tc);
        end
        sm_maps = 0.5 * (sm_maps_rest{1} + sm_maps_rest{2});
        Phi_sub = sm_maps * Phi_all;
    else
        sm_maps = [];
        Phi_sub = [];
    end
end

%% Save
if isnan(slurm_id)
    filename = 'results/MDTB_estimated_betas';
else
    if fine_or_not
        filename = sprintf('results/MDTB_estimated_betas_FCs_fine_smoothed_sub%02d', slurm_id);
    else
        filename = sprintf('results/MDTB_estimated_betas_FCs_smoothed_sub%02d', slurm_id);
    end
end

save(filename, 'beta_*_vals', 'FC_*_vals', 'task_list', 'sm_maps', 'Phi_sub');
save(sprintf('results/MDTB_FCs_sub%02d', slurm_id), 'FC_ROIs_vals', '-v7.3');
save(sprintf('results/MDTB_FCs_rest_sub%02d', slurm_id), 'FC_rest_ROIs_vals', '-v7.3');