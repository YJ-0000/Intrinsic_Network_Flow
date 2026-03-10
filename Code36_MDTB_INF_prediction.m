clear; clc;

startup;

%% Settings
target_dim = 27;
REST_num   = 'REST1';
TR         = 1;
t_step     = 1;

slurm_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if ~isnan(slurm_id)
    rng('default');
    rng(77 * slurm_id + 123);
    fprintf('SLURM Mode: Seed set to %d\n', 77 * slurm_id + 123);
end

%% Load paths and group-level results
path_info = load('secure_info/path_info.mat');
mdtb_data_path        = path_info.mdtb_preproc_data_path;
mdtb_data_unprocessed = path_info.mdtb_bids_data_path;

g_ica_files = dir(sprintf('results/INF_G_lev_%s_ALL_%03d_MIGP_results_*.mat', REST_num, target_dim));
load_results = load(fullfile(g_ica_files(end).folder, g_ica_files(end).name));
lambda      = load_results.lambda;
source_maps = load_results.source_maps';
inv_source  = load_results.inv_source;
Phi_all     = load_results.Phi_all;

%% Find subjects
sub_dirs       = dir(fullfile(mdtb_data_path, 'sub-*'));
sub_dirs(~[sub_dirs.isdir]) = [];
sub_unproc_dirs = dir(fullfile(mdtb_data_unprocessed, 'sub-*'));

if ~isnan(slurm_id)
    sub_dirs        = sub_dirs(slurm_id);
    sub_unproc_dirs = sub_unproc_dirs(slurm_id);
end

%% Prediction method names
method_names = {'Null', 'Rest_STF', 'Task_STF', 'Task_TF', 'Task_SF', 'Task_A'};
num_methods  = length(method_names);

% Accumulate results across subjects/runs
R2_all    = [];  % (num_methods x total_runs)
B_abs_all = [];
task_list = cell(length(sub_dirs), 32);

%% Main loop
for nsub = 1:length(sub_dirs)
    fprintf('\n>> Subject %d / %d\n', nsub, length(sub_dirs));

    cifti_files = dir(fullfile(sub_dirs(nsub).folder, sub_dirs(nsub).name, ...
        '**', 'func', '*desc-8mmSmoothedDenoised_bold.dtseries.nii'));
    event_files = dir(fullfile(sub_unproc_dirs(nsub).folder, sub_unproc_dirs(nsub).name, ...
        '**', 'func', '*events.tsv'));

    if length(cifti_files) <= 32, continue; end

    %% Rest runs → subject spatial maps and rest temporal fingerprint
    rest_files = cifti_files(33:end);
    sm_maps_rest = cell(1, 2);
    for nrest = 1:2
        data = cifti_read(fullfile(rest_files(nrest).folder, rest_files(nrest).name));
        data_norm = normalize(data.cdata')';
        data_norm(isnan(data_norm)) = 0;
        IC_tc = inv_source * data_norm;
        sm_maps_rest{nrest} = data_norm * pinv(IC_tc);
        if nrest == 1
            data_rest1 = data_norm;
        else
            data_rest2 = data_norm;
        end
    end

    sm_maps = 0.5 * (sm_maps_rest{1} + sm_maps_rest{2});
    Phi_sub = sm_maps * Phi_all;

    % Rest temporal fingerprint
    X_rest = [data_rest1(:, 1+t_step:end), data_rest2(:, 1+t_step:end)];
    Y_rest = [data_rest1(:, 1:end-t_step), data_rest2(:, 1:end-t_step)];
    D_rest = computeDMcoefficients([], Phi_sub, [], X_rest, Y_rest);

    % Rest engagement levels
    B_tc = pinv(Phi_sub) * Y_rest;
    B_abs_all = [B_abs_all, mean(abs(B_tc), 2)]; %#ok<AGROW>

    %% Load all task runs
    num_runs = length(event_files);
    data_task_sub   = cell(1, num_runs);
    hrf_weights_sub = cell(1, num_runs);

    for nrun = 1:num_runs
        fprintf('  Loading run %d / %d\n', nrun, num_runs);

        data = cifti_read(fullfile(cifti_files(nrun).folder, cifti_files(nrun).name));
        data_norm = normalize(data.cdata')';
        data_norm(isnan(data_norm)) = 0;
        data_task_sub{nrun} = data_norm;

        % Events
        event_table = readtable(fullfile(event_files(nrun).folder, event_files(nrun).name), ...
            'FileType', 'text', 'Delimiter', '\t');
        trial_types = event_table.taskName;
        type_list = sort(unique(trial_types));
        type_list(strcmp(type_list, 'instruct')) = [];

        num_frames = size(data_norm, 2);
        total_time = TR * num_frames;
        hrf = zeros(num_frames, length(type_list));
        for nt = 1:length(type_list)
            mask_task = strcmp(trial_types, type_list{nt});
            conv_reg = createTaskRegressor(TR, total_time, event_table.onset(mask_task), ...
                event_table.duration(mask_task), 16);
            hrf(:, nt) = conv_reg(1:num_frames);
        end
        hrf_weights_sub{nrun} = hrf;
        task_list{nsub, nrun} = type_list;
    end

    %% Prediction — leave-one-run-out
    for nrun = 1:num_runs
        fprintf('  Prediction run %d / %d\n', nrun, num_runs);
        tic

        data_test = data_task_sub{nrun};
        [voxel_dim, time_dim] = size(data_test);

        % Build training data from other runs
        X_train = nan(voxel_dim, (num_runs-1)*(time_dim-t_step));
        Y_train = nan(voxel_dim, (num_runs-1)*(time_dim-t_step));
        sm_maps_tasks = nan(voxel_dim, target_dim, num_runs-1);
        count = 0;
        for nrun_inner = 1:num_runs
            if nrun_inner == nrun, continue; end
            count = count + 1;
            seg = (count-1)*(time_dim-t_step)+1 : count*(time_dim-t_step);
            X_train(:, seg) = data_task_sub{nrun_inner}(:, 1+t_step:end);
            Y_train(:, seg) = data_task_sub{nrun_inner}(:, 1:end-t_step);

            d_inner = data_task_sub{nrun_inner};
            sm_maps_tasks(:,:,count) = d_inner * pinv(inv_source * d_inner);
        end

        sm_maps_task_gen = mean(sm_maps_tasks, 3);
        Phi_sub_task_gen = sm_maps_task_gen * Phi_all;

        % Task temporal fingerprints
        D_task_onlyTF = computeDMcoefficients([], Phi_sub, [], X_train, Y_train);
        D_task_STF    = computeDMcoefficients([], Phi_sub_task_gen, [], X_train, Y_train);

        % Task-general transition matrix
        A_task_gen = (inv_source * X_train) / (inv_source * Y_train);

        % Precompute test quantities
        Y_test   = data_test(:, 1:end-t_step);
        Y_target = data_test(:, 1+t_step:end);
        SS_tot   = mean((Y_target - mean(Y_target, 2)).^2, 2);
        valid_vox = ~isinf(1 ./ SS_tot);

        % --- Run all prediction methods ---
        R2_run = zeros(num_methods, 1);

        % 1) Null model (scalar autoregression)
        a_null = real(sum(dot(Y_test, Y_target, 1)) / sum(dot(Y_test, Y_test, 1)));
        R2_run(1) = compute_R2(Y_target, a_null * Y_test, SS_tot, valid_vox);

        % 2) Rest STF (rest spatial + rest temporal fingerprint)
        R2_run(2) = predict_INF(Y_target, Y_test, Phi_sub, D_rest, SS_tot, valid_vox);

        % 3) Task STF (task spatial + task temporal fingerprint)
        R2_run(3) = predict_INF(Y_target, Y_test, Phi_sub_task_gen, D_task_STF, SS_tot, valid_vox);

        % 4) Task TF only (rest spatial + task temporal fingerprint)
        R2_run(4) = predict_INF(Y_target, Y_test, Phi_sub, D_task_onlyTF, SS_tot, valid_vox);

        % 5) Task SF only (task spatial + rest temporal fingerprint)
        R2_run(5) = predict_INF(Y_target, Y_test, Phi_sub_task_gen, D_rest, SS_tot, valid_vox);

        % 6) Task-general A matrix
        IC_Y = inv_source * Y_test;
        X_pred = source_maps * (A_task_gen * IC_Y);
        X_resid = Y_target - X_pred;
        X_pred = real(D_task_onlyTF(1)) * X_resid + X_pred;
        R2_run(6) = compute_R2(Y_target, X_pred, SS_tot, valid_vox);

        R2_all = [R2_all, R2_run]; %#ok<AGROW>

        for m = 1:num_methods
            fprintf('    %s: R2 = %.4f\n', method_names{m}, R2_run(m));
        end
        toc
    end
end

%% Summary
fprintf('\n========== PREDICTION SUMMARY ==========\n');
for m = 1:num_methods
    vals = R2_all(m, :);
    fprintf('%-12s: mean R2 = %.4f +/- %.4f\n', method_names{m}, mean(vals), std(vals));
end

%% Save
if isnan(slurm_id)
    filename = 'results/MDTB_BOLD_prediction';
else
    filename = sprintf('results/MDTB_BOLD_prediction_sub%02d', slurm_id);
end
save(filename, 'R2_all', 'B_abs_all', 'task_list', 'method_names');

%% ======================== Helper Functions ========================

function R2 = predict_INF(Y_target, Y_test, Phi, D, SS_tot, valid_vox)
    % INF-based prediction: Phi * (D .* pinv(Phi)*Y) + residual correction
    B = pinv(Phi) * Y_test;
    X_pred = real(Phi * (B .* D(2:end)));
    X_resid = Y_target - X_pred;
    X_pred = real(D(1)) * X_resid + X_pred;
    R2 = compute_R2(Y_target, X_pred, SS_tot, valid_vox);
end

function R2 = compute_R2(Y_target, X_pred, SS_tot, valid_vox)
    SS_res = mean((Y_target - X_pred).^2, 2);
    R2 = 1 - mean(SS_res(valid_vox) ./ SS_tot(valid_vox), 'omitnan');
end