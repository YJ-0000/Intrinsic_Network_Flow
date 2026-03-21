clear; clc;

startup;

%% Settings
target_dim = 27;
REST_num   = 'REST1';
TR         = 2.47;
num_runs   = 5;
is_perm    = false;

flow_include = [1,3,5,7,9,11,13,15,17,19,21,23,27];

% Condition labels and display order
cond_labels = {'AV (noise) 0 SNR', 'AV (No noise)', 'AV (noise) -10 SNR', ...
               'AV (noise) -5 SNR', 'AV (noise) +5 SNR', 'A only (No noise)', ...
               'Null', 'V only'};
cond_display_order = [7, 6, 2, 5, 1, 4, 3, 8];
cond_classify      = [2,6,8,7];  % AV, A, V, Null

slurm_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if ~isnan(slurm_id)
    rng('default');
    rng(77 * slurm_id + 123);
    fprintf('SLURM Mode: Seed set to %d\n', 77 * slurm_id + 123);
end

%% Load paths and group-level results
path_info = load('secure_info/path_info.mat');
as_data_path        = path_info.audiovisual_preproc_data_path;
as_data_unprocessed = path_info.audiovisual_bids_data_path;

g_ica_result_files = dir(sprintf('results/INF_G_lev_%s_ALL_%03d_MIGP_results_*.mat', REST_num, target_dim));
load_results = load(fullfile(g_ica_result_files(end).folder, g_ica_result_files(end).name));
lambda      = load_results.lambda;
source_maps = load_results.source_maps';
inv_source  = load_results.inv_source;
Phi_all     = load_results.Phi_all;

%% Find subjects
sub_dirs       = dir(fullfile(as_data_path, 'sub-*'));
sub_dirs(~[sub_dirs.isdir]) = [];
sub_unproc_dirs = dir(fullfile(as_data_unprocessed, 'sub-*'));
num_subjects = length(sub_dirs);

%% Preallocate
beta_ampl_vals = nan(num_subjects, num_runs, 8, target_dim);
beta_sin_vals  = nan(num_subjects, num_runs, 8, target_dim);
beta_cos_vals  = nan(num_subjects, num_runs, 8, target_dim);
beta_ic_vals   = nan(num_subjects, num_runs, 8, target_dim);

%% Main estimation loop
for nsub = 1:num_subjects
    fprintf('\n>> Subject %d / %d\n', nsub, num_subjects);

    cifti_files_denoised = dir(fullfile(sub_dirs(nsub).folder, sub_dirs(nsub).name, ...
        'func', '*desc-8mmSmoothedDenoised_bold.dtseries.nii'));
    cifti_files_smoothed = dir(fullfile(sub_dirs(nsub).folder, sub_dirs(nsub).name, ...
        '**', 'func', '*desc-8mmSmoothed_bold.dtseries.nii'));
    event_files = dir(fullfile(sub_unproc_dirs(nsub).folder, sub_unproc_dirs(nsub).name, ...
        'func', '*events.tsv'));
    confounds_files = dir(fullfile(sub_dirs(nsub).folder, sub_dirs(nsub).name, ...
        '**', 'func', '*confounds_timeseries.tsv'));

    for nrun = 1:num_runs
        tic
        fprintf('  Run %d / %d\n', nrun, num_runs);

        % Load and normalize denoised data
        data = cifti_read(fullfile(cifti_files_denoised(nrun).folder, cifti_files_denoised(nrun).name));
        data_norm = normalize(data.cdata')';
        data_norm(isnan(data_norm)) = 0;
        num_frames = size(data_norm, 2);

        % Load events
        event_table = readtable(fullfile(event_files(nrun).folder, event_files(nrun).name), ...
            'FileType', 'text', 'Delimiter', '\t');
        trial_types = event_table.trial_type;

        % Optional permutation of condition labels
        if is_perm
            selected_type = {'AV_Quiet', 'NULL', 'V_only_Quiet_InScanner', 'A_Quiet'};
            perm_order = selected_type(randperm(length(selected_type)));
            for nt = 1:length(selected_type)
                trial_types(strcmp(trial_types, selected_type{nt})) = perm_order(nt);
            end
        end

        type_list = sort(unique(trial_types));
        num_type  = length(type_list);
        total_time = TR * num_frames;

        % Build HRF-convolved design matrix
        hrf_weights = zeros(num_frames, num_type);
        for nt = 1:num_type
            onsets    = event_table.onset(strcmp(trial_types, type_list{nt}));
            durations = 2.3 * ones(size(onsets));
            conv_reg  = createTaskRegressor(TR, total_time, onsets, durations, 16);
            hrf_weights(:, nt) = conv_reg(1:num_frames);
        end
        X = [ones(num_frames, 1), hrf_weights];

        % Subject spatial maps and modes
        IC_tc   = inv_source * data_norm;
        sm_maps = data_norm * pinv(IC_tc);
        Phi_sub = sm_maps * Phi_all;

        % Compute temporal coefficients via sliding DMD
        D = computeDMcoefficients(data_norm, Phi_sub);
        f = 3;
        dd = zeros(target_dim, f);
        for k = 1:f
            dd(:, k) = D(2:end) .^ (k - ceil(f/2));
        end

        B_tc = nan(target_dim, num_frames);
        for t = 1:floor(f/2)
            B_tc(:, t)       = pinv(Phi_sub) * data_norm(:, t);
            B_tc(:, end-t+1) = pinv(Phi_sub) * data_norm(:, end-t+1);
        end
        for t = 1+floor(f/2):num_frames-floor(f/2)
            t_start = t - floor(f/2);
            B_tc(:, t) = compute_B_from_Y_PHI_D(data_norm(:, t_start:t_start+f-1), Phi_sub, dd);
        end

        % INF betas via ReML
        b = REML_spm(X, abs(B_tc)');
        beta_ampl_vals(nsub, nrun, :, :) = b(2:end, :);

        b = REML_spm(X, cos(angle(B_tc))');
        beta_cos_vals(nsub, nrun, :, :) = b(2:end, :);

        b = REML_spm(X, sin(angle(B_tc))');
        beta_sin_vals(nsub, nrun, :, :) = b(2:end, :);

        % ICA betas via ReML
        % IC_sub_tc = pinv(sm_maps) * data_norm;
        b = REML_spm(X, IC_tc');
        beta_ic_vals(nsub, nrun, :, :) = b(2:end, :);

        % Activation betas (smoothed data + motion confounds)
        data_smooth = cifti_read(fullfile(cifti_files_smoothed(nrun).folder, cifti_files_smoothed(nrun).name));
        data_smooth = data_smooth.cdata;
        data_smooth = 100 * data_smooth / mean(data_smooth, 'all');

        confounds_table = readtable(fullfile(confounds_files(nrun).folder, confounds_files(nrun).name), ...
            'FileType', 'text', 'Delimiter', '\t');
        motion_cols = contains(confounds_table.Properties.VariableNames, {'trans','rot'}, 'IgnoreCase', true) & ...
                     ~contains(confounds_table.Properties.VariableNames, 'power2', 'IgnoreCase', true);
        motion_confounds = confounds_table{:, motion_cols};
        motion_confounds(isnan(motion_confounds)) = 0;

        toc
    end
end

%% LOSO Classification
beta_act_mean = squeeze(mean(beta_act_vals, 2));
num_cond = length(cond_classify);

feature_configs = {
    'Phase (cos+sin)',     flow_include, {'cos','sin'};
    'Amplitude',           flow_include, {'ampl'};
    'IC activations',      1:target_dim, {'ic'};
    'Phase + Amplitude',   flow_include, {'cos','sin','ampl'};
};

fprintf('\n========== CLASSIFICATION (LOSO-CV) ==========\n');

for i_cfg = 1:size(feature_configs, 1)
    cfg_name     = feature_configs{i_cfg, 1};
    cfg_flows    = feature_configs{i_cfg, 2};
    cfg_features = feature_configs{i_cfg, 3};

    % Build feature matrix: (num_subjects * num_cond) x num_features
    X_all = [];
    for flow_num = cfg_flows
        for i_f = 1:length(cfg_features)
            feat_type = cfg_features{i_f};
            switch feat_type
                case 'ampl', src = beta_ampl_vals;
                case 'cos',  src = beta_cos_vals;
                case 'sin',  src = beta_sin_vals;
                case 'ic',   src = beta_ic_vals;
            end
            mean_beta = squeeze(mean(src(:,:,:,flow_num), 2));
            X_col = reshape(mean_beta(:, cond_classify)', [], 1);
            X_all = [X_all, X_col]; %#ok<AGROW>
        end
    end

    Y_all  = repmat((1:num_cond)', [num_subjects, 1]);
    groups = reshape(repmat((1:num_subjects)', [1, num_cond])', [], 1);

    accuracy_CV = zeros(num_subjects, 1);
    Y_pred_all  = zeros(size(Y_all));

    for nsub = 1:num_subjects
        test_mask  = (groups == nsub);
        train_mask = ~test_mask;

        X_train = X_all(train_mask, :);
        X_test  = X_all(test_mask, :);
        Y_train = Y_all(train_mask);
        Y_test  = Y_all(test_mask);

        Mdl    = fitcecoc(X_train, Y_train, ...
                    'Learners', templateSVM('KernelFunction', 'linear'), ...
                    'Coding', 'onevsall');
        Y_pred = predict(Mdl, X_test);

        accuracy_CV(nsub) = mean(Y_test == Y_pred);
        Y_pred_all(test_mask) = Y_pred;
    end

    fprintf('%-25s  Acc: %.2f%% +/- %.2f%%\n', cfg_name, ...
        mean(accuracy_CV)*100, std(accuracy_CV)*100);

    figure;
    confusionchart(Y_all, Y_pred_all, ...
        'Title', sprintf('LOSO: %s', cfg_name), ...
        'RowSummary', 'row-normalized', ...
        'ColumnSummary', 'column-normalized');
end

%% Display activation maps
mkdir('Speech_results');
save_dir = fullfile(pwd, 'Speech_results');

labels = cifti_read('atlas/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii');
lh_surface_file = 'atlas/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii';
rh_surface_file = 'atlas/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii';

fig = figure('Position', [0, 0, 554, 416]);
for n_task = cond_classify
    % Mean beta map
    act_vec = mean(squeeze(beta_act_mean(:, n_task, :)))';
    snapshot = cifti_struct_create_from_template(labels, act_vec, 'dscalar');
    max_scale = max(abs(act_vec));
    display_cifti_cortex(fig, snapshot, lh_surface_file, rh_surface_file, [], -max_scale, max_scale, []);
    print(fig, fullfile(save_dir, sprintf('Act_MeanBeta_%s.jpg', cond_labels{n_task})), '-djpeg', '-r300');

    % T-statistic map
    [~, ~, ~, st] = ttest(squeeze(beta_act_mean(:, n_task, :)));
    act_vec = st.tstat';
    snapshot = cifti_struct_create_from_template(labels, act_vec, 'dscalar');
    max_scale = max(abs(act_vec));
    display_cifti_cortex(fig, snapshot, lh_surface_file, rh_surface_file, [], -max_scale, max_scale, []);
    print(fig, fullfile(save_dir, sprintf('Act_Tval_%s.jpg', cond_labels{n_task})), '-djpeg', '-r300');
end
close(fig);

%% ======================== Helper Functions ========================

function plot_bar_with_sem(data, labels, idx)
    num_cond = length(idx);
    means = zeros(1, num_cond);
    sems  = zeros(1, num_cond);
    for n = 1:num_cond
        means(n) = mean(data(:, idx(n)));
        sems(n)  = 1.96 * std(data(:, idx(n))) / sqrt(size(data, 1));
    end

    figure;
    bar(means); hold on;
    errorbar(1:num_cond, means, sems, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    set(gca, 'XTickLabel', labels(idx), 'FontName', 'Times New Roman', 'FontSize', 12);
    ylabel('Beta');
    hold off;
end