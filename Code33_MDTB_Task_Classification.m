clear; clc;

%% Settings
target_dim   = 27;
num_subjects = 24;
num_sessions = 32;

% Propagating mode indices (one per conjugate pair)
flow_include = [1,3,5,7,9,11,13,15,17,19,21,23,27];

%% Load per-subject results
result_files = dir('results/MDTB_estimated_betas_FCs_smoothed_sub*');

beta_ampl_vals = cell(num_subjects, num_sessions);
beta_sin_vals  = cell(num_subjects, num_sessions);
beta_cos_vals  = cell(num_subjects, num_sessions);
beta_ic_vals   = cell(num_subjects, num_sessions);
task_list      = cell(num_subjects, num_sessions);

for nsub = 1:num_subjects
    res = load(fullfile(result_files(nsub).folder, result_files(nsub).name));
    beta_ampl_vals(nsub, :) = res.beta_ampl_vals;
    beta_cos_vals(nsub, :)  = res.beta_cos_vals;
    beta_sin_vals(nsub, :)  = res.beta_sin_vals;
    beta_ic_vals(nsub, :)   = res.beta_ic_vals;
    task_list(nsub, :)      = res.task_list;
end

%% Build unified task type list
all_tasks = vertcat(task_list{1, :});
task_types = sort(unique(all_tasks));

% Remove duplicates, instruction, and merge movie subtypes
task_types(contains(task_types, '2'))         = [];
task_types(contains(task_types, 'instruct'))  = [];
task_types(contains(task_types, 'nBackPic'))  = [];
task_types(contains(task_types, 'Movie'))     = [];
task_types{end+1} = 'Movie';
task_types = sort(task_types);

num_tasks = length(task_types);

%% Aggregate betas per subject x task (average across sessions)
beta_sources = {beta_ampl_vals, beta_cos_vals, beta_sin_vals, beta_ic_vals};
beta_names   = {'ampl', 'cos', 'sin', 'ic'};

num_rows = num_subjects * num_tasks;
sub_idx    = zeros(num_rows, 1);
task_labels = cell(num_rows, 1);
beta_agg    = struct();

for i_src = 1:length(beta_sources)
    beta_vals = beta_sources{i_src};
    beta_all  = nan(num_rows, target_dim);

    for nsub = 1:num_subjects
        for ntask = 1:num_tasks
            row = (nsub-1) * num_tasks + ntask;

            % Collect betas matching this task across all sessions
            betas_collected = [];
            for nses = 1:num_sessions
                if isempty(beta_vals{nsub, nses}), continue; end
                task_match = contains(task_list{nsub, nses}, task_types{ntask});
                if any(task_match)
                    betas_collected = [betas_collected; beta_vals{nsub, nses}(task_match, :)]; %#ok<AGROW>
                end
            end

            beta_all(row, :) = mean(betas_collected, 1);
            sub_idx(row)     = nsub;
            task_labels{row} = task_types{ntask};
        end
    end

    beta_agg.(beta_names{i_src}) = beta_all;
end

%% Leave-one-subject-out classification
feature_sets = {
    [beta_agg.cos(:, flow_include), beta_agg.sin(:, flow_include)],  'Phase (cos+sin)';
    beta_agg.ampl(:, flow_include),                                   'Amplitude';
    beta_agg.ic,                                                      'IN activations';
    [beta_agg.cos(:, flow_include), beta_agg.sin(:, flow_include), ...
     beta_agg.ampl(:, flow_include)],                                 'Phase + Amplitude';
};

fprintf('\n========== TASK DECODING (LOSO-CV) ==========\n');
fprintf('%-25s  Mean Accuracy\n', 'Feature Set');
fprintf('%s\n', repmat('-', 1, 45));

for i_feat = 1:size(feature_sets, 1)
    X_all = feature_sets{i_feat, 1};
    feat_name = feature_sets{i_feat, 2};

    accuracy_CV = zeros(num_subjects, 1);
    for nsub = 1:num_subjects
        test_mask  = (sub_idx == nsub);
        train_mask = ~test_mask;

        Mdl    = fitcecoc(X_all(train_mask, :), task_labels(train_mask));
        Y_pred = predict(Mdl, X_all(test_mask, :));
        Y_test = task_labels(test_mask);

        accuracy_CV(nsub) = mean(strcmp(Y_test, Y_pred));
    end

    fprintf('%-25s  %.3f +/- %.3f\n', feat_name, mean(accuracy_CV), std(accuracy_CV));
end