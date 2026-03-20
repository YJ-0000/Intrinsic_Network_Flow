clear; clc;

%% Load per-subject prediction results
result_files = dir('results/MDTB_BOLD_prediction_sub*');
num_subjects = length(result_files);

% Load first file to get method names and count
res0 = load(fullfile(result_files(1).folder, result_files(1).name));
method_names = fieldnames(res0);
num_run  = size(res0.R2_all,2);

% Collect mean R2 per subject per method
R2_INF_per_subject = nan(num_subjects, num_run);
R2_null_per_subject = nan(num_subjects, num_run);

for nsub = 1:num_subjects
    res = load(fullfile(result_files(nsub).folder, result_files(nsub).name));
    % R2_all is (num_methods x num_runs) — average across runs per subject
    R2_INF_per_subject(nsub, :) = res.R2_all(2,:);
    R2_null_per_subject(nsub, :) = res.R2_all(1,:);
end

%% Summary
fprintf('\n========== MDTB BOLD Prediction Summary (%d subjects) ==========\n', num_subjects);
fprintf('%-15s  %8s  %8s  %8s\n', 'Method', 'Mean R2', 'Std', '95% CI');
fprintf('%s\n', repmat('-', 1, 50));

inf_vals = mean(R2_INF_per_subject,2);
mu  = mean(inf_vals);
sd  = std(inf_vals);
ci  = 1.96 * sd / sqrt(num_subjects);
fprintf('%-15s  %8.4f  %8.4f  %8.4f\n', 'R2 (INF)', mu, sd, ci);

null_vals = mean(R2_null_per_subject,2);
mu  = mean(null_vals);
sd  = std(null_vals);
ci  = 1.96 * sd / sqrt(num_subjects);
fprintf('%-15s  %8.4f  %8.4f  %8.4f\n', 'Null', mu, sd, ci);


%% Pairwise comparisons against null
fprintf('\n--- Paired t-tests INF vs. Null ---\n');

[~, p, ~, st] = ttest(inf_vals, null_vals);
fprintf('%-15s vs Null: t=%.4f, p=%.6f\n', 'R2 (INF)', st.tstat, p);
