clear; clc;

mkdir('figures');

%% Settings
rest_num = 'REST1';
% rest_num = 'REST2';
target_dim = 27;
TRtarget = 0.72;
fdr_threshold = 0.05;

%% Load group-level results (eigenvalues)
g_ica_result_files = dir(sprintf('results/INF_G_lev_%s_ALL_027_MIGP_results_*.mat', rest_num));
load_results = load(fullfile(g_ica_result_files(end).folder, g_ica_result_files(end).name));
lambda = load_results.lambda;

abs_DM    = abs(lambda);
period_DM = 2*pi*TRtarget ./ angle(lambda);

%% Select non-redundant modes (one per conjugate pair + aperiodic modes)
thres_period = 1e10;
dup_idx = find(period_DM <= thres_period);           % conjugate pairs
select_idx = dup_idx(1:2:end);                        % one per pair
select_idx = [select_idx; find(period_DM > thres_period)];  % + aperiodic
num_conjugate = length(dup_idx) / 2;

%% Load temporal fingerprints
result_files = dir(sprintf('results/Temporal_Fingerprints_%s_ALL_027_results_20*.mat', rest_num));
load_result = load(fullfile(result_files(end).folder, result_files(end).name));

sub_ids = load_result.sub_ids;

engagement_level_list = load_result.engagement_level_list(:, select_idx);
persistent_rate_list  = load_result.persistent_rate_list(:, select_idx);
progression_rate_list = load_result.progression_rate_list(:, select_idx);

% Double amplitude for conjugate pair modes
engagement_level_list(:, 1:num_conjugate) = 2 * engagement_level_list(:, 1:num_conjugate);

% Sort modes by mean engagement (descending)
[~, sort_idx] = sort(mean(engagement_level_list), 'descend');
engagement_level_list = engagement_level_list(:, sort_idx);
persistent_rate_list  = persistent_rate_list(:, sort_idx);
progression_rate_list = progression_rate_list(:, sort_idx);

%% Load behavioral and demographic data
behav_factor = readtable('./data/scores_04_REST1_2_intersect.csv');

path_info = load('secure_info/path_info.mat');
gene_data_table      = readtable(path_info.gene_data_path, 'VariableNamingRule', 'preserve');
behav_data_table     = readtable(path_info.behav_data_path, 'VariableNamingRule', 'preserve');
freesurfer_data_table = readtable(path_info.freesurfer_data_path, 'VariableNamingRule', 'preserve');

% Filter tables to included subjects
gene_data_table       = gene_data_table(ismember(gene_data_table.Subject, sub_ids), :);
behav_data_table      = behav_data_table(ismember(behav_data_table.Subject, sub_ids), :);
behav_factor          = behav_factor(ismember(behav_factor.Subject, sub_ids), :);
freesurfer_data_table = freesurfer_data_table(ismember(freesurfer_data_table.Subject, sub_ids), :);

% Align subjects across all tables
[~, ia] = intersect(sub_ids, behav_factor.Subject);
engagement_level_list = engagement_level_list(ia, :);
persistent_rate_list  = persistent_rate_list(ia, :);
progression_rate_list = progression_rate_list(ia, :);
sub_ids = sub_ids(ia);

gene_data_table       = sortrows(gene_data_table, 'Subject');
behav_data_table      = sortrows(behav_data_table, 'Subject');
behav_factor          = sortrows(behav_factor, 'Subject');
freesurfer_data_table = sortrows(freesurfer_data_table, 'Subject');

%% Build confound matrix
age        = gene_data_table.Age_in_Yrs;
sex        = strcmp(behav_data_table.Gender, 'F');
ICV        = freesurfer_data_table.FS_IntraCranial_Vol;
TGMV       = freesurfer_data_table.FS_Total_GM_Vol;
handedness = gene_data_table.Handedness;

X_rest = [ones(size(age)), age, sex, ICV.^(1/3), TGMV.^(1/3), handedness];
X_rest = X_rest(ia);

% Remove subjects with NaN confounds
nan_mask = any(isnan(X_rest), 2);
X_rest(nan_mask, :)                = [];
engagement_level_list(nan_mask, :) = [];
persistent_rate_list(nan_mask, :)  = [];
progression_rate_list(nan_mask, :) = [];
behav_factor(nan_mask, :)          = [];
gene_data_table(nan_mask, :)       = [];
behav_data_table(nan_mask, :)      = [];

FA_values = table2array(behav_factor(:, 2:end));
num_FA = size(FA_values, 2);
num_DM = size(engagement_level_list, 2);

fprintf('Analysis with %d subjects, %d modes, %d behavioral factors\n', size(FA_values,1), num_DM, num_FA);

%% Partial correlation analysis (confound regression via GLM)
fprintf('Running partial correlation analysis...\n');

metric_names = {'engagement', 'persistent', 'progression'};
metric_data  = {engagement_level_list, persistent_rate_list, progression_rate_list};

r_values = zeros(num_DM, num_FA, 3);
t_values = zeros(num_DM, num_FA, 3);
p_values = zeros(num_DM, num_FA, 3);

for i_metric = 1:3
    for i_DM = 1:num_DM
        for i_FA = 1:num_FA
            [r, t, p] = run_glm_correlation(X_rest, metric_data{i_metric}(:, i_DM), FA_values(:, i_FA));
            r_values(i_DM, i_FA, i_metric) = r;
            t_values(i_DM, i_FA, i_metric) = t;
            p_values(i_DM, i_FA, i_metric) = p;
        end
    end
end
fprintf('Correlation analysis complete.\n');

%% FDR correction
% Remove aperiodic mode (last) from significance testing
p_for_fdr = p_values(1:end-1, :, :);

% Joint FDR across all three metrics
all_p = p_for_fdr(:);
fdr_q_all = mafdr(all_p, 'BHFDR', true);
q_values = reshape(fdr_q_all, num_DM-1, num_FA, 3);

% Separate FDR for engagement only
fdr_q_engagement_alone = mafdr(reshape(p_for_fdr(:,:,1),[],1), 'BHFDR', true);
q_engagement_alone = reshape(fdr_q_engagement_alone, num_DM-1, num_FA);

fprintf('FDR correction complete.\n');

%% Print significant results
fprintf('\n===============================================================\n');
fprintf('           SIGNIFICANT CORRELATIONS (q < %.2f)\n', fdr_threshold);
fprintf('===============================================================\n');

titles = {'Amplitude', 'Persistence', 'Speed'};
for i_metric = 1:3
    print_results(titles{i_metric}, ...
        r_values(:,:,i_metric), t_values(:,:,i_metric), ...
        p_values(1:end-1,:,i_metric), q_values(:,:,i_metric), fdr_threshold);
end

%% Heatmap — 1 x 3 layout
col_names = {'Mental Health', 'Cognition', 'Processing Speed', 'Substance Use'};
DM_names = arrayfun(@(i) sprintf('INF %02d', i), 1:num_conjugate, 'UniformOutput', false);

% Coolwarm colormap
nColors = 256;
positions = [0, 0.5, 1];
colors = [0.2298, 0.2980, 0.7530; 1, 1, 1; 0.7530, 0.1504, 0.0980];
x_interp = linspace(0, 1, nColors);
coolwarm = interp1(positions, colors, x_interp, 'linear');
save('results/colormap_coolwarm.mat', 'coolwarm');

% Common color scale
t_oscillatory = t_values(1:num_conjugate, :, :);
t_max = max(abs(t_oscillatory(:)));

figure('Position', [0, 0, 1000, 400]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for i_metric = 1:3
    nexttile;

    t_mat = t_values(1:num_conjugate, :, i_metric);

    % Use engagement-alone FDR for amplitude, joint FDR for others
    if i_metric == 1
        q_mat = q_engagement_alone(1:num_conjugate, :);
    else
        q_mat = q_values(1:num_conjugate, :, i_metric);
    end

    imagesc(t_mat);
    clim([-t_max, t_max]);
    colormap(coolwarm);
    colorbar;

    ax = gca;
    ax.FontName = 'Times New Roman';
    ax.FontSize = 14;
    ax.XTick = 1:numel(col_names);
    ax.XTickLabel = col_names;
    ax.YTick = 1:numel(DM_names);
    ax.YTickLabel = DM_names;
    title(titles{i_metric}, 'FontSize', 18, 'FontName', 'Times New Roman');

    % Overlay significance stars
    hold on;
    for i = 1:size(q_mat, 1)
        for j = 1:size(q_mat, 2)
            if q_mat(i,j) < 0.001
                star_str = '***';
            elseif q_mat(i,j) < 0.01
                star_str = '**';
            elseif q_mat(i,j) < 0.05
                star_str = '*';
            else
                continue;
            end
            text(j, i + 0.1, star_str, 'HorizontalAlignment', 'center', ...
                'FontSize', 15, 'Color', 'k');
        end
    end
    hold off;
end

print(gcf, fullfile('figures', sprintf('Behavior_sig_%s.jpg', rest_num)), '-djpeg', '-r600');

save(sprintf('results/Cognitive_relevance_%s.mat', rest_num), 't_values', 'q_values', 'q_engagement_alone');

%% ======================== Helper Functions ========================

function print_results(label, r_vals, t_vals, p_vals, q_vals, threshold)
    [rows, cols] = find(q_vals < threshold);
    if isempty(rows)
        fprintf('\n--- No significant results for %s ---\n', label);
    else
        fprintf('\n--- Significant Results for %s ---\n', label);
        fprintf('%-8s %-8s %-12s %-12s %-12s %-12s\n', 'DM_Idx', 'FA_Idx', 'r-value', 't-value', 'p-value', 'q-value');
        for i = 1:length(rows)
            fprintf('%-8d %-8d %-12.4f %-12.4f %-12.4e %-12.4e\n', ...
                rows(i), cols(i), r_vals(rows(i),cols(i)), t_vals(rows(i),cols(i)), ...
                p_vals(rows(i),cols(i)), q_vals(rows(i),cols(i)));
        end
    end
end

function [r, t, p] = run_glm_correlation(X_confound, y_metric, x_behav)
    X_full = [X_confound, x_behav];
    [~, ~, stats] = glmfit(X_full, y_metric, 'normal', 'constant', 'off');
    t = stats.t(end);
    p = stats.p(end);
    df_residual = stats.dfe;
    r = t / sqrt(t^2 + df_residual);
end