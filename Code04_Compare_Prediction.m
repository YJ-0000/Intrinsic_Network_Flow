%% 
clear; clc;

mkdir('figures');

sample_set = 'discovery';
% sample_set = 'replication';

%% Load G-ICA-DMD results
target_dim_list = 10:10:100;
num_subjects = 105;
predict_time_window_list = [1,2,4,8];

n_dims    = length(target_dim_list);
n_windows = length(predict_time_window_list);

% Group-level methods: DR, noSF (GroupSF), noTF (GroupTF)
all_R2_DR   = zeros(n_dims, num_subjects, n_windows);
all_R2_noSF = zeros(n_dims, num_subjects, n_windows);
all_R2_noTF = zeros(n_dims, num_subjects, n_windows);
all_R2_null = zeros(n_dims, num_subjects, n_windows);

mean_R2_DR   = zeros(n_windows, n_dims);
mean_R2_noSF = zeros(n_windows, n_dims);
mean_R2_noTF = zeros(n_windows, n_dims);
mean_R2_null = zeros(n_windows, n_dims);

std_R2_DR   = zeros(n_windows, n_dims);
std_R2_noSF = zeros(n_windows, n_dims);
std_R2_noTF = zeros(n_windows, n_dims);
std_R2_null = zeros(n_windows, n_dims);

for i_dim = 1:n_dims
    target_dim = target_dim_list(i_dim);

    if strcmp(sample_set, 'discovery')
        g_ica_result_files = dir(sprintf('results/disc_INF_G_lev_%03d_MIGP_results_*.mat', target_dim));
    elseif strcmp(sample_set, 'replication')
        g_ica_result_files = dir(sprintf('results/repl_INF_G_lev_%03d_MIGP_results_*.mat', target_dim));
    else
        error('Undefined sample set!!');
    end
    load_results = load(fullfile(g_ica_result_files(end).folder, g_ica_result_files(end).name));

    % DR (full spatiotemporal fingerprinting)
    R2_temp = load_results.R2_DM_DR_array_list;
    mean_R2_DR(:, i_dim)   = squeeze(mean(R2_temp(:,:,1,:), [1,2]));
    std_R2_DR(:, i_dim)    = squeeze(std(squeeze(mean(R2_temp(:,:,1,:), 2)), 0, 1));
    all_R2_DR(i_dim, :, :) = squeeze(mean(R2_temp(:,:,1,:), 2));

    % noSF (GroupSF — no spatial fingerprinting)
    R2_temp = load_results.R2_DM_GroupSF_array_list;
    mean_R2_noSF(:, i_dim)   = squeeze(mean(R2_temp(:,:,1,:), [1,2]));
    std_R2_noSF(:, i_dim)    = squeeze(std(squeeze(mean(R2_temp(:,:,1,:), 2)), 0, 1));
    all_R2_noSF(i_dim, :, :) = squeeze(mean(R2_temp(:,:,1,:), 2));

    % noTF (GroupTF — no temporal fingerprinting)
    R2_temp = load_results.R2_DM_GroupTF_array_list;
    mean_R2_noTF(:, i_dim)   = squeeze(mean(R2_temp(:,:,1,:), [1,2]));
    std_R2_noTF(:, i_dim)    = squeeze(std(squeeze(mean(R2_temp(:,:,1,:), 2)), 0, 1));
    all_R2_noTF(i_dim, :, :) = squeeze(mean(R2_temp(:,:,1,:), 2));

    % Null
    R2_temp = load_results.R2_null_array_list;
    mean_R2_null(:, i_dim)   = squeeze(mean(R2_temp(:,:,1,:), [1,2]));
    std_R2_null(:, i_dim)    = squeeze(std(squeeze(mean(R2_temp(:,:,1,:), 2)), 0, 1));
    all_R2_null(i_dim, :, :) = squeeze(mean(R2_temp(:,:,1,:), 2));
end

%% Load subject-wise ICA-DMD results
all_R2_sub_ica  = zeros(n_dims, num_subjects, n_windows);
mean_R2_sub_ica = zeros(n_windows, n_dims);
std_R2_sub_ica  = zeros(n_windows, n_dims);

for i_dim = 1:n_dims
    target_dim = target_dim_list(i_dim);

    if strcmp(sample_set, 'discovery')
        sub_result_files = dir(sprintf('results/disc_subject_wise_ica%03d_dmd_results_*.mat', target_dim));
    elseif strcmp(sample_set, 'replication')
        sub_result_files = dir(sprintf('results/repl_subject_wise_ica%03d_dmd_results_*.mat', target_dim));
    else
        error('Undefined sample set!!');
    end
    load_results = load(fullfile(sub_result_files(end).folder, sub_result_files(end).name));

    R2_temp = load_results.R2_DM_array_list;
    mean_R2_sub_ica(:, i_dim)   = squeeze(mean(R2_temp(:,:,1,:), [1,2]));
    std_R2_sub_ica(:, i_dim)    = squeeze(std(squeeze(mean(R2_temp(:,:,1,:), 2)), 0, 1));
    all_R2_sub_ica(i_dim, :, :) = squeeze(mean(R2_temp(:,:,1,:), 2));
end

%% DR vs Null — t-tests
fprintf('\n=== DR vs Null ===\n');
for i_dim = 1:n_dims
    for i_win = 1:1  % 1-step ahead only
        x1 = squeeze(all_R2_DR(i_dim, :, i_win));
        x2 = squeeze(all_R2_null(i_dim, :, i_win));
        [~, p] = ttest(x1, x2);
        fprintf('dim=%d, win=%d: p = %.12f\n', target_dim_list(i_dim), predict_time_window_list(i_win), p);
    end
end

%% Ablation: DR vs noSF vs noTF — paired t-tests & heatmaps
% DR vs noSF
ablation_pairs = {
    'DR',   'noSF',  all_R2_DR,   all_R2_noSF;
    'DR',   'noTF',  all_R2_DR,   all_R2_noTF;
    'noSF', 'noTF',  all_R2_noSF, all_R2_noTF;
};

% Diverging colormap: blue -> white -> red
ncol = 256; half = ncol/2;
blue_cmap = [linspace(0,1,half)' linspace(0,1,half)' ones(half,1)];
red_cmap  = [ones(half,1) linspace(1,0,half)' linspace(1,0,half)'];
divCmap   = [blue_cmap; red_cmap];

for ip = 1:size(ablation_pairs, 1)
    name1 = ablation_pairs{ip, 1};
    name2 = ablation_pairs{ip, 2};
    data1 = ablation_pairs{ip, 3};
    data2 = ablation_pairs{ip, 4};

    tmat = zeros(n_dims, n_windows);
    pmat = zeros(n_dims, n_windows);

    for i_dim = 1:n_dims
        for i_win = 1:n_windows
            x1 = squeeze(data1(i_dim, :, i_win));
            x2 = squeeze(data2(i_dim, :, i_win));
            [~, p, ~, stats] = ttest(x1, x2);
            tmat(i_dim, i_win) = stats.tstat;
            pmat(i_dim, i_win) = p;
        end
    end

    maxT = max(abs(tmat(:)));

    figure('Position', [100 100 900 600]);
    tiledlayout(2, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    % t-value heatmap
    ax1 = nexttile(1);
    imagesc(tmat.');
    ax1.Colormap = divCmap;
    caxis([-maxT, maxT]);
    colorbar;
    set(ax1, 'XTick', 1:n_dims, 'XTickLabel', target_dim_list, ...
             'YTick', 1:n_windows, 'YTickLabel', predict_time_window_list, 'FontSize', 12);
    title(sprintf('T-value: %s vs. %s', name1, name2), 'FontSize', 22, 'FontName', 'Times New Roman');
    xlabel('Number of ICs (Q)', 'FontSize', 18, 'FontName', 'Times New Roman');
    ylabel('Prediction horizon (TR)', 'FontSize', 18, 'FontName', 'Times New Roman');

    hold(ax1, 'on');
    for i_win = 1:n_windows
        for i_dim = 1:n_dims
            text(i_dim, i_win, sprintf('%.2f', tmat(i_dim, i_win)), ...
                'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'k');
        end
    end
    hold(ax1, 'off');

    % p-value heatmap
    ax2 = nexttile(2);
    imagesc(pmat.');
    ax2.Colormap = copper;
    caxis([0, 1]);
    colorbar;
    set(ax2, 'XTick', 1:n_dims, 'XTickLabel', target_dim_list, ...
             'YTick', 1:n_windows, 'YTickLabel', predict_time_window_list, 'FontSize', 12);
    title(sprintf('p-value: %s vs. %s', name1, name2), 'FontSize', 22, 'FontName', 'Times New Roman');
    xlabel('Number of ICs (Q)', 'FontSize', 18, 'FontName', 'Times New Roman');
    ylabel('Prediction horizon (TR)', 'FontSize', 18, 'FontName', 'Times New Roman');

    hold(ax2, 'on');
    for i_win = 1:n_windows
        for i_dim = 1:n_dims
            text(i_dim, i_win, sprintf('%.4f', pmat(i_dim, i_win)), ...
                'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'w');
        end
    end
    hold(ax2, 'off');

    print(gcf, sprintf('figures/ablation_%s_vs_%s.png', name1, name2), '-dpng', '-r600');
end

%% DR vs subject-wise ICA — paired t-tests & heatmap
tmat = zeros(n_dims, n_windows);
pmat = zeros(n_dims, n_windows);

for i_dim = 1:n_dims
    for i_win = 1:n_windows
        x1 = squeeze(all_R2_DR(i_dim, :, i_win));
        x2 = squeeze(all_R2_sub_ica(i_dim, :, i_win));
        [~, p, ~, stats] = ttest(x1, x2);
        tmat(i_dim, i_win) = stats.tstat;
        pmat(i_dim, i_win) = p;
    end
end

maxT = max(abs(tmat(:)));

figure('Position', [100 100 900 600]);
tiledlayout(2, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact');

ax1 = nexttile(1);
imagesc(tmat.');
ax1.Colormap = divCmap;
caxis([-maxT, maxT]);
colorbar;
set(ax1, 'XTick', 1:n_dims, 'XTickLabel', target_dim_list, ...
         'YTick', 1:n_windows, 'YTickLabel', predict_time_window_list, 'FontSize', 12);
title('T-value: DR vs. subject-wise ICA', 'FontSize', 22, 'FontName', 'Times New Roman');
xlabel('Number of ICs (Q)', 'FontSize', 18, 'FontName', 'Times New Roman');
ylabel('Prediction horizon (TR)', 'FontSize', 18, 'FontName', 'Times New Roman');

hold(ax1, 'on');
for i_win = 1:n_windows
    for i_dim = 1:n_dims
        text(i_dim, i_win, sprintf('%.2f', tmat(i_dim, i_win)), ...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'k');
    end
end
hold(ax1, 'off');

ax2 = nexttile(2);
imagesc(pmat.');
ax2.Colormap = copper;
caxis([0, 1]);
colorbar;
set(ax2, 'XTick', 1:n_dims, 'XTickLabel', target_dim_list, ...
         'YTick', 1:n_windows, 'YTickLabel', predict_time_window_list, 'FontSize', 12);
title('p-value: DR vs. subject-wise ICA', 'FontSize', 22, 'FontName', 'Times New Roman');
xlabel('Number of ICs (Q)', 'FontSize', 18, 'FontName', 'Times New Roman');
ylabel('Prediction horizon (TR)', 'FontSize', 18, 'FontName', 'Times New Roman');

hold(ax2, 'on');
for i_win = 1:n_windows
    for i_dim = 1:n_dims
        text(i_dim, i_win, sprintf('%.4f', pmat(i_dim, i_win)), ...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'w');
    end
end
hold(ax2, 'off');

print(gcf, 'figures/comparison_DR_vs_subwise_ICA.png', '-dpng', '-r600');

%% Comparing best method across Q — 1-step ahead
fprintf('\n=== Best method per Q (1-step ahead) ===\n');
methods_all  = {all_R2_DR, all_R2_noSF, all_R2_noTF, all_R2_sub_ica};
method_names = {'DR', 'noSF', 'noTF', 'subject-wise ICA'};

for i_dim = 1:n_dims
    temp_means = zeros(1, length(methods_all));
    for i_m = 1:length(methods_all)
        temp_means(i_m) = mean(squeeze(methods_all{i_m}(i_dim, :, 1)));
    end
    [~, idx_best] = max(temp_means);
    fprintf('Q=%03d — Best: %s (mean R2=%.4f)\n', target_dim_list(i_dim), method_names{idx_best}, temp_means(idx_best));
end

%% Paired t-tests across Q for DR — heatmap
p_mat_Q = zeros(n_dims, n_dims);
for i_dim = 1:n_dims
    for j_dim = 1:n_dims
        x_i = squeeze(all_R2_DR(i_dim, :, 1));
        x_j = squeeze(all_R2_DR(j_dim, :, 1));
        [~, p_mat_Q(i_dim, j_dim)] = ttest(x_i, x_j);
    end
end

% Mask lower triangle
mask_tri = tril(true(n_dims), -1);
p_mat_Q(mask_tri) = NaN;

figure('Position', [100 100 800 600]);
h = heatmap(string(target_dim_list), string(target_dim_list), p_mat_Q, ...
    'ColorLimits', [0 1], 'CellLabelFormat', '%.4f');
h.Title    = sprintf('p-values (1 TR-ahead prediction, DR)');
h.XLabel   = 'Number of ICs (Q)';
h.YLabel   = 'Number of ICs (Q)';
h.FontSize = 18;
h.FontName = 'Times New Roman';

print(gcf, 'figures/DR_Q_comparison.png', '-dpng', '-r600');

%% Bar plot: DR vs noSF vs noTF vs subject-wise ICA
figure;
x_labels = {'predict 1s ahead', 'predict 2s ahead', 'predict 4s ahead', 'predict 8s ahead'};
grey_shades = [0.8 0.8 0.8; 0.6 0.6 0.6; 0.4 0.4 0.4; 0.2 0.2 0.2];

for i_pred = 1:4
    subplot(4, 1, i_pred);

    data = [mean_R2_DR(i_pred,:)', mean_R2_noSF(i_pred,:)', mean_R2_noTF(i_pred,:)', mean_R2_sub_ica(i_pred,:)'];
    h = bar(target_dim_list, data);
    hold on;

    numGroups = size(data, 1);
    for s = 1:4
        h(s).FaceColor = 'flat';
        h(s).CData = repmat(grey_shades(s,:), numGroups, 1);
    end

    [~, linearIdx] = max(data(:));
    [maxGroup, maxSeries] = ind2sub(size(data), linearIdx);
    h(maxSeries).CData(maxGroup, :) = [1 0 0];

    title(x_labels{i_pred});
    ylabel('Mean R^2');
    legend({'DR', 'noSF', 'noTF', 'subject-wise ICA'}, 'Location', 'southeast');
    ylim([0, inf]);
    hold off;
end

print(gcf, 'figures/bar_comparison_all_methods.png', '-dpng', '-r600');

%% Interpolation to find optimal Q per method
x  = 10:10:100;
xq = 10:1:100;

fprintf('\n=== Optimal Q (interpolated, 1-step ahead) ===\n');
interp_methods = {mean_R2_DR, mean_R2_noSF, mean_R2_noTF, mean_R2_sub_ica};
interp_names   = {'DR', 'noSF', 'noTF', 'subject-wise ICA'};

for i_pred = 1:4
    fprintf('predict %ds ahead: ', i_pred);
    for i_m = 1:length(interp_methods)
        yq = interp1(x, interp_methods{i_m}(i_pred, :), xq, 'spline');
        [~, idx] = max(yq);
        fprintf('%s max at Q=%d, ', interp_names{i_m}, xq(idx));
    end
    fprintf('\n');
end

%% Fine-grained Q search (DR only)
target_dim_fine_list = 20:30;
n_dims_fine = length(target_dim_fine_list);

all_R2_DR_fine  = zeros(n_dims_fine, num_subjects, n_windows);
mean_R2_DR_fine = zeros(n_windows, n_dims_fine);
std_R2_DR_fine  = zeros(n_windows, n_dims_fine);

for i_dim = 1:n_dims_fine
    target_dim = target_dim_fine_list(i_dim);

    if strcmp(sample_set, 'discovery')
        g_ica_result_files = dir(sprintf('results/disc_INF_G_lev_%03d_MIGP_results_*.mat', target_dim));
    elseif strcmp(sample_set, 'replication')
        g_ica_result_files = dir(sprintf('results/repl_INF_G_lev_%03d_MIGP_results_*.mat', target_dim));
    else
        error('Undefined sample set!!');
    end
    load_results = load(fullfile(g_ica_result_files(end).folder, g_ica_result_files(end).name));

    R2_temp = load_results.R2_DM_DR_array_list;
    mean_R2_DR_fine(:, i_dim)       = squeeze(mean(R2_temp(:,:,1,:), [1,2]));
    std_R2_DR_fine(:, i_dim)        = squeeze(std(squeeze(mean(R2_temp(:,:,1,:), 2)), 0, 1));
    all_R2_DR_fine(i_dim, :, :)     = squeeze(mean(R2_temp(:,:,1,:), 2));
end

%% Fine Q — paired t-tests heatmap
p_mat_fine = zeros(n_dims_fine, n_dims_fine);
for i_dim = 1:n_dims_fine
    for j_dim = 1:n_dims_fine
        x_i = squeeze(all_R2_DR_fine(i_dim, :, 1));
        x_j = squeeze(all_R2_DR_fine(j_dim, :, 1));
        [~, p_mat_fine(i_dim, j_dim)] = ttest(x_i, x_j);
    end
end

mask_tri_fine = tril(true(n_dims_fine), -1);
p_mat_fine(mask_tri_fine) = NaN;

figure('Position', [100 100 800 600]);
h = heatmap(string(target_dim_fine_list), string(target_dim_fine_list), p_mat_fine, ...
    'ColorLimits', [0 1], 'CellLabelFormat', '%.4f');
h.Title    = sprintf('p-values (1 TR-ahead, DR, fine Q)');
h.XLabel   = 'Number of ICs (Q)';
h.YLabel   = 'Number of ICs (Q)';
h.FontSize = 18;
h.FontName = 'Times New Roman';

print(gcf, 'figures/DR_Q_comparison_fine.png', '-dpng', '-r600');

%% Fine Q — bar plot (DR only)
figure;
for i_pred = 1:4
    subplot(4, 1, i_pred);

    data = mean_R2_DR_fine(i_pred, :)';
    h = bar(target_dim_fine_list, data);
    hold on;

    h.FaceColor = 'flat';
    h.CData = repmat([0.6 0.6 0.6], length(data), 1);

    [~, idx_max] = max(data);
    h.CData(idx_max, :) = [1 0 0];

    title(x_labels{i_pred});
    ylabel('Mean R^2');
    ylim([0, inf]);
    hold off;
end

print(gcf, 'figures/bar_DR_fine_Q.png', '-dpng', '-r600');