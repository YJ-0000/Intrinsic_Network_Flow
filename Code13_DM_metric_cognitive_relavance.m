clear; clc;

%% load DM metrics
result_files = dir('results/DM_metric_*.mat');
load_result = load(fullfile(result_files(end).folder,result_files(end).name));
sub_ids = load_result.sub_ids(load_result.sub_idx_include);

target_dim = 27;
TRtarget = 0.72;
g_ica_result_files = dir(['results/loop_gica_',num2str(target_dim,'%03d'),'_dmd_results_normalized_*.mat']);
load_results = load(fullfile(g_ica_result_files(end).folder,g_ica_result_files(end).name));
lambda = load_results.lambda;

abs_DM = abs(lambda);
period_DM = 2*pi*TRtarget ./ angle(lambda);

thres_period = 1e10;
dup_idx = find(period_DM <= thres_period);
select_idx = dup_idx(1:2:end);
select_idx = [select_idx; find(period_DM > thres_period)];

engagement_level_list = load_result.engagement_level_list;
persistent_rate_list = load_result.persistent_rate_list;
progression_rate_list = load_result.progression_rate_list;

engagement_level_list_tl_cc = squeeze(engagement_level_list(2,:,select_idx));
persistent_rate_list_tl_cc = squeeze(persistent_rate_list(2,:,select_idx));
progression_rate_list_tl_cc = squeeze(progression_rate_list(2,:,select_idx));
engagement_level_list_dr = squeeze(engagement_level_list(1,:,select_idx));
persistent_rate_list_dr = squeeze(persistent_rate_list(1,:,select_idx));
progression_rate_list_dr = squeeze(progression_rate_list(1,:,select_idx));

engagement_level_list_tl_cc(:,1:length(dup_idx)/2) = 2 * engagement_level_list_tl_cc(:,1:length(dup_idx)/2);
engagement_level_list_dr(:,1:length(dup_idx)/2) = 2 * engagement_level_list_dr(:,1:length(dup_idx)/2);

engagement_level_list_tl_cc = engagement_level_list_tl_cc ./ sum(engagement_level_list_tl_cc,2);
engagement_level_list_dr = engagement_level_list_dr ./ sum(engagement_level_list_dr,2);

[~,sort_idx] = sort(mean(engagement_level_list_dr),'descend');
engagement_level_list_tl_cc = engagement_level_list_tl_cc(:,sort_idx);
persistent_rate_list_tl_cc = persistent_rate_list_tl_cc(:,sort_idx);
progression_rate_list_tl_cc = progression_rate_list_tl_cc(:,sort_idx);
engagement_level_list_dr = engagement_level_list_dr(:,sort_idx);
persistent_rate_list_dr = persistent_rate_list_dr(:,sort_idx);
progression_rate_list_dr = progression_rate_list_dr(:,sort_idx);

%% load behavior
behav_factor = readtable('./results/scores_04_selected.csv');
load secure_info/path_info;
gene_data_table = readtable(gene_data_path,'VariableNamingRule','preserve');
behav_data_table = readtable(behav_data_path,'VariableNamingRule','preserve');
freesurfer_data_table = readtable(freesurfer_data_path,'VariableNamingRule','preserve');
for nrow = size(gene_data_table,1):-1:1
    if ~any(sub_ids==gene_data_table(nrow,'Subject').Variables)
        gene_data_table(nrow,:) = [];
    end
end
for nrow = size(behav_data_table,1):-1:1
    if ~any(sub_ids==behav_data_table(nrow,'Subject').Variables)
        behav_data_table(nrow,:) = [];
    end
end
for nrow = size(behav_factor,1):-1:1
    if ~any(sub_ids==behav_factor(nrow,'Subject').Variables)
        behav_factor(nrow,:) = [];
    end
end

for nrow = size(freesurfer_data_table,1):-1:1
    if ~any(sub_ids==freesurfer_data_table(nrow,'Subject').Variables)
        freesurfer_data_table(nrow,:) = [];
    end
end
gene_data_table = sortrows(gene_data_table, 'Subject');
behav_data_table = sortrows(behav_data_table, 'Subject');
behav_factor = sortrows(behav_factor, 'Subject');
freesurfer_data_table = sortrows(freesurfer_data_table, 'Subject');

age = gene_data_table.Age_in_Yrs;
sex = strcmp(behav_data_table.Gender,'F');
bmi = gene_data_table.BMI;
race1 = strcmp(gene_data_table.Race,'Am. Indian/Alaskan Nat.');
race2 = strcmp(gene_data_table.Race,'Asian/Nat. Hawaiian/Othr Pacific Is.');
race3 = strcmp(gene_data_table.Race,'Black or African Am.');
race4 = strcmp(gene_data_table.Race,'More than one');
race5 = strcmp(gene_data_table.Race,'Unknown or Not Reported');
race6 = strcmp(gene_data_table.Race,'White');
ICV = freesurfer_data_table.FS_IntraCranial_Vol;
TGMV = freesurfer_data_table.FS_Total_GM_Vol;
handedness = gene_data_table.Handedness;
edu_year = gene_data_table.SSAGA_Educ;

X_rest = [ones(size(age,1),1),age,sex, ICV.^(1/3),TGMV.^(1/3),handedness];
% X_rest = [ones(size(age,1),1),age,sex, ICV.^(1/3),TGMV.^(1/3)];
% X_rest = [ones(size(age,1),1),age,sex];
% X_rest = [ones(size(age,1),1)];


nan_sub = sum(isnan(X_rest),2)>0;
X_rest(nan_sub,:) = [];
gene_data_table(nan_sub,:) = [];
engagement_level_list_tl_cc(nan_sub,:) = [];
persistent_rate_list_tl_cc(nan_sub,:) = [];
progression_rate_list_tl_cc(nan_sub,:) = [];
engagement_level_list_dr(nan_sub,:) = [];
persistent_rate_list_dr(nan_sub,:) = [];
progression_rate_list_dr(nan_sub,:) = [];
behav_data_table(nan_sub,:) = [];
behav_factor(nan_sub,:) = [];

%%
num_FA = size(behav_factor,2)-1;
num_DM = size(progression_rate_list_dr,2);
t_values_engagement_tl_cc = zeros(num_DM, num_FA);
t_values_persistent_tl_cc = zeros(num_DM, num_FA);
t_values_progression_tl_cc = zeros(num_DM, num_FA);
p_values_engagement_tl_cc = zeros(num_DM, num_FA);
p_values_persistent_tl_cc = zeros(num_DM, num_FA);
p_values_progression_tl_cc = zeros(num_DM, num_FA);

r_values_engagement_tl_cc = zeros(num_DM, num_FA);
r_values_persistent_tl_cc = zeros(num_DM, num_FA);
r_values_progression_tl_cc = zeros(num_DM, num_FA);

t_values_engagement_dr = zeros(num_DM, num_FA);
t_values_persistent_dr = zeros(num_DM, num_FA);
t_values_progression_dr = zeros(num_DM, num_FA);
p_values_engagement_dr = zeros(num_DM, num_FA);
p_values_persistent_dr = zeros(num_DM, num_FA);
p_values_progression_dr = zeros(num_DM, num_FA);

r_values_engagement_dr = zeros(num_DM, num_FA);
r_values_persistent_dr = zeros(num_DM, num_FA);
r_values_progression_dr = zeros(num_DM, num_FA);

FA_values = table2array(behav_factor(:,2:end));

%% Outlier Detection and Removal
fprintf('Starting outlier detection...\n');

% Get the number of subjects before removal
num_subjects_before = size(FA_values, 1);

% Create a mask to identify subjects who are outliers in ANY of the relevant columns.
% isoutlier uses the median absolute deviation (MAD) method by default, which is robust.
total_outlier_mask = false(num_subjects_before, 1);

% 1. Check for outliers in FA_values
for i = 1:num_FA
    total_outlier_mask = total_outlier_mask | isoutlier(FA_values(:, i));
end

% 2. Check for outliers in all DM metric lists
all_metrics = {
    engagement_level_list_tl_cc, persistent_rate_list_tl_cc, progression_rate_list_tl_cc, ...
    engagement_level_list_dr, persistent_rate_list_dr, progression_rate_list_dr
};

for i = 1:length(all_metrics)
    current_metric_list = all_metrics{i};
    for j = 1:size(current_metric_list, 2)
        total_outlier_mask = total_outlier_mask | isoutlier(current_metric_list(:, j));
    end
end

total_outlier_mask = false(1,size(FA_values,1));
% Remove subjects identified as outliers
FA_values_filt = FA_values(~total_outlier_mask, :);
X_rest_filt = X_rest(~total_outlier_mask, :);
engagement_level_list_tl_cc_filt = engagement_level_list_tl_cc(~total_outlier_mask, :);
persistent_rate_list_tl_cc_filt = persistent_rate_list_tl_cc(~total_outlier_mask, :);
progression_rate_list_tl_cc_filt = progression_rate_list_tl_cc(~total_outlier_mask, :);
engagement_level_list_dr_filt = engagement_level_list_dr(~total_outlier_mask, :);
persistent_rate_list_dr_filt = persistent_rate_list_dr(~total_outlier_mask, :);
progression_rate_list_dr_filt = progression_rate_list_dr(~total_outlier_mask, :);

num_subjects_after = size(FA_values_filt, 1);
fprintf('%d outliers were detected and removed.\n', num_subjects_before - num_subjects_after);
fprintf('Analysis will proceed with %d subjects.\n\n', num_subjects_after);


%% Perform partial correlation analysis using glmfit
fprintf('Starting partial correlation analysis with confound regression...\n');



for i_DM = 1:num_DM
    % Provide progress feedback in the command window
    if mod(i_DM, 5) == 0
        fprintf('Processing DM %d/%d...\n', i_DM, num_DM);
    end
    
    for i_FA = 1:num_FA
        current_FA_values = FA_values_filt(:, i_FA);

        % --- TL-CC Condition ---
        [r, t, p] = run_glm_correlation(X_rest_filt, engagement_level_list_tl_cc_filt(:, i_DM), current_FA_values);
        r_values_engagement_tl_cc(i_DM, i_FA) = r;
        t_values_engagement_tl_cc(i_DM, i_FA) = t;
        p_values_engagement_tl_cc(i_DM, i_FA) = p;

        [r, t, p] = run_glm_correlation(X_rest_filt, persistent_rate_list_tl_cc_filt(:, i_DM), current_FA_values);
        r_values_persistent_tl_cc(i_DM, i_FA) = r;
        t_values_persistent_tl_cc(i_DM, i_FA) = t;
        p_values_persistent_tl_cc(i_DM, i_FA) = p;

        [r, t, p] = run_glm_correlation(X_rest_filt, progression_rate_list_tl_cc_filt(:, i_DM), current_FA_values);
        r_values_progression_tl_cc(i_DM, i_FA) = r;
        t_values_progression_tl_cc(i_DM, i_FA) = t;
        p_values_progression_tl_cc(i_DM, i_FA) = p;

        % --- DR Condition ---
        [r, t, p] = run_glm_correlation(X_rest_filt, engagement_level_list_dr_filt(:, i_DM), current_FA_values);
        r_values_engagement_dr(i_DM, i_FA) = r;
        t_values_engagement_dr(i_DM, i_FA) = t;
        p_values_engagement_dr(i_DM, i_FA) = p;

        [r, t, p] = run_glm_correlation(X_rest_filt, persistent_rate_list_dr_filt(:, i_DM), current_FA_values);
        r_values_persistent_dr(i_DM, i_FA) = r;
        t_values_persistent_dr(i_DM, i_FA) = t;
        p_values_persistent_dr(i_DM, i_FA) = p;

        [r, t, p] = run_glm_correlation(X_rest_filt, progression_rate_list_dr_filt(:, i_DM), current_FA_values);
        r_values_progression_dr(i_DM, i_FA) = r;
        t_values_progression_dr(i_DM, i_FA) = t;
        p_values_progression_dr(i_DM, i_FA) = p;
    end
end
fprintf('Correlation analysis finished.\n\n');

%% Calculate FDR corrected p-values (q-values)
% This requires the Bioinformatics Toolbox for the mafdr function.
% We will use the Benjamini-Hochberg method ('BHFDR').

p_values_progression_tl_cc(end,:) = [];
p_values_progression_dr(end,:) = [];

% --- TL-CC Condition ---
all_p_tl_cc = [p_values_engagement_tl_cc(:); p_values_persistent_tl_cc(:); p_values_progression_tl_cc(:)];
fdr_q_tl_cc = mafdr(all_p_tl_cc, 'BHFDR', true);

num_tests_per_metric = num_DM * num_FA;
q_values_engagement_tl_cc = reshape(fdr_q_tl_cc(1:num_tests_per_metric), num_DM, num_FA);
q_values_persistent_tl_cc = reshape(fdr_q_tl_cc(num_tests_per_metric+1:2*num_tests_per_metric), num_DM, num_FA);
q_values_progression_tl_cc = reshape(fdr_q_tl_cc(2*num_tests_per_metric+1:end), num_DM-1, num_FA);

% --- DR Condition ---
all_p_dr = [p_values_engagement_dr(:); p_values_persistent_dr(:); p_values_progression_dr(:)];
fdr_q_dr = mafdr(all_p_dr, 'BHFDR', true);

q_values_engagement_dr = reshape(fdr_q_dr(1:num_tests_per_metric), num_DM, num_FA);
q_values_persistent_dr = reshape(fdr_q_dr(num_tests_per_metric+1:2*num_tests_per_metric), num_DM, num_FA);
q_values_progression_dr = reshape(fdr_q_dr(2*num_tests_per_metric+1:end), num_DM-1, num_FA);

fprintf('FDR correction finished.\n\n');

%% Print significant results
threshold = 0.05; % Significance threshold for q-values

fprintf('===============================================================\n');
fprintf('           SIGNIFICANT CORRELATIONS (q < %.2f)        \n', threshold);
fprintf('===============================================================\n');

% Print for TL-CC condition
print_results('TL-CC Engagement', r_values_engagement_tl_cc, t_values_engagement_tl_cc, p_values_engagement_tl_cc, q_values_engagement_tl_cc, threshold);
print_results('TL-CC Persistence', r_values_persistent_tl_cc, t_values_persistent_tl_cc, p_values_persistent_tl_cc, q_values_persistent_tl_cc, threshold);
print_results('TL-CC Progression', r_values_progression_tl_cc, t_values_progression_tl_cc, p_values_progression_tl_cc, q_values_progression_tl_cc, threshold);

fprintf('\n-------------------------------------------------------------\n');

% Print for DR condition
print_results('DR Engagement', r_values_engagement_dr, t_values_engagement_dr, p_values_engagement_dr, q_values_engagement_dr, threshold);
print_results('DR Persistence', r_values_persistent_dr, t_values_persistent_dr, p_values_persistent_dr, q_values_persistent_dr, threshold);
print_results('DR Progression', r_values_progression_dr, t_values_progression_dr, p_values_progression_dr, q_values_progression_dr, threshold);

fprintf('\n===============================================================\n');
fprintf('                       ANALYSIS COMPLETE                       \n');
fprintf('===============================================================\n');

% Save all workspace variables to a results file
results_filename = ['results/correlation_analysis_DM', num2str(target_dim), '_outlier_removed_', datestr(now, 'yyyymmdd_HHMMSS'), '.mat'];
save(results_filename);
fprintf('All workspace variables saved to %s\n', results_filename);

%% Heatmap
t_all = {t_values_engagement_tl_cc, t_values_persistent_tl_cc,t_values_progression_tl_cc; 
    t_values_engagement_dr,t_values_persistent_dr,t_values_progression_dr;
    };

q_all = {q_values_engagement_tl_cc, q_values_persistent_tl_cc,q_values_progression_tl_cc; 
    q_values_engagement_dr,q_values_persistent_dr,q_values_progression_dr;
    };

titles = {'Engagement level (normalized)','Persistence rate', 'Progression rate'};

temp_max = max(abs(cell2mat(t_all)),[],'all');
minValue = -temp_max;
maxValue = temp_max;

col_names = {'Mental Health', 'Cognition', 'Processing Speed', 'Substance Use'};
DM_names = cell(1, 14); 
for i = 1:14
    DM_names{i} = sprintf('DM%02d', i);
end

figure; set(gcf, 'Position', [0, 0, 1200, 800]);
plot_idx = 0;
for jj = 1:2
    for ii = 1:3
        t_values_mat = t_all{jj,ii};
        q_FDR_mat = q_all{jj,ii};
        
        plot_idx = plot_idx + 1;
        subplot(2,3,plot_idx);
        % Display the matrix as an image
        if ii ~= 3
            imagesc(t_values_mat);
        else
            imagesc(t_values_mat(1:length(dup_idx)/2,:));
        end
        caxis([minValue maxValue]);

        % Positions (normalized between 0 and 1)
        positions = [0, 0.5, 1];
        % Define the number of colors in the colormap
        nColors = 256;
        % Corresponding RGB colors
        % Blue: [0.2298, 0.2980, 0.7530]
        % White: [1, 1, 1]
        % Red: [0.7530, 0.2314, 0.0980]
        colors = [
            0.2298, 0.2980, 0.7530;    % Blue
            1.0000, 1.0000, 1.0000;    % White
            0.7530, 0.1504, 0.0980     % Red
        ];

        % Create a matrix for interpolation
        % Each row represents R, G, B channels respectively
        R = colors(:,1);
        G = colors(:,2);
        B = colors(:,3);

        % Define the query points for interpolation
        x = linspace(0, 1, nColors);

        % Interpolate each color channel
        R_interp = interp1(positions, R, x, 'linear');
        G_interp = interp1(positions, G, x, 'linear');
        B_interp = interp1(positions, B, x, 'linear');

        % Combine the interpolated channels into a colormap
        coolwarm = [R_interp', G_interp', B_interp'];

        % Apply the custom 'coolwarm' colormap
        colormap(coolwarm);

        save results/colormap_coolwarm coolwarm

        % Add a colorbar to the figure
        colorbar;

        % Get the current axes handle
        ax = gca;

        % Set the font name and size for the axes
        ax.FontName = 'Times New Roman';
        ax.FontSize = 14;

        % Define X-axis ticks and labels
        ax.XTick = 1:numel(col_names);
        ax.XTickLabel = col_names;

        % Optionally, rotate X-axis labels for better readability
        % Uncomment the next line if you want to rotate the labels by 45 degrees
        % ax.XTickLabelRotation = 45;

        % Define Y-axis ticks and labels
        ax.YTick = 1:numel(DM_names);
        ax.YTickLabel = DM_names;
        
        % figure title
        title(titles{ii});

        % Hold the current plot to add text annotations
        hold on;

        y_offset = 0.1;
        for i = 1:size(q_FDR_mat, 1)
            for j = 1:size(q_FDR_mat, 2)
                if q_FDR_mat(i,j) < 0.001
                    text(j, i + y_offset, '***', ...
                        'HorizontalAlignment', 'center', ...
                        'FontSize', 15, ...                    
                        'Color', 'k');                         % Set text color to red
                elseif q_FDR_mat(i,j) < 0.01
                    text(j, i + y_offset, '**', ...
                        'HorizontalAlignment', 'center', ...
                        'FontSize', 15, ...                    
                        'Color', 'k');                         % Set text color to red
                elseif q_FDR_mat(i,j) < 0.05
                    text(j, i + y_offset, '*', ...
                        'HorizontalAlignment', 'center', ...
                        'FontSize', 15, ...                    
                        'Color', 'k');                         % Set text color to red
                end
            end
        end

        % Release the hold on the current plot
        hold off;

        % Optional: Adjust the figure's overall properties for better aesthetics
        set(gcf, 'Color', 'w'); % Set figure background to white
    end
end
print(gcf, fullfile('figures','Behavior_sig.jpg'), '-djpeg', '-r600');

%% =================================================================
%      "TL-CC" vs "DR" Set Comparison
%  =================================================================
fprintf('\n\n===== Detailed Comparison of "TL-CC" vs "DR" sets by Rate Type =====\n');
fprintf('This analysis individually compares t-value distributions for engagement, persistence, and progression.\n');

figure('Name', 'Rate-specific Distribution of Absolute t-values', 'Color', 'w', 'Position', [100, 100, 700, 900]);

%%% 1. Engagement Rate Comparison
% -----------------------------------------------------------------
fprintf('\n--- 1. Comparison for ENGAGEMENT Rates ---\n');

t_tl_cc_eng = abs(t_values_engagement_tl_cc(:));
t_dr_eng = abs(t_values_engagement_dr(:));

[p_val_eng, h_eng] = ranksum(t_tl_cc_eng, t_dr_eng);

% Print
fprintf('Mean absolute t-value for TL-CC Engagement: %.4f\n', mean(t_tl_cc_eng));
fprintf('Mean absolute t-value for DR Engagement:     %.4f\n', mean(t_dr_eng));
fprintf('Wilcoxon rank-sum test p-value:              %.4f\n', p_val_eng);

if h_eng
    fprintf('Conclusion: The difference is statistically significant. ');
    if mean(t_tl_cc_eng) > mean(t_dr_eng)
        fprintf('The "TL-CC" set shows stronger effects for Engagement rates.\n');
    else
        fprintf('The "DR" set shows stronger effects for Engagement rates.\n');
    end
else
    fprintf('Conclusion: The difference is not statistically significant for Engagement rates.\n');
end

% Visualization
subplot(3, 1, 1);
hold on;
histogram(t_tl_cc_eng, 'DisplayName', 'TL-CC', 'Normalization', 'pdf', 'FaceAlpha', 0.6);
histogram(t_dr_eng, 'DisplayName', 'DR', 'Normalization', 'pdf', 'FaceAlpha', 0.6);
title('Engagement Rates: Absolute t-value Distribution', 'FontSize', 14);
xlabel('Absolute t-value');
ylabel('Probability Density');
legend('show');
grid on;
hold off;


%%% 2. Persistent Rate Comparison
% -----------------------------------------------------------------
fprintf('\n--- 2. Comparison for PERSISTENT Rates ---\n');

t_tl_cc_per = abs(t_values_persistent_tl_cc(:));
t_dr_per = abs(t_values_persistent_dr(:));

[p_val_per, h_per] = ranksum(t_tl_cc_per, t_dr_per);

% Print
fprintf('Mean absolute t-value for TL-CC Persistence: %.4f\n', mean(t_tl_cc_per));
fprintf('Mean absolute t-value for DR Persistence:     %.4f\n', mean(t_dr_per));
fprintf('Wilcoxon rank-sum test p-value:               %.4f\n', p_val_per);

if h_per
    fprintf('Conclusion: The difference is statistically significant. ');
    if mean(t_tl_cc_per) > mean(t_dr_per)
        fprintf('The "TL-CC" set shows stronger effects for Persistent rates.\n');
    else
        fprintf('The "DR" set shows stronger effects for Persistent rates.\n');
    end
else
    fprintf('Conclusion: The difference is not statistically significant for Persistent rates.\n');
end

% Visualization
subplot(3, 1, 2);
hold on;
histogram(t_tl_cc_per, 'DisplayName', 'TL-CC', 'Normalization', 'pdf', 'FaceAlpha', 0.6);
histogram(t_dr_per, 'DisplayName', 'DR', 'Normalization', 'pdf', 'FaceAlpha', 0.6);
title('Persistent Rates: Absolute t-value Distribution', 'FontSize', 14);
xlabel('Absolute t-value');
ylabel('Probability Density');
legend('show');
grid on;
hold off;


%%% 3. Progression Rate Comparison
% -----------------------------------------------------------------
fprintf('\n--- 3. Comparison for PROGRESSION Rates ---\n');

t_tl_cc_prog = abs(t_values_progression_tl_cc(:));
t_dr_prog = abs(t_values_progression_dr(:));

[p_val_prog, h_prog] = ranksum(t_tl_cc_prog, t_dr_prog);

% Print
fprintf('Mean absolute t-value for TL-CC Progression: %.4f\n', mean(t_tl_cc_prog));
fprintf('Mean absolute t-value for DR Progression:     %.4f\n', mean(t_dr_prog));
fprintf('Wilcoxon rank-sum test p-value:               %.4f\n', p_val_prog);

if h_prog
    fprintf('Conclusion: The difference is statistically significant. ');
    if mean(t_tl_cc_prog) > mean(t_dr_prog)
        fprintf('The "TL-CC" set shows stronger effects for Progression rates.\n');
    else
        fprintf('The "DR" set shows stronger effects for Progression rates.\n');
    end
else
    fprintf('Conclusion: The difference is not statistically significant for Progression rates.\n');
end

% Visualization
subplot(3, 1, 3);
hold on;
histogram(t_tl_cc_prog, 'DisplayName', 'TL-CC', 'Normalization', 'pdf', 'FaceAlpha', 0.6);
histogram(t_dr_prog, 'DisplayName', 'DR', 'Normalization', 'pdf', 'FaceAlpha', 0.6);
title('Progression Rates: Absolute t-value Distribution', 'FontSize', 14);
xlabel('Absolute t-value');
ylabel('Probability Density');
legend('show');
grid on;
hold off;

sgtitle('Rate-Specific Comparison of Effect Sizes (t-values)', 'FontSize', 16, 'FontWeight', 'bold');
fprintf('\n=================================================================\n');

%% Chi-sqaured test & Fisher exact test
n_total = 3 * numel(q_values_engagement_dr) - 4;
n_sig_dr = sum(q_values_engagement_dr < 0.05,'all') ...
    + sum(q_values_persistent_dr < 0.05,'all') ...
    + sum(q_values_progression_dr < 0.05,'all');
n_sig_tl_cc = sum(q_values_engagement_tl_cc < 0.05,'all') ...
    + sum(q_values_persistent_tl_cc < 0.05,'all') ...
    + sum(q_values_progression_tl_cc < 0.05,'all');
compareFeatureSets(n_sig_tl_cc,n_total,n_sig_dr,n_total);

%% ============================================================
%  Scatter plots for significant (q <= 0.05) partial correlations
%  - Draws residual-vs-residual scatter
%  - Adds linear fit, t/p values, and asterisks on each subplot
%  - Separate figures for TL-CC and DR
% ============================================================

% Preparation: Behavioral factor names (excluding Subject)
if istable(behav_factor)
    FA_names = behav_factor.Properties.VariableNames(2:end);
else
    % Fallback
    FA_names = arrayfun(@(k)sprintf('FA_%02d',k), 1:size(FA_values_filt,2), 'uni', 0);
end

% DM names (reuse if DM_names already created above, otherwise generate)
if ~exist('DM_names','var') || isempty(DM_names)
    DM_names = arrayfun(@(k)sprintf('DM%02d',k), 1:size(engagement_level_list_dr_filt,2), 'uni', 0);
end




% -------- TL-CC --------
Y_blocks_tl = {engagement_level_list_tl_cc_filt, ...
               persistent_rate_list_tl_cc_filt, ...
               progression_rate_list_tl_cc_filt};

Q_blocks_tl = {q_values_engagement_tl_cc, ...
               q_values_persistent_tl_cc, ...
               q_values_progression_tl_cc};

T_blocks_tl = {t_values_engagement_tl_cc, ...
               t_values_persistent_tl_cc, ...
               t_values_progression_tl_cc};

P_blocks_tl = {p_values_engagement_tl_cc, ...
               p_values_persistent_tl_cc, ...
               p_values_progression_tl_cc};

make_condition_figure('TL-CC-based construction', Y_blocks_tl, Q_blocks_tl, T_blocks_tl, P_blocks_tl, ...
                      X_rest_filt, FA_values_filt, DM_names, col_names);

% -------- DR --------
Y_blocks_dr = {engagement_level_list_dr_filt, ...
               persistent_rate_list_dr_filt, ...
               progression_rate_list_dr_filt};

Q_blocks_dr = {q_values_engagement_dr, ...
               q_values_persistent_dr, ...
               q_values_progression_dr};

T_blocks_dr = {t_values_engagement_dr, ...
               t_values_persistent_dr, ...
               t_values_progression_dr};

P_blocks_dr = {p_values_engagement_dr, ...
               p_values_persistent_dr, ...
               p_values_progression_dr};

make_condition_figure('DR-based construction', Y_blocks_dr, Q_blocks_dr, T_blocks_dr, P_blocks_dr, ...
                      X_rest_filt, FA_values_filt, DM_names, col_names);


% Utility for scatter plotting: DM vs Behavior(residual)
function plot_partial_scatter(Xc, x_raw, y_raw, tval, pval, qval, xlab, ylab, ttl)
    % Xc : confounds
    % x_raw : DM metric (no residual)
    % y_raw : behavior variable (residualize with confounds)

    % Compute residuals for behavior only (include intercept)
    Xc_full = Xc;
    by = Xc_full \ y_raw;
    y_res = y_raw - Xc_full*by;

    % Scatter plot: x=DM metric (raw), y=behavior residual
    scatter(x_raw, y_res, 14, 'filled'); hold on; grid on;

    % Linear fit between DM (x) and behavior residual (y)
    mdl = fitlm(x_raw, y_res, 'Intercept', true);
    xr = linspace(min(x_raw), max(x_raw), 200)';
    yr = predict(mdl, xr);
    plot(xr, yr, 'LineWidth', 2);

    % Labels/title
    xlabel(xlab);
    ylabel([ylab ' (residual)']);

    % Title with stats
    title(sprintf('%s  |  t=%.3f, q=%.3g %s', ttl, tval, qval, sigstar(qval)));

    % Zero reference lines
%     yline(0, ':'); xline(0, ':');
    hold off;
end


% Helper function for grouped processing by condition
function make_condition_figure(cond_label, Y_blocks, Q_blocks, T_blocks, P_blocks, ...
                               Xc, X_behav, DM_names, FA_names)
    % Y_blocks: {engagement, persistent, progression} matrices (Nsub x Ndm)
    % Q/T/P_blocks: corresponding q/t/p matrices
    metric_labels = {'Engagement Level','Persistence rate','Progression rate'};

    % Collect indices of significant (DM, FA, metric) pairs
    sig_list = {};   % each element: struct('metric',m,'iDM',i,'iFA',j)
    for m = 1:numel(Y_blocks)
        qmat = Q_blocks{m};
        if isempty(qmat), continue; end
        [rows, cols] = find(qmat <= 0.05);
        for k = 1:numel(rows)
            sig_list{end+1} = struct('metric', m, 'iDM', rows(k), 'iFA', cols(k)); %#ok<AGROW>
        end
    end

    if isempty(sig_list)
        figure('Color','w','Name',[cond_label ' (no q<=0.05)']);
        annotation('textbox',[.1 .45 .8 .1],'String',['No significant pairs (q > 0.05) for ' cond_label],...
                   'EdgeColor','none','HorizontalAlignment','center','FontSize',14);
        return;
    end

    % Determine subplot layout
    N = numel(sig_list);
    ncol = min(4, N);
    nrow = ceil(N / ncol);

    figure('Color','w','Name',[cond_label ' | significant partial correlations']);
    tl = tiledlayout(nrow, ncol, 'TileSpacing','compact', 'Padding','compact');
    title(tl, sprintf('%s: q < 0.05 (partial residual plots)', cond_label), ...
          'FontWeight','bold','FontSize',14);

    % Plot each significant pair
    for idx = 1:N
        S = sig_list{idx};
        m = S.metric; iDM = S.iDM; jFA = S.iFA;

        % Boundary check (especially for progression DM length)
        Ymat = Y_blocks{m};
        if iDM < 1 || iDM > size(Ymat,2), continue; end

        % Extract data
        y_vec = X_behav(:, jFA);      % behavior
        x_vec = Ymat(:, iDM);         % DM metric

        % Statistics
        tval = T_blocks{m}(iDM, jFA);
        pval = P_blocks{m}(iDM, jFA);
        qval = Q_blocks{m}(iDM, jFA);

        nexttile; %#ok<*UNRCH>
        plot_partial_scatter(Xc, x_vec, y_vec, tval, pval, qval, ...
            sprintf('%s (%s)', metric_labels{m}, DM_names{iDM}), FA_names{jFA}, ...
            sprintf('%s vs %s', DM_names{iDM}, FA_names{jFA}));
    end
end

% Function to assign significance symbols
function a = sigstar(q)
    if q < 0.001
        a = '***';
    elseif q < 0.01
        a = '**';
    elseif q < 0.05
        a = '*';
    end
end


%% Helper function for printing results
function Z = steiger_z_test(r1, r2, r12, n)
    R_sq = r1^2 + r2^2;
    f = (1 - r12) / (2 * (1 - R_sq));
    if f > 1, f = 1; end
    h = (1 - f * R_sq) / (1 - R_sq);
    Z = (r1 - r2) * sqrt((n - 3) * (1 + r12) / (2 * (1 - R_sq) * h^2));
end

function print_results(label, r_vals, t_vals, p_vals, q_vals, threshold)
    [rows, cols] = find(q_vals < threshold);
    if isempty(rows)
        fprintf('\n--- No significant results for %s ---\n', label);
    else
        fprintf('\n--- Significant Results for %s ---\n', label);
        fprintf('%-8s %-8s %-12s %-12s %-12s %-12s\n', 'DM_Idx', 'FA_Idx', 'r-value', 't-value', 'p-value', 'q-value');
        for i = 1:length(rows)
            r = rows(i);
            c = cols(i);
            fprintf('%-8d %-8d %-12.4f %-12.4f %-12.4e %-12.4e\n', r, c, r_vals(r,c), t_vals(r,c), p_vals(r,c), q_vals(r,c));
        end
    end
end

% Helper function to perform regression and extract values
function [r, t, p] = run_glm_correlation(X_confound, y_metric, x_behav)
    % Combine confounds and the behavioral variable of interest
    X_full = [X_confound, x_behav];
    
    % Use glmfit to perform ordinary least squares regression
    % The 'normal' distribution makes it equivalent to standard linear regression.
    [~, ~, stats] = glmfit(X_full, y_metric, 'normal','constant','off');
    
    % Extract results for the behavior variable (the last one in the model)
    t = stats.t(end);
    p = stats.p(end);
    
    % Calculate partial correlation coefficient 'r' from 't'
    % r = t / sqrt(t^2 + df_residual)
    df_residual = stats.dfe; % Degrees of freedom for error
    r = t / sqrt(t^2 + df_residual);
end