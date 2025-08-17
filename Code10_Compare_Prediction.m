clear;clc;

current_path = pwd;
mkdir('figures');

% sample_set = 'discovery';
sample_set = 'replication';

%% load G-ICA-DMD results
target_dim_list = 10:10:100;
num_subjects = 50;
predict_time_window_list = [1,2,4,8];

mean_R2_g_ica_DR = [];
mean_R2_g_ica_TL_cov = [];
std_R2_g_ica_DR = [];
std_R2_g_ica_TL_cov = [];
mean_R2_null = [];
std_R2_null = [];

all_R2_g_ica_DR = zeros(length(target_dim_list),num_subjects,length(predict_time_window_list));
all_R2_g_ica_TL_cov = zeros(length(target_dim_list),num_subjects,length(predict_time_window_list));
all_R2_null = zeros(length(target_dim_list),num_subjects,length(predict_time_window_list));

i_dim = 0;
for target_dim = target_dim_list
    i_dim = i_dim + 1;
    
    if strcmp(sample_set,'discovery')
        g_ica_result_files = dir(['results/loop_gica_',num2str(target_dim,'%03d'),'_dmd_results_normalized_*.mat']);
    elseif strcmp(sample_set, 'replication')
        g_ica_result_files = dir(['results/repl_gica_',num2str(target_dim,'%03d'),'_dmd_results_normalized_*.mat']);
    else
        error('Undefined sample set!!');
    end
    load_results = load(fullfile(g_ica_result_files(end).folder,g_ica_result_files(end).name));
    
    R2_DM_DR_array_list = load_results.R2_DM_DR_array_list;        
    mean_R2_g_ica_DR = [mean_R2_g_ica_DR, squeeze(mean(R2_DM_DR_array_list(:,:,1,:),[1,2]))]; %#ok<AGROW>
    std_R2_g_ica_DR = [std_R2_g_ica_DR, squeeze(std(squeeze(mean(R2_DM_DR_array_list(:,:,1,:),2)),0,1))']; %#ok<AGROW>
    all_R2_g_ica_DR(i_dim,:,:) = squeeze(mean(R2_DM_DR_array_list(:,:,1,:),2));
    R2_DM_TL_cov_array_list = load_results.R2_DM_TL_cov_array_list;        
    mean_R2_g_ica_TL_cov = [mean_R2_g_ica_TL_cov, squeeze(mean(R2_DM_TL_cov_array_list(:,:,1,:),[1,2]))]; %#ok<AGROW>
    std_R2_g_ica_TL_cov = [std_R2_g_ica_TL_cov, squeeze(std(squeeze(mean(R2_DM_TL_cov_array_list(:,:,1,:),2)),0,1))']; %#ok<AGROW>
    all_R2_g_ica_TL_cov(i_dim,:,:) = squeeze(mean(R2_DM_TL_cov_array_list(:,:,1,:),2));
    
    R2_DM_null_array_list = load_results.R2_null_array_list; 
    mean_R2_null = [mean_R2_null, squeeze(mean(R2_DM_null_array_list(:,:,1,:),[1,2]))]; %#ok<AGROW>
    std_R2_null = [std_R2_null, squeeze(std(squeeze(mean(R2_DM_null_array_list(:,:,1,:),2)),0,1))']; %#ok<AGROW>
    all_R2_null(i_dim,:,:) = squeeze(mean(R2_DM_null_array_list(:,:,1,:),2));
end

%% Comparing with Null results
p_TL_CC_null = zeros(length(target_dim_list),length(predict_time_window_list));
p_DR_null = zeros(length(target_dim_list),length(predict_time_window_list));

for n_dim = 1:length(target_dim_list)
    for n_ph = 1:length(predict_time_window_list)
        % TL-CC vs Null
        R2_temp1 = squeeze(all_R2_g_ica_TL_cov(n_dim,:,n_ph));
        R2_temp2 = squeeze(all_R2_null(n_dim,:,n_ph));
        
        [~,p,~,~] = ttest(R2_temp1,R2_temp2);
        p_TL_CC_null(n_dim,n_ph) = p;
        fprintf('TL-CC vs Null | dim=%d, win=%d: p = %.12f\n', target_dim_list(n_dim), predict_time_window_list(n_ph), p);

        % DR vs Null
        R2_temp1 = squeeze(all_R2_g_ica_DR(n_dim,:,n_ph));
        R2_temp2 = squeeze(all_R2_null(n_dim,:,n_ph));
        
        [~,p,~,~] = ttest(R2_temp1,R2_temp2);
        p_DR_null(n_dim,n_ph) = p;
        fprintf('DR vs Null    | dim=%d, win=%d: p = %.12f\n', target_dim_list(n_dim), predict_time_window_list(n_ph), p);
    end
end

%% Paired t-test & heatmaps with overlaid values
% Assumes target_dim_list, predict_time_window_list, all_R2_g_ica_DR,
% and all_R2_g_ica_TL_cov are already in the workspace.

n_dims    = length(target_dim_list);
n_windows = length(predict_time_window_list);

% Preallocate t‑statistic and p‑value matrices
t_DR_TL = zeros(n_dims, n_windows);
p_DR_TL = zeros(n_dims, n_windows);

% Compute paired t-test: DR vs TL_cov
for i_dim = 1:n_dims
    for i_win = 1:n_windows
        x_DR = squeeze(all_R2_g_ica_DR   (i_dim, :, i_win));
        x_TL = squeeze(all_R2_g_ica_TL_cov(i_dim, :, i_win));
        [~, p, ~, stats] = ttest(x_DR, x_TL);
        t_DR_TL(i_dim, i_win) = stats.tstat;
        p_DR_TL(i_dim, i_win) = p;
    end
end

%%% Plot heatmaps: t-value (red–white–blue) & p-value (grayscale),  
% with X = target_dim, Y = predict_window, and values overlaid
figure('Position',[100 100 900 600]);
tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');

% Create diverging colormap: blue → white → red
ncol = 256;
half = ncol/2;
blue = [linspace(0,1,half)' linspace(0,1,half)' ones(half,1)];   % blue→white
red  = [ones(half,1)    linspace(1,0,half)' linspace(1,0,half)']; % white→red
divCmap = [blue; red];

maxT = max(abs(t_DR_TL(:)));  % symmetric color scaling limit

% --- t-value heatmap ---
ax1 = nexttile(1);
imagesc(-t_DR_TL.');            % transpose so X=dim, Y=window
ax1.Colormap = divCmap;        
caxis([-maxT, maxT]);          % center at zero
colorbar;
set(ax1, ...
    'XTick', 1:n_dims,    'XTickLabel', target_dim_list, ...
    'YTick', 1:n_windows, 'YTickLabel', predict_time_window_list);
set(ax1, 'FontSize', 12);
title('T-value: TL-CC vs. DR', 'FontSize', 22, 'FontName','Times New Roman');
xlabel('Number of ICs (Q)', 'FontSize', 18, 'FontName','Times New Roman');
ylabel('Prediction horizon (TR)', 'FontSize', 18, 'FontName','Times New Roman');

% Overlay t-values on heatmap
hold(ax1, 'on');
for i_win = 1:n_windows
    for i_dim = 1:n_dims
        val = t_DR_TL(i_dim, i_win);
        text(i_dim, i_win, sprintf('%.2f', val), ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 12, ...
            'Color', 'k');
    end
end
hold(ax1, 'off');

% --- p-value heatmap ---
ax2 = nexttile(2);
imagesc(p_DR_TL.');            % transpose to match axes
ax2.Colormap = copper;           
caxis([0,1]);                  % p in [0,1]
colorbar;
set(ax2, ...
    'XTick', 1:n_dims,    'XTickLabel', target_dim_list, ...
    'YTick', 1:n_windows, 'YTickLabel', predict_time_window_list);
set(ax2, 'FontSize', 12);
title('p-value: TL-CC vs. DR', 'FontSize', 22, 'FontName','Times New Roman');
xlabel('Number of ICs (Q)', 'FontSize', 18, 'FontName','Times New Roman');
ylabel('Prediction horizon (TR)', 'FontSize', 18, 'FontName','Times New Roman');

% Overlay p-values on heatmap
hold(ax2, 'on');
for i_win = 1:n_windows
    for i_dim = 1:n_dims
        val = p_DR_TL(i_dim, i_win);
        text(i_dim, i_win, sprintf('%.4f', val), ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 12, ...
            'Color', 'w');
    end
end
hold(ax2, 'off');

fig = gcf;
print(fig, 'figures/comparison_G_ICA_TL_vs_DR.png', '-dpng', '-r600');

%%


all_R2_sub_ica = zeros(length(target_dim_list),num_subjects,length(predict_time_window_list));
all_R2_sub_exact = zeros(length(target_dim_list),num_subjects,length(predict_time_window_list));

mean_R2_sub_ica = [];
mean_R2_sub_exact = [];
std_R2_sub_ica = [];
std_R2_sub_exact = [];

i_dim = 0;
for target_dim = 10:10:100
    i_dim = i_dim + 1;
    g_ica_result_files = dir(['results/subject_wise_ica',num2str(target_dim,'%03d'),'_dmd_results_normalized_*.mat']);
    load_results = load(fullfile(g_ica_result_files(end).folder,g_ica_result_files(end).name));
    
    R2_DM_array_list = load_results.R2_DM_array_list;        
    mean_R2_sub_ica = [mean_R2_sub_ica, squeeze(mean(R2_DM_array_list(:,:,1,:),[1,2]))]; %#ok<AGROW>
    std_R2_sub_ica = [std_R2_sub_ica, squeeze(std(squeeze(mean(R2_DM_array_list(:,:,1,:),2)),0,1))']; %#ok<AGROW>
    all_R2_sub_ica(i_dim,:,:) = squeeze(mean(R2_DM_array_list(:,:,1,:),2));
    
    if strcmp(sample_set,'discovery')
        g_ica_result_files = dir(['results/subject_wise_exact_',num2str(target_dim,'%03d'),'_dmd_results_normalized_*.mat']);
    elseif strcmp(sample_set, 'replication')
        g_ica_result_files = dir(['results/repl_subject_wise_exact_',num2str(target_dim,'%03d'),'_dmd_results_normalized_*.mat']);
    else
        error('Undefined sample set!!');
    end
    load_results = load(fullfile(g_ica_result_files(end).folder,g_ica_result_files(end).name));
    
    R2_DM_array_list = load_results.R2_DM_array_list;        
    mean_R2_sub_exact = [mean_R2_sub_exact, squeeze(mean(R2_DM_array_list(:,:,1,:),[1,2]))]; %#ok<AGROW>
    std_R2_sub_exact = [std_R2_sub_exact, squeeze(std(squeeze(mean(R2_DM_array_list(:,:,1,:),2)),0,1))']; %#ok<AGROW>
    all_R2_sub_exact(i_dim,:,:) = squeeze(mean(R2_DM_array_list(:,:,1,:),2));
end

%%
all_R2_sub_ica_filtered = all_R2_sub_ica;
% all_R2_sub_ica_filtered(all_R2_sub_ica_filtered<-1) = nan;
fprintf('.. subject-wise ICA-based DMD results %d elements were excluded. \n',sum(isnan(all_R2_sub_ica_filtered),'all'));
all_R2_sub_exact_filtered = all_R2_sub_exact;
% all_R2_sub_exact_filtered(all_R2_sub_exact_filtered<-1) = nan;
fprintf('.. subject-wise SVD-based DMD results %d elements were excluded.\n',sum(isnan(all_R2_sub_exact_filtered),'all'));
mean_R2_sub_ica_filtered = squeeze(mean(all_R2_sub_ica_filtered,2,'omitnan'));
mean_R2_sub_exact_filtered = squeeze(mean(all_R2_sub_exact_filtered,2,'omitnan'));
std_R2_sub_ica_filtered = squeeze(std(all_R2_sub_ica_filtered,0,2,'omitnan'));
CI_95_R2_sub_ica_filtered = std_R2_sub_ica_filtered * 1.96/sqrt(50);
std_R2_sub_exact_filtered = squeeze(std(all_R2_sub_exact_filtered,0,2,'omitnan'));
CI_95_R2_sub_exact_filtered = std_R2_sub_exact_filtered * 1.96/sqrt(50);

%% Comparing best and second-best methods
methods = { all_R2_g_ica_DR, all_R2_g_ica_TL_cov, all_R2_sub_exact, all_R2_sub_ica };
method_names = {'G-ICA-DMD with DR', 'G-ICA-DMD with TL-CC', 'sub-wise SVD-based DMD', 'sub-wise ICA-based DMD'};
pred_horizon_list = [1,2,4,8];

for n_ph = 1:length(pred_horizon_list)
    fprintf('%d-ahead prediction ... \n',pred_horizon_list(n_ph));
    temp_R2_list = zeros(10,4);
    
    for n_method = 1:length(methods)
        temp_all_R2 = methods{n_method};
        temp_R2_list(:,n_method) = squeeze(mean(temp_all_R2(:,:,n_ph),2));
    end
        
    for n_dm = 1:length(target_dim_list)
        [~,idx] = maxk(temp_R2_list(n_dm,:),2);
        
        if idx(1) == 1 || idx(1) == 2
            [~,idx2] = max(temp_R2_list(n_dm,[3,4]));
            idx2 = idx2 + 2;
        else
            [~,idx2] = max(temp_R2_list(n_dm,[1,2]));
        end
        
        fprintf('.. at Q=%03d -- Best: %s, 2nd-best (in other): %s \n',target_dim_list(n_dm),method_names{idx(1)},method_names{idx2});
        
        best_R2 = methods{idx(1)};
        second_best_R2 = methods{idx2};
        
        [~,p,~,st] = ttest(squeeze(best_R2(n_dm,:,n_ph)),squeeze(second_best_R2(n_dm,:,n_ph)));
        fprintf('.... t = %.4f, p = %.4f \n',st.tstat,p);
    end

end


%% Define methods and pairings (exclude DR vs TL)
pairs = [1 3; 1 4; 2 3; 2 4; 3 4];  % indices into methods

n_dims    = length(target_dim_list);
n_windows = length(predict_time_window_list);

% Precompute diverging colormap for t-values
ncol = 256; half = ncol/2;
blue = [linspace(0,1,half)' linspace(0,1,half)' ones(half,1)];   % blue→white
red  = [ones(half,1)    linspace(1,0,half)' linspace(1,0,half)']; % white→red
divCmap = [blue; red];

for ip = 1:size(pairs,1)
    i1 = pairs(ip,1);
    i2 = pairs(ip,2);
    
    % Allocate t- and p-matrices
    tmat = zeros(n_dims, n_windows);
    pmat = zeros(n_dims, n_windows);
    
    % Compute paired t-test for this pair
    for i_dim = 1:n_dims
        for i_win = 1:n_windows
            x1 = squeeze( methods{i1}(i_dim, :, i_win) );
            x2 = squeeze( methods{i2}(i_dim, :, i_win) );
            [~, p, ~, stats] = ttest(x1, x2);
            tmat(i_dim, i_win) = stats.tstat;
            pmat(i_dim, i_win) = p;
        end
    end
    
    % Determine symmetric color scale for t
    maxT = max(abs(tmat(:)));
    
    % Create new figure for this comparison
    figure('Position',[100 100 900 600]);
    tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
    
    % --- t-value heatmap ---
    ax1 = nexttile(1);
    imagesc(tmat.');                % transpose so X=dim, Y=window
    ax1.Colormap = divCmap;
    caxis([-maxT, maxT]);
    colorbar;
    set(ax1, 'FontSize', 12);
    title(sprintf('T-value: %s vs. %s', method_names{i1}, method_names{i2}), ...
          'FontSize', 22, 'FontName', 'Times New Roman');
    xlabel('Number of ICs (Q)',       'FontSize', 18, 'FontName', 'Times New Roman');
    ylabel('Prediction horizon (TR)',     'FontSize', 18, 'FontName', 'Times New Roman');
    set(ax1, ...
        'XTick', 1:n_dims,    'XTickLabel', target_dim_list, ...
        'YTick', 1:n_windows, 'YTickLabel', predict_time_window_list);
    
    hold(ax1, 'on');
    for i_win = 1:n_windows
        for i_dim = 1:n_dims
            text(i_dim, i_win, sprintf('%.2f', tmat(i_dim, i_win)), ...
                'HorizontalAlignment','center','FontSize',12,'Color','k');
        end
    end
    hold(ax1, 'off');
    
    % --- p-value heatmap ---
    ax2 = nexttile(2);
    imagesc(pmat.');                % transpose so X=dim, Y=window
    ax2.Colormap = copper;
    caxis([0,1]);
    colorbar;
    set(ax2, 'FontSize', 12);
    title(sprintf('p-value: %s vs. %s', method_names{i1}, method_names{i2}), ...
          'FontSize', 22, 'FontName', 'Times New Roman');
    xlabel('Number of ICs (Q)',       'FontSize', 18, 'FontName', 'Times New Roman');
    ylabel('Prediction horizon (TR)',     'FontSize', 18, 'FontName', 'Times New Roman');
    set(ax2, ...
        'XTick', 1:n_dims,    'XTickLabel', target_dim_list, ...
        'YTick', 1:n_windows, 'YTickLabel', predict_time_window_list);
    
    hold(ax2, 'on');
    for i_win = 1:n_windows
        for i_dim = 1:n_dims
            text(i_dim, i_win, sprintf('%.4f', pmat(i_dim, i_win)), ...
                'HorizontalAlignment','center','FontSize',12,'Color','w');
        end
    end
    hold(ax2, 'off');
    
    fig = gcf;
    print(fig, sprintf('figures/comparison_G_ICA_%s_vs_%s.png', method_names{i1}, method_names{i2}), '-dpng', '-r600');
end

%% Paired t-tests across all target-dimension pairs for each prediction window
% Using heatmap instead of imagesc
% Assumes target_dim_list, predict_time_window_list, and all_R2_g_ica_DR
% are already in the workspace.

n_dims    = length(target_dim_list);
n_windows = length(predict_time_window_list);

for i_win = 1
    % Preallocate p‑value matrix for this window
    p_mat = zeros(n_dims, n_dims);
    
    % Compute paired t‑test for every pair of target dimensions
    for i_dim = 1:n_dims
        for j_dim = 1:n_dims
            x_i = squeeze(all_R2_g_ica_DR(i_dim, :, i_win));  
            x_j = squeeze(all_R2_g_ica_DR(j_dim, :, i_win));  
            [~, p_mat(i_dim,j_dim)] = ttest(x_i, x_j);
        end
    end
    
    % Mask out the lower triangle (i > j) with NaN
    mask = tril(true(n_dims), -1);
    p_mat(mask) = NaN;
    
    % Plot p‑value heatmap for this prediction window
    figure('Position',[100 100 800 600]);
    h = heatmap( ...
        string(target_dim_list), ...   % X-axis labels
        string(target_dim_list), ...   % Y-axis labels
        p_mat, ...                     % data matrix
        'ColorLimits', [0 1], ...      % p-value range
        'CellLabelFormat','%.4f' ...   % format cell labels
    );
    
    % Customize appearance
    h.Title = sprintf('p-values (%d TR-ahead prediction)', ...
                     predict_time_window_list(i_win));
    h.XLabel = 'Number of ICs (Q)';
    h.YLabel = 'Number of ICs (Q)';
    h.FontSize = 18;                  % axis label and tick font size
    h.FontName = 'Times New Roman';
    
    fig = gcf;
    print(fig, 'figures/G_ICA_DR_Q_compare.png', '-dpng', '-r600');
end


%%
figure;
x_labels = {'predict 1s ahead','predict 2s ahead','predict 4s ahead','predict 8s ahead'};

% Define grey tones for the three series
grey_shades = [0.8 0.8 0.8;  % light grey for DR
               0.6 0.6 0.6;  % medium grey for TL-COV
               0.4 0.4 0.4; % dark grey for sub-wise
               0.2 0.2 0.2]; % dark grey for sub-wise

for i_pred = 1:4
    subplot(4,1,i_pred);
    
    % Combine the DR, TL-COV, and sub-wise data into one matrix
    data = [mean_R2_g_ica_DR(i_pred,:)', mean_R2_g_ica_TL_cov(i_pred,:)', mean_R2_sub_ica(i_pred,:)', mean_R2_sub_exact(i_pred,:)'];
    
    % Create grouped bar chart
    h = bar(10:10:100,data);
    hold on;
    
    % Apply flat face color and assign grey tones via CData
    numGroups = size(data,1);
    for s = 1:3
        h(s).FaceColor = 'flat';
        % replicate the grey shade for each bar in this series
        h(s).CData = repmat(grey_shades(s,:), numGroups, 1);
    end
    
    % Find the overall maximum value and its indices
    [~, linearIdx] = max(data(:));
    [maxGroup, maxSeries] = ind2sub(size(data), linearIdx);
    
    % Highlight the max-value bar in red
    h(maxSeries).CData(maxGroup, :) = [1 0 0];
    
    % Add title, labels, legend
    title(x_labels{i_pred});
    ylabel('Mean R^2');
    legend({'DR','TL-COV','sub-wise (ica)','sub-wise (exact)'}, 'Location', 'southeast');
    ylim([0, inf]);
    
    hold off;
end

%%
% Define original x‐positions and query points
x  = 10:10:100;     % original points at 10,20,…,100
xq = 10:1:100;      % integer points for interpolation

% Loop over each prediction horizon
for i_pred = 1:4
    % Extract the three series for this horizon
    y1 = mean_R2_g_ica_DR(i_pred, :);
    y2 = mean_R2_g_ica_TL_cov(i_pred, :);
    y3 = mean_R2_sub_exact(i_pred, :);
    
    % Perform cubic (spline) interpolation at integer points
    y1q = interp1(x, y1, xq, 'spline');
    y2q = interp1(x, y2, xq, 'spline');
    y3q = interp1(x, y3, xq, 'spline');
    
    % Find the integer x where each interpolated curve attains its maximum
    [~, idx1] = max(y1q);
    [~, idx2] = max(y2q);
    [~, idx3] = max(y3q);
    
    x_max1 = xq(idx1);
    x_max2 = xq(idx2);
    x_max3 = xq(idx3);
    
    % Print out the predictions
    fprintf('predict %ds ahead: DR max at %d, TL-CC max at %d, sub-wise SVD-based max at %d\n', ...
            i_pred, x_max1, x_max2, x_max3);
end


%%
mean_R2_g_ica_DR_fine = [];
std_R2_g_ica_DR_fine = [];
mean_R2_g_ica_TL_cov_fine = [];
target_dim_fine_list = 20:30;

all_R2_g_ica_DR_fine = zeros(length(target_dim_fine_list),num_subjects,length(predict_time_window_list));
i_dim = 0;
for target_dim = target_dim_fine_list
    i_dim = i_dim + 1;
    if strcmp(sample_set,'discovery')
        g_ica_result_files = dir(['results/loop_gica_',num2str(target_dim,'%03d'),'_dmd_results_normalized_*.mat']);
    elseif strcmp(sample_set, 'replication')
        g_ica_result_files = dir(['results/repl_gica_',num2str(target_dim,'%03d'),'_dmd_results_normalized_*.mat']);
    else
        error('Undefined sample set!!');
    end
    load_results = load(fullfile(g_ica_result_files(end).folder,g_ica_result_files(end).name));
    
    R2_DM_DR_array_list = load_results.R2_DM_DR_array_list;        
    mean_R2_g_ica_DR_fine = [mean_R2_g_ica_DR_fine, squeeze(mean(R2_DM_DR_array_list(:,:,1,:),[1,2]))]; %#ok<AGROW>
    std_R2_g_ica_DR_fine = [std_R2_g_ica_DR_fine, squeeze(std(squeeze(mean(R2_DM_DR_array_list(:,:,1,:),2)),0,1))']; %#ok<AGROW>
    all_R2_g_ica_DR_fine(i_dim,:,:) = squeeze(mean(R2_DM_DR_array_list(:,:,1,:),2));
    R2_DM_TL_cov_array_list = load_results.R2_DM_TL_cov_array_list;        
    mean_R2_g_ica_TL_cov_fine = [mean_R2_g_ica_TL_cov_fine, squeeze(mean(R2_DM_TL_cov_array_list(:,:,1,:),[1,2]))]; %#ok<AGROW>
end

%% Paired t-tests across all target-dimension pairs for each prediction window
% Using heatmap instead of imagesc
% Assumes target_dim_list, predict_time_window_list, and all_R2_g_ica_DR
% are already in the workspace.
n_dims    = length(target_dim_fine_list);
n_windows = length(predict_time_window_list);

for i_win = 1 %:n_windows
    % Preallocate p‑value matrix for this window
    p_mat = zeros(n_dims, n_dims);
    
    % Compute paired t‑test for every pair of target dimensions
    for i_dim = 1:n_dims
        for j_dim = 1:n_dims
            x_i = squeeze(all_R2_g_ica_DR_fine(i_dim, :, i_win));  
            x_j = squeeze(all_R2_g_ica_DR_fine(j_dim, :, i_win));  
            [~, p_mat(i_dim,j_dim)] = ttest(x_i, x_j);
        end
    end
    
    % Mask out the lower triangle (i > j) with NaN
    mask = tril(true(n_dims), -1);
    p_mat(mask) = NaN;
    
    % Plot p‑value heatmap for this prediction window
    figure('Position',[100 100 800 600]);
    h = heatmap( ...
        string(target_dim_fine_list), ...   % X-axis labels
        string(target_dim_fine_list), ...   % Y-axis labels
        p_mat, ...                     % data matrix
        'ColorLimits', [0 1], ...      % p-value range
        'CellLabelFormat','%.4f' ...   % format cell labels
    );
    
    % Customize appearance
    h.Title = sprintf('p-values (%d TR-ahead prediction)', ...
                     predict_time_window_list(i_win));
    h.XLabel = 'Number of ICs (Q)';
    h.YLabel = 'Number of ICs (Q)';
    h.FontSize = 18;                  % axis label and tick font size
    h.FontName = 'Times New Roman';
    
    fig = gcf;
    print(fig, 'figures/G_ICA_DR_Q_compare_fine.png', '-dpng', '-r600');
end

%%
figure;
x_labels = {'predict 1s ahead','predict 2s ahead','predict 4s ahead','predict 8s ahead'};

for i_pred = 1:4
    subplot(4,1,i_pred);
    
    % Combine the DR and TL-COV data into one matrix
    data = [mean_R2_g_ica_DR_fine(i_pred,:)', mean_R2_g_ica_TL_cov_fine(i_pred,:)'];
    
    % Create grouped bar chart
    h = bar(20:30,data);
    hold on;
    
    % Define grey tones for the two series
    grey_shades = [0.7 0.7 0.7;  % lighter grey for series 1
                   0.3 0.3 0.3]; % darker grey for series 2
    
    % Apply flat face color and assign grey tones via CData
    numGroups = size(data,1);
    for s = 1:2
        h(s).FaceColor = 'flat';
        % replicate the grey color for each bar in this series
        h(s).CData = repmat(grey_shades(s,:), numGroups, 1);
    end
    
    % Find the overall maximum value and its indices
    [~, linearIdx] = max(data(:));
    [maxGroup, maxSeries] = ind2sub(size(data), linearIdx);
    
    % Highlight the max-value bar in red
    h(maxSeries).CData(maxGroup, :) = [1 0 0];
    
    % Add title, labels, legend
    title(x_labels{i_pred});
    ylabel('Mean R^2');
    legend({'DR','TL-COV'}, 'Location', 'southeast');
    ylim([0, inf]);
    
    hold off;
end

%% load GIG-ICA results
target_dim_list = 10:10:30;
num_subjects = 50;
predict_time_window_list = [1,2,4,8];

mean_R2_gig_ica_DR = [];
mean_R2_gig_ica_TL_cov = [];
std_R2_gig_ica_DR = [];
std_R2_gig_ica_TL_cov = [];

all_R2_gig_ica_DR = zeros(length(target_dim_list),num_subjects,length(predict_time_window_list));
all_R2_gig_ica_TL_cov = zeros(length(target_dim_list),num_subjects,length(predict_time_window_list));

i_dim = 0;
for target_dim = target_dim_list
    i_dim = i_dim + 1;
    
    if strcmp(sample_set,'discovery')
        gig_ica_result_files = dir(['results/loop_gig_ica_',num2str(target_dim,'%03d'),'_dmd_results_normalized_*.mat']);
    elseif strcmp(sample_set, 'replication')
        continue
    else
        error('Undefined sample set!!');
    end
    load_results = load(fullfile(gig_ica_result_files(end).folder,gig_ica_result_files(end).name));
    
    R2_DM_DR_array_list = load_results.R2_DM_DR_array_list;        
    mean_R2_gig_ica_DR = [mean_R2_gig_ica_DR, squeeze(mean(R2_DM_DR_array_list(:,:,1,:),[1,2]))]; %#ok<AGROW>
    std_R2_gig_ica_DR = [std_R2_gig_ica_DR, squeeze(std(squeeze(mean(R2_DM_DR_array_list(:,:,1,:),2)),0,1))']; %#ok<AGROW>
    all_R2_gig_ica_DR(i_dim,:,:) = squeeze(mean(R2_DM_DR_array_list(:,:,1,:),2));
    R2_DM_TL_cov_array_list = load_results.R2_DM_TL_cov_array_list;        
    mean_R2_gig_ica_TL_cov = [mean_R2_gig_ica_TL_cov, squeeze(mean(R2_DM_TL_cov_array_list(:,:,1,:),[1,2]))]; %#ok<AGROW>
    std_R2_gig_ica_TL_cov = [std_R2_gig_ica_TL_cov, squeeze(std(squeeze(mean(R2_DM_TL_cov_array_list(:,:,1,:),2)),0,1))']; %#ok<AGROW>
    all_R2_gig_ica_TL_cov(i_dim,:,:) = squeeze(mean(R2_DM_TL_cov_array_list(:,:,1,:),2));
end

%% Statistical test (paired t-test) between standard and GIG
t_stats_stan_gig_tl_cc = zeros(length(predict_time_window_list),length(target_dim_list));
t_stats_stan_gig_dr = zeros(length(predict_time_window_list),length(target_dim_list));
p_stats_stan_gig_tl_cc = zeros(length(predict_time_window_list),length(target_dim_list));
p_stats_stan_gig_dr = zeros(length(predict_time_window_list),length(target_dim_list));

for n_ph = 1:length(predict_time_window_list)
    for n_dim = 1:length(target_dim_list)
        pred_hori = predict_time_window_list(n_ph);
        target_dim = target_dim_list(n_dim);
        
        [~,p,~,st] = ttest(squeeze(all_R2_gig_ica_TL_cov(n_dim,:,n_ph)), ...
                            squeeze(all_R2_g_ica_TL_cov(n_dim,:,n_ph)));
        t_stats_stan_gig_tl_cc(n_ph,n_dim) = st.tstat;
        p_stats_stan_gig_tl_cc(n_ph,n_dim) = p;
                        
        [~,p,~,st] = ttest(squeeze(all_R2_gig_ica_DR(n_dim,:,n_ph)), ...
                            squeeze(all_R2_g_ica_DR(n_dim,:,n_ph)));
        t_stats_stan_gig_dr(n_ph,n_dim) = st.tstat;
        p_stats_stan_gig_dr(n_ph,n_dim) = p;
    end
end
