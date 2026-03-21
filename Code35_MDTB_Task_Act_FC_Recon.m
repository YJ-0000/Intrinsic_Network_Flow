clear; clc;

mkdir('MDTB_results');
save_dir = fullfile(pwd, 'MDTB_results');

%% Settings
target_dim   = 27;
REST_num     = 'REST1';
num_subjects = 24;
N_voxel      = 91282;
N_cortex     = 59412;
N_ROIs_cortex = 360;

%% ====================================================================
%  SECTION 1: Load task activation betas and compute task-averaged maps
%  ====================================================================
result_files = dir('results/MDTB_task_act_sub*');

beta_act_vals = cell(num_subjects, 32);
task_list     = cell(num_subjects, 32);

for nsub = 1:num_subjects
    res = load(fullfile(result_files(nsub).folder, result_files(nsub).name));
    beta_act_vals(nsub, :) = res.beta_act_vals;
    task_list(nsub, :)     = res.task_list;
end

% Build unified task names (merge Movie subtypes, remove duplicates)
task_all_ref = vertcat(task_list{1, :});
task_names = sort(unique(task_all_ref));
task_names(contains(task_names, '2'))         = [];
task_names(contains(task_names, 'nBackPic'))  = [];
task_names(contains(task_names, 'Movie'))     = [];
task_names{end+1} = 'Movie';
task_names = sort(task_names);
nn = length(task_names);

% Average betas across sessions per subject x task
beta_act_mean           = nan(num_subjects, N_voxel, nn);

for nsub = 1:num_subjects
    task_all_sub = vertcat(task_list{nsub, :});
    beta_all_sub = horzcat(beta_act_vals{nsub, :});

    for n_task = 1:nn
        match_idx = contains(task_all_sub, task_names{n_task});
        temp = beta_all_sub(:, match_idx);
        beta_act_mean(nsub, :, n_task) = mean(temp, 2);

        
    end
end


%% ====================================================================
%  SECTION 2: Load rest FC and compute group-average
%  ====================================================================
result_files = dir('results/MDTB_FCs_rest_sub*');
N_ROIs = 718;

FC_rest_sub_all = zeros(num_subjects, N_ROIs, N_ROIs, 2);
nan_idx = true(num_subjects,1);
for nsub = 1:num_subjects
    res = load(fullfile(result_files(nsub).folder, result_files(nsub).name));
    for nrun = 1:2
        if ~isempty(res.FC_rest_ROIs_vals{1, nrun})
            nan_idx(nsub) = false;
            FC_rest_sub_all(nsub, :, :, nrun) = res.FC_rest_ROIs_vals{1, nrun};
        end
    end
end

FC_rest_sub_all(nan_idx,:,:,:) = [];

% Fisher-z average across runs, then across subjects
FC_rest_sub_all = tanh(mean(atanh(FC_rest_sub_all), 4));
Group_FC_rest   = squeeze(tanh(mean(atanh(FC_rest_sub_all), 1)));

%% ====================================================================
%  SECTION 3: INF-based rest FC reconstruction
%  ====================================================================

% Load parcellation
labels_cortical    = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
labels_subcortical = cifti_read('atlas/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii');

label_data = round(labels_cortical.cdata);
try
    subcort = round(labels_subcortical.cdata);
    label_data = [label_data; subcort(length(label_data)+1:end)];
catch
end

label_idx_list = sort(unique(label_data));
label_idx_list(label_idx_list == 0) = [];
N_ROIs = numel(label_idx_list);

% Load group-level ICA-DMD results
g_ica_files = dir(sprintf('results/INF_G_lev_%s_ALL_%03d_MIGP_results_*.mat', REST_num, target_dim));
load_results = load(fullfile(g_ica_files(end).folder, g_ica_files(end).name));
lambda      = load_results.lambda;
source_maps = load_results.source_maps';
inv_source  = load_results.inv_source;
Phi_orig    = load_results.Phi_orig_DR;

% Load temporal fingerprints
tf_files = dir(sprintf('results/Temporal_Fingerprints_%s_ALL_027_results_20*.mat', REST_num));
tf_result = load(fullfile(tf_files(end).folder, tf_files(end).name));
engagement_mean       = mean(tf_result.engagement_level_list)';
progression_rate_list = tf_result.progression_rate_list;

% Reconstruct rest FC
mode_select = [1:24, 26:27];
FC_Recon_rest = cal_FC_from_INF(Phi_orig, [], engagement_mean, mode_select, label_data);
FC_Recon_rest_cortex = FC_Recon_rest(1:N_ROIs_cortex, 1:N_ROIs_cortex);
Group_FC_rest_cortex = Group_FC_rest(1:N_ROIs_cortex, 1:N_ROIs_cortex);

tri_idx = find(tril(true(N_ROIs_cortex), -1));
r_rest_fc = corr(Group_FC_rest_cortex(tri_idx), FC_Recon_rest_cortex(tri_idx));
fprintf('Rest FC reconstruction (cortex): R = %.4f\n', r_rest_fc);

%% Rest FC plots
fig = figure('Position', [100,100,900,700]);
imagesc(Group_FC_rest_cortex, [-0.83, 1]); colorbar; axis square;
title('Real Rest FC (cortex)');
print(fig, fullfile(save_dir, 'FC_REST_real.jpg'), '-djpeg', '-r300');

fig = figure('Position', [100,100,900,700]);
imagesc(FC_Recon_rest_cortex, [-0.83, 1]); colorbar; axis square;
title('Reconstructed Rest FC (cortex)');
print(fig, fullfile(save_dir, 'FC_REST_recon.jpg'), '-djpeg', '-r300');

%% ====================================================================
%  SECTION 4: Load task FC and compute group-average
%  ====================================================================
result_files = dir('results/MDTB_FCs_sub*');
FCs_vals = cell(num_subjects, 32);
for nsub = 1:num_subjects
    res = load(fullfile(result_files(nsub).folder, result_files(nsub).name));
    FCs_vals(nsub, :) = res.FC_ROIs_vals;
end

FC_mean = nan(num_subjects, N_ROIs, N_ROIs, nn);
for nsub = 1:num_subjects
    FC_all_sub = cat(3, FCs_vals{nsub, :});
    task_all_sub = vertcat(task_list{nsub, :});
    for n_task = 1:nn
        match_idx = contains(task_all_sub, task_names{n_task});
        FC_mean(nsub, :, :, n_task) = tanh(mean(atanh(FC_all_sub(:,:,match_idx)), 3));
    end
end
Group_FC_task = tanh(squeeze(mean(atanh(FC_mean), 1)));

%% ====================================================================
%  SECTION 5: Task activation & FC reconstruction via INF
%  ====================================================================

Task_act_real_tval  = zeros(N_voxel, nn);
Task_act_recon_tval = zeros(N_voxel, nn);
recon_beta_all      = zeros(target_dim, nn);

act_corr    = zeros(1, nn);
R_rest_all  = zeros(1, nn);
R_restR_all = zeros(1, nn);
R_taskR_all = zeros(1, nn);

for n_task = 1:nn
    % Real t-stat map
    [~,~,~,st] = ttest(beta_act_mean(:,:,n_task));
    act_vec = st.tstat';
    Task_act_real_tval(:, n_task) = act_vec;

    mask = ~isnan(act_vec);
    act_vec_clean = act_vec; act_vec_clean(~mask) = 0;

    % Project onto INF modes
    b1 = pinv(Phi_orig(mask,:)) * act_vec_clean(mask);
    b1 = engagement_mean .* (b1 ./ abs(b1));
    recon_beta_all(:, n_task) = b1;

    % Reconstruct activation
    new_act_vec = real(Phi_orig * b1);
    Task_act_recon_tval(:, n_task) = new_act_vec;
    act_corr(n_task) = corr(act_vec_clean(mask), new_act_vec(mask));
    fprintf('(Act) %s: R = %.4f\n', task_names{n_task}, act_corr(n_task));

    % Reconstruct task FC via analytic formula
    % FC_recon_deriv = cal_FC_from_INF(Phi_orig(1:N_cortex,:), lambda, engagement_mean, ...
    %     mode_select, label_data(1:N_cortex), angle(b1));
    % FC_recon_deriv = FC_recon_deriv(1:N_ROIs_cortex, 1:N_ROIs_cortex);

    % Compare with real task FC
    temp_FC_task = Group_FC_task(1:N_ROIs_cortex, 1:N_ROIs_cortex, n_task);
    
    R_rest_all(n_task)  = corr(Group_FC_rest_cortex(tri_idx), temp_FC_task(tri_idx));
    R_restR_all(n_task) = corr(FC_Recon_rest_cortex(tri_idx), temp_FC_task(tri_idx));
    
    fprintf('(FC) %s: rest-real R=%.4f, rest-recon R=%.4f\n', ...
        task_names{n_task}, R_rest_all(n_task), R_restR_all(n_task));
end

fprintf('\n=== Summary ===\n');
fprintf('Act reconstruction:        mean R = %.4f +/- %.4f\n', mean(act_corr), std(act_corr));
fprintf('Real rest FC with taks FC: mean R = %.4f +/- %.4f\n', mean(R_rest_all), std(R_rest_all));
fprintf('Recon FC with taks FC:     mean R = %.4f +/- %.4f\n', mean(R_restR_all), std(R_restR_all));

%% PCC seed FC reconstruction (visualization)
labels = cifti_read('atlas/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii');
lh_surface_file = 'atlas/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii';
rh_surface_file = 'atlas/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii';

PCC_label = [35, 35+180];
pcc_FC = mean(FC_Recon_rest(1:N_ROIs_cortex, PCC_label), 2);
pcc_vec = zeros(N_voxel, 1);
for nr = 1:N_ROIs_cortex
    pcc_vec(label_data == label_idx_list(nr)) = pcc_FC(nr);
end

snapshot = cifti_struct_create_from_template(labels, pcc_vec, 'dscalar');
max_scale = max(abs(pcc_vec));

fig = figure('Position', [0, 0, 554, 416]);
display_cifti_cortex(fig, snapshot, lh_surface_file, rh_surface_file, [], -max_scale, max_scale, parula, 'flat');
print(fig, fullfile(save_dir, 'PCC_FC_Recon.jpg'), '-djpeg', '-r300');
close(fig);