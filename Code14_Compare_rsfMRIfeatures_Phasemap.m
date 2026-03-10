clear; clc;

%% Settings
target_dim = 27;
TRtarget = 0.72;
max_flow = 4;
n_permutations = 1000;

cortex_idx = 1:59412;
left_idx   = 1:29696;
right_idx  = 29697:59412;

% Phase correction per flow (flow 4 gets pi/4 rotation)
phase_signs = ones(1, max_flow);
phase_signs(4) = exp(1i * pi / 4);

%% Load spatial modes
g_ica_result_files = dir('results/INF_G_lev_REST1_ALL_027_MIGP_results_*.mat');
load_results = load(fullfile(g_ica_result_files(end).folder, g_ica_result_files(end).name));

if isfield(load_results, 'Phi_orig')
    Phi_orig = -load_results.Phi_orig;
elseif isfield(load_results, 'Phi_orig_DR')
    Phi_orig = -load_results.Phi_orig_DR;
else
    Phi_orig = -load_results.Phi_orig_DL;
end

%% Load rsfMRI features
load_results = load('results/seed_based_conn.mat');
conn_maps = struct( ...
    'DMN', load_results.seed_conn_PCC, ...
    'CEN', load_results.seed_conn_SMG, ...
    'SN',  load_results.seed_conn_aIns_R);

load_results = load('results/MLI.mat');
MLI = load_results.MLI;

feature_names = {'DMN', 'CEN', 'SN', 'MLI'};
feature_maps  = {conn_maps.DMN, conn_maps.CEN, conn_maps.SN, MLI};

%% Load atlas for spin test
ciftiData = cifti_read('atlas/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii');
lh_surface_file = 'atlas/S1200.L.sphere.32k_fs_LR.surf.gii';
rh_surface_file = 'atlas/S1200.R.sphere.32k_fs_LR.surf.gii';

lh_vertlist = ciftiData.diminfo{1,1}.models{1,1}.vertlist + 1;
rh_vertlist = ciftiData.diminfo{1,1}.models{1,2}.vertlist + 1;
lh_indices  = 1:ciftiData.diminfo{1,1}.models{1,1}.count;
rh_start    = ciftiData.diminfo{1,1}.models{1,2}.start;
rh_indices  = rh_start:(rh_start + ciftiData.diminfo{1,1}.models{1,2}.count - 1);

%% Section 1: Spin-test corrected correlations
fprintf('\n========== SPIN-TEST CORRECTED CORRELATIONS ==========\n');

for i = 1:max_flow
    flow_idx = 2*i - 1;
    fprintf('\n# ----- Flow #%02d ----- #\n', i);

    mode_data = Phi_orig(cortex_idx, flow_idx) * phase_signs(i);
    r_amp = abs(mode_data);
    theta = angle(mode_data);
    data_cos = -r_amp .* cos(theta);
    data_sin = -r_amp .* sin(theta);

    % Generate spin permutations
    lh_cos = zeros(ciftiData.diminfo{1,1}.models{1,1}.numvert, 1);
    lh_cos(lh_vertlist) = data_cos(lh_indices);
    rh_cos = zeros(ciftiData.diminfo{1,1}.models{1,2}.numvert, 1);
    rh_cos(rh_vertlist) = data_cos(rh_indices);

    lh_sin = zeros(ciftiData.diminfo{1,1}.models{1,1}.numvert, 1);
    lh_sin(lh_vertlist) = data_sin(lh_indices);
    rh_sin = zeros(ciftiData.diminfo{1,1}.models{1,2}.numvert, 1);
    rh_sin(rh_vertlist) = data_sin(rh_indices);

    y_rand = spin_permutations( ...
        {[lh_cos, lh_sin], [rh_cos, rh_sin]}, ...
        {lh_surface_file, rh_surface_file}, ...
        n_permutations, 'random_state', 0);

    cos_rotated = squeeze([y_rand{1}(lh_vertlist,1,:); y_rand{2}(rh_vertlist,1,:)]);
    sin_rotated = squeeze([y_rand{1}(lh_vertlist,2,:); y_rand{2}(rh_vertlist,2,:)]);

    % Test each feature
    for f = 1:length(feature_names)
        feat = feature_maps{f};
        r1 = corr(feat, data_cos);
        r2 = corr(feat, data_sin);
        p1 = mean(abs(corr(feat, cos_rotated)) > abs(r1));
        p2 = mean(abs(corr(feat, sin_rotated)) > abs(r2));
        fprintf('%s: R_cos = %.5f (p=%.5f), R_sin = %.5f (p=%.5f)\n', ...
            feature_names{f}, r1, p1, r2, p2);
    end

    writematrix([conn_maps.DMN, conn_maps.CEN, conn_maps.SN, MLI, data_cos, data_sin], ...
        sprintf('rsfmri_features_INF%d.txt', flow_idx), 'Delimiter', ',');
end

%% Section 2: Correlations without spin test (+ lateralized MLI)
fprintf('\n========== CORRELATIONS (NO SPIN TEST) ==========\n');

for i = 1:max_flow
    flow_idx = 2*i - 1;
    fprintf('\n# ----- Flow #%02d ----- #\n', i);

    mode_data = Phi_orig(cortex_idx, flow_idx) * phase_signs(i);
    r_amp = abs(mode_data);
    theta = angle(mode_data);
    data_cos = -r_amp .* cos(theta);
    data_sin = -r_amp .* sin(theta);

    for f = 1:length(feature_names)
        feat = feature_maps{f};
        fprintf('%s: R_cos = %.5f, R_sin = %.5f\n', ...
            feature_names{f}, corr(feat, data_cos), corr(feat, data_sin));
    end

    % Lateralized MLI
    fprintf('MLI (L): R_cos = %.5f, R_sin = %.5f\n', ...
        corr(MLI(left_idx), data_cos(left_idx)), corr(MLI(left_idx), data_sin(left_idx)));
    fprintf('MLI (R): R_cos = %.5f, R_sin = %.5f\n', ...
        corr(MLI(right_idx), data_cos(right_idx)), corr(MLI(right_idx), data_sin(right_idx)));
end

%% Section 3: Phase sweep — max |correlation| across theta offsets
fprintf('\n========== PHASE SWEEP (MAX CORRELATION) ==========\n');

theta_step = 0:0.1:pi;

for i = 1:max_flow
    flow_idx = 2*i - 1;
    fprintf('\n# ----- Flow #%02d ----- #\n', i);

    r_amp = abs(Phi_orig(cortex_idx, flow_idx));
    theta = angle(Phi_orig(cortex_idx, flow_idx));

    % Sweep phase offsets
    theta_all = theta + theta_step;  % (voxels x offsets)
    proj_all  = -r_amp .* cos(theta_all);

    for f = 1:length(feature_names)
        r_sweep = corr(feature_maps{f}, proj_all);
        fprintf('%s: max|R| = %.5f\n', feature_names{f}, max(abs(r_sweep)));
    end
end

%% Section 4: Temporal evolution — network correlation over oscillation cycle
fprintf('\n========== TEMPORAL EVOLUTION ==========\n');

theta_cycle = 0:0.1:4*pi;
network_names = {'DMN', 'CEN', 'SN'};
network_maps  = {conn_maps.DMN, conn_maps.CEN, conn_maps.SN};

figure('Position', [100 100 1200 800]);
for i = 1:max_flow
    flow_idx = 2*i - 1;
    target = Phi_orig(cortex_idx, flow_idx);

    % Reconstruct activation over one full oscillation cycle
    act_timecourse = 2 * real(target .* exp(1i * theta_cycle));

    subplot(2, 2, i);
    r_networks = zeros(length(theta_cycle), length(network_names));
    for n = 1:length(network_names)
        r_networks(:, n) = corr(network_maps{n}, act_timecourse)';
    end
    plot(theta_cycle, r_networks, 'LineWidth', 1.5);
    legend(network_names, 'Location', 'best');
    xlabel('Phase (rad)');
    ylabel('Correlation');
    title(sprintf('Flow #%02d', i));
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
end

print(gcf, 'figures/temporal_evolution_networks.jpg', '-djpeg', '-r300');