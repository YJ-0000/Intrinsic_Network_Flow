clear; clc;

%% Settings
target_dim = 27;
TRtarget = 0.72;
flow_idx = 1;
g1_idx = 1;
g2_idx = 2;
is_test = true;
n_permutations = 1000;

%% Load group-level spatial modes
g_ica_result_files = dir('results/INF_G_lev_REST1_ALL_027_MIGP_results_*.mat');
load_results = load(fullfile(g_ica_result_files(end).folder, g_ica_result_files(end).name));

if isfield(load_results, 'Phi_orig')
    Phi_orig = -load_results.Phi_orig;
elseif isfield(load_results, 'Phi_orig_DR')
    Phi_orig = -load_results.Phi_orig_DR;
else
    Phi_orig = -load_results.Phi_orig_DL;
end

%% Define brain regions
% Each region: name, voxel indices, gradient file, cos/sin assignment for g1/g2,
%              test method ('spin' for cortex, 'moran' for subcortical),
%              model indices for Moran test [L models; R models] (if applicable)
regions = struct( ...
    'name',        {'cortex',     'cerebellum',  'thalamus',    'hippocampus', 'striatum',              'amygdala',    'brainstem'}, ...
    'target_idx',  {1:59412,      65290:83142,   88747:91282,   84561:86119,   [59413:59687,63807:65289,86677:88746], 59688:60334, 60335:63806}, ...
    'grad_file',   {'results/gradient_rsfmri_cortex', 'results/gradient_rsfmri_cerebellum', 'results/gradient_rsfmri_thalamus', ...
                    'results/gradient_rsfmri_hippocampus', 'results/gradient_rsfmri_striatum', 'results/gradient_rsfmri_amygdala', ...
                    'results/gradient_rsfmri_brainstem'}, ...
    'g1_component', {'cos',  'cos',  'cos',  'sin',  'cos',  'sin',  'cos'}, ...
    'g2_component', {'sin',  'sin',  'sin',  'cos',  'sin',  'cos',  'sin'}, ...
    'phase_sign',   {exp(1i*0), exp(1i*0), exp(1i*0), exp(1i*0), exp(1i*0), exp(1i*0), exp(1i*0)}, ...
    'test_method',  {'spin', 'moran', 'moran', 'moran', 'moran', 'moran', 'moran'}, ...
    'model_idx_L',  {[],  [10],  [20],  [14],  [3,8,18],  [5],  [7]}, ...
    'model_idx_R',  {[],  [11],  [21],  [15],  [4,9,19],  [6],  []} ...
);

% Special phase for cortex if flow_idx == 7
if flow_idx == 7
    regions(1).phase_sign = exp(1i * pi / 4);
end

%% Load atlas (needed for spin/moran tests)
if is_test
    ciftiData = cifti_read('atlas/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii');
    lh_surface_file = 'atlas/S1200.L.sphere.32k_fs_LR.surf.gii';
    rh_surface_file = 'atlas/S1200.R.sphere.32k_fs_LR.surf.gii';
    models = ciftiData.diminfo{1}.models;
    sform  = ciftiData.diminfo{1,1}.vol.sform;
end

%% Run gradient correlation for each region
for ir = 1:length(regions)
    reg = regions(ir);
    fprintf('\n=== %s ===\n', upper(reg.name));

    % Load gradients
    load_results = load(reg.grad_file);
    grad = load_results.gm.gradients{1};

    % Decompose mode into amplitude and phase
    mode_data = Phi_orig(reg.target_idx, flow_idx) * reg.phase_sign;
    r_amp = abs(mode_data);
    theta = angle(mode_data);
    data_cos = -r_amp .* cos(theta);
    data_sin = -r_amp .* sin(theta);

    % Assign components to gradients
    if strcmp(reg.g1_component, 'cos')
        g1_data = data_cos;  g2_data = data_sin;
        g1_label = 'cos';    g2_label = 'sin';
    else
        g1_data = data_sin;  g2_data = data_cos;
        g1_label = 'sin';    g2_label = 'cos';
    end

    r1 = corr(grad(:, g1_idx), g1_data);
    r2 = corr(grad(:, g2_idx), g2_data);
    fprintf('R1 (%s) = %.5f, R2 (%s) = %.5f\n', g1_label, r1, g2_label, r2);

    % Save gradient data
    writematrix([grad(:, g1_idx), grad(:, g2_idx), data_cos, data_sin], ...
        sprintf('Grad_%s.txt', reg.name), 'Delimiter', ',');

    % Permutation test
    if ~is_test, continue; end

    if strcmp(reg.test_method, 'spin')
        % Spin test for cortex
        lh_vertlist = ciftiData.diminfo{1,1}.models{1,1}.vertlist + 1;
        rh_vertlist = ciftiData.diminfo{1,1}.models{1,2}.vertlist + 1;
        lh_indices  = 1:ciftiData.diminfo{1,1}.models{1,1}.count;
        rh_start    = ciftiData.diminfo{1,1}.models{1,2}.start;
        rh_indices  = rh_start:(rh_start + ciftiData.diminfo{1,1}.models{1,2}.count - 1);

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

        cos_rotated = squeeze([y_rand{1}(lh_vertlist, 1, :); y_rand{2}(rh_vertlist, 1, :)]);
        sin_rotated = squeeze([y_rand{1}(lh_vertlist, 2, :); y_rand{2}(rh_vertlist, 2, :)]);

    else
        % Moran test for subcortical regions
        [MEM_L, MEM_R, n_L] = build_moran_mems(models, sform, reg.model_idx_L, reg.model_idx_R);

        if isempty(reg.model_idx_R)
            % Single structure (e.g., brainstem)
            cos_rotated = moran_randomization(data_cos, MEM_L, n_permutations);
            sin_rotated = moran_randomization(data_sin, MEM_L, n_permutations);
        else
            cos_L = moran_randomization(data_cos(1:n_L), MEM_L, n_permutations);
            cos_R = moran_randomization(data_cos(n_L+1:end), MEM_R, n_permutations);
            cos_rotated = [squeeze(cos_L); squeeze(cos_R)];

            sin_L = moran_randomization(data_sin(1:n_L), MEM_L, n_permutations);
            sin_R = moran_randomization(data_sin(n_L+1:end), MEM_R, n_permutations);
            sin_rotated = [squeeze(sin_L); squeeze(sin_R)];
        end
    end

    % Assign rotated data to gradients (same mapping as original)
    if strcmp(reg.g1_component, 'cos')
        r1_rand = corr(grad(:, g1_idx), cos_rotated);
        r2_rand = corr(grad(:, g2_idx), sin_rotated);
    else
        r1_rand = corr(grad(:, g1_idx), sin_rotated);
        r2_rand = corr(grad(:, g2_idx), cos_rotated);
    end

    p1 = mean(abs(r1_rand) > abs(r1));
    p2 = mean(abs(r2_rand) > abs(r2));
    fprintf('Permutation test: R1 (%s) = %.5f, p = %.5f, R2 (%s) = %.5f, p = %.5f\n', ...
        g1_label, r1, p1, g2_label, r2, p2);
end

%% ======================== Helper Functions ========================

function [MEM_L, MEM_R, n_L] = build_moran_mems(models, sform, model_idx_L, model_idx_R)
    % Build MNI coordinates and MEMs for left and right structures
    MNI_coords_L = [];
    for i = 1:length(model_idx_L)
        m = models{1, model_idx_L(i)};
        MNI_coords_L = [MNI_coords_L, sform * [m.voxlist; ones(1, m.count)]]; %#ok<AGROW>
    end
    n_L = size(MNI_coords_L, 2);
    MEM_L = coords_to_mem(MNI_coords_L);

    MEM_R = [];
    if ~isempty(model_idx_R)
        MNI_coords_R = [];
        for i = 1:length(model_idx_R)
            m = models{1, model_idx_R(i)};
            MNI_coords_R = [MNI_coords_R, sform * [m.voxlist; ones(1, m.count)]]; %#ok<AGROW>
        end
        MEM_R = coords_to_mem(MNI_coords_R);
    end
end

function MEM = coords_to_mem(MNI_coords)
    D = pdist2(MNI_coords', MNI_coords');
    W = 1 ./ D;
    W(1:size(W,1)+1:end) = 1;  % set diagonal to 1
    MEM = compute_mem(W);
end