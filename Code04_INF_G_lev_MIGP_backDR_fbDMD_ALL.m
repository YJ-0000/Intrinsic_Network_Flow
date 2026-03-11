clear; clc;

startup;

%% Get env variables
slurm_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));

%% Run options

method_subject_proj_list = {'DR'};

max_dim = 100;
target_dim_list = 27;
if isnan(slurm_id) || slurm_id == 1
    REST_num = 'REST1';
elseif slurm_id == 2
    REST_num = 'REST2';
end

fprintf('Processing ses-%s data in HCP...\n', REST_num);

%% Load file list
file_list_data = load('secure_info/hcp_rest_file_list.mat');

if strcmp(REST_num, 'REST1')
    file_list_all = file_list_data.hcp_rest1_file_list;
    sub_ids       = file_list_data.sub_ids_rest1;
elseif strcmp(REST_num, 'REST2')
    file_list_all = file_list_data.hcp_rest2_file_list;
    sub_ids       = file_list_data.sub_ids_rest2;
else
    error('Unrecognized session name: %s', REST_num);
end

num_subjects = size(file_list_all, 1);
fprintf('Total subjects for %s: %d\n', REST_num, num_subjects);

%% Load one file to get dimensions
data = cifti_read(file_list_all{1, 1});
data = data.cdata;
voxel_num    = size(data, 1);
temporal_dim = size(data, 2);

%% Group PCA (MIGP)
file_path_list = reshape(file_list_all', [], 1);  % interleave LR, RL per subject

rng('default');
rng(42);
data_concat_max_dim = MIGP_forHCP(file_path_list, temporal_dim, max_dim, 'none');

latent = sqrt(sum(data_concat_max_dim.^2));
data_concat_max_dim_whiten = data_concat_max_dim ./ latent;

%% ICA + fbDMD (loop over target_dim)
for target_dim = target_dim_list

    data_concat_reduced = data_concat_max_dim_whiten(:, 1:target_dim);

    %% ICA — Group-level ICA
    data_concat_reduced_removeNAN = data_concat_reduced;
    mask_vec = sum(isnan(data_concat_reduced), 2) == 0;
    data_concat_reduced_removeNAN(~mask_vec, :) = [];
    [~, W, ~, source_maps] = icatb_icaAlgorithm('infomax', data_concat_reduced_removeNAN');

    source_maps_oridim = zeros(size(data_concat_reduced))';
    source_maps_oridim(:, mask_vec) = source_maps;

    %% Back reconstruction
    inv_source = pinv(source_maps_oridim');

    data_concat_prev = zeros(target_dim, (temporal_dim-1)*2*num_subjects, 'single');
    data_concat_next = zeros(target_dim, (temporal_dim-1)*2*num_subjects, 'single');

    tcount = 1;
    for nfolder = 1:num_subjects
        fprintf('Back-reconstruction >> %d / %d\n', nfolder, num_subjects);
        for ndir = 1:2
            fprintf('  Loading %s\n', file_list_all{nfolder, ndir});
            data = cifti_read(file_list_all{nfolder, ndir});
            data_normalized = normalize(data.cdata')';

            time_course_sub = (inv_source * data_normalized)';

            data_concat_prev(:, tcount:tcount+(temporal_dim-2)) = time_course_sub(1:end-1, :)';
            data_concat_next(:, tcount:tcount+(temporal_dim-2)) = time_course_sub(2:end,   :)';
            tcount = tcount + (temporal_dim - 1);
        end
    end

    X = data_concat_next; clear data_concat_next
    Y = data_concat_prev; clear data_concat_prev

    %% fbDMD
    TRtarget = 0.72;

    A1 = X*Y';  A2 = Y*Y';
    A_f = A1 * pinv(A2);
    B1 = Y*X';  B2 = X*X';
    A_b = B1 * pinv(B2);
    A = real((A_f / A_b)^0.5);

    [Phi_all, D] = eig(A);
    lambda = diag(D);

    abs_DM    = abs(lambda);
    period_DM = 2*pi*TRtarget ./ angle(lambda);

    %% DM in original space
    temp_v = pinv(Y) * Phi_all;
    Z = A2 / num_subjects;
    Phi_orig_DR = source_maps' * Phi_all;

    %% Reconstruct exact mode
    Phi_orig_exact = zeros(voxel_num, size(Phi_all, 2));

    tcount = 1;
    for nfolder = 1:num_subjects
        fprintf('Exact mode reconstruction >> %d / %d\n', nfolder, num_subjects);
        for ndir = 1:2
            fprintf('  Loading %s\n', file_list_all{nfolder, ndir});
            data = cifti_read(file_list_all{nfolder, ndir});
            data_normalized = normalize(data.cdata')';

            temp_prod = data_normalized(:, 2:end) * temp_v(tcount:tcount+(temporal_dim-2), :);
            tcount = tcount + (temporal_dim - 1);
            Phi_orig_exact = Phi_orig_exact + temp_prod;
        end
    end

    %% Save
    dtStr = datestr(now, 'yyyymmdd_HHMMSS');
    filename = sprintf('results/INF_G_lev_%s_ALL_%03d_MIGP_results_%s.mat', REST_num, target_dim, dtStr);

    save(filename, 'lambda', 'Phi_all', 'Phi_orig_DR', 'Phi_orig_exact', 'D', 'W', 'A', 'Z', ...
        'source_maps', 'inv_source', 'latent', 'data_concat_max_dim');
end