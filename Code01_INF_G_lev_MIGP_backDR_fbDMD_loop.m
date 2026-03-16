clear; clc;

startup;

%% Get env variables
slurm_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));

%% Run options

method_subject_proj_list = {'DR','GroupSF','GroupTF','GroupSTF'};

max_dim = 100;
if isnan(slurm_id) || slurm_id == 1 || slurm_id == 3
    target_dim_list = 10:10:100;
elseif slurm_id == 2 || slurm_id == 4
    target_dim_list = 21:1:29;
    pause(60);
end

if isnan(slurm_id) || slurm_id == 1 || slurm_id == 2
    sample_set = 'discovery';
elseif slurm_id == 3 || slurm_id == 4
    sample_set = 'replication';
end

fitting_time_window_list  = [1];
predict_time_window_list  = [1,2,4,8];

fprintf("sample_set=%s, max_dim=%d, target_dim_list=[%s], fitting_time_window_list=[%s], predict_time_window_list=[%s]\n", ...
    sample_set, max_dim, num2str(target_dim_list), num2str(fitting_time_window_list), num2str(predict_time_window_list));

%% Load file list and subject split
% hcp_rest1_file_list: (num_all_subjects x 2) cell array of file paths
%   column 1 = REST1_LR,  column 2 = REST1_RL
data_load = load('results/split_subjects_105.mat');
file_list_data = load('secure_info/hcp_rest_file_list.mat', 'hcp_rest1_file_list','sub_ids_rest1');
hcp_rest1_file_list = file_list_data.hcp_rest1_file_list;

sub_ids = file_list_data.sub_ids_rest1;

if strcmp(sample_set, 'discovery')
    sub_ids_train = data_load.sub_ids_set1_explore;
    sub_ids_test  = data_load.sub_ids_set1_test;
elseif strcmp(sample_set, 'replication')
    sub_ids_train = data_load.sub_ids_set2_explore;
    sub_ids_test  = data_load.sub_ids_set2_test;
else
    error('Undefined sample set!!');
end

[~, idx_train] = intersect(sub_ids, sub_ids_train);
[~, idx_test]  = intersect(sub_ids, sub_ids_test);

file_list_train = hcp_rest1_file_list(idx_train, :);  % (num_train x 2)
file_list_test  = hcp_rest1_file_list(idx_test, :);    % (num_test  x 2)

num_subjects_train = size(file_list_train, 1);
num_subjects_test  = size(file_list_test, 1);

%% Print demographic info
data_load = load('secure_info/path_info.mat');
behav_data_path = data_load.behav_data_path;   % adjust if stored elsewhere
gene_data_path  = data_load.gene_data_path;    % adjust if stored elsewhere

sub_ids_sets = {sub_ids_train, sub_ids_test};
set_labels   = {'training', 'test'};
for ii = 1:2
    gene_data_table  = readtable(gene_data_path,  'VariableNamingRule', 'preserve');
    behav_data_table = readtable(behav_data_path, 'VariableNamingRule', 'preserve');

    gene_data_table  = gene_data_table(ismember(gene_data_table.Subject,  sub_ids_sets{ii}), :);
    behav_data_table = behav_data_table(ismember(behav_data_table.Subject, sub_ids_sets{ii}), :);

    gene_data_table  = sortrows(gene_data_table,  'Subject');
    behav_data_table = sortrows(behav_data_table, 'Subject');

    ages    = gene_data_table.Age_in_Yrs;
    genders = behav_data_table.Gender;

    fprintf('Total number of %s subjects: %d, Female=%d, mean age=%.2f, std=%.2f\n', ...
        set_labels{ii}, length(ages), sum(strcmp(genders,'F')), mean(ages), std(ages));
end

%% Load one file to get dimensions
data = cifti_read(file_list_train{1, 1});
data = data.cdata;
voxel_num    = size(data, 1);
temporal_dim = size(data, 2);

%% Group PCA (MIGP) — build file path list from both runs
file_path_list = reshape(file_list_train', [], 1);  % interleave LR, RL per subject

rng('default');
rng(42);
data_concat_max_dim = MIGP_forHCP(file_path_list, temporal_dim, max_dim, 'none');

data_concat_max_dim_norm   = sqrt(sum(data_concat_max_dim.^2));
data_concat_max_dim_whiten = data_concat_max_dim ./ data_concat_max_dim_norm;

%% ICA + fbDMD + Prediction (loop over target_dim)
for target_dim = target_dim_list

    data_concat_reduced = data_concat_max_dim_whiten(:, 1:target_dim);

    %% ICA — Group-level ICA
    data_concat_reduced_removeNAN = data_concat_reduced;
    mask_vec = sum(isnan(data_concat_reduced), 2) == 0;
    data_concat_reduced_removeNAN(~mask_vec, :) = [];
    [~, ~, ~, source_maps] = icatb_icaAlgorithm('infomax', data_concat_reduced_removeNAN');

    source_maps_oridim = zeros(size(data_concat_reduced))';
    source_maps_oridim(:, mask_vec) = source_maps;

    %% Back reconstruction
    inv_source = pinv(source_maps_oridim');

    data_concat_prev = zeros(target_dim, (temporal_dim-1)*2*num_subjects_train, 'single');
    data_concat_next = zeros(target_dim, (temporal_dim-1)*2*num_subjects_train, 'single');

    tcount = 1;
    for nfolder = 1:num_subjects_train
        fprintf('Back-reconstruction >> %d / %d\n', nfolder, num_subjects_train);
        for ndir = 1:2
            fprintf('  Loading %s\n', file_list_train{nfolder, ndir});
            data = cifti_read(file_list_train{nfolder, ndir});
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
    Z = A2 / num_subjects_train;
    Phi_orig_DR = source_maps' * Phi_all;

    %% Reconstruct exact mode & Spatiotemporal fingerprinting from training set
    Phi_orig_exact    = zeros(voxel_num, size(Phi_all, 2));
    D_training_noSF   = nan(num_subjects_train, 2, target_dim+1);
    D_training_SF     = nan(num_subjects_train, 2, target_dim+1);

    tcount = 1;
    for nfolder = 1:num_subjects_train
        fprintf('Fingerprinting (train) >> %d / %d\n', nfolder, num_subjects_train);
        for ndir = 1:2
            fprintf('  Loading %s\n', file_list_train{nfolder, ndir});
            data = cifti_read(file_list_train{nfolder, ndir});
            data_normalized = normalize(data.cdata')';

            % Exact mode accumulation
            temp_prod = data_normalized(:, 2:end) * temp_v(tcount:tcount+(temporal_dim-2), :);
            tcount = tcount + (temporal_dim - 1);
            Phi_orig_exact = Phi_orig_exact + temp_prod;

            % Subject-level spatial fingerprint
            time_course_sub = inv_source * data_normalized;
            temp_v1 = pinv(time_course_sub) * Phi_all;
            Phi_orig_sub = data_normalized * temp_v1;

            [Dsub, ~] = computeDMcoefficients(data_normalized, Phi_orig_sub);
            D_training_SF(nfolder, ndir, :) = Dsub;

            [Dsub, ~] = computeDMcoefficients(data_normalized, Phi_orig_DR);
            D_training_noSF(nfolder, ndir, :) = Dsub;
        end
    end

    D_group_SF   = squeeze(mean(D_training_SF, [1,2]));
    D_group_SF(1) = [];

    D_group_noSF   = squeeze(mean(D_training_noSF, [1,2]));
    D_group_noSF(1) = [];

    %% Test prediction
    for i_proj = 1:length(method_subject_proj_list)
        method_subject_proj = method_subject_proj_list{i_proj};

        R2_DM_array_list              = zeros(num_subjects_test, 2, length(fitting_time_window_list), length(predict_time_window_list));
        R2_null_array_list            = zeros(num_subjects_test, 2, length(fitting_time_window_list), length(predict_time_window_list));
        R2_DM_array_cortex_list       = zeros(num_subjects_test, 2, length(fitting_time_window_list), length(predict_time_window_list));
        R2_null_array_cortex_list     = zeros(num_subjects_test, 2, length(fitting_time_window_list), length(predict_time_window_list));
        R2_DM_array_subcortical_list  = zeros(num_subjects_test, 2, length(fitting_time_window_list), length(predict_time_window_list));
        R2_null_array_subcortical_list= zeros(num_subjects_test, 2, length(fitting_time_window_list), length(predict_time_window_list));

        for nfolder = 1:num_subjects_test
            tic
            fprintf('Test prediction [%s] >> %d / %d\n', method_subject_proj, nfolder, num_subjects_test);

            for ndir = 1:2
                fprintf('  Loading %s\n', file_list_test{nfolder, ndir});
                data = cifti_read(file_list_test{nfolder, ndir});
                data_normalized = normalize(data.cdata')';

                time_course_sub = inv_source * data_normalized(:, 1:720);

                % Subject projection
                if strcmp(method_subject_proj, 'DR') || strcmp(method_subject_proj, 'GroupTF')
                    tv = pinv(time_course_sub) * Phi_all;
                    Phi_orig_sub = data_normalized(:, 1:720) * tv;
                elseif strcmp(method_subject_proj, 'TL-cov')
                    tv = (time_course_sub(:, 1:end-1)' / Z) * Phi_all;
                    Phi_orig_sub = data_normalized(:, 2:720) * tv;
                elseif strcmp(method_subject_proj, 'GroupSF') || strcmp(method_subject_proj, 'GroupSTF')
                    Phi_orig_sub = Phi_orig_DR;
                else
                    error('Not defined projection method: %s', method_subject_proj);
                end

                % Group temporal coefficients
                if ismember(method_subject_proj, {'DR','TL-cov','GroupSF'})
                    D_group_use = [];
                elseif strcmp(method_subject_proj, 'GroupTF')
                    D_group_use = D_group_SF;
                elseif strcmp(method_subject_proj, 'GroupSTF')
                    D_group_use = D_group_noSF;
                else
                    error('Not defined projection method: %s', method_subject_proj);
                end

                [R2_DM_array, R2_DM_array_cortex, R2_DM_array_subcortical, ...
                 R2_null_array, R2_null_array_cortex, R2_null_array_subcortical] ...
                    = benchmark_using_DMs(Phi_orig_sub, data_normalized, 1, 1, 0.6, ...
                        fitting_time_window_list, predict_time_window_list, D_group_use, false);

                R2_DM_array_list(nfolder, ndir, :, :)              = R2_DM_array;
                R2_null_array_list(nfolder, ndir, :, :)            = R2_null_array;
                R2_DM_array_cortex_list(nfolder, ndir, :, :)       = R2_DM_array_cortex;
                R2_null_array_cortex_list(nfolder, ndir, :, :)     = R2_null_array_cortex;
                R2_DM_array_subcortical_list(nfolder, ndir, :, :)  = R2_DM_array_subcortical;
                R2_null_array_subcortical_list(nfolder, ndir, :, :)= R2_null_array_subcortical;
            end
            toc
        end

        % Store results per method
        switch method_subject_proj
            case 'DR'
                R2_DM_DR_array_list              = R2_DM_array_list;
                R2_DM_DR_array_cortex_list       = R2_DM_array_cortex_list;
                R2_DM_DR_array_subcortical_list  = R2_DM_array_subcortical_list;
            case 'TL-cov'
                R2_DM_TL_cov_array_list              = R2_DM_array_list;
                R2_DM_TL_cov_array_cortex_list       = R2_DM_array_cortex_list;
                R2_DM_TL_cov_array_subcortical_list  = R2_DM_array_subcortical_list;
            case 'GroupSF'
                R2_DM_GroupSF_array_list              = R2_DM_array_list;
                R2_DM_GroupSF_array_cortex_list       = R2_DM_array_cortex_list;
                R2_DM_GroupSF_array_subcortical_list  = R2_DM_array_subcortical_list;
            case 'GroupTF'
                R2_DM_GroupTF_array_list              = R2_DM_array_list;
                R2_DM_GroupTF_array_cortex_list       = R2_DM_array_cortex_list;
                R2_DM_GroupTF_array_subcortical_list  = R2_DM_array_subcortical_list;
            case 'GroupSTF'
                R2_DM_GroupSTF_array_list              = R2_DM_array_list;
                R2_DM_GroupSTF_array_cortex_list       = R2_DM_array_cortex_list;
                R2_DM_GroupSTF_array_subcortical_list  = R2_DM_array_subcortical_list;
        end
    end

    %% Save results
    dtStr = datestr(now, 'yyyymmdd_HHMMSS');

    if strcmp(sample_set, 'discovery')
        filename = sprintf('results/disc_INF_G_lev_%03d_MIGP_results_%s.mat', target_dim, dtStr);
    elseif strcmp(sample_set, 'replication')
        filename = sprintf('results/repl_INF_G_lev_%03d_MIGP_results_%s.mat', target_dim, dtStr);
    else
        error('Undefined sample set!!');
    end
    save(filename, 'lambda', 'Phi_all', 'Phi_orig_DR', 'Phi_orig_exact', 'D', 'A', 'Z', ...
        'source_maps', 'inv_source', 'R2*list', '*time_window_list', 'D_group_*');
end