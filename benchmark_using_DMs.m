function [R2_DM_array,R2_DM_array_cortex,R2_DM_array_subcortical, R2_null_array,R2_null_array_cortex,R2_null_array_subcortical, R2_lin_array, fitting_time_window_list, predict_time_window_list] = benchmark_using_DMs(...
    DM_vectors, data, source_maps, inv_source_maps, training_ratio, ...
    fitting_time_window_list, predict_time_window_list, D_group, is_source_fit, cortex_dim)
% BENCHMARK_USING_DMS Performance benchmarking using Dynamic Mode Decomposition (DMD)
%
%
%   Inputs:
%     DM_vectors             - DMD modes
%     data                   - Observation data matrix (space × time)
%     source_maps            - Mapping from latent to observation space
%     inv_source_maps        - Mapping from observation to latent space
%     training_ratio         - Ratio of data for training (default 0.6)
%     fitting_time_window_list  - Fitting horizons (default [1, 2, 4])
%     predict_time_window_list  - Prediction horizons (default [1, 2, 4, 8])
%
%   Outputs:
%     R2_DM_array           - R² scores using DMD-based predictions
%     R2_null_array         - R² scores of null (persistence) model
%     R2_lin_array          - R² scores of linear auto-regressive model
%     fitting_time_window_list
%     predict_time_window_list

% Set default parameters if not provided
if nargin < 5
    training_ratio = 0.6;
end
if nargin < 6
    fitting_time_window_list = [1,2,4,8,16];
end
if nargin < 7
    predict_time_window_list = [1, 2, 4, 8];
end
if nargin < 8
    D_group = [];
end
if nargin < 9
    is_source_fit = false;
end
if nargin < 10
    cortex_dim = min(59412,size(data,1)-1);
end

% Split data into training and test sets
[spatial_dim, time_len] = size(data);
split_idx = round(time_len * training_ratio);
training_data = data(:, 1:split_idx);
test_data     = data(:, split_idx+1:end);


% Compute DMD eigenvalues (D) from training data
if ~is_source_fit
    if isempty(D_group)
        [D, ~] = computeDMcoefficients(training_data, DM_vectors);
    else
        B_seg_array = pinv(DM_vectors) * training_data(:,1:end-1);
        x_obs_pred = real(DM_vectors * (B_seg_array .* D_group));
        resid_training = training_data(:,1:end-1) - real(DM_vectors * B_seg_array);
        resid_pred_training = training_data(:,2:end) - x_obs_pred;
        c0 = sum(dot(resid_training, resid_training, 1));
        b0 = sum(dot(resid_training, resid_pred_training, 1));
        D = [real(b0 / c0); D_group];
    end
    
else
    % Compute linear mapping A_linear in latent space
        source_training = inv_source_maps * training_data;
        [D, ~] = computeDMcoefficients(source_training, DM_vectors);
        resid_training = training_data - source_maps * source_training;
        c0 = sum(dot(resid_training(:, 1:end-1), resid_training(:, 1:end-1), 1));
        b0 = sum(dot(resid_training(:, 1:end-1), resid_training(:, 2:end), 1));
        D(1) = real(b0 / c0);
    
end
% Compute null model coefficient a_null
c0 = sum(dot(training_data(:, 1:end-1), training_data(:, 1:end-1), 1));
b0 = sum(dot(training_data(:, 1:end-1), training_data(:, 2:end), 1));
a_null = real(b0 / c0);

% Initialize result arrays with NaN
n_fit = length(fitting_time_window_list);
n_pred = length(predict_time_window_list);
R2_DM_array   = nan(n_fit, n_pred);
R2_DM_array_cortex = nan(n_fit, n_pred);
R2_DM_array_subcortical = nan(n_fit, n_pred);
R2_null_array = nan(n_fit, n_pred);
R2_null_array_cortex = nan(n_fit, n_pred);
R2_null_array_subcortical = nan(n_fit, n_pred);
R2_lin_array = []; % For compatiblity

% Prepare for sliding-window evaluation
N_test = size(test_data, 2);
max_f = max(fitting_time_window_list);
max_p = max(predict_time_window_list);
N_t   = N_test - max_f - max_p + 1;

for i_fit = 1:n_fit
    f = fitting_time_window_list(i_fit);

    % Precompute eigenvalue powers up to fitting window f
    num_modes = numel(D) - 1;
    dd = zeros(num_modes, f);
    for k = 1:f
        dd(:, k) = D(2:end) .^ (k-1);
    end
    
    % Allocate storage for predictions and actuals
    all_actual      = zeros(spatial_dim, length(predict_time_window_list), N_t);
    all_pred        = zeros(spatial_dim, length(predict_time_window_list), N_t);
    all_pred_null   = zeros(spatial_dim, length(predict_time_window_list), N_t);

    t_count = 0;
    
    if f == 1
        if ~is_source_fit
            B_seg_array = pinv(DM_vectors) * test_data(:, max_f-f+1:(N_test - max_p));
        else
            B_seg_array = pinv(DM_vectors) * inv_source_maps * test_data(:, max_f-f+1:(N_test - max_p));
        end
    else    
        B_seg_array = zeros(size(DM_vectors,2),N_test - max_p - max_f);
        for t = max_f:(N_test - max_p)
            t_count = t_count + 1;

            % Extract recent f time steps and compute DMD mode amplitudes
            if ~is_source_fit
                B_seg = compute_B_from_Y_PHI_D(test_data(:, t-f+1:t), DM_vectors, dd);
            else
                Y_seg = inv_source_maps * test_data(:, t-f+1:t);
                B_seg = compute_B_from_Y_PHI_D(Y_seg, DM_vectors, dd);
            end
            B_seg = B_seg .* dd(:, end);  % scale by last time-step amplitudes
            B_seg_array(:,t_count) = B_seg;
        end
    end
    
    for i_pred = 1:n_pred
        p = predict_time_window_list(i_pred);
        N_seg = N_test - f - p + 1;

        if N_seg < 1
            continue;
        end

        % Predict p steps ahead in latent space
        lam_p = D(2:end) .^ p;
        x_latent_pred = real(DM_vectors * (B_seg_array .* lam_p));

        % Reconstruct observations and add residual trend
        x_obs_pred = x_latent_pred;
        residual_pred = (real(D(1)).^ p) * (test_data(:, max_f:(N_test - max_p)) - real(DM_vectors * B_seg_array));

        % Store results
        all_pred(:, i_pred, :)        = x_obs_pred + residual_pred;
        all_actual(:, i_pred, :)      = test_data(:, max_f+p:(N_test - max_p)+p);
        all_pred_null(:, i_pred, :)   = (a_null^p) * test_data(:, max_f:(N_test - max_p));
    end
    % Compute R² scores for each model
    SS_res      = mean((all_actual - all_pred)      .^2, 3);
    SS_res_null = mean((all_actual - all_pred_null) .^2, 3);
    SS_tot      = mean((all_actual - mean(all_actual,3)).^2, 3);

    R2_DM_array(i_fit,   :) = 1 - mean(SS_res      ./ SS_tot,1);
    R2_null_array(i_fit, :) = 1 - mean(SS_res_null ./ SS_tot,1);
    R2_DM_array_cortex(i_fit,   :) = 1 - mean(SS_res(1:cortex_dim,:)      ./ SS_tot(1:cortex_dim,:),1);
    R2_null_array_cortex(i_fit, :) = 1 - mean(SS_res_null(1:cortex_dim,:) ./ SS_tot(1:cortex_dim,:),1);
    R2_DM_array_subcortical(i_fit,   :) = 1 - mean(SS_res(cortex_dim+1:end,:)      ./ SS_tot(cortex_dim+1:end,:),1);
    R2_null_array_subcortical(i_fit, :) = 1 - mean(SS_res_null(cortex_dim+1:end,:) ./ SS_tot(cortex_dim+1:end,:),1);
    
end
end
