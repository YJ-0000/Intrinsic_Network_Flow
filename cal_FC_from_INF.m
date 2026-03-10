function FC_matrix = cal_FC_from_INF(phi, D, B_mean, target_phi_idx, parcels, phase)
% CAL_FC_FROM_INF
% Computes a functional connectivity (FC) matrix from INF modes.
%
% The FC is based on closed-form Pearson correlation between signals modeled as
% a weighted sum of cosine components:
%
%   x_i(t) = sum_k a_{i,k} cos(omega_k t + theta_{i,k})
%
% Under long-time averaging and orthogonality of different frequencies,
% Pearson correlation becomes:
%
%   rho_ij = ( sum_k a_{i,k} a_{j,k} cos(theta_{i,k} - theta_{j,k}) ) /
%            sqrt( (sum_k a_{i,k}^2) (sum_k a_{j,k}^2) )
%
% If 'phase' is provided, a short-window / derivative-sign approximation
% is used instead (correlation dominated by instantaneous slope alignment).
%
% INPUTS
%   phi             : Complex spatial modes (nodes x modes)
%   D               : Eigenvalues or dynamical coefficients per mode
%   B_mean          : Mean modal amplitudes
%   target_phi_idx  : Optional index selection of modes
%   parcels         : Optional parcel labels for ROI averaging
%   phase           : Optional phase offsets (enables derivative-based mode)
%
% OUTPUT
%   FC_matrix       : ROI-to-ROI functional connectivity matrix

    % -------------------------------------------------------------
    % Select only target modes if indices are provided
    % -------------------------------------------------------------
    if ~(nargin < 4 || isempty(target_phi_idx))
        phi = phi(:,target_phi_idx);
        if ~isempty(D)
            D = D(target_phi_idx);
        end
        B_mean = B_mean(target_phi_idx);
        if nargin >= 6 && ~isempty(phase)
            phase = phase(target_phi_idx);
        end
    end

    % -------------------------------------------------------------
    % Parcel averaging:
    % If parcels are provided, average phi within each ROI.
    % Otherwise, each node is treated as one ROI.
    % -------------------------------------------------------------
    if nargin < 5 || isempty(parcels)
        Phi_parcel = phi;
        Num_ROIs = size(phi,1);
    else
        ROI_list = sort(unique(parcels));
        ROI_list(ROI_list == 0) = [];

        Num_ROIs = length(ROI_list);
        
        Phi_parcel = zeros(Num_ROIs,size(phi,2));
        for n_roi = 1:Num_ROIs
            Phi_parcel(n_roi,:) = ...
                mean(phi(parcels==ROI_list(n_roi),:));
        end 
    end
    
    % =============================================================
    % CASE 1: No phase provided
    % Long-time averaged Pearson correlation based on cosine phase differences
    % =============================================================
    if nargin < 6 || isempty(phase)
        B_mean_temp = 2 * B_mean(1:2:end);
        Phi_parcel_temp = Phi_parcel(:,1:2:end);
        Num_INF = length(B_mean_temp);
        abs_prod_phi = zeros(Num_ROIs,Num_ROIs,Num_INF);
        cos_phi = zeros(Num_ROIs,Num_ROIs,Num_INF);
        abs_phi1_squared = zeros(Num_ROIs,Num_ROIs,Num_INF);
        abs_phi2_squared = zeros(Num_ROIs,Num_ROIs,Num_INF);
        for n_inf = 1:Num_INF
            abs_phi_single = B_mean_temp(n_inf) * abs(Phi_parcel_temp(:,n_inf));
            ang_phi_single = angle(Phi_parcel_temp(:,n_inf));
            cos_phi(:,:,n_inf) = cos(ang_phi_single - ang_phi_single');
            abs_prod_phi(:,:,n_inf) = abs_phi_single .* abs_phi_single';
            abs_phi1_squared(:,:,n_inf) = repmat(abs_phi_single.^2,[1,Num_ROIs]);
            abs_phi2_squared(:,:,n_inf) = repmat(abs_phi_single.^2,[1,Num_ROIs])';
        end
        
        numer_phi = sum(abs_prod_phi .* cos_phi,3);
        denom_phi = sqrt(sum(abs_phi1_squared,3) .* sum(abs_phi2_squared,3));

        FC_matrix = numer_phi ./ denom_phi;
    
    % =============================================================
    % CASE 2: Phase provided
    % Short-window / derivative-sign approximation of correlation
    % Correlation sign is determined by instantaneous slope alignment
    % =============================================================
    else
        B_mean_temp = 2 * B_mean(1:2:end);
        D_temp = D(1:2:end);
        Phi_parcel_temp = Phi_parcel(:,1:2:end);
        phase_temp = phase(1:2:end);
        Num_INF = length(B_mean_temp);
        
        abs_prod_phi = zeros(Num_ROIs,Num_ROIs,Num_INF);
        abs_prod_minus_sin_phi = zeros(Num_ROIs,Num_INF);
        for n_inf = 1:Num_INF
            abs_phi_single = B_mean_temp(n_inf) * angle(D_temp(n_inf)) * abs(Phi_parcel_temp(:,n_inf));
            ang_phi_single = angle(Phi_parcel_temp(:,n_inf));
            minus_sin_phi_single = - sin(ang_phi_single + phase_temp(n_inf));

            abs_prod_minus_sin_phi(:,n_inf) = ...
                abs_phi_single .* minus_sin_phi_single;

            abs_phi_single = B_mean_temp(n_inf) * abs(Phi_parcel_temp(:,n_inf));
            
            abs_prod_phi(:,:,n_inf) = abs_phi_single .* abs_phi_single';
        end
        deriv_phi = sum(abs_prod_minus_sin_phi,2);

        deriv_prod_phi = deriv_phi .* deriv_phi';
        numer_phi = sum(abs_prod_phi,3);

        FC_matrix = numer_phi .* sign(deriv_prod_phi);
    end

    

end