function [D, B_mean,a0_vertex] = computeDMcoefficients(X, U, B, X_seg, Y_seg, w)
% computeDMcoefficients: Estimates evolution coefficients using (Weighted) Least Squares.
%
%   INPUT:
%     X       - (P x T) Original fMRI time-series
%     U       - (P x M) Intrinsic Network Flow (INF) spatial modes
%     B       - (M x T-1) Mode amplitudes (optional, projected VY if empty)
%     X_seg   - (P x T') Target time-points (optional,X_next)
%     Y_seg   - (P x T') Source time-points (optional,X_prev)
%     w       - (1 x T') Temporal weights (optional; e.g., task block with HRF convolution)
%
%   OUTPUT:
%     D       - ((M+1) x 1) Estimated evolution coefficients [a0; d1; d2; ...; dM]
%     B_mean  - (M x 1) Weighted mean engagement level for each mode
%     a0_vertex - (1 x 1) Scalar vertex coefficient for residual dynamics
    
    if nargin < 4 || isempty(X_seg) || isempty(Y_seg)
        X_seg = X(:,2:end);
        Y_seg = X(:,1:end-1);
    end
    if nargin < 6
        w = [];
    else
        w = reshape(w, 1, []);
        w_ = w(1:end-1);
    end
    
    V=pinv(U);
    
    num_DMs = size(U, 2);

    % initialize
    Z_temp = zeros(num_DMs+1, num_DMs+1);
    W_temp = zeros(num_DMs+1, 1);

    % Projection for mode amplitude   
    if nargin < 3 || isempty(B)
        VY      = V * Y_seg;      % (M x T')
    else
        VY      = B;
    end
    resid_Y = Y_seg - U*(V*Y_seg);       % (R x T')
    
    if isempty(w)
        % Z_temp(1,1)
        Z_temp(1,1) = sum(dot(resid_Y, resid_Y, 1));
    
        % Z_temp(1,2:end) 및 (2:end,1)
        for j = 1:num_DMs
            val = sum(dot(U(:,j)' * resid_Y, VY(j,:), 1));
            Z_temp(1, j+1) = val;
            Z_temp(j+1, 1) = val;
        end
    
        % Z_temp(2:end,2:end)
        for i = 1:num_DMs
            for j = 1:num_DMs
                Z_temp(i+1, j+1) = (U(:,i)' * U(:,j)) * sum(dot(VY(i,:), VY(j,:), 1));
            end
        end
    
        % W_temp
        W_temp(1) = sum(dot(resid_Y, X_seg, 1));
        for k = 1:num_DMs
            W_temp(k+1) = sum(dot(U(:,k) * VY(k,:), X_seg, 1));
        end
    
        % Least-square solution
        if size(U,2) < size(X_seg,1)
            D       = Z_temp \ W_temp;            
            if nargout > 2            
                resid_X = X_seg - U*(V*X_seg);  
                c1 = (dot(resid_Y, resid_Y, 2));
                b1 = (dot(resid_Y,resid_X,2));
                a0_vertex = real(b1./c1);
            else
                a0_vertex = [];
            
            end
            
            
        else
            D       = Z_temp(2:end,2:end) \ W_temp(2:end);
            D       = [0;D];
            
            a0_vertex = [];
        end
        B_mean  = mean(abs(VY), 2);
    else
        % --- System Matrix Z Construction (Weighted) ---
        % Z(1,1): Interaction of residual components weighted by w
        Z_temp(1,1) = sum(dot(resid_Y, resid_Y, 1) .* w_);
    
        % Z(1, 2:end) & Z(2:end, 1): Interactions between residuals and INF modes
        for j = 1:num_DMs
            % Correlation weighted by temporal importance
            val = sum(dot(U(:,j)' * resid_Y, VY(j,:), 1) .* w_);
            Z_temp(1, j+1) = val;
            Z_temp(j+1, 1) = val;
        end
    
        % Z(2:end, 2:end): Pairwise interactions between INF modes
        for i = 1:num_DMs
            for j = 1:num_DMs
                % Product of spatial correlation and weighted temporal correlation
                Z_temp(i+1, j+1) = (U(:,i)' * U(:,j)) * sum(dot(VY(i,:), VY(j,:), 1) .* w_);
            end
        end
    
        % --- Target Vector W Construction (Weighted) ---
        % W(1): Weighted projection of target signal onto residual space
        W_temp(1) = sum(dot(resid_Y, X_seg, 1) .* w_);
        for k = 1:num_DMs
            W_temp(k+1) = sum(dot(U(:,k) * VY(k,:), X_seg, 1) .* w_);
        end
    
        % --- Solving for Evolution Coefficients ---
        if size(U, 2) < size(X_seg, 1)
            % Standard WLS solution
            D = Z_temp \ W_temp;
            
            if nargout > 2
                % Estimate local vertex coefficient (a0) based on residuals
                resid_X = X_seg - U * (V * X_seg);  
                c1 = sum(dot(resid_Y, resid_Y, 1) .* w_);
                b1 = sum(dot(resid_Y, resid_X, 1) .* w_);
                a0_vertex = real(b1 / c1);
            else
                a0_vertex = [];
            end
        else
            % Subspace-only solution if P is small
            D = Z_temp(2:end, 2:end) \ W_temp(2:end);
            D = [0; D];
            a0_vertex = [];
        end
    
        % Calculate the weighted mean engagement level for each mode
        % This represents the characteristic "Amplitude" of the flow during the weighted period
        B_mean = sum(abs(VY) .* w_, 2) ./ sum(w_);
    end
end
