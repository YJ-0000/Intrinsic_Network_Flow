function [D, B_mean,a0_vertex] = computeDMcoefficients(X, U, B, X_seg, Y_seg)
%COMPUTESEGMENTCOEFFICIENTS 한 구간에 대해 Z·D = W 를 풀고
%                       V*Y_seg의 절댓값 평균/분산을 구함
%
%   INPUT:
%     X               - (P x T) fMRI time-series
%     U               - (R x M) DM vectors
%     B               - (P x (T-1)) mode amplitude (optional)
%   OUTPUT:
%     D       - ((M+1) x 1)   evolution coefficient
%     B_mean  - (M x 1)       Mean engagement level
    
    if nargin < 4
        X_seg = X(:,2:end);
        Y_seg = X(:,1:end-1);
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
    
end
