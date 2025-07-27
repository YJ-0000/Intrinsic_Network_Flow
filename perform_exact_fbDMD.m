function [Phi_all,lambda] = perform_exact_fbDMD(X,Y,r)
    % compute forward DMD
    [U, S1, V1] = svd(Y, 'econ'); 
    U_r = U(:, 1:r); % truncate to rank-r
    S_r = S1(1:r, 1:r);
    V_r = V1(:, 1:r);
    f_Atilde = U_r' * X * V_r / S_r;

    % compute backward DMD
    [U, S, V] = svd(X, 'econ');
    U_r = U(:, 1:r); % truncate to rank-r
    S_r = S(1:r, 1:r);
    V_r = V(:, 1:r);
    b_Atilde = U_r' * Y * V_r / S_r;

    % estimate forward/backward DMD
    Atilde = real((f_Atilde / b_Atilde) ^ 0.5);


    [Phi_all,D] = eig(Atilde);
    lambda = diag(D);
    
    Phi_all = X * (V1(:, 1:r) / S1(1:r, 1:r)) * Phi_all;
end