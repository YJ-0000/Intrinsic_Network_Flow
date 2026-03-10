function [smac, theta] = SpaceMAC(V1, V2)
% SpaceMAC  Subspace correlation via principal angles between two subspaces.
%
%   [smac, theta] = SpaceMAC(V1, V2)
%
%   Compares two modal subspaces by computing the principal angles between
%   them. Each subspace is specified by a set of mode vectors (columns of
%   V1 and V2). Complex conjugate pairs are automatically converted to
%   real bases {Re, Im}, so this function naturally handles DMD modes.
%
%   The SpaceMAC value is defined as:
%       SpaceMAC = 1 - sin(theta_max)
%   where theta_max is the largest principal angle between the two
%   subspaces (Vinze, Allemang & Phillips, 2021).
%
%   SpaceMAC = 1  -> subspaces are identical
%   SpaceMAC = 0  -> subspaces are maximally separated
%
%   INPUTS:
%       V1 - (N x p) matrix whose columns span the first subspace.
%            Can be real or complex. Complex columns are expanded to
%            real bases {Re, Im} automatically.
%       V2 - (N x q) matrix whose columns span the second subspace.
%            Same convention as V1.
%
%   OUTPUTS:
%       smac  - scalar in [0, 1], the SpaceMAC correlation value.
%       theta - vector of all principal angles (radians), sorted in
%               ascending order. theta(end) is the largest angle used
%               for SpaceMAC.
%
%   EXAMPLES:
%       % 1) Compare two DMD conjugate pairs
%       phi1 = randn(10,1) + 1i*randn(10,1);
%       phi2 = randn(10,1) + 1i*randn(10,1);
%       smac = SpaceMAC([phi1, conj(phi1)], [phi2, conj(phi2)])
%
%       % 2) Identical subspace -> SpaceMAC = 1
%       phi = randn(10,1) + 1i*randn(10,1);
%       R = randn(2);  % arbitrary mixing
%       V1 = [phi, conj(phi)];
%       V2 = [phi, conj(phi)] * R;  % same span, different basis
%       smac = SpaceMAC(V1, V2)
%
%       % 3) Real modes (e.g. from symmetric FE model)
%       v1 = randn(10,1); v2 = randn(10,1);
%       v3 = randn(10,1); v4 = randn(10,1);
%       smac = SpaceMAC([v1, v2], [v3, v4])
%
%   REFERENCE:
%       Vinze, P.M., Allemang, R.J., Phillips, A.W. (2021). "Developing
%       a Correlation Criterion (SpaceMAC) for Repeated and Pseudo-repeated
%       Modes." Topics in Modal Analysis & Testing, Vol. 8, Springer.

% --- Input validation ---
if size(V1, 1) ~= size(V2, 1)
    error('SpaceMAC:dimMismatch', ...
        'V1 and V2 must have the same number of rows (DOFs).');
end

% --- Convert complex columns to real basis {Re, Im} ---
V1 = complex_to_real_basis(V1);
V2 = complex_to_real_basis(V2);

% --- Orthonormalize each subspace via QR with column pivoting ---
[Q1, R1, E1] = qr(V1, 0);
tol1 = max(size(V1)) * eps(R1(1,1));
r1 = sum(abs(diag(R1)) > tol1);
Q1 = Q1(:, 1:r1);

[Q2, R2, E2] = qr(V2, 0);
tol2 = max(size(V2)) * eps(R2(1,1));
r2 = sum(abs(diag(R2)) > tol2);
Q2 = Q2(:, 1:r2);

% --- Principal angles via SVD of Q1' * Q2 ---
sigma = svd(Q1' * Q2);

% Clamp singular values to [0, 1] for numerical safety
sigma = min(max(sigma, 0), 1);

% Principal angles (ascending order: smallest angle first)
theta = sort(acos(sigma), 'ascend');

% --- SpaceMAC = 1 - sin(largest principal angle) ---
theta_max = theta(end);
smac = 1 - sin(theta_max);

end


function Vr = complex_to_real_basis(V)
% Convert complex columns to real basis pairs {Re, Im}.
% Real columns are kept as-is. Near-zero imaginary parts are dropped.

    tol = eps(norm(V, 'fro')) * max(size(V));
    Vr = zeros(size(V, 1), 0);

    for j = 1:size(V, 2)
        col = V(:, j);
        if norm(imag(col)) > tol
            % Complex column: expand to Re and Im
            Vr = [Vr, real(col), imag(col)]; %#ok<AGROW>
        else
            % Real column: keep as-is
            Vr = [Vr, real(col)]; %#ok<AGROW>
        end
    end
end