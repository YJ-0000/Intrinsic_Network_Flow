function S2mac = S2MAC(ui, vj, vk)
% S2MAC  Subspace of order 2 Modal Assurance Criterion for complex mode vectors
%
%   S2mac = S2MAC(ui, vj, vk)
%
%   Computes the correlation between a modal vector ui and the 2D subspace
%   spanned by modal vectors vj and vk, as defined in:
%
%   D'Ambrogio & Fregolent, "Higher-order MAC for the correlation of close
%   and multiple modes", Mech. Syst. Signal Process., 17(3), 599-610, 2003.
%
%   Definition (Eq. 1):
%       S2MAC_{i,jk} = max_{a,b} |ui^H (a*vj + b*vk)|^2 / (||ui||^2 * ||a*vj + b*vk||^2)
%
%   For real mode vectors, the closed-form expression (Eq. 4) is used.
%   For complex mode vectors, the maximum is found via eigenvalue decomposition
%   of the projection matrix onto the subspace spanned by vj and vk.
%
%   INPUTS:
%       ui  - (N x 1) modal vector (real or complex)
%       vj  - (N x 1) modal vector spanning the subspace (real or complex)
%       vk  - (N x 1) modal vector spanning the subspace (real or complex)
%
%   OUTPUT:
%       S2mac - scalar in [0, 1], the S2MAC correlation coefficient
%
%   NOTES:
%       - All vectors are automatically normalized to unit length.
%       - If vj == vk (or vj and vk are linearly dependent), S2MAC reduces to MAC.
%       - Works for both real and complex mode vectors.
%
%   EXAMPLES:
%       % Two orthogonal vectors spanning a subspace
%       ui = [1; 1; 0; 0];
%       vj = [1; 0; 0; 0];
%       vk = [0; 1; 0; 0];
%       s = S2MAC(ui, vj, vk)   % returns 1.0
%
%       % Reduces to standard MAC when vj == vk
%       s = S2MAC(ui, vj, vj)   % same as MAC(ui, vj)

% Input validation
if ~isequal(size(ui), size(vj)) || ~isequal(size(ui), size(vk))
    error('S2MAC:dimMismatch', 'All input vectors must have the same size.');
end

% Ensure column vectors
ui = ui(:);
vj = vj(:);
vk = vk(:);

% Normalize to unit vectors
ui = ui / norm(ui);
vj = vj / norm(vj);
vk = vk / norm(vk);

% Check if vj and vk are linearly dependent (reduces to standard MAC)
cross_val = abs(vj' * vk);
if cross_val > 1 - 1e-12
    % vj and vk are essentially the same direction -> standard MAC
    S2mac = abs(ui' * vj)^2;
    return;
end

% --- Use projection approach (works for both real and complex vectors) ---
%
% The key insight: S2MAC is the squared norm of the projection of ui onto
% the subspace spanned by {vj, vk}.
%
% Let V = [vj, vk]. The orthogonal projector onto the column space of V is:
%   P = V * (V^H * V)^{-1} * V^H
%
% Then:  S2MAC = ||P * ui||^2 / ||ui||^2
%
% Since ui is already unit-normalized:
%   S2MAC = ui^H * P * ui = ui^H * V * (V^H * V)^{-1} * V^H * ui
%
% This is equivalent to the definition (Eq. 1) because maximizing
% |ui^H * w|^2 / ||w||^2 over w in span{vj,vk} gives the squared length
% of the projection of ui onto that subspace.

V = [vj, vk];
G = V' * V;          % Gram matrix (2x2, Hermitian positive definite)
c = V' * ui;          % projection coefficients (2x1)

% S2mac = real(c' * (G \ c));  % ui^H * V * inv(V^H*V) * V^H * ui

v_a = V * pinv(V) * ui;
v_a = v_a / norm(v_a);
S2mac = abs(ui' * v_a)^2;


% Clamp to [0, 1] to handle numerical noise
S2mac = max(0, min(1, S2mac));

end