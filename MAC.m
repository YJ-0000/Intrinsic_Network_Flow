function mac = MAC(psi1, psi2)
%MAC Compute Modal Assurance Criterion (MAC) between two complex mode shapes.
%
%   mac = MAC(phi, psi)
%
%   Inputs:
%     psi1 : Complex mode vector (Nx1 or 1xN). Can be real or complex.
%     psi2 : Complex mode vector (Nx1 or 1xN). Can be real or complex.
%
%   Output:
%     mac : Scalar MAC value in [0, 1] (within numerical tolerance).
%
%   Definition (complex-valued):
%     MAC(phi,psi) = |psi1^H * psi2|^2 / ((psi1^H * psi1) * (psi2^H * psi2))
%   where ^H denotes conjugate transpose.

    % Ensure column vectors
    psi1 = psi1(:);
    psi2 = psi2(:);

    % Basic validation
    if length(psi1) ~= length(psi2)
        error('MAC:DimensionMismatch', 'phi and psi must have the same length.');
    end

    % Compute numerator and denominator using conjugate inner products
    num = abs(psi1' * psi2)^2;          % |phi^H psi|^2
    den = (psi1' * psi1) * (psi2' * psi2);% (phi^H phi)(psi^H psi)

    % Handle degenerate cases (zero vectors)
    if den == 0
        warning('MAC:ZeroVector', 'One or both input vectors have zero norm. Returning NaN.');
        mac = NaN;
        return;
    end

    mac = real(num / den); % Should be real; take real part for numerical robustness

    % Clamp to [0,1] to avoid tiny numerical overshoots
    mac = max(0, min(1, mac));
end
