function B = compute_B_from_Y_PHI_D(Y, PHI, D)
    
    P = (PHI' * PHI) .* conj(D * D');
    
    q = conj(diag(D*Y'*PHI));
    
    B = P \ q;
    
end