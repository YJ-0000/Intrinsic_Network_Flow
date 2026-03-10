function beta = REML_spm(X,Y)

        xKXs = spm_sp('Set', X);
        trRV         = spm_SpUtil('trRV',xKXs);

        Y_ = Y(:,var(Y)>1e-6);

        res = spm_sp('r',xKXs,Y_); 
        ResSS   = sum(res.^2); 

        q   = spdiags(sqrt(trRV./ResSS'),0,size(Y_,2),size(Y_,2));
        Y_ = Y_ * q;

        % Intialization
        YY = Y_ * Y_' / size(Y,2);
        nscan = size(Y,1);

        Q = spm_Ce('ar',nscan,0.2);
        
        % 1) Whitening matrix W = V^(-1/2)
        [V,h,Ph,F,Fa,Fc] = spm_reml(YY,X,Q);
        V = V*nscan/trace(V);
        W = spm_sqrtm(spm_inv(V));

        %-Design space and projector matrix [pseudoinverse] for WLS
        %--------------------------------------------------------------------------
        xKXs        = spm_sp('Set',W*X);    % KWX

        pKX         = spm_sp('x-',xKXs);
        
        %-Whiten/Weight data and remove filter confounds
        %======================================================================
        KWY          = W*Y;

        %-Weighted Least Squares estimation
        %======================================================================
        beta         = pKX*KWY;                     %-Parameter estimates
end