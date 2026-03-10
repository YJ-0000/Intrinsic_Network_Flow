function beta = REML_1stLV_spm(X,Y,RT,cutoff_highpass)

        % Get DCT basis for high-pass filter
        K_struct.RT     = RT;
        K_struct.row    = 1:size(X,1); 
        K_struct.HParam = cutoff_highpass;

        K_out = spm_filter(K_struct);

        DCT_basis = K_out.X0;

        xKXs = spm_sp('Set', spm_filter(K_out,X));
        trRV         = spm_SpUtil('trRV',xKXs);

        Y_ = Y(:,var(Y)>1e-6);
        % Y_ = Y;

        KWY      = spm_filter(K_struct,Y_);

        res = spm_sp('r',xKXs,KWY); 
        ResSS   = sum(res.^2); 

        q   = spdiags(sqrt(trRV./ResSS'),0,size(Y_,2),size(Y_,2));
        Y_ = Y_ * q;
        
        % Intialization
        YY = Y_ * Y_'/size(Y_,2);
        nscan = size(Y,1);

        Q = spm_Ce('ar',nscan,0.2);

        % design space for ReML (with confounds in filter)
        %------------------------------------------------------------------
        X_ = [X,DCT_basis];
        
        % 1) Whitening matrix W = V^(-1/2)
        [V,h,Ph,F,Fa,Fc] = spm_reml(YY,X_,Q);
        V = V*nscan/trace(V);
        W = spm_sqrtm(spm_inv(V));

        %-Design space and projector matrix [pseudoinverse] for WLS
        %--------------------------------------------------------------------------
        xKXs        = spm_sp('Set',spm_filter(K_out,W*X));    % KWX

        pKX         = spm_sp('x-',xKXs);
        
        %-Whiten/Weight data and remove filter confounds
        %======================================================================
        KWY          = spm_filter(K_out,W*Y);

        %-Weighted Least Squares estimation
        %======================================================================
        beta         = pKX*KWY;                     %-Parameter estimates
end