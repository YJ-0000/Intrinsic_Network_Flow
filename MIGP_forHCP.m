function W = MIGP_forHCP(file_path_list, t, n, is_prev)
    % MIGP_forHCP: Apply MIGP algorithm to a list of HCP data files.
    %
    % INPUT:
    %   file_path_list - Cell array containing paths to HCP data files.
    %   t              - Number of time points to consider for eigen-decomposition.
    %   n              - Number of strongest spatial eigenvectors to retain.
    %   is_prev        - Boolean flag: 
    %                    true  = use all but the last time point,
    %                    false = use all but the first time point.
    %
    % OUTPUT:
    %   W - Weighted spatial eigenvectors (v x n matrix).

    s = length(file_path_list);    % Total number of subjects/files.
    r = randperm(s);               % Randomize the order of subjects.
    
    % Load the first (randomly chosen) subject data.
    W = HCPread(file_path_list{r(1)}, is_prev);
    fprintf('.. initilized \n');

    for i = 2:s
        fprintf('.. MIGP processing (%d/%d) \n',i,s);
        % Concatenate W with the next subject data.
        W = [W; HCPread(file_path_list{r(i)}, is_prev)]; %#ok<AGROW>

        % Compute the top "temporal" eigenvectors of W.
        [U, ~] = eigs(double(W * W'), t * 2 - 1);

        % Project W onto the eigenvectors to obtain weighted spatial eigenvectors.
        W = U' * W;
    end

    % Keep only the top `n` strongest spatial eigenvectors.
    W = W(1:n, :)';
end

function data = HCPread(file_path, is_prev)
    % HCPread: Read and normalize HCP CIFTI data.
    %
    % INPUT:
    %   file_path - Path to the HCP CIFTI file.
    %   is_prev   - Boolean flag:
    %               true  = remove the last time point,
    %               false = remove the first time point.
    %
    % OUTPUT:
    %   data - Normalized time-series matrix (T x V).

    cifti = cifti_read(file_path);
    data = normalize(cifti.cdata'); % Normalize data into (T x V) matrix.
    
    % Adjust time points based on `is_prev`.
    if is_prev
        data = data(1:end-1, :);   % Exclude the last time point.
    else
        data = data(2:end, :);     % Exclude the first time point.
    end
    
    data = data - mean(data,2);
end
