function Z = tSNE_plot(X, Y, varargin)
% tSNE_plot : Run t-SNE on feature matrix X and visualize with labels Y
%
% Usage:
%   Z = tSNE_plot(X, Y)
%   Z = tSNE_plot(X, Y, 'Perplexity', 30, 'NumDimensions', 2)
%
% Inputs
%   X : [N x D] feature matrix (double/single)
%   Y : [N x 1] cell / string / categorical labels
%
% Optional Name-Value Parameters
%   'Perplexity'    : default = automatically selected
%   'NumDimensions' : default = 2 (2D or 3D supported)
%
% Output
%   Z : [N x NumDimensions] t-SNE embedding coordinates

% -----------------------------
% 1. Parse input arguments
% -----------------------------
p = inputParser;
addRequired(p, 'X', @(x)isnumeric(x));
addRequired(p, 'Y');
addParameter(p, 'Perplexity', [], @(x)isempty(x) || (isscalar(x) && x>0));
addParameter(p, 'NumDimensions', 2, @(x)isscalar(x) && (x==2 || x==3));
addParameter(p, 'FileName', 2, @(x)isempty(x) || ischar(x));
parse(p, X, Y, varargin{:});

perp   = p.Results.Perplexity;
numDim = p.Results.NumDimensions;
FileName = p.Results.FileName;

% -----------------------------
% 2. Data validation and cleanup
% -----------------------------
assert(size(X,1) == numel(Y), ...
    'Number of rows in X must match the number of labels in Y.');

X = double(X);

% Convert labels to categorical type
if iscell(Y)
    Ycat = categorical(string(Y));
elseif isstring(Y)
    Ycat = categorical(Y);
elseif iscategorical(Y)
    Ycat = Y;
else
    error('Y must be cell, string, or categorical.');
end

% Remove samples containing NaN or Inf
bad = any(~isfinite(X), 2);
if any(bad)
    warning('Removing %d samples containing NaN or Inf.', nnz(bad));
    X(bad,:) = [];
    Ycat(bad) = [];
end

% Automatically select perplexity if not specified
N = size(X,1);
if isempty(perp)
    perp = min(30, max(5, floor((N-1)/3)));
end

% -----------------------------
% 3. Run t-SNE
%    → Standardization is handled internally by tsne
% -----------------------------
rng(0); % For reproducibility
Z = tsne(X, ...
    'NumDimensions', numDim, ...
    'Perplexity', perp, ...
    'Standardize', true, ...
    'Verbose', 1);

% -----------------------------
% 4. Visualization
% -----------------------------
figure('Color','w','Position',[100,100,900,600]); hold on; grid on;

classes = categories(Ycat);
C = turbo(numel(classes));   % Color map for each class

if numDim == 2
    for i = 1:numel(classes)
        idx = (Ycat == classes{i});
        scatter(Z(idx,1), Z(idx,2), 25, ...
            'filled', ...
            'MarkerFaceColor', C(i,:), ...
            'DisplayName', classes{i});
    end
    xlabel('t-SNE 1');
    ylabel('t-SNE 2');

else   % 3D visualization
    for i = 1:numel(classes)
        idx = (Ycat == classes{i});
        scatter3(Z(idx,1), Z(idx,2), Z(idx,3), 25, ...
            'filled', ...
            'MarkerFaceColor', C(i,:), ...
            'DisplayName', classes{i});
    end
    xlabel('t-SNE 1');
    ylabel('t-SNE 2');
    zlabel('t-SNE 3');
    view(3);
end

title(sprintf('t-SNE (Perplexity = %d)', perp));
legend('Location','bestoutside');
axis tight;

if ~isempty(FileName)
    print(gcf, FileName, '-djpeg', '-r300');
end

end
