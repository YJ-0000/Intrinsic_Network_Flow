function [p_chi2, chi2_stat, p_fisher, tbl] = compareFeatureSets(m, t1, n, t2)
%compareFeatureSets Statistically compares the association strength of two feature sets.
%
%   [p_chi2, chi2_stat, p_fisher, tbl] = compareFeatureSets(m, t1, n, t2)
%   performs a Chi-Squared test and Fisher's Exact test based on the number
%   of significant associations and total features in two sets.
%
%   INPUTS:
%       m         - Number of significant features in Set 1
%       t1        - Total number of features in Set 1
%       n         - Number of significant features in Set 2
%       t2        - Total number of features in Set 2
%
%   OUTPUTS:
%       p_chi2    - p-value from the Chi-Squared test
%       chi2_stat - The Chi-Squared statistic (χ²)
%       p_fisher  - p-value from Fisher's Exact Test
%       tbl       - The 2x2 contingency table used for the analysis
%
%   EXAMPLE USAGE:
%       [p_c, chi2, p_f, tbl] = compareFeatureSets(50, 1000, 80, 900);

% --- 1. Input Validation ---
if m > t1 || n > t2 || m < 0 || n < 0 || t1 < 0 || t2 < 0
    error('Invalid input. Significant counts cannot exceed total counts.');
end

% --- 2. Create Contingency Table ---
% Rows: Groups (Set 1, Set 2) | Columns: Outcome (Significant, Not Significant)
tbl = [m, t1 - m;  % Row 1: Counts for Set 1 (significant, not significant)
       n, t2 - n]; % Row 2: Counts for Set 2 (significant, not significant)

% --- 3. Run Fisher's Exact Test (Direct Method) ---
[~, p_fisher] = fishertest(tbl);

% --- 4. Run Chi-Squared Test (Data Reconstruction for crosstab) ---
group_var   = [ones(t1, 1); 2 * ones(t2, 1)];
outcome_var = [ones(m, 1); zeros(t1 - m, 1); ones(n, 1); zeros(t2 - n, 1)];
[~, chi2_stat, p_chi2] = crosstab(group_var, outcome_var);

% --- 5. Display Results ---
disp('-------------------------------------------');
disp('      Feature Set Association Comparison');
disp('-------------------------------------------');
disp('Contingency Table:');
fprintf('                 Significant   Not Significant\n');
fprintf('Feature Set 1:    %6d        %6d\n', tbl(1,1), tbl(1,2));
fprintf('Feature Set 2:    %6d        %6d\n', tbl(2,1), tbl(2,2));
disp('-------------------------------------------');
disp('Statistical Test Results:');
fprintf('Chi-Squared (χ²) Test: statistic = %.4f, p-value = %.4f\n', chi2_stat, p_chi2);
fprintf('Fisher''s Exact Test:                          p-value = %.4f\n', p_fisher);
disp('-------------------------------------------');

% --- Interpretation Guide ---
if p_fisher < 0.05
    fprintf('Conclusion: The difference between the feature sets is statistically significant (p < 0.05).\n');
else
    fprintf('Conclusion: The difference between the feature sets is not statistically significant (p ≥ 0.05).\n');
end
disp(' '); % Add a newline for readability

end