clear; clc; close all;

mkdir('figures');

%% Settings
rest_num = 'REST1';
% rest_num = 'REST2';
target_dim = 27;
TRtarget = 0.72;

%% Load temporal fingerprints
result_files = dir(sprintf('results/Temporal_Fingerprints_%s_ALL_027_results_20*.mat', rest_num));
load_result = load(fullfile(result_files(end).folder, result_files(end).name));

sub_ids = load_result.sub_ids;
B = load_result.engagement_level_list;

%% Load eigenvalues and select oscillatory modes
g_ica_result_files = dir(sprintf('results/INF_G_lev_%s_ALL_027_MIGP_results_*.mat', rest_num));
load_results = load(fullfile(g_ica_result_files(end).folder, g_ica_result_files(end).name));
lambda = load_results.lambda;

period_DM = 2*pi*TRtarget ./ angle(lambda);

% Keep only oscillatory modes, select one per conjugate pair
oscillatory_idx = find(period_DM <= 1e6);
select_idx = oscillatory_idx(1:2:end);
num_modes = length(select_idx);

B = B(:, select_idx);

%% Load demographic data
path_info = load('secure_info/path_info.mat');
gene_data_table       = readtable(path_info.gene_data_path, 'VariableNamingRule', 'preserve');
behav_data_table      = readtable(path_info.behav_data_path, 'VariableNamingRule', 'preserve');
freesurfer_data_table = readtable(path_info.freesurfer_data_path, 'VariableNamingRule', 'preserve');

% Filter to included subjects
gene_data_table       = gene_data_table(ismember(gene_data_table.Subject, sub_ids), :);
behav_data_table      = behav_data_table(ismember(behav_data_table.Subject, sub_ids), :);
freesurfer_data_table = freesurfer_data_table(ismember(freesurfer_data_table.Subject, sub_ids), :);

% Remove subjects without zygosity info
has_zygosity = ~cellfun(@isempty, gene_data_table.ZygositySR) & ...
               ~ismissing(gene_data_table.ZygositySR);
gene_data_table       = gene_data_table(has_zygosity, :);
behav_data_table      = behav_data_table(has_zygosity, :);
freesurfer_data_table = freesurfer_data_table(has_zygosity, :);
B = B(has_zygosity, :);

gene_data_table       = sortrows(gene_data_table, 'Subject');
behav_data_table      = sortrows(behav_data_table, 'Subject');
freesurfer_data_table = sortrows(freesurfer_data_table, 'Subject');

%% Build confound matrix
age        = gene_data_table.Age_in_Yrs;
sex        = strcmp(behav_data_table.Gender, 'F');
bmi        = gene_data_table.BMI;
ICV        = freesurfer_data_table.FS_IntraCranial_Vol;
TGMV       = freesurfer_data_table.FS_Total_GM_Vol;
handedness = gene_data_table.Handedness;
edu_year   = gene_data_table.SSAGA_Educ;

race_categories = {'Am. Indian/Alaskan Nat.', ...
                   'Asian/Nat. Hawaiian/Othr Pacific Is.', ...
                   'Black or African Am.', ...
                   'More than one', ...
                   'Unknown or Not Reported', ...
                   'White'};
race_dummies = zeros(length(age), length(race_categories));
for i = 1:length(race_categories)
    race_dummies(:, i) = strcmp(gene_data_table.Race, race_categories{i});
end

X = [ones(size(age)), age, sex, age.*sex, age.^2, (age.^2).*sex, ...
     bmi, ICV.^(1/3), TGMV.^(1/3), race_dummies, handedness, edu_year];

% Remove NaN subjects
nan_mask = any(isnan(X), 2);
X(nan_mask, :)              = [];
gene_data_table(nan_mask,:) = [];
B(nan_mask, :)              = [];

%% Prepare zygosity table for APACE
new_table = gene_data_table(:, {'Subject', 'Mother_ID', 'Father_ID', 'ZygositySR'});
new_table = renamevars(new_table, 'ZygositySR', 'Zygosity');

%% Regress out confounds
Y = B - X * (pinv(X) * B);

%% Run APACE (ACE model)
writetable(new_table, 'data/zygosityT.csv');

ACEfit_Par.Model     = 'ACE';
ACEfit_Par.P_nm      = Y';
ACEfit_Par.InfMx     = fullfile(pwd, 'data', 'zygosityT.csv');
ACEfit_Par.ResDir     = fullfile('.', 'ACE_model_results_amp2');
mkdir(ACEfit_Par.ResDir);
ACEfit_Par.Subset    = [];
ACEfit_Par.Pmask     = '';
ACEfit_Par.Dsnmtx    = '';
ACEfit_Par.Nlz       = 1;
ACEfit_Par.AggNlz    = 0;
ACEfit_Par.NoImg     = 0;
ACEfit_Par.alpha_CFT = [];
ACEfit_Par.nPerm     = 5000;
ACEfit_Par.nBoot     = 1000;
ACEfit_Par.nParallel = 1;

% Run APACE pipeline
ACEfit_Par = PrepData(ACEfit_Par);
ACEfit_Par = ACEfit(ACEfit_Par);
PrepParallel(ACEfit_Par);

% Permutation inference
if ACEfit_Par.nPerm > 0
    load(fullfile(ACEfit_Par.ResDir, 'ACEfit_Par.mat'));
    ACEfit_Perm_Parallel(ACEfit_Par, 1);

    load(fullfile(ACEfit_Par.ResDir, 'ACEfit_Par.mat'));
    ACEfit_Perm_Parallel_Results(ACEfit_Par);

    load(fullfile(ACEfit_Par.ResDir, 'ACEfit_Par.mat'));
    ACEfit_Results(ACEfit_Par);
end

% Bootstrap inference
if ACEfit_Par.nBoot > 0
    load(fullfile(ACEfit_Par.ResDir, 'ACEfit_Par.mat'));
    ACEfit_Boot_Parallel(ACEfit_Par, 1);

    load(fullfile(ACEfit_Par.ResDir, 'ACEfit_Par.mat'));
    ACEfit_Boot_Parallel_Results(ACEfit_Par);

    load(fullfile(ACEfit_Par.ResDir, 'ACEfit_Par.mat'));
    Boot_CIs(ACEfit_Par);
end

% Aggregate heritability
load(fullfile(ACEfit_Par.ResDir, 'ACEfit_Par.mat'));
AgHe_Method(ACEfit_Par);

% Summary
APACEsummary(ACEfit_Par, 'ResultSummary');

%% Violin plot — cosine distance by kinship group
mkdir('figures');

% Compute cosine distances
Y_transposed = Y;  % (num_subjects x num_modes)
cosine_distances = pdist(Y_transposed, 'cosine');

% Generate all pairs
num_subs = size(Y, 1);
pairs = nchoosek(1:num_subs, 2);

zygosity1 = new_table.Zygosity(pairs(:,1));
zygosity2 = new_table.Zygosity(pairs(:,2));
mother1 = new_table.Mother_ID(pairs(:,1));
mother2 = new_table.Mother_ID(pairs(:,2));
father1 = new_table.Father_ID(pairs(:,1));
father2 = new_table.Father_ID(pairs(:,2));

% Assign kinship groups
groupLabels = repmat("Unrelated", size(pairs,1), 1);

isMZ = strcmp(zygosity1, 'MZ') & strcmp(zygosity2, 'MZ') & ...
       (mother1 == mother2) & (father1 == father2);
groupLabels(isMZ) = "MZ";

isDZ = strcmp(zygosity1, 'NotMZ') & strcmp(zygosity2, 'NotMZ') & ...
       (mother1 == mother2) & (father1 == father2);
groupLabels(isDZ) = "DZ";

isSibling = ~isMZ & ~isDZ & (mother1 == mother2 | father1 == father2);
groupLabels(isSibling) = "Sibling";

distanceTable = table(groupLabels, cosine_distances', 'VariableNames', {'Group', 'CosineDistance'});

% Plot
groupNames = {'MZ', 'DZ', 'Sibling', 'Unrelated'};
[~, group_idx] = ismember(distanceTable.Group, groupNames);

colors = [0.45, 0.80, 0.69; ...
          0.98, 0.40, 0.35; ...
          0.55, 0.60, 0.79; ...
          0.90, 0.70, 0.30];

fig = figure('Name', 'Cosine Distance by Kinship', 'NumberTitle', 'off');
daviolinplot(distanceTable.CosineDistance, 'groups', group_idx, ...
    'color', colors, 'xtlabels', groupNames, 'violin', 'full');
ylim([-1, 3]);
set(gca, 'FontName', 'Serif', 'FontSize', 24);
print(fig, 'figures/Herit_violin_amp', '-djpeg', '-r300');

%% Pairwise group comparisons
comparison_pairs = {'MZ', 'DZ'; 'DZ', 'Sibling'; 'DZ', 'Unrelated'; 'Sibling', 'Unrelated'};

for i = 1:size(comparison_pairs, 1)
    g1 = distanceTable.CosineDistance(strcmp(distanceTable.Group, comparison_pairs{i,1}));
    g2 = distanceTable.CosineDistance(strcmp(distanceTable.Group, comparison_pairs{i,2}));
    [~, p, ~, stats] = ttest2(g1, g2);
    d = (mean(g1) - mean(g2)) / stats.sd;
    fprintf('%s vs. %s: t=%.4f, Cohen''s d=%.4f, p=%.8f (one-tailed: %.8f)\n', ...
        comparison_pairs{i,1}, comparison_pairs{i,2}, stats.tstat, d, p, p/2);
end