clear; clc;

current_path = pwd;

data_load = load('secure_info\path_info.mat');
HCP_denoised_path = data_load.HCP_denoised_path;
HCP_preprocessed_rest_path = data_load.HCP_preprocessed_rest_path;
behav_data_path = data_load.behav_data_path;
gene_data_path = data_load.gene_data_path;
%%
folder_denoised_list = dir([HCP_denoised_path,filesep,'*3T_rfMRI*REST1*']);
folder_denoised_list = folder_denoised_list([folder_denoised_list.isdir]);

sub_ids_foler = zeros(size(folder_denoised_list));
for n_fol = 1:length(folder_denoised_list)
    fol_denoised_name = folder_denoised_list(n_fol).name;
    aa = split(fol_denoised_name,'_');
    sub_ids_foler(n_fol) = str2double(aa{1});
end

sub_ids = unique(sub_ids_foler);
num_rest_session = zeros(size(sub_ids));
for n_sub = 1:length(sub_ids)
    num_rest_session(n_sub) = sum(sub_ids_foler==sub_ids(n_sub));
end

num_subjects = length(sub_ids);

%%
len_time = 1200;

data_folders = folder_denoised_list;

does_not_data_exist = false(length(sub_ids),1);
for nfolder = 1:num_subjects
    fprintf('>> %d / %d \n',nfolder,num_subjects);
    
    subname = split(data_folders(nfolder).name,'_');
    REST_num = subname{4};
    subname = subname{1};
    
    try
        datafile_path_LR = fullfile(data_folders(nfolder).folder, ...
                            data_folders(nfolder).name, ...
                            subname,'MNINonLinear','Results', ...
                            ['rfMRI_',REST_num,'_LR'], ...
                            ['rfMRI_',REST_num,'_LR_Atlas_hp2000_clean.dtseries.nii']);
        fprintf(['Loading data ... sub-',subname,' rfMRI LR \n']);
        data = cifti_read(datafile_path_LR);
        data = data.cdata;
        assert(size(data,2) == len_time);

        datafile_path_RL = fullfile(data_folders(nfolder).folder, ...
                            data_folders(nfolder).name, ...
                            subname,'MNINonLinear','Results', ...
                            ['rfMRI_',REST_num,'_RL'], ...
                            ['rfMRI_',REST_num,'_RL_Atlas_hp2000_clean.dtseries.nii']);
        fprintf(['Loading data ... sub-',subname,' rfMRI RL \n']);
        data = cifti_read(datafile_path_RL);
        data = data.cdata;
        assert(size(data,2) == len_time);
    catch
        does_not_data_exist(nfolder) = true;
        fprintf('*** <Data don''t exist>%s (%d / %d) ***\n',subname,nfolder,num_subjects);
    end
end

%% Identifying problematic subjects
cd(current_path);
behav_data_table = readtable(behav_data_path,'VariableNamingRule','preserve');
for nrow = size(behav_data_table,1):-1:1
    if ~any(sub_ids==behav_data_table(nrow,'Subject').Variables)
        behav_data_table(nrow,:) = [];
    end
end
behav_data_table = sortrows(behav_data_table, 'Subject');

does_have_MMSE = ~strcmp(behav_data_table.MMSE_Compl,'true');
MMSE_thres = 26;
is_cognitive_impaired = behav_data_table.MMSE_Score <= MMSE_thres;

HCP_preproc_dir = dir(fullfile(HCP_preprocessed_rest_path,'*rfMRI*'));

sub_ids_foler = zeros(size(HCP_preproc_dir));
movement_thres = 0.15;
is_excluded_due_movement = false(length(sub_ids),2);
for n_fol = 1:length(HCP_preproc_dir)
    fol_prerpoc_name = HCP_preproc_dir(n_fol).name;
    aa = split(fol_prerpoc_name,'_');
    sub_idx = find(sub_ids == str2double(aa{1}));
    
    if isempty(sub_idx)
        continue;
    end
    
    fprintf('Checking %s whether mean relative movement exceeds %0.4f mm \n',aa{1},movement_thres)
    
    REST_1_or_2 = aa{4};
    rest_idx = str2double(REST_1_or_2(end));
    try
        movement_rel_RMS_LR_path = fullfile(HCP_preprocessed_rest_path,HCP_preproc_dir(n_fol).name,aa{1},...
            'MNINonLinear','Results',['rfMRI_',REST_1_or_2,'_LR'],'Movement_RelativeRMS_mean.txt');
        movement_rel_RMS_LR = readmatrix(movement_rel_RMS_LR_path);
        movement_rel_RMS_RL_path = fullfile(HCP_preprocessed_rest_path,HCP_preproc_dir(n_fol).name,aa{1},...
            'MNINonLinear','Results',['rfMRI_',REST_1_or_2,'_RL'],'Movement_RelativeRMS_mean.txt');
        movement_rel_RMS_RL = readmatrix(movement_rel_RMS_RL_path);
        is_excluded_due_movement(sub_idx,rest_idx) = ...
            any(movement_rel_RMS_LR > movement_thres) || any(movement_rel_RMS_RL > movement_thres);
    catch
        is_excluded_due_movement(sub_idx,rest_idx) = true;
    end
end

RL_processing_errors_subs = [103010;113417;116423;120010;121719;127226;130114;143830;169040;185038;189652;202820;204218;329844;385046;401422;462139;469961;644246;688569;723141;908860;943862;969476;971160];
is_RL_processing_errors = false(length(sub_ids),1);
for n_rl = 1:length(RL_processing_errors_subs)
    is_RL_processing_errors(sub_ids==RL_processing_errors_subs(n_rl)) = true;
end

%%

sub_idx_exclude = does_not_data_exist | is_cognitive_impaired | is_excluded_due_movement(:,1) | is_RL_processing_errors;
sub_idx_include = ~sub_idx_exclude;

%%
gene_data_table = readtable(gene_data_path,'VariableNamingRule','preserve');
for nrow = size(gene_data_table,1):-1:1
    if ~any(sub_ids==gene_data_table(nrow,'Subject').Variables)
        gene_data_table(nrow,:) = [];
    end
end
gene_data_table = sortrows(gene_data_table, 'Subject');

%%
rng('default');
rng(111);

sub_ids_include = sub_ids(sub_idx_include);

sub_ids_include = sub_ids_include(randperm(length(sub_ids_include)));

sub_ids_set1_explore = sort(sub_ids_include(1:50));
sub_ids_set1_test = sort(sub_ids_include(51:100));
sub_ids_set2_explore = sort(sub_ids_include(101:150));
sub_ids_set2_test = sort(sub_ids_include(151:200));
ages = gene_data_table.Age_in_Yrs;
genders = behav_data_table.Gender;

[~,IA,~] = intersect(sub_ids,sub_ids_set1_explore);
num_female = sum(strcmp(genders(IA),'F'));
mean_age = mean(ages(IA));
std_age = std(ages(IA));
fprintf('Total number of subjects (Set1--exploration): %d, Female=%d, mean age=%0.2f, std=%0.2f \n', ...
    length(IA),num_female,mean_age,std_age);

[~,IA,~] = intersect(sub_ids,sub_ids_set1_test);
num_female = sum(strcmp(genders(IA),'F'));
mean_age = mean(ages(IA));
std_age = std(ages(IA));
fprintf('Total number of subjects (Set1--test): %d, Female=%d, mean age=%0.2f, std=%0.2f \n', ...
    length(IA),num_female,mean_age,std_age);

[~,IA,~] = intersect(sub_ids,sub_ids_set2_explore);
num_female = sum(strcmp(genders(IA),'F'));
mean_age = mean(ages(IA));
std_age = std(ages(IA));
fprintf('Total number of subjects (Set2--exploration): %d, Female=%d, mean age=%0.2f, std=%0.2f \n', ...
    length(IA),num_female,mean_age,std_age);

[~,IA,~] = intersect(sub_ids,sub_ids_set2_test);
num_female = sum(strcmp(genders(IA),'F'));
mean_age = mean(ages(IA));
std_age = std(ages(IA));
fprintf('Total number of subjects (Set2--test): %d, Female=%d, mean age=%0.2f, std=%0.2f \n', ...
    length(IA),num_female,mean_age,std_age);
%%
save results/split_subjects sub_idx_exclude sub_idx_include does_not_data_exist is_cognitive_impaired is_excluded_due_movement is_RL_processing_errors sub_ids data_folders sub_ids_set*

%%
sub_ids_include_net = sort(sub_ids_include(1:200));

writematrix(sub_ids_include_net,'results/sub_ids_include.txt');