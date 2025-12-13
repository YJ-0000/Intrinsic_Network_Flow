clear; clc;

startup;

current_path = pwd;

%% Get env variables
slurm_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));

%% Run options

method_subject_proj_list = {'DR','TL-cov'}; % DR or TL-cov

max_dim = 100;
if isnan(slurm_id) || slurm_id == 1 || slurm_id == 3 
    target_dim_list = 10:10:100; 
elseif slurm_id == 2 || slurm_id == 4
    target_dim_list = 21:1:29; 
end

if isnan(slurm_id) || slurm_id == 1 || slurm_id == 2
    sample_set = 'discovery';
elseif slurm_id == 3 || slurm_id == 4
    sample_set = 'replication';
end

fitting_time_window_list =  [1];
predict_time_window_list = [1,2,4,8];

fprintf("sample_set=%s, max_dim=%d, target_dim_list=[%s], fitting_time_window_list=[%s], predict_time_window_list=[%s]\n", sample_set, max_dim, num2str(target_dim_list), num2str(fitting_time_window_list), num2str(predict_time_window_list));

%%

path_info = load('secure_info/path_info.mat');
HCP_data_path = path_info.HCP_denoised_path;
data_load = load('results/split_subjects_105.mat');

data_folders = dir(HCP_data_path);
isFolder = [data_folders.isdir];
names = {data_folders.name};
mask = isFolder & ~ismember(names, {'.','..'});
data_folders = data_folders(mask);

sub_ids = data_load.sub_ids;
if strcmp(sample_set,'discovery')
    sub_ids_set_explore = data_load.sub_ids_set1_explore;
elseif strcmp(sample_set, 'replication')
    sub_ids_set_explore = data_load.sub_ids_set2_explore;
else
    error('Undefined sample set!!');
end
[~,IA,~] = intersect(sub_ids,sub_ids_set_explore);
data_folders_train = data_folders(IA);

sub_ids = data_load.sub_ids;
if strcmp(sample_set,'discovery')
    sub_ids_set_test = data_load.sub_ids_set1_test;
elseif strcmp(sample_set, 'replication')
    sub_ids_set_test = data_load.sub_ids_set2_test;
else
    error('Undefined sample set!!');
end
[~,IA,~] = intersect(sub_ids,sub_ids_set_test);
data_folders_test = data_folders(IA);

%%
behav_data_path = path_info.behav_data_path;
gene_data_path = path_info.gene_data_path;
sub_ids_sets = {sub_ids_set_explore,sub_ids_set_test};
for ii = 1:2
    gene_data_table = readtable(gene_data_path,'VariableNamingRule','preserve');
    behav_data_table = readtable(behav_data_path,'VariableNamingRule','preserve');
    for nrow = size(gene_data_table,1):-1:1
        if ~any(sub_ids_sets{ii}==gene_data_table(nrow,'Subject').Variables)
            gene_data_table(nrow,:) = [];
        end
    end
    for nrow = size(behav_data_table,1):-1:1
        if ~any(sub_ids_sets{ii}==behav_data_table(nrow,'Subject').Variables)
            behav_data_table(nrow,:) = [];
        end
    end
    gene_data_table = sortrows(gene_data_table, 'Subject');
    behav_data_table = sortrows(behav_data_table, 'Subject');
    
    ages = gene_data_table.Age_in_Yrs;
    genders = behav_data_table.Gender;
    
    num_subjects = length(ages);
    num_female = sum(strcmp(genders,'F'));
    mean_age = mean(ages);
    std_age = std(ages);
    if ii == 1
        fprintf('Total number in discovery sample subjects: %d, Female=%d, mean age=%0.2f, std=%0.2f \n', ...
            num_subjects,num_female,mean_age,std_age);
    elseif ii == 2
        fprintf('Total number of test samplesubjects: %d, Female=%d, mean age=%0.2f, std=%0.2f \n', ...
            num_subjects,num_female,mean_age,std_age);
    end
end

%% Loading one file

num_data_used = 1;

comp_est_accum = zeros(num_data_used*2,1);

nfolder = 1;

subname = split(data_folders_train(nfolder).name,'_');
REST_num = 'REST1';
subname = subname{1};

datafile_path_LR = fullfile(data_folders_train(nfolder).folder, subname,...
                    ['s6_rfMRI_',REST_num,'_LR_Atlas_hp2000_clean.dtseries.nii']);
fprintf(['Loading data ... sub-',subname,' rfMRI LR \n']);
data = cifti_read(datafile_path_LR);
data = data.cdata;
voxel_num = size(data,1);
temporal_dim = size(data,2);

%% Group PCA

num_subjects = 50;

file_path_list = cell(2*num_subjects,1);
tcount = 1;
for nfolder = 1:num_subjects
    fprintf('>> %d / %d \n',nfolder,num_subjects);

    subname = split(data_folders_train(nfolder).name,'_');
    subname = subname{1};

    datafile_path_LR = fullfile(HCP_data_path, ...
                        subname,...
                        ['s6_rfMRI_',REST_num,'_LR_Atlas_hp2000_clean.dtseries.nii']);
    
    datafile_path_RL = fullfile(data_folders_train(nfolder).folder, ...
                        subname,...
                        ['s6_rfMRI_',REST_num,'_RL_Atlas_hp2000_clean.dtseries.nii']);
    

    file_path_list{tcount} = datafile_path_LR;
    tcount = tcount + 1;
    file_path_list{tcount} = datafile_path_RL;
    tcount = tcount + 1;
end

rng('default');
rng(42);
data_concat_max_dim = MIGP_forHCP(file_path_list, temporal_dim, max_dim, 'none');

data_concat_max_dim_norm = sqrt(sum(data_concat_max_dim.^2));

data_concat_max_dim_whiten = data_concat_max_dim ./ data_concat_max_dim_norm;

%% ICA - data concatanation
for target_dim = target_dim_list
    

    dim_reduced_1 = round(1.5*target_dim);
    dim_reduced_2 = target_dim;
    
    data_concat_reduced = data_concat_max_dim_whiten(:,1:target_dim);
    
    
    %% ICA - Group-level ICA
    data_concat_reduced_removeNAN = data_concat_reduced;
    mask_vec = sum(isnan(data_concat_reduced),2)==0;
    data_concat_reduced_removeNAN(sum(isnan(data_concat_reduced),2)>0,:) = [];
    [icaAlgo, W, W_inv, source_maps] = icatb_icaAlgorithm('infomax', data_concat_reduced_removeNAN');
    
    source_maps_oridim = zeros(size(data_concat_reduced))';
    source_maps_oridim(:,mask_vec) = source_maps;
    

    %% Back reconstruction

    inv_source = pinv(source_maps_oridim');

    data_concat_prev = zeros(dim_reduced_2,(temporal_dim-1)*2*num_subjects,'single');
    data_concat_next = zeros(dim_reduced_2,(temporal_dim-1)*2*num_subjects,'single');

    tcount = 1;
    ncount = 1;
    for nfolder = 1:num_subjects
        fprintf('>> %d / %d \n',nfolder,num_subjects);

        subname = split(data_folders_train(nfolder).name,'_');
        subname = subname{1};

        datafile_path_LR = fullfile(data_folders_train(nfolder).folder, ...
                            subname, ...
                            ['s6_rfMRI_',REST_num,'_LR_Atlas_hp2000_clean.dtseries.nii']);
        fprintf(['Loading data ... sub-',subname,' rfMRI LR \n']);
        data = cifti_read(datafile_path_LR);
        data = data.cdata;

        fprintf(['Normalizing data ... sub-',subname,' rfMRI LR \n']);
        data_normalized = normalize(data')';

        time_course_sub = (inv_source * data_normalized)';

        data_concat_prev(:,tcount:tcount+(temporal_dim-1)-1) = time_course_sub(1:end-1,:)';
        data_concat_next(:,tcount:tcount+(temporal_dim-1)-1) = time_course_sub(2:end,:)';
        tcount = tcount + (temporal_dim-1);
        ncount = ncount + dim_reduced_1;
        fprintf(['Time course extracted ... sub-',subname,' rfMRI LR \n']);

        datafile_path_RL = fullfile(data_folders_train(nfolder).folder, ...
                            subname, ...
                            ['s6_rfMRI_',REST_num,'_RL_Atlas_hp2000_clean.dtseries.nii']);
        fprintf(['Loading data ... sub-',subname,' rfMRI RL \n']);
        data = cifti_read(datafile_path_RL);
        data = data.cdata;

        fprintf(['Normalizing data ... sub-',subname,' rfMRI LR \n']);
        data_normalized = normalize(data')';

        time_course_sub = (inv_source * data_normalized)';

        data_concat_prev(:,tcount:tcount+(temporal_dim-1)-1) = time_course_sub(1:end-1,:)';
        data_concat_next(:,tcount:tcount+(temporal_dim-1)-1) = time_course_sub(2:end,:)';
        tcount = tcount + (temporal_dim-1);
        ncount = ncount + dim_reduced_1;
        fprintf(['Time course extracted ... sub-',subname,' rfMRI RL \n']);
    end

    X = data_concat_next; clear data_concat_next
    Y = data_concat_prev; clear data_concat_prev

    %% fbDMD
    TRtarget = 0.72;
    time_len = 1200;

    A1 = X*Y'; A2 = Y*Y';
    A_f = A1 * pinv(A2);
    B1 = Y*X'; B2 = X*X';
    A_b = B1 * pinv(B2);
    A = (A_f/A_b)^0.5;
    A = real(A);


    [Phi_all,D] = eig(A);
    lambda = diag(D);

    abs_DM = abs(lambda);
    period_DM = 2*pi*TRtarget ./ angle(lambda);

    %% DM in original space

    [U,S,V] = svd(Y,'econ');
    temp_v = pinv(Y) * Phi_all;
    
    Z = A2/num_subjects;

    Phi_orig_DL = source_maps' * Phi_all;

    %% reconstruct exact mode
    Phi_orig_exact = zeros(voxel_num,size(Phi_all,2));
    tcount = 1;
    for nfolder = 1:num_subjects
        fprintf('>> %d / %d \n',nfolder,num_subjects);

        subname = split(data_folders_train(nfolder).name,'_');
        subname = subname{1};

        datafile_path_LR = fullfile(data_folders_train(nfolder).folder, ...
                            subname, ...
                            ['s6_rfMRI_',REST_num,'_LR_Atlas_hp2000_clean.dtseries.nii']);
        fprintf(['Loading data ... sub-',subname,' rfMRI LR \n']);
        data = cifti_read(datafile_path_LR);
        data = data.cdata;

        fprintf(['Normalizing data ... sub-',subname,' rfMRI LR \n']);
        data_normalized = normalize(data')';

        temp_prod = data_normalized(:,2:end) * temp_v(tcount:tcount+(temporal_dim-1)-1,:);
        tcount = tcount + (temporal_dim-1);
        Phi_orig_exact = Phi_orig_exact + temp_prod;

        datafile_path_RL = fullfile(data_folders_train(nfolder).folder, ...
                            subname, ...
                            ['s6_rfMRI_',REST_num,'_RL_Atlas_hp2000_clean.dtseries.nii']);
        fprintf(['Loading data ... sub-',subname,' rfMRI RL \n']);
        data = cifti_read(datafile_path_RL);
        data = data.cdata;

        fprintf(['Normalizing data ... sub-',subname,' rfMRI LR \n']);
        data_normalized = normalize(data')';

        temp_prod = data_normalized(:,2:end) * temp_v(tcount:tcount+(temporal_dim-1)-1,:);
        tcount = tcount + (temporal_dim-1);
        Phi_orig_exact = Phi_orig_exact + temp_prod;

        fprintf(['Time course extracted ... sub-',subname,' rfMRI RL \n']);
    end



    %% Test prediction
    
    for i_proj = 1:length(method_subject_proj_list)
        method_subject_proj = method_subject_proj_list{i_proj};        
        
        R2_DM_array_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));
        R2_lin_array_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));
        R2_null_array_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));
        R2_DM_array_cortex_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));
        R2_null_array_cortex_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));
        R2_DM_array_subcortical_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));
        R2_null_array_subcortical_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));

        for nfolder = 1:num_subjects
            tic
            fprintf('>> %d / %d \n',nfolder,num_subjects);

            subname = split(data_folders_test(nfolder).name,'_');
            subname = subname{1};

            datafile_path_LR = fullfile(data_folders_test(nfolder).folder, ...
                                subname, ...
                                ['s6_rfMRI_',REST_num,'_LR_Atlas_hp2000_clean.dtseries.nii']);
            fprintf(['Loading data ... sub-',subname,' rfMRI LR \n']);
            data = cifti_read(datafile_path_LR);
            data = data.cdata;

            fprintf(['Normalizing data ... sub-',subname,' rfMRI LR \n']);
            data_normalized = normalize(data')';

            time_course_sub = (inv_source * data_normalized(:,1:720));        
            if strcmp(method_subject_proj, 'DR')
                temp_v = pinv(time_course_sub(:,2:end)) * Phi_all;
            elseif strcmp(method_subject_proj, 'TL-cov')
                temp_v = (time_course_sub(:,1:end-1)' / (Z)) * Phi_all;
            else
                error('Not defined projection methods!!');
            end
            Phi_orig_sub = data_normalized(:,2:720) * temp_v;
            [R2_DM_array,R2_DM_array_cortex,R2_DM_array_subcortical, R2_null_array,R2_null_array_cortex,R2_null_array_subcortical] ...
            =benchmark_using_DMs(Phi_orig_sub,data_normalized,1,1,0.6,fitting_time_window_list,predict_time_window_list,false);

            R2_DM_array_list(nfolder,1,:,:) = R2_DM_array;
            R2_null_array_list(nfolder,1,:,:) = R2_null_array;
    %         R2_lin_array_list(nfolder,1,:,:) = R2_lin_array;
            R2_DM_array_cortex_list(nfolder,1,:,:) = R2_DM_array_cortex;
            R2_null_array_cortex_list(nfolder,1,:,:) = R2_null_array_cortex;
            R2_DM_array_subcortical_list(nfolder,1,:,:) = R2_DM_array_subcortical;
            R2_null_array_subcortical_list(nfolder,1,:,:) = R2_null_array_subcortical;

            fprintf(['Prediction complete ... sub-',subname,' rfMRI LR \n']);

            datafile_path_RL = fullfile(data_folders_test(nfolder).folder, ...
                                subname, ...
                                ['s6_rfMRI_',REST_num,'_RL_Atlas_hp2000_clean.dtseries.nii']);
            fprintf(['Loading data ... sub-',subname,' rfMRI RL \n']);
            data = cifti_read(datafile_path_RL);
            data = data.cdata;

            fprintf(['Normalizing data ... sub-',subname,' rfMRI LR \n']);
            data_normalized = normalize(data')';

            time_course_sub = (inv_source * data_normalized(:,1:720));        
            if strcmp(method_subject_proj, 'DR')
                temp_v = pinv(time_course_sub(:,2:end)) * Phi_all;
            elseif strcmp(method_subject_proj, 'TL-cov')
                temp_v = (time_course_sub(:,1:end-1)' / (Z)) * Phi_all;
            else
                error('Not defined projection methods!!');
            end
            Phi_orig_sub = data_normalized(:,2:720) * temp_v;
            [R2_DM_array,R2_DM_array_cortex,R2_DM_array_subcortical, R2_null_array,R2_null_array_cortex,R2_null_array_subcortical] ...
            =benchmark_using_DMs(Phi_orig_sub,data_normalized,1,1,0.6,fitting_time_window_list,predict_time_window_list,false);

            R2_DM_array_list(nfolder,2,:,:) = R2_DM_array;
            R2_null_array_list(nfolder,2,:,:) = R2_null_array;
    %         R2_lin_array_list(nfolder,2,:,:) = R2_lin_array;
            R2_DM_array_cortex_list(nfolder,2,:,:) = R2_DM_array_cortex;
            R2_null_array_cortex_list(nfolder,2,:,:) = R2_null_array_cortex;
            R2_DM_array_subcortical_list(nfolder,2,:,:) = R2_DM_array_subcortical;
            R2_null_array_subcortical_list(nfolder,2,:,:) = R2_null_array_subcortical;

            fprintf(['Prediction complete ... sub-',subname,' rfMRI RL \n']);
            toc
        end
        
        if i_proj == 1
            R2_DM_DR_array_list = R2_DM_array_list;
            R2_DM_DR_array_cortex_list = R2_DM_array_cortex_list;
            R2_DM_DR_array_subcortical_list = R2_DM_array_subcortical_list;
        elseif i_proj == 2
            R2_DM_TL_cov_array_list = R2_DM_array_list;
            R2_DM_TL_cov_array_cortex_list = R2_DM_array_cortex_list;
            R2_DM_TL_cov_array_subcortical_list = R2_DM_array_subcortical_list;
        end

        %% display
        % figure;
        % bar( [...
        % squeeze(mean(R2_DM_array_list(:,:,1,:),[1,2])), ...
        % squeeze(mean(R2_lin_array_list(:,:,1,:),[1,2])),...
        % squeeze(mean(R2_null_array_list(:,:,1,:),[1,2])) ...
        % ]');
        % legend({'predict 1s ahead','predict 2s ahead','predict 4s ahead','predict 8s ahead'});
        % clear R2_DM_array_list R2_DM_array_cortex_list R2_DM_array_subcortical_list
    end
    %% 
    % Get the current date and time formatted as YYYYMMDD_HHMMSS
    dtStr = datestr(now, 'yyyymmdd_HHMMSS');

    % Create a filename using the datetime string
    if strcmp(sample_set,'discovery')
        filename = ['results/disc_INF_G_lev_', num2str(target_dim,'%03d'), '_MIGP_results_' dtStr '.mat'];
    elseif strcmp(sample_set, 'replication')
        filename = ['results/repl_INF_G_lev_', num2str(target_dim,'%03d'), '_MIGP_results_' dtStr '.mat'];
    else
        error('Undefined sample set!!');
    end
    save(filename, 'lambda', 'Phi_all', 'Phi_orig_DL', 'Phi_orig_exact', 'D', 'W', 'A', 'Z', 'source_maps', 'inv_source','R2*list','*time_window_list');
end
