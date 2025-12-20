clear; clc;

startup;

current_path = pwd;

%% Get env variables
slurm_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));

%% Run options

method_subject_proj_list = {'DR','TL-cov'}; % DR or TL-cov

max_dim = 100;
target_dim_list = 27;
if isnan(slurm_id) || slurm_id == 1 
    % target_dim_list = 10:10:100; 
    REST_num = 'REST1';
elseif slurm_id == 2 
    % target_dim_list = 21:1:29; 
    REST_num = 'REST2';
end

fprintf('Processing ses-%s data in HCP... \n',REST_num);

%%

data_load = load('secure_info/path_info.mat');
HCP_data_path = data_load.HCP_denoised_path;
data_load = load('results/split_subjects.mat');
sub_idx_include = data_load.sub_idx_include;
is_cognitive_impaired = data_load.is_cognitive_impaired;
is_excluded_due_movement = data_load.is_excluded_due_movement;
is_RL_processing_errors = data_load.is_RL_processing_errors;

data_folders = dir(HCP_data_path);
isFolder = [data_folders.isdir];
names = {data_folders.name};
mask = isFolder & ~ismember(names, {'.','..'});
data_folders = data_folders(mask);

fprintf('Checking all data.. \n\n')
if strcmp(REST_num,'REST1')
    data_folders_train = data_folders(sub_idx_include);
elseif strcmp(REST_num,'REST2')
    does_not_data_exist = true(length(data_folders),1);
    for nsub = 1:length(data_folders)
        fprintf('.. checking %03d/%03d \n',nsub,length(data_folders));
        ciftifile_LR = fullfile(data_folders(nsub).folder,...
                data_folders(nsub).name,...
                '*REST2_LR*.nii');
        ciftifile_LR = dir(ciftifile_LR);
        ciftifile_RL = fullfile(data_folders(nsub).folder,...
                data_folders(nsub).name,...
                '*REST2_RL*.nii');
        ciftifile_RL = dir(ciftifile_RL);
        if length(ciftifile_LR) == 1 && length(ciftifile_RL) == 1
            data_LR = cifti_read(fullfile(ciftifile_LR(1).folder,ciftifile_LR(1).name));
            if size(data_LR.cdata,2) == 1200
                data_RL = cifti_read(fullfile(ciftifile_RL(1).folder,ciftifile_RL(1).name));
                if size(data_RL.cdata,2) == 1200
                    does_not_data_exist(nsub) = false;
                end
            end
            
        end
    end

    sub_idx_exclude = does_not_data_exist | is_cognitive_impaired | is_excluded_due_movement(:,2) | is_RL_processing_errors;
    sub_idx_include = ~sub_idx_exclude;
    data_folders_train = data_folders(sub_idx_include);
else
    error('Unrecognized session name.')
end


%% Loading one file

num_data_used = 1;

comp_est_accum = zeros(num_data_used*2,1);

nfolder = 1;

subname = split(data_folders_train(nfolder).name,'_');
% REST_num = 'REST1';
subname = subname{1};

datafile_path_LR = fullfile(data_folders_train(nfolder).folder, subname,...
                    ['s6_rfMRI_',REST_num,'_LR_Atlas_hp2000_clean.dtseries.nii']);
fprintf(['Loading data ... sub-',subname,' rfMRI LR \n']);
data = cifti_read(datafile_path_LR);
data = data.cdata;
voxel_num = size(data,1);
temporal_dim = size(data,2);

%% Group PCA

num_subjects = length(data_folders_train);

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
    
    assert(exist(datafile_path_LR,'file') && exist(datafile_path_RL,'file'),sprintf('No file for sub-%s',subname));

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



    %% 
    % Get the current date and time formatted as YYYYMMDD_HHMMSS
    dtStr = datestr(now, 'yyyymmdd_HHMMSS');

    % Create a filename using the datetime string
    filename = ['results/INF_G_lev_',REST_num,'_ALL_', num2str(target_dim,'%03d'), '_MIGP_results_' dtStr '.mat'];
    
    save(filename, 'lambda', 'Phi_all', 'Phi_orig_DL', 'Phi_orig_exact', 'D', 'W', 'A', 'Z', 'source_maps', 'inv_source');
end
