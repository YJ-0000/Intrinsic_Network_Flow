clear; clc;

% gift; close all;

current_path = pwd;

%% Run options

method_subject_proj_list = {'DR','TL-cov'}; % DR or TL-cov

target_dim_list = 10:10:100; 
% target_dim_list = 21:1:29; 

sample_set = 'discovery';
% sample_set = 'replication';

fitting_time_window_list =  [1];
predict_time_window_list = [1,2,4,8];

%%

data_load = load('secure_info\path_info.mat');
HCP_data_path = data_load.HCP_denoised_path;
data_load = load('results/split_subjects.mat');

data_folders = data_load.data_folders;
sub_ids = data_load.sub_ids;
if strcmp(sample_set,'discovery')
    sub_ids_set_explore = data_load.sub_ids_set1_explore;
elseif strcmp(sample_set, 'replication')
    sub_ids_set_explore = data_load.sub_ids_set2_explore;
else
    error('Undefined sample set!!');
end
[~,IA,~] = intersect(sub_ids,sub_ids_set_explore);
data_folders = data_folders(IA);

sub_ids = data_load.sub_ids;
if strcmp(sample_set,'discovery')
    sub_ids_set_test = data_load.sub_ids_set1_test;
elseif strcmp(sample_set, 'replication')
    sub_ids_set_test = data_load.sub_ids_set2_test;
else
    error('Undefined sample set!!');
end
[~,IA,~] = intersect(sub_ids,sub_ids_set_test);
data_folders_test = data_load.data_folders;
data_folders_test = data_folders_test(IA);

%% Loading one file

num_data_used = 1;

comp_est_accum = zeros(num_data_used*2,1);

nfolder = 1;

subname = split(data_folders(nfolder).name,'_');
REST_num = subname{4};
subname = subname{1};

datafile_path_LR = fullfile(data_folders(nfolder).folder, ...
                    data_folders(nfolder).name, ...
                    subname,'MNINonLinear','Results', ...
                    ['rfMRI_',REST_num,'_LR'], ...
                    ['s6_rfMRI_',REST_num,'_LR_Atlas_hp2000_clean.dtseries.nii']);
fprintf(['Loading data ... sub-',subname,' rfMRI LR \n']);
data = cifti_read(datafile_path_LR);
data = data.cdata;
voxel_num = size(data,1);
temporal_dim = size(data,2);

%% ICA - data concatanation


for target_dim = target_dim_list

    num_subjects = 50;
    
    dim_reduced_1 = round(1.5*target_dim);
    dim_reduced_2 = target_dim;

    data_concat = zeros(voxel_num,dim_reduced_1*2*num_subjects,'single');

    tcount = 1;
    for nfolder = 1:num_subjects
        fprintf('>> %d / %d \n',nfolder,num_subjects);

        subname = split(data_folders(nfolder).name,'_');
        REST_num = subname{4};
        subname = subname{1};

        datafile_path_LR = fullfile(data_folders(nfolder).folder, ...
                            data_folders(nfolder).name, ...
                            subname,'MNINonLinear','Results', ...
                            ['rfMRI_',REST_num,'_LR'], ...
                            ['s6_rfMRI_',REST_num,'_LR_Atlas_hp2000_clean.dtseries.nii']);
        fprintf(['Loading data ... sub-',subname,' rfMRI LR \n']);
        data = cifti_read(datafile_path_LR);
        data = data.cdata;

        fprintf(['Normalizing data ... sub-',subname,' rfMRI LR \n']);
        data_normalized = normalize(data')';
    %     data_normalized_centered = data_normalized - mean(data_normalized);

        fprintf(['Reducing data ... sub-',subname,' rfMRI LR \n']);
        [~,score,lambda_k] = pca(data_normalized,'NumComponents',dim_reduced_1,'Centered',true);
        data_reduced = score*diag(1 ./ sqrt(lambda_k(1:dim_reduced_1) + eps));

        data_concat(:,tcount:tcount+dim_reduced_1-1) = data_reduced;
        tcount = tcount + dim_reduced_1;
        fprintf(['Data saved ... sub-',subname,' rfMRI LR \n']);

        datafile_path_RL = fullfile(data_folders(nfolder).folder, ...
                            data_folders(nfolder).name, ...
                            subname,'MNINonLinear','Results', ...
                            ['rfMRI_',REST_num,'_RL'], ...
                            ['s6_rfMRI_',REST_num,'_RL_Atlas_hp2000_clean.dtseries.nii']);
        fprintf(['Loading data ... sub-',subname,' rfMRI RL \n']);
        data = cifti_read(datafile_path_RL);
        data = data.cdata;

        fprintf(['Normalizing data ... sub-',subname,' rfMRI LR \n']);
        data_normalized = normalize(data')';
    %     data_normalized_centered = data_normalized - mean(data_normalized);

         fprintf(['Reducing data ... sub-',subname,' rfMRI RL \n']);
        [~,score,lambda_k] = pca(data_normalized,'NumComponents',dim_reduced_1,'Centered',true);
        data_reduced = score*diag(1 ./ sqrt(lambda_k(1:dim_reduced_1) + eps));

        data_concat(:,tcount:tcount+dim_reduced_1-1) = data_reduced;
        tcount = tcount + dim_reduced_1;
        fprintf(['Data saved ... sub-',subname,' rfMRI RL \n']);
    end

    %% ICA - Group-level reduction
    [G_group,score,lambda_k] = pca(data_concat,'NumComponents',dim_reduced_2,'Centered',true);
    data_concat_reduced = score*diag(1 ./ sqrt(lambda_k(1:dim_reduced_2) + eps));
    G_group_whitening = G_group * diag(1 ./ sqrt(lambda_k(1:dim_reduced_2) + eps));
    lambda_k_group = lambda_k;
    %% ICA - Group-level ICA
    data_concat_reduced_removeNAN = data_concat_reduced;
    mask_vec = sum(isnan(data_concat_reduced),2)==0;
    data_concat_reduced_removeNAN(sum(isnan(data_concat_reduced),2)>0,:) = [];
    [icaAlgo, W, W_inv, source_maps] = icatb_icaAlgorithm('infomax', data_concat_reduced_removeNAN');

    %% save source maps
    labels = cifti_read('atlas/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii');
    source_maps_oridim = zeros(size(data_concat_reduced))';
    source_maps_oridim(:,mask_vec) = source_maps;
    mkdir('source_maps');
    for k = 1:target_dim
        source_map_cifti = cifti_struct_create_from_template(labels,source_maps_oridim(k,:)', 'dscalar');
        cifti_write(source_map_cifti, ['source_maps/source_',num2str(k),'.dscalar.nii']);
    end

    %% Back reconstruction

    inv_source = pinv(source_maps_oridim');

    data_concat_prev = zeros(dim_reduced_2,(temporal_dim-1)*2*num_subjects,'single');
    data_concat_next = zeros(dim_reduced_2,(temporal_dim-1)*2*num_subjects,'single');

    tcount = 1;
    ncount = 1;
    for nfolder = 1:num_subjects
        fprintf('>> %d / %d \n',nfolder,num_subjects);

        subname = split(data_folders(nfolder).name,'_');
        REST_num = subname{4};
        subname = subname{1};

        datafile_path_LR = fullfile(data_folders(nfolder).folder, ...
                            data_folders(nfolder).name, ...
                            subname,'MNINonLinear','Results', ...
                            ['rfMRI_',REST_num,'_LR'], ...
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

        datafile_path_RL = fullfile(data_folders(nfolder).folder, ...
                            data_folders(nfolder).name, ...
                            subname,'MNINonLinear','Results', ...
                            ['rfMRI_',REST_num,'_RL'], ...
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
    
    YY = A2/num_subjects;

    %% reconstruct exact mode
    Phi_orig = zeros(voxel_num,size(Phi_all,2));
    tcount = 1;
    for nfolder = 1:num_subjects
        fprintf('>> %d / %d \n',nfolder,num_subjects);

        subname = split(data_folders(nfolder).name,'_');
        REST_num = subname{4};
        subname = subname{1};

        datafile_path_LR = fullfile(data_folders(nfolder).folder, ...
                            data_folders(nfolder).name, ...
                            subname,'MNINonLinear','Results', ...
                            ['rfMRI_',REST_num,'_LR'], ...
                            ['s6_rfMRI_',REST_num,'_LR_Atlas_hp2000_clean.dtseries.nii']);
        fprintf(['Loading data ... sub-',subname,' rfMRI LR \n']);
        data = cifti_read(datafile_path_LR);
        data = data.cdata;

        fprintf(['Normalizing data ... sub-',subname,' rfMRI LR \n']);
        data_normalized = normalize(data')';
        data_normalized = data_normalized - mean(data_normalized);

        temp_prod = data_normalized(:,2:end) * temp_v(tcount:tcount+(temporal_dim-1)-1,:);
        tcount = tcount + (temporal_dim-1);
        Phi_orig = Phi_orig + temp_prod;

        datafile_path_RL = fullfile(data_folders(nfolder).folder, ...
                            data_folders(nfolder).name, ...
                            subname,'MNINonLinear','Results', ...
                            ['rfMRI_',REST_num,'_RL'], ...
                            ['s6_rfMRI_',REST_num,'_RL_Atlas_hp2000_clean.dtseries.nii']);
        fprintf(['Loading data ... sub-',subname,' rfMRI RL \n']);
        data = cifti_read(datafile_path_RL);
        data = data.cdata;

        fprintf(['Normalizing data ... sub-',subname,' rfMRI LR \n']);
        data_normalized = normalize(data')';
        data_normalized = data_normalized - mean(data_normalized);

        temp_prod = data_normalized(:,2:end) * temp_v(tcount:tcount+(temporal_dim-1)-1,:);
        tcount = tcount + (temporal_dim-1);
        Phi_orig = Phi_orig + temp_prod;

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
            REST_num = subname{4};
            subname = subname{1};

            datafile_path_LR = fullfile(data_folders_test(nfolder).folder, ...
                                data_folders_test(nfolder).name, ...
                                subname,'MNINonLinear','Results', ...
                                ['rfMRI_',REST_num,'_LR'], ...
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
                temp_v = (time_course_sub(:,1:end-1)' / (YY)) * Phi_all;
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
                                data_folders_test(nfolder).name, ...
                                subname,'MNINonLinear','Results', ...
                                ['rfMRI_',REST_num,'_RL'], ...
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
                temp_v = (time_course_sub(:,1:end-1)' / (YY)) * Phi_all;
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
        figure;
        bar( [...
        squeeze(mean(R2_DM_array_list(:,:,1,:),[1,2])), ...
        squeeze(mean(R2_lin_array_list(:,:,1,:),[1,2])),...
        squeeze(mean(R2_null_array_list(:,:,1,:),[1,2])) ...
        ]');
        legend({'predict 1s ahead','predict 2s ahead','predict 4s ahead','predict 8s ahead'});
        clear R2_DM_array_list R2_DM_array_cortex_list R2_DM_array_subcortical_list
    end
    %% 
    % Get the current date and time formatted as YYYYMMDD_HHMMSS
    dtStr = datestr(now, 'yyyymmdd_HHMMSS');

    % Create a filename using the datetime string
    if strcmp(sample_set,'discovery')
        filename = ['results/loop_gica_', num2str(target_dim,'%03d'), '_dmd_results_normalized_' dtStr '.mat'];
    elseif strcmp(sample_set, 'replication')
        filename = ['results/repl_gica_', num2str(target_dim,'%03d'), '_dmd_results_normalized_' dtStr '.mat'];
    else
        error('Undefined sample set!!');
    end
    save(filename, 'lambda', 'Phi_all', 'Phi_orig', 'D', 'W', 'A', 'source_maps', 'G_group','lambda_k_group','inv_source','R2*list','*time_window_list');
end
