clear; clc;

% gift; close all;

current_path = pwd;

%% Run options

method_subject_proj_list = {'DR','TL-cov'}; % DR or TL-cov

target_dim_list = 10:10:100; 
% target_dim_list = 21:1:29; 

fitting_time_window_list =  [1];
predict_time_window_list = [1,2,4,8];

%%
data_load = load('secure_info\path_info.mat');
HCP_data_path = data_load.HCP_denoised_path;
data_load = load('results/split_subjects.mat');
data_folders = data_load.data_folders;

sub_ids = data_load.sub_ids;
sub_ids_set1_explore = data_load.sub_ids_set1_explore;
[~,IA,~] = intersect(sub_ids,sub_ids_set1_explore);
data_folders = data_folders(IA);

sub_ids = data_load.sub_ids;
sub_ids_set1_test = data_load.sub_ids_set1_test;
[~,IA,~] = intersect(sub_ids,sub_ids_set1_test);
data_folders_test = data_load.data_folders;
data_folders_test = data_folders_test(IA);

%%

num_subjects = length(data_folders_test);

for target_dim = target_dim_list

    %% Test prediction
    TRtarget = 0.72;

    time_len = 1200;

    R2_DM_array_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));
    R2_lin_array_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));
    R2_null_array_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));
    R2_DM_array_cortex_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));
    R2_null_array_cortex_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));
    R2_DM_array_subcortical_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));
    R2_null_array_subcortical_list = zeros(num_subjects,2,length(fitting_time_window_list),length(predict_time_window_list));

    for nfolder = 1:num_subjects
        fprintf('>> %d / %d \n',nfolder,num_subjects);
        tic

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
        data_normalized_centered = data_normalized - mean(data_normalized);

        fprintf(['Performing single subj ICA ... sub-',subname,' rfMRI LR \n']);

        [~,score,lambda_k] = pca(data_normalized_centered(:,1:720),'NumComponents',target_dim,'Centered',true);
        score_whitened = score * diag(1 ./ sqrt(lambda_k(1:target_dim) + eps));
        [~, ~, ~, source_maps] = icatb_icaAlgorithm('infomax', score_whitened', {'verbose', 'off'});

        inv_source = pinv(source_maps');

        time_course_sub = (inv_source * data_normalized(:,1:720));

        [Phi_all,~] = performDMD(time_course_sub(:,2:end),time_course_sub(:,1:end-1),time_len,TRtarget);

    %     [R2_DM_array, R2_null_array, R2_lin_array,~, ~] ...
    %     = benchmark_using_DMs(Phi_all,data_normalized,source_maps',inv_source);

        temp_v = pinv(time_course_sub(:,1:end-1)) * Phi_all;
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
        data_normalized_centered = data_normalized - mean(data_normalized);

        [~,score,lambda_k] = pca(data_normalized_centered(:,1:720),'NumComponents',target_dim,'Centered',true);
        score_whitened = score * diag(1 ./ sqrt(lambda_k(1:target_dim) + eps));
        [~, ~, ~, source_maps] = icatb_icaAlgorithm('infomax', score_whitened', {'verbose', 'off'});

        inv_source = pinv(source_maps');

        time_course_sub = (inv_source * data_normalized(:,1:720));

        [Phi_all,~] = performDMD(time_course_sub(:,2:end),time_course_sub(:,1:end-1),time_len,TRtarget);

    %     [R2_DM_array, R2_null_array, R2_lin_array,~, ~] ...
    %     = benchmark_using_DMs(Phi_all,data_normalized,source_maps',inv_source);

        temp_v = pinv(time_course_sub(:,1:end-1)) * Phi_all;
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


    %% display
    figure;
    bar( [...
    squeeze(mean(R2_DM_array_list(:,:,1,:),[1,2])), ...
    squeeze(mean(R2_lin_array_list(:,:,1,:),[1,2])),...
    squeeze(mean(R2_null_array_list(:,:,1,:),[1,2])) ...
    ]');
    legend({'predict 1s ahead','predict 2s ahead','predict 4s ahead','predict 8s ahead'});

    %% 
    % Get the current date and time formatted as YYYYMMDD_HHMMSS
    dtStr = datestr(now, 'yyyymmdd_HHMMSS');

    % Create a filename using the datetime string
    filename = ['results/subject_wise_ica',num2str(target_dim,'%03d'),'_dmd_results_normalized_' dtStr '.mat'];
    save(filename, 'R2*array_list','*time_window_list');
end

%%
function [Phi_all,lambda] = performDMD(X,Y,time_len,TRtarget)
    A1 = X*Y'; A2 = Y*Y';
    A_f = A1 * pinv(A2);
    B1 = Y*X'; B2 = X*X';
    A_b = B1 * pinv(B2);
    A = (A_f/A_b)^0.5;
    A = real(A);


    [Phi_all,D] = eig(A);
    lambda = diag(D);
    
%     idx_exclude = abs(lambda) < 1e-4 | 2*pi*TRtarget ./ abs(angle(lambda)) > time_len * TRtarget;
%     lambda(idx_exclude) = [];
%     Phi_all(:,idx_exclude) = [];
end
