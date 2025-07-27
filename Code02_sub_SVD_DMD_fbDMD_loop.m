clear; clc;

% gift; close all;

current_path = pwd;

%% Run options

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
%         data_normalized_centered = data_normalized - mean(data_normalized);

        [Phi_orig_sub,~] = perform_exact_fbDMD(data_normalized(:,2:720),data_normalized(:,1:720-1),target_dim);

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
%         data_normalized_centered = data_normalized - mean(data_normalized);
        
        [Phi_orig_sub,~] = perform_exact_fbDMD(data_normalized(:,2:720),data_normalized(:,1:720-1),target_dim);

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
    filename = ['results/subject_wise_exact_',num2str(target_dim,'%03d'),'_dmd_results_normalized_' dtStr '.mat'];
    save(filename, 'R2*array_list','*time_window_list');
end

