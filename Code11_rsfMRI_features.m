clear; clc;

current_path = pwd;

sample_set = 'discovery';
% sample_set = 'replication';

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
left_cortex_dim = data.diminfo{1, 1}.models{1, 2}.start - 1;
cortex_dim = data.diminfo{1, 1}.models{1, 3}.start - 1;
data = data.cdata;

spatial_dim = size(data,1);
temporal_dim = size(data,2);
num_subjects = 50;

%%

labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

PCC_label = [35,35+180];
aIns_R_label = 111;
SMG_label = [148,148+180];
FEF_label = [10,10+180];
TPJ_label = [139,139+180];

PCC_idx = false(cortex_dim,1);
PCC_idx(labels.cdata == PCC_label(1) | labels.cdata == PCC_label(2)) = true;
% PCC_idx = [PCC_idx;false(voxel_num-cortex_dim,1)];

aIns_R_idx = false(cortex_dim,1);
aIns_R_idx(labels.cdata == aIns_R_label) = true;
% aIns_R_idx = [aIns_R_idx;false(voxel_num-cortex_dim,1)];

SMG_idx = false(cortex_dim,1);
SMG_idx(labels.cdata == SMG_label(1) | labels.cdata == SMG_label(2)) = true;
% SMG_idx = [SMG_idx;false(voxel_num-cortex_dim,1)];

FEF_idx = false(cortex_dim,1);
FEF_idx(labels.cdata == FEF_label(1) | labels.cdata == FEF_label(2)) = true;

TPJ_idx = false(cortex_dim,1);
TPJ_idx(labels.cdata == TPJ_label(1) | labels.cdata == TPJ_label(2)) = true;

seed_conn_PCC = zeros(cortex_dim,1,'single');
seed_conn_aIns_R = zeros(cortex_dim,1,'single');
seed_conn_SMG = zeros(cortex_dim,1,'single');
seed_conn_FEF = zeros(cortex_dim,1,'single');
seed_conn_TPJ = zeros(cortex_dim,1,'single');

for nfolder = 1:num_subjects
    fprintf('>> %d / %d \n',nfolder,num_subjects); tic;
    
    subname = split(data_folders(nfolder).name,'_');
    REST_num = subname{4};
    subname = subname{1};
    
    %%%% LR
    datafile_path_LR = fullfile(data_folders(nfolder).folder, ...
                        data_folders(nfolder).name, ...
                        subname,'MNINonLinear','Results', ...
                        ['rfMRI_',REST_num,'_LR'], ...
                        ['s6_rfMRI_',REST_num,'_LR_Atlas_hp2000_clean.dtseries.nii']);
    fprintf(['Loading data ... sub-',subname,' rfMRI LR \n']);
    data = cifti_read(datafile_path_LR);
    data = data.cdata;
    
    data_filtered = bandpass(data(1:cortex_dim,:)',[0.01,0.1],1/0.72)';
    
    %%% GSR
    GS = mean(data_filtered(1:cortex_dim,:));
    GS_mat = [ones(temporal_dim,1),GS'];
    data_filtered_gsr = ((eye(temporal_dim) - GS_mat * pinv(GS_mat)) * data_filtered(1:cortex_dim,:)')';
    
    fprintf(['GSR ... sub-',subname,' rfMRI LR \n']);
    
    fprintf(['Calculating connectivity ... sub-',subname,' rfMRI LR \n']);
    seed_conn_PCC = seed_conn_PCC + corr(mean(data_filtered_gsr(PCC_idx,:))',data_filtered_gsr(1:cortex_dim,:)')';
    seed_conn_aIns_R = seed_conn_aIns_R + corr(mean(data_filtered_gsr(aIns_R_idx,:))',data_filtered_gsr(1:cortex_dim,:)')';
    seed_conn_SMG = seed_conn_SMG + corr(mean(data_filtered_gsr(SMG_idx,:))',data_filtered_gsr(1:cortex_dim,:)')';
    seed_conn_FEF = seed_conn_FEF + corr(mean(data_filtered_gsr(FEF_idx,:))',data_filtered_gsr(1:cortex_dim,:)')';
    seed_conn_TPJ = seed_conn_TPJ + corr(mean(data_filtered_gsr(TPJ_idx,:))',data_filtered_gsr(1:cortex_dim,:)')';

    fprintf(['Connectivity saved ... sub-',subname,' rfMRI LR \n']);
    
    %%%% RL
    datafile_path_RL = fullfile(data_folders(nfolder).folder, ...
                        data_folders(nfolder).name, ...
                        subname,'MNINonLinear','Results', ...
                        ['rfMRI_',REST_num,'_RL'], ...
                        ['s6_rfMRI_',REST_num,'_RL_Atlas_hp2000_clean.dtseries.nii']);
    fprintf(['Loading data ... sub-',subname,' rfMRI RL \n']);
    data = cifti_read(datafile_path_RL);
    data = data.cdata;
    
    data_filtered = bandpass(data(1:cortex_dim,:)',[0.01,0.1],1/0.72)';
    
    %%% GSR
    GS = mean(data_filtered(1:cortex_dim,:));
    GS_mat = [ones(temporal_dim,1),GS'];
    data_filtered_gsr = ((eye(temporal_dim) - GS_mat * pinv(GS_mat)) * data_filtered(1:cortex_dim,:)')';
    
    fprintf(['GSR ... sub-',subname,' rfMRI RL \n']);
    
    fprintf(['Calculating connectivity ... sub-',subname,' rfMRI RL \n']);
    seed_conn_PCC = seed_conn_PCC + corr(mean(data_filtered_gsr(PCC_idx,:))',data_filtered_gsr(1:cortex_dim,:)')';
    seed_conn_aIns_R = seed_conn_aIns_R + corr(mean(data_filtered_gsr(aIns_R_idx,:))',data_filtered_gsr(1:cortex_dim,:)')';
    seed_conn_SMG = seed_conn_SMG + corr(mean(data_filtered_gsr(SMG_idx,:))',data_filtered_gsr(1:cortex_dim,:)')';
    seed_conn_FEF = seed_conn_FEF + corr(mean(data_filtered_gsr(FEF_idx,:))',data_filtered_gsr(1:cortex_dim,:)')';
    seed_conn_TPJ = seed_conn_TPJ + corr(mean(data_filtered_gsr(TPJ_idx,:))',data_filtered_gsr(1:cortex_dim,:)')';

    fprintf(['Connectivity saved ... sub-',subname,' rfMRI RL \n']); toc
end
seed_conn_PCC = seed_conn_PCC ./ (2 * num_subjects);
seed_conn_aIns_R = seed_conn_aIns_R ./ (2 * num_subjects);
seed_conn_SMG = seed_conn_SMG ./ (2 * num_subjects);
seed_conn_FEF = seed_conn_FEF ./ (2 * num_subjects);
seed_conn_TPJ = seed_conn_TPJ ./ (2 * num_subjects);

%%
save results/seed_based_conn seed_conn*

%% Dispaly seed-based fc
lh_surface_file = 'atlas/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii';
rh_surface_file = 'atlas/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii';

save_dir = 'rsfMRI_features';
mkdir(save_dir);

figure; fig = gcf;

snapshot_cifti = cifti_struct_create_from_template(labels,seed_conn_PCC, 'dscalar');
display_cifti_cortex(fig, snapshot_cifti, lh_surface_file, rh_surface_file, [], [], [], []);
fname = fullfile(save_dir, 'seed_conn_PCC.jpg');
print(fig, fname, '-djpeg', '-r300');

snapshot_cifti = cifti_struct_create_from_template(labels,seed_conn_aIns_R, 'dscalar');
display_cifti_cortex(fig, snapshot_cifti, lh_surface_file, rh_surface_file, [], [], [], []);
fname = fullfile(save_dir, 'seed_conn_aIns_R.jpg');
print(fig, fname, '-djpeg', '-r300');

snapshot_cifti = cifti_struct_create_from_template(labels,seed_conn_SMG, 'dscalar');
display_cifti_cortex(fig, snapshot_cifti, lh_surface_file, rh_surface_file, [], [], [], []);
fname = fullfile(save_dir, 'seed_conn_SMG.jpg');
print(fig, fname, '-djpeg', '-r300');

snapshot_cifti = cifti_struct_create_from_template(labels,seed_conn_FEF, 'dscalar');
display_cifti_cortex(fig, snapshot_cifti, lh_surface_file, rh_surface_file, [], [], [], []);
fname = fullfile(save_dir, 'seed_conn_FEFjpg');
print(fig, fname, '-djpeg', '-r300');


snapshot_cifti = cifti_struct_create_from_template(labels,seed_conn_TPJ, 'dscalar');
display_cifti_cortex(fig, snapshot_cifti, lh_surface_file, rh_surface_file, [], [], [], []);
fname = fullfile(save_dir, 'seed_conn_TPJ.jpg');
print(fig, fname, '-djpeg', '-r300');

close(fig);

%% MLI

MLI_all = zeros(num_subjects,2,cortex_dim);

frame_dt = 0.72;
time_window_length = 30;
time_window_interval = 10;

for nfolder = 1:num_subjects
    fprintf('>> %d / %d \n',nfolder,num_subjects); tic;
    
    subname = split(data_folders(nfolder).name,'_');
    REST_num = subname{4};
    subname = subname{1};
    
    %%%% LR
    datafile_path_LR = fullfile(data_folders(nfolder).folder, ...
                        data_folders(nfolder).name, ...
                        subname,'MNINonLinear','Results', ...
                        ['rfMRI_',REST_num,'_LR'], ...
                        ['s6_rfMRI_',REST_num,'_LR_Atlas_hp2000_clean.dtseries.nii']);
    fprintf(['Loading data ... sub-',subname,' rfMRI LR \n']);
    data = cifti_read(datafile_path_LR);
    data = data.cdata;
    
    data_filtered = bandpass(data(1:cortex_dim,:)',[0.01,0.1],1/0.72)';
    
    GS_L = mean(data_filtered(1:left_cortex_dim,:));
    GS_R = mean(data_filtered(left_cortex_dim+1:cortex_dim,:));
    
    fprintf(['Calculating MLI ... sub-',subname,' rfMRI LR \n']);
    
    MLI = zeros(cortex_dim,1);
    
    current_frame = 1;
    current_time = current_frame*frame_dt;
    accum_num = 0;
    while current_time+time_window_length < frame_dt*size(data_filtered,2)
        accum_num = accum_num + 1;

        current_GS_L = GS_L(current_frame:round(current_frame+time_window_length/frame_dt));
        current_GS_R = GS_R(current_frame:round(current_frame+time_window_length/frame_dt));

        current_time_series_matrix = data_filtered(:,current_frame:round(current_frame+time_window_length/frame_dt));

        r_L = corr(current_GS_L',current_time_series_matrix');
        r_R = corr(current_GS_R',current_time_series_matrix');

        MLI = MLI + (r_L-r_R)';

        current_frame = current_frame + round(time_window_interval/frame_dt);
        current_time = current_frame*frame_dt;
    end

    MLI = MLI / accum_num;

    MLI_all(nfolder,1,:) = MLI;

    fprintf(['MLI saved ... sub-',subname,' rfMRI LR \n']);
    
    %%%% RL
    datafile_path_RL = fullfile(data_folders(nfolder).folder, ...
                        data_folders(nfolder).name, ...
                        subname,'MNINonLinear','Results', ...
                        ['rfMRI_',REST_num,'_RL'], ...
                        ['s6_rfMRI_',REST_num,'_RL_Atlas_hp2000_clean.dtseries.nii']);
    fprintf(['Loading data ... sub-',subname,' rfMRI RL \n']);
    data = cifti_read(datafile_path_RL);
    data = data.cdata;
    
    data_filtered = bandpass(data(1:cortex_dim,:)',[0.01,0.1],1/0.72)';
    
    GS_L = mean(data_filtered(1:left_cortex_dim,:));
    GS_R = mean(data_filtered(left_cortex_dim+1:cortex_dim,:));
    
    fprintf(['Calculating MLI ... sub-',subname,' rfMRI LR \n']);
    
    MLI = zeros(cortex_dim,1);
    
    current_frame = 1;
    current_time = current_frame*frame_dt;
    accum_num = 0;
    while current_time+time_window_length < frame_dt*size(data_filtered,2)
        accum_num = accum_num + 1;

        current_GS_L = GS_L(current_frame:round(current_frame+time_window_length/frame_dt));
        current_GS_R = GS_R(current_frame:round(current_frame+time_window_length/frame_dt));

        current_time_series_matrix = data_filtered(:,current_frame:round(current_frame+time_window_length/frame_dt));

        r_L = corr(current_GS_L',current_time_series_matrix');
        r_R = corr(current_GS_R',current_time_series_matrix');

        MLI = MLI + (r_L-r_R)';

        current_frame = current_frame + round(time_window_interval/frame_dt);
        current_time = current_frame*frame_dt;
    end

    MLI = MLI / accum_num;

    MLI_all(nfolder,2,:) = MLI;

    fprintf(['MLI saved ... sub-',subname,' rfMRI LR \n']); toc
end

%%
save results/MLI MLI_all

%% Dispaly MLI
labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

lh_surface_file = 'atlas/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii';
rh_surface_file = 'atlas/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii';

save_dir = 'rsfMRI_features';
mkdir(save_dir);

mean_MLI = squeeze(mean(MLI_all,[1,2]));

figure; fig = gcf;

max_abs_val = max(abs(mean_MLI));

snapshot_cifti = cifti_struct_create_from_template(labels, mean_MLI, 'dscalar');
display_cifti_cortex(fig, snapshot_cifti, lh_surface_file, rh_surface_file, [], -max_abs_val, max_abs_val, []);
fname = fullfile(save_dir, 'MLI.jpg');
print(fig, fname, '-djpeg', '-r300');

close(fig);

%%
num_subjects = 50;

conn_mat = zeros(cortex_dim,cortex_dim,'single');

for nfolder = 1:num_subjects
    fprintf('>> %d / %d \n',nfolder,num_subjects); tic;
    
    subname = split(data_folders(nfolder).name,'_');
    REST_num = subname{4};
    subname = subname{1};
    
    %%%% LR
    datafile_path_LR = fullfile(data_folders(nfolder).folder, ...
                        data_folders(nfolder).name, ...
                        subname,'MNINonLinear','Results', ...
                        ['rfMRI_',REST_num,'_LR'], ...
                        ['s6_rfMRI_',REST_num,'_LR_Atlas_hp2000_clean.dtseries.nii']);
    fprintf(['Loading data ... sub-',subname,' rfMRI LR \n']);
    data = cifti_read(datafile_path_LR);
    data = data.cdata;
    
    data_filtered = bandpass(data(1:cortex_dim,:)',[0.01,0.1],1/0.72)';

    fprintf(['Calculating connectivity ... sub-',subname,' rfMRI LR \n']);
    conn_mat = conn_mat + corr(data_filtered(1:cortex_dim,:)');
    
    %%%% RL
    datafile_path_RL = fullfile(data_folders(nfolder).folder, ...
                        data_folders(nfolder).name, ...
                        subname,'MNINonLinear','Results', ...
                        ['rfMRI_',REST_num,'_RL'], ...
                        ['s6_rfMRI_',REST_num,'_RL_Atlas_hp2000_clean.dtseries.nii']);
    fprintf(['Loading data ... sub-',subname,' rfMRI RL \n']);
    data = cifti_read(datafile_path_RL);
    data = data.cdata;
    
    data_filtered = bandpass(data(1:cortex_dim,:)',[0.01,0.1],1/0.72)';
    
    fprintf(['Calculating connectivity ... sub-',subname,' rfMRI RL \n']);
    conn_mat = conn_mat + corr(data_filtered(1:cortex_dim,:)');

    fprintf(['Connectivity saved ... sub-',subname,' rfMRI RL \n']); toc
end

conn_mat = conn_mat / (2 * num_subjects);
