clear; clc;

% gift; close all;

current_path = pwd;

%% Run options

sample_set = 'discovery';
% sample_set = 'replication';

target_dim = 27;
TRtarget = 0.72;
cortex_dim = 59412;

%% load G-ICA results
g_ica_result_files = dir(['results/loop_gica_',num2str(target_dim,'%03d'),'_dmd_results_normalized_*.mat']);
load_results = load(fullfile(g_ica_result_files(end).folder,g_ica_result_files(end).name));

source_maps = load_results.source_maps;
inv_source = load_results.inv_source;
Phi_all = load_results.Phi_all;
Phi_orig = load_results.Phi_orig;
lambda = load_results.lambda;

abs_DM = abs(lambda);
period_DM = 2*pi*TRtarget ./ angle(lambda);

thres_period = 1e10;
Phi_orig_complex = Phi_orig(:,period_DM <= thres_period);
Phi_orig_real = Phi_orig(:,period_DM > thres_period);

%% source display

labels = cifti_read('atlas/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii');

lh_surface_file = 'atlas/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii';
rh_surface_file = 'atlas/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii';

save_dir = 'source_maps';
mkdir(save_dir);


source_maps_norm = zscore(source_maps');
thres = 2;

structure_label_list = {'hippocampus','striatum','thalamus','amygdala','cerebellum','brainstem'};
for n_source = 1:target_dim
    figure; fig = gcf;
%     source_maps_norm_thresholded = source_maps_norm .* (abs(source_maps_norm) > thres);

    max_abs_val = max(abs(source_maps(n_source,:)));
    
    snapshot_cifti = cifti_struct_create_from_template(labels,source_maps(n_source,:)', 'dscalar');
    display_cifti_cortex(fig, snapshot_cifti, lh_surface_file, rh_surface_file, [], -max_abs_val, max_abs_val, []);
    fname = fullfile(save_dir, sprintf('IC%02d_map_cortex.jpg',n_source));
    print(fig, fname, '-djpeg', '-r300');
    
    for i_struc = 1:length(structure_label_list)
        structure_label = structure_label_list{i_struc};
        plot_subcortex(fig, snapshot_cifti, [], [], -max_abs_val, max_abs_val, structure_label);
        fname = fullfile(save_dir, sprintf('IC%02d_map_%s.jpg', n_source, structure_label));
        print(fig, fname, '-djpeg', '-r300');
    end
    close(fig);
end

%% source labels
source_labels = cell(target_dim,2);
n = 1;
source_labels{1,1} = 'CEN'; source_labels{n,2} = 'control';
n = n + 1;
source_labels{n,1} = 'SM L'; source_labels{n,2} = 'SM';
n = n + 1;
source_labels{n,1} = '?'; source_labels{n,2} = '?';
n = n + 1;
source_labels{n,1} = 'cereb R'; source_labels{n,2} = 'cereb';
n = n + 1;
source_labels{n,1} = '?'; source_labels{n,2} = '?';
n = n + 1;
source_labels{n,1} = 'Auditory'; source_labels{n,2} = 'Auditory';
n = n + 1;
source_labels{n,1} = 'CEN'; source_labels{n,2} = 'control';
n = n + 1;
source_labels{n,1} = 'Vis'; source_labels{n,2} = 'Vis';
n = n + 1;
source_labels{n,1} = 'DMN L'; source_labels{n,2} = 'DMN';
n = n + 1;
source_labels{n,1} = 'CEN'; source_labels{n,2} = 'control';




%% load rsfMRI features
rsfMRI_features = load('results\seed_based_conn.mat');
seed_conn_PCC = rsfMRI_features.seed_conn_PCC;
seed_conn_aIns_R = rsfMRI_features.seed_conn_aIns_R;
seed_conn_SMG = rsfMRI_features.seed_conn_SMG;
seed_conn_FEF = rsfMRI_features.seed_conn_FEF;
seed_conn_TPJ = rsfMRI_features.seed_conn_TPJ;
rsfMRI_features = load('results\MLI.mat');
MLI_all = rsfMRI_features.MLI_all;
MLI = squeeze(mean(MLI_all,[1,2]));

%%
frame_dt = 0.5;
for pair_num = 1:4%size(Phi_orig_complex,2)/2-1
        
    DM_conjugate1_num = 2*(pair_num-1)+1;
    DM_conjugate2_num = 2*pair_num;
    
    lambda_conjugate1 = lambda(DM_conjugate1_num);
    lambda_conjugate2 = lambda(DM_conjugate2_num);
    
    if pair_num == 1
        ref_t = 0;
    elseif pair_num == 2
        ref_t = -3;
    elseif pair_num == 3
        ref_t = -1.5;
    elseif pair_num == 4
        ref_t = -6;
    else
        ref_t = 0;
    end
    
    frame_length = ceil(TRtarget * 2*pi / abs(angle(lambda_conjugate1))  / frame_dt);
    corr_DMN = zeros(1,frame_length);
    corr_SN = zeros(1,frame_length);
    corr_CEN = zeros(1,frame_length);
    corr_FEF = zeros(1,frame_length);
    corr_TPJ = zeros(1,frame_length);
    corr_MLI = zeros(1,frame_length);
    
    for frame = 1:frame_length
        current_time = frame_dt * (frame-1);
        
        source_snapshot = real( ...
            (lambda_conjugate1^((frame*frame_dt+ref_t)/TRtarget)) * Phi_orig_complex(:,DM_conjugate1_num) + ...
            (lambda_conjugate2^((frame*frame_dt+ref_t)/TRtarget)) * Phi_orig_complex(:,DM_conjugate2_num) ...
        );
    
        corr_DMN(frame) = corr(source_snapshot(1:cortex_dim),seed_conn_PCC);
        corr_SN(frame) = corr(source_snapshot(1:cortex_dim),seed_conn_aIns_R);
        corr_CEN(frame) = corr(source_snapshot(1:cortex_dim),seed_conn_SMG);
        corr_FEF(frame) = corr(source_snapshot(1:cortex_dim),seed_conn_FEF);
        corr_TPJ(frame) = corr(source_snapshot(1:cortex_dim),seed_conn_TPJ);
        corr_MLI(frame) = corr(source_snapshot(1:cortex_dim),MLI);
    end
    
    figure; hold on;
    plot(frame_dt*(1:frame_length),corr_DMN,'LineWidth',2);
    plot(frame_dt*(1:frame_length),corr_SN,'LineWidth',2);
    plot(frame_dt*(1:frame_length),corr_CEN,'LineWidth',2);
    plot(frame_dt*(1:frame_length),corr_FEF,'LineWidth',2);
    plot(frame_dt*(1:frame_length),corr_TPJ,'LineWidth',2);
    plot(frame_dt*(1:frame_length),corr_MLI,'LineWidth',2);
    yline(0.7,'r--');
    ylim([-1,1]);
    hold off;
    legend('DMN','SN','CEN','DAN','VAN','MLI');
    
end

