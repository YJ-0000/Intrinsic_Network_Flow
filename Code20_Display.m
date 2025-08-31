clear;clc;

target_dim = 27;
TRtarget = 0.72;
g_ica_result_files = dir(['results/loop_gica_',num2str(target_dim,'%03d'),'_dmd_results_normalized_*.mat']);
load_results = load(fullfile(g_ica_result_files(end).folder,g_ica_result_files(end).name));
lambda = load_results.lambda;
source_maps = load_results.source_maps;
Phi_all = load_results.Phi_all;
Phi_orig = load_results.Phi_orig;

abs_DM = abs(lambda);
period_DM = 2*pi*TRtarget ./ angle(lambda);

thres_period = 1e10;
Phi_orig_complex = Phi_orig(:,period_DM <= thres_period);
Phi_orig_real = Phi_orig(:,period_DM > thres_period);
Phi_all_complex = Phi_all(:,period_DM <= thres_period);
Phi_all_real = Phi_all(:,period_DM > thres_period);

labels = cifti_read('atlas/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii');

mkdir('G_ICA_DMD_video_HCP_REST_fbDMD_REST1_exact_mode');
save_dir = [pwd filesep 'G_ICA_DMD_video_HCP_REST_fbDMD_REST1_exact_mode'];

lh_surface_file = 'atlas/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii';
rh_surface_file = 'atlas/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii';

%% Evolution of oscillating DMs
for pair_num = 1:10%size(Phi_orig_complex,2)/2-1
    
    if pair_num <= 5
        frame_dt = 0.5;
    elseif pair_num <= 10  
        frame_dt = 1; 
    else
        frame_dt = 2;
    end
    
    DM_conjugate1_num = 2*(pair_num-1)+1;
    DM_conjugate2_num = 2*pair_num;
    
    lambda_conjugate1 = lambda(DM_conjugate1_num);
    lambda_conjugate2 = lambda(DM_conjugate2_num);
    
    
    ref_t = 0;
    
    % 1. 비디오 작성 객체 생성 (MPEG-4 형식)
    v = VideoWriter([save_dir, filesep, 'DM_pair', num2str(pair_num),'_cortex'], 'MPEG-4');
    v.FrameRate = 10; % Set frame rate
    open(v);
    
    figure;
    fig = gcf;
    
    min_scale = -2*max(abs(Phi_orig_complex(:,DM_conjugate1_num)));
    max_scale = 2*max(abs(Phi_orig_complex(:,DM_conjugate1_num)));
    
    frame_length = ceil(TRtarget * 2*pi / abs(angle(lambda_conjugate1))  / frame_dt);
    
    for frame = 1:frame_length
        
        current_time = frame_dt * (frame-1);
        
        source_snapshot = real( ...
            (lambda_conjugate1^((frame*frame_dt+ref_t)/TRtarget)) * Phi_orig_complex(:,DM_conjugate1_num) + ...
            (lambda_conjugate2^((frame*frame_dt+ref_t)/TRtarget)) * Phi_orig_complex(:,DM_conjugate2_num) ...
        );
    
        activation_snapshot = source_snapshot;
        
        activation_snapshot_cifti = cifti_struct_create_from_template(labels,activation_snapshot, 'dscalar');
        
        % display_cifti_cortex를 이용하여 figure handle에 플롯 생성
        display_cifti_cortex(fig, activation_snapshot_cifti, lh_surface_file, rh_surface_file, [], min_scale, max_scale);
        
        % 좌상단에 시간 표시 (annotation 사용, figure 기준 좌표)
        annotation_handle = annotation('textbox', [0.01, 0.94, 0.2, 0.05], ...
            'String', sprintf('t = %.1f s', current_time), ...
            'EdgeColor', 'none', ...
            'FontSize', 14, ...
            'FontWeight', 'bold', ...
            'Color', 'k', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top');
        
%         plot_subcortex(cifti_data, [], [], min_scale, max_scale, 'thalamus')
        
        % 프레임 캡처 후 동영상에 기록
        writeVideo(v, getframe(fig));
    end
    
    close(v);
    close(fig);
    
    figure; fig = gcf; set(fig, 'Position', [0, 0, 800, 800]);
    structure_label_list = {'cerebellum_flat','hippocampus','striatum','thalamus','amygdala','cerebellum','brainstem'};
    
    for i_struc = 1:length(structure_label_list)
        structure_label = structure_label_list{i_struc};
        
        % 1. 비디오 작성 객체 생성 (MPEG-4 형식)
        v = VideoWriter([save_dir, filesep, 'DM_pair', num2str(pair_num),'_',structure_label], 'MPEG-4');
        v.FrameRate = 10; % Set frame rate
        open(v);        
        
         for frame = 1:frame_length
            current_time = frame_dt * (frame-1);

            source_snapshot = real( ...
                (lambda_conjugate1^((frame*frame_dt+ref_t)/TRtarget)) * Phi_orig_complex(:,DM_conjugate1_num) + ...
                (lambda_conjugate2^((frame*frame_dt+ref_t)/TRtarget)) * Phi_orig_complex(:,DM_conjugate2_num) ...
            );

            activation_snapshot = source_snapshot;

            activation_snapshot_cifti = cifti_struct_create_from_template(labels,activation_snapshot, 'dscalar');
            if ~strcmp(structure_label,'cerebellum_flat')
                plot_subcortex(fig, activation_snapshot_cifti, [], [], min_scale, max_scale, structure_label);
            else
                plot_cifti_on_cereb_flatmap(activation_snapshot_cifti, 'FigHandle',fig,'CLim',[min_scale, max_scale]);
            end
            % 프레임 캡처 후 동영상에 기록
            writeVideo(v, getframe(fig));
         end
        close(v);
    end
    
    close(fig);
end

%% Evolution of stationary DMs
structure_label_list = {'hippocampus','striatum','thalamus','amygdala','cerebellum','brainstem'};

for dm_num = 1:size(Phi_orig_real,2)
    
    figure; fig = gcf;
    activation_snapshot_cifti = cifti_struct_create_from_template(labels,real(Phi_orig_real(:,dm_num)), 'dscalar');
    
    min_scale = -2*max(abs(Phi_orig(:,DM_conjugate1_num)));
    max_scale = 2*max(abs(Phi_orig(:,DM_conjugate1_num)));
    
    display_cifti_cortex(fig, activation_snapshot_cifti, lh_surface_file, rh_surface_file, [], min_scale, max_scale, []);
    
    fname = fullfile(save_dir, sprintf('Statinary_DM_Map_%03d_Cortex.jpg', dm_num));
    print(fig, fname, '-djpeg', '-r300');
    
    for i_struc = 1:length(structure_label_list)
        structure_label = structure_label_list{i_struc};
        plot_subcortex(fig, activation_snapshot_cifti, [], [], min_scale, max_scale, structure_label);
        fname = fullfile(save_dir, sprintf('Statinary_DM_Map_%03d_%s.jpg', dm_num, structure_label));
        print(fig, fname, '-djpeg', '-r300');
    end
end

%% Evolution in source space

mkdir('G_ICA_DMD_video_HCP_REST_fbDMD_REST1_source_space');
save_dir = [pwd filesep 'G_ICA_DMD_video_HCP_REST_fbDMD_REST1_source_space'];

norm_source_maps = sqrt(mean(source_maps.^2,2));
norm_coeff = 1./norm_source_maps;

for pair_num = 1:6%size(Phi_all_complex,2)/2-1
    
    if pair_num <= 5
        frame_dt = 0.5;
    else
        frame_dt = 1;
    end
    
    DM_conjugate1_num = 2*(pair_num-1)+1;
    DM_conjugate2_num = 2*pair_num;
    
    lambda_conjugate1 = lambda(DM_conjugate1_num);
    lambda_conjugate2 = lambda(DM_conjugate2_num);
    
    
    ref_t = 0;
    
    % 1. 비디오 작성 객체 생성 (MPEG-4 형식)
    v = VideoWriter([save_dir, filesep, 'DM_pair', num2str(pair_num)], 'MPEG-4');
    v.FrameRate = 10; % Set frame rate
    open(v);
    
    figure;
    fig = gcf;
    
    min_scale = -2*max(abs(diag(norm_coeff)\Phi_all_complex(:,DM_conjugate1_num)));
    max_scale = 2*max(abs(diag(norm_coeff)\Phi_all_complex(:,DM_conjugate1_num)));
    
    frame_length = ceil(TRtarget * 2*pi / abs(angle(lambda_conjugate1))  / frame_dt);
    
    for frame = 1:frame_length
        
        current_time = frame_dt * (frame-1);
        
        source_snapshot = real( ...
            (lambda_conjugate1^((frame*frame_dt+ref_t)/TRtarget)) * Phi_all_complex(:,DM_conjugate1_num) + ...
            (lambda_conjugate2^((frame*frame_dt+ref_t)/TRtarget)) * Phi_all_complex(:,DM_conjugate2_num) ...
        );
        
        figure(fig); clf(fig);
        bar(diag(norm_coeff)\source_snapshot,'k');
        ylim([min_scale,max_scale]);
        
        % 좌상단에 시간 표시 (annotation 사용, figure 기준 좌표)
        annotation_handle = annotation('textbox', [0.01, 0.94, 0.2, 0.05], ...
            'String', sprintf('t = %.1f s', current_time), ...
            'EdgeColor', 'none', ...
            'FontSize', 14, ...
            'FontWeight', 'bold', ...
            'Color', 'k', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top');
        
%         plot_subcortex(cifti_data, [], [], min_scale, max_scale, 'thalamus')
        
        % 프레임 캡처 후 동영상에 기록
        writeVideo(v, getframe(fig));
    end
    
    close(v);
    close(fig);
    
end

