clear;clc;

target_dim = 27;
TRtarget = 0.72;
g_ica_result_files = dir('results/INF_G_lev_REST1_ALL_027_MIGP_results_*.mat');

load_results = load(fullfile(g_ica_result_files(end).folder,g_ica_result_files(end).name));
lambda = load_results.lambda;
source_maps = load_results.source_maps;
Phi_all = load_results.Phi_all;
if isfield(load_results, 'Phi_orig')
    Phi_orig = -load_results.Phi_orig;
elseif isfield(load_results, 'Phi_orig_DR')
    Phi_orig = -load_results.Phi_orig_DR;
else
    Phi_orig = -load_results.Phi_orig_DL;
end
abs_DM = abs(lambda);
period_DM = 2*pi*TRtarget ./ angle(lambda);

thres_period = 1e10;
Phi_orig_complex = Phi_orig(:,period_DM <= thres_period);
Phi_orig_real = Phi_orig(:,period_DM > thres_period);
Phi_all_complex = Phi_all(:,period_DM <= thres_period);
Phi_all_real = Phi_all(:,period_DM > thres_period);
lambda_complex = lambda(period_DM < thres_period);

labels = cifti_read('atlas/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii');

mkdir('G_ICA_DMD_video_HCP_REST_fbDMD_REST1_ALL_DR_mode_whole');
save_dir = [pwd filesep 'G_ICA_DMD_video_HCP_REST_fbDMD_REST1_ALL_DR_mode_whole'];

lh_surface_file = 'atlas/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii';
rh_surface_file = 'atlas/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii';

%% Evolution of oscillating DMs
is_cortex_plot = true;
structure_label_list = {'cortex','hippocampus','amygdala','thalamus','striatum','brainstem','cerebellum_flat'};
structure_names_display = {'Cortex','Hippocampus','Amygdala','Thalamus','Striatum','Brainstem','Cerebellum'};

% idx_sort = [1,2,3,4,5,6,7,11,8,9,10,12,13,14];
idx_sort = 1:13;

for pair_num = 13%4:size(Phi_orig_complex,2)/2
    
    save_dir_dm = fullfile(save_dir,['DM_pair_',num2str(pair_num,'%02d')]);
    mkdir(save_dir_dm);
    
    if pair_num <= 5
        frame_dt = 0.5;
    elseif pair_num <=10
        frame_dt = 2;
    else 
        frame_dt = 4;
    end
    
    DM_conjugate1_num = 2*(idx_sort(pair_num)-1)+1;
    DM_conjugate2_num = 2*idx_sort(pair_num);
    
    lambda_conjugate1 = lambda_complex(DM_conjugate1_num);
    lambda_conjugate2 = lambda_complex(DM_conjugate2_num);
    
    frame_length = ceil(TRtarget * 2*pi / abs(angle(lambda_conjugate1))  / frame_dt);
    
    % if pair_num == 1
    %     ref_t = 0;
    % elseif pair_num == 2
    %     ref_t = 0;
    % elseif pair_num == 3
    %     ref_t = 0;
    % elseif pair_num == 4
    %     ref_t = 0;
    % else
    %     ref_t = 0;
    % end

    ref_t = 0;

    if pair_num == 1
        phase_sign = -1;
    elseif pair_num == 2
        phase_sign = -1;
    elseif  pair_num == 3
        phase_sign = 1;
    elseif  pair_num == 4
        phase_sign = exp(1i * pi / 4);
    else 
        phase_sign  = 1;
    end
    
    figure('Position', [0, 0, 554, 416]);
    fig = gcf;
    
    min_scale = -2*max(abs(Phi_orig_complex(:,DM_conjugate1_num)));
    max_scale = 2*max(abs(Phi_orig_complex(:,DM_conjugate1_num)));
    
    fprintf('Scale -- Max: %0.4f, Min: %0.4f \n',max_scale,min_scale);
    
    if is_cortex_plot
        for frame = 1:frame_length

            current_time = frame_dt * (frame-1);

            source_snapshot = real( ...
                (lambda_conjugate1^((frame*frame_dt+ref_t)/TRtarget)) * phase_sign * Phi_orig_complex(:,DM_conjugate1_num) + ...
                (lambda_conjugate2^((frame*frame_dt+ref_t)/TRtarget)) * phase_sign * Phi_orig_complex(:,DM_conjugate2_num) ...
            );

            activation_snapshot = source_snapshot;

            activation_snapshot_cifti = cifti_struct_create_from_template(labels,activation_snapshot, 'dscalar');

            % display_cifti_cortex를 이용하여 figure handle에 플롯 생성
            display_cifti_cortex(fig, activation_snapshot_cifti, lh_surface_file, rh_surface_file, [], min_scale, max_scale);

            % 좌상단에 시간 표시 (annotation 사용, figure 기준 좌표)
            annotation_handle = annotation('textbox', [0.01, 0.94, 0.2, 0.05], ...
                'String', sprintf('t = %.1f s', current_time), ...
                'EdgeColor', 'none', ...
                'FontSize', 12, ...
                'FontWeight', 'bold', ...
                'Color', 'k', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top');

            
            fname = fullfile(save_dir_dm, sprintf('DM%02d_map_%s_%04d.jpg', pair_num, 'cortex',frame));
            print(fig, fname, '-djpeg', '-r300');
        end
    end
    close(fig);
    
    for i_struc = 2:length(structure_label_list)
        figure; fig = gcf; set(fig, 'Position', [0, 0, 800, 800]);
        structure_label = structure_label_list{i_struc};
        
         for frame = 1:frame_length
            current_time = frame_dt * (frame-1);

            source_snapshot = real( ...
                (lambda_conjugate1^((frame*frame_dt+ref_t)/TRtarget)) * phase_sign * Phi_orig_complex(:,DM_conjugate1_num) + ...
                (lambda_conjugate2^((frame*frame_dt+ref_t)/TRtarget)) * phase_sign * Phi_orig_complex(:,DM_conjugate2_num) ...
            );

            activation_snapshot = source_snapshot;

            activation_snapshot_cifti = cifti_struct_create_from_template(labels,activation_snapshot, 'dscalar');
            
            if ~strcmp(structure_label,'cerebellum_flat')
                plot_subcortex(fig, activation_snapshot_cifti, [], [], min_scale, max_scale, structure_label);
            else
                plot_cifti_on_cereb_flatmap(activation_snapshot_cifti, 'FigHandle',fig,'CLim',[min_scale, max_scale]);
            end
            
            fname = fullfile(save_dir_dm, sprintf('DM%02d_map_%s_%04d.jpg', pair_num, structure_label,frame));
            print(fig, fname, '-djpeg', '-r300');
         end
         close(fig);
    end
    
    %%
    positions = [
        0.00, 0.3, 0.7, 0.7;  % Position for image 1
        0.00, 0.00, 0.25, 0.3;  % Position for image 2
        0.22, 0.00, 0.25, 0.3;  % Position for image 3
        0.45, 0.00, 0.25, 0.3;  % Position for image 4
        0.7, 0.7, 0.25, 0.3;  % Position for image 5
        0.7, 0.40, 0.25, 0.3;  % Position for image 6
        0.62, 0, 0.4, 0.4;  % Position for image 7
    ];
    title_positions = [
        0.51, 0.9;
        0.5, 0.9;
        0.5, 0.9;
        0.5, 0.9;
        0.5, 0.9;
        0.5, 0.9;
        0.5, 0.85;
    ];
    
    v = VideoWriter([save_dir, filesep, 'DM', num2str(pair_num,'%02d')], 'Motion JPEG AVI');
    v.FrameRate = 10; % Set frame rate
    open(v);
%%
    figure('Position', [100, 100, 900, 700]); fig = gcf;
    for frame = 1:frame_length
        clf(fig);
        % Loop through the seven images
        for n_roi = [7,1:6]
            structure_label = structure_label_list{n_roi};
            % Construct the filename for each image (assuming names as 'image1.jpg' to 'image7.jpg')
            filename = sprintf('DM%02d_map_%s_%04d.jpg', pair_num, structure_label,frame);

            % Load the image
            img = imread(fullfile(save_dir_dm,filename));

            % Create axes with specified position
            ax = axes('InnerPosition', positions(n_roi, :));

            % Display the image within the axes
            imshow(img, 'Parent', ax);

            % Optionally, turn off the axis for the image
            axis(ax, 'off');

            text(title_positions(n_roi,1), title_positions(n_roi,2), structure_names_display{n_roi}, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold', 'Parent', ax);
        end
        
        %%% colorbar
        img = imread(fullfile(save_dir,'colorbar.png'));

        % Create axes with specified position
        ax = axes('InnerPosition', [0.56, 0.4, 0.25, 0.25]);
        imshow(img, 'Parent', ax);
        axis(ax, 'off');
        text(4.5, 0.9, num2str(max_scale,'%.1f'), 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Parent', ax);
        text(4.5, 0.1, num2str(min_scale,'%.1f'), 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                 'VerticalAlignment', 'top', 'FontSize', 10, 'Parent', ax);
        
        
        set(gca, 'color', 'none'); set(gcf, 'color', 'w');
        writeVideo(v, getframe(fig));
    end
    close(v);
    close(fig);

end

%%

for pair_num = 1:size(Phi_orig_real,2)
    
    save_dir_dm = fullfile(save_dir,['DM_pair_',num2str(size(Phi_orig_complex,2)/2+pair_num,'%02d')]);
    mkdir(save_dir_dm);
    
    DM_conjugate1_num = pair_num;
    DM_conjugate2_num = pair_num;
    
    figure;
    fig = gcf;
    
    min_scale = -max(abs(Phi_orig_real(:,pair_num)));
    max_scale = max(abs(Phi_orig_real(:,pair_num)));
    
    fprintf('Scale -- Max: %0.4f, Min: %0.4f \n',max_scale,min_scale);
    
    frame_length = 1;
    
    if is_cortex_plot
        for frame = 1:frame_length

            current_time = frame_dt * (frame-1);

            activation_snapshot_cifti =  cifti_struct_create_from_template(labels,real(Phi_orig_real(:,pair_num)), 'dscalar');

            % display_cifti_cortex를 이용하여 figure handle에 플롯 생성
            display_cifti_cortex(fig, activation_snapshot_cifti, lh_surface_file, rh_surface_file, [], min_scale, max_scale);

            
            fname = fullfile(save_dir_dm, sprintf('DM%02d_map_%s_%04d.jpg', pair_num, 'cortex',frame));
            print(fig, fname, '-djpeg', '-r300');
        end
    end
    close(fig);
    
    for i_struc = 2:length(structure_label_list)
        figure; fig = gcf; set(fig, 'Position', [0, 0, 800, 800]);
        structure_label = structure_label_list{i_struc};
        
         for frame = 1:frame_length
             
            if ~strcmp(structure_label,'cerebellum_flat')
                plot_subcortex(fig, activation_snapshot_cifti, [], [], min_scale, max_scale, structure_label);
            else
                plot_cifti_on_cereb_flatmap(activation_snapshot_cifti, 'FigHandle',fig,'CLim',[min_scale, max_scale]);
            end
            
            fname = fullfile(save_dir_dm, sprintf('DM%02d_map_%s_%04d.jpg', pair_num, structure_label,frame));
            print(fig, fname, '-djpeg', '-r300');
         end
         close(fig);
    end
    
    %%
    positions = [
        0.00, 0.3, 0.7, 0.7;  % Position for image 1
        0.00, 0.00, 0.25, 0.3;  % Position for image 2
        0.22, 0.00, 0.25, 0.3;  % Position for image 3
        0.45, 0.00, 0.25, 0.3;  % Position for image 4
        0.7, 0.7, 0.25, 0.3;  % Position for image 5
        0.7, 0.40, 0.25, 0.3;  % Position for image 6
        0.62, 0, 0.4, 0.4;  % Position for image 7
    ];
    title_positions = [
        0.51, 0.9;
        0.5, 0.9;
        0.5, 0.9;
        0.5, 0.9;
        0.5, 0.9;
        0.5, 0.9;
        0.5, 0.85;
    ];
    
%%
    figure('Position', [100, 100, 900, 700]); fig = gcf;
    for frame = 1:frame_length
        clf(fig);
        % Loop through the seven images
        for n_roi = [7,1:6]
            structure_label = structure_label_list{n_roi};
            % Construct the filename for each image (assuming names as 'image1.jpg' to 'image7.jpg')
            filename = sprintf('DM%02d_map_%s_%04d.jpg', pair_num, structure_label,frame);

            % Load the image
            img = imread(fullfile(save_dir_dm,filename));

            % Create axes with specified position
            ax = axes('InnerPosition', positions(n_roi, :));

            % Display the image within the axes
            imshow(img, 'Parent', ax);

            % Optionally, turn off the axis for the image
            axis(ax, 'off');

            text(title_positions(n_roi,1), title_positions(n_roi,2), structure_names_display{n_roi}, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold', 'Parent', ax);
        end
        
        %%% colorbar
        img = imread(fullfile(save_dir,'colorbar.png'));

        % Create axes with specified position
        ax = axes('InnerPosition', [0.56, 0.4, 0.25, 0.25]);
        imshow(img, 'Parent', ax);
        axis(ax, 'off');
        text(2.8, 0.9, num2str(max_scale,'%.1f'), 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Parent', ax);
        text(2.8, 0.1, num2str(min_scale,'%.1f'), 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                 'VerticalAlignment', 'top', 'FontSize', 10, 'Parent', ax);
        
        
        set(gca, 'color', 'none'); set(gcf, 'color', 'w');
        print(fig, fullfile(save_dir,sprintf('DM%02d_map.jpg', size(Phi_orig_complex,2)/2+pair_num)), '-djpeg', '-r300');
    end
    close(fig);

end



