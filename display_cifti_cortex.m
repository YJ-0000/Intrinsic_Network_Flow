function display_cifti_cortex(fig_handle, ciftiData, lh_surface_file, rh_surface_file, output_plot_path, min_scale, max_scale, color_map, interp)
% display_cifti_cortex plots BOLD activation 데이터를 HCP cortex surface에 overlay하여 시각화합니다.
%
%   display_cifti_cortex(fig_handle, ciftiData, lh_surface_file, rh_surface_file, output_plot_path, min_scale, max_scale, color_map)
%
%   입력 인자:
%       fig_handle       - figure handle (빈 값이면 새 figure 생성)
%       ciftiData        - BOLD activation 값이 저장된 cifti 파일 경로 또는 cifti 구조체
%       lh_surface_file  - 왼쪽 hemisphere의 surface를 담은 gifti 파일 경로
%       rh_surface_file  - 오른쪽 hemisphere의 surface를 담은 gifti 파일 경로
%       output_plot_path - (옵션) 저장할 플롯 파일 경로 (비어있으면 저장하지 않음)
%       min_scale        - activation 값의 최소 스케일
%       max_scale        - activation 값의 최대 스케일
%       color_map        - (옵션) 사용자 지정 colormap. 입력하지 않으면 기본 colormap 생성.
%
% 예시:
%   fig = figure;
%   display_cifti_cortex(fig, 'myCIFTI.dtseries.nii', 'lh.surface.gii', 'rh.surface.gii', 'output.png', -2, 2, []);

    %% 0. Figure handle 처리
    if nargin < 1 || isempty(fig_handle) || ~ishandle(fig_handle)
        fig_handle = figure;
    else
        figure(fig_handle);
        clf(fig_handle);  % 기존 내용을 지웁니다.
    end

    %% 1. CIFTI 파일 읽기 (cifti_read 사용)
    if ischar(ciftiData) || isstring(ciftiData)
        if ~isfile(ciftiData)
            error('The specified file does not exist: %s', ciftiData);
        end
        try
            ciftiData = cifti_read(ciftiData);
        catch ME
            error('Error reading CIFTI file: %s', ME.message);
        end
    end
    % ciftiData 구조체의 cdata 필드에 BOLD activation 값이 있다고 가정
    data = ciftiData.cdata;
    
    %% 2. Colormap 생성 (입력 인자가 없으면 기본 colormap 생성)
    if nargin < 8 || isempty(color_map)
        n = 128; % 각 구간의 포인트 수 (전체 256)
        % Red to Grey Transition
        R_red_to_grey = linspace(1, 0.5, n);
        G_red_to_grey = linspace(0, 0.5, n);
        B_red_to_grey = linspace(0, 0.5, n);
        % Grey to Blue Transition
        R_grey_to_blue = linspace(0.5, 0, n);
        G_grey_to_blue = linspace(0.5, 0, n);
        B_grey_to_blue = linspace(0.5, 1, n);
        % Concatenate
        R = [R_red_to_grey, R_grey_to_blue];
        G = [G_red_to_grey, G_grey_to_blue];
        B = [B_red_to_grey, B_grey_to_blue];
        cmap = [R', G', B'];
        color_map = flipud(cmap);
    end

    if nargin < 9
        interp = 'interp';
    end

    %% 3. Cortex 부분 추출
    % 좌/우 모델의 정보 추출 (예: HCP cifti 구조체의 diminfo 이용)
    lh_indices = 1:ciftiData.diminfo{1, 1}.models{1, 1}.count;
    rh_indices = ciftiData.diminfo{1, 1}.models{1, 2}.start:(...
                  ciftiData.diminfo{1, 1}.models{1, 2}.start + ciftiData.diminfo{1, 1}.models{1, 2}.count - 1);
    
    lh_vertlist = ciftiData.diminfo{1, 1}.models{1, 1}.vertlist + 1;
    rh_vertlist = ciftiData.diminfo{1, 1}.models{1, 2}.vertlist + 1;
    
    lh_activation = zeros(ciftiData.diminfo{1, 1}.models{1, 1}.numvert, 1);
    lh_activation(lh_vertlist) = data(lh_indices);
    
    rh_activation = zeros(ciftiData.diminfo{1, 1}.models{1, 2}.numvert, 1);
    rh_activation(rh_vertlist) = data(rh_indices);
    
    plot_max_value = max([max(abs(data(lh_indices))),max(abs(data(rh_indices)))]);
    plot_min_value = -plot_max_value;

    %% 4. gifti 파일을 이용하여 surface 정보 읽기
    lh_surface = gifti(lh_surface_file);
    rh_surface = gifti(rh_surface_file);

    %% 5. 네 가지 view로 plot 그리기
    % figure handle에 plot을 그립니다.
    figure(fig_handle);
    
    zoomFactor = 3;  % cortex 확대 배율

    % 5-1. 왼쪽 hemisphere, 왼쪽 view
    subplot(2,2,1, 'Parent', fig_handle);
    patch('Vertices', lh_surface.vertices, 'Faces', lh_surface.faces, ...
          'FaceVertexCData', lh_activation, 'FaceColor', interp, 'EdgeColor', 'none');
    axis equal; axis off;
    view([-90, 0]);  % 왼쪽 view
    camzoom(zoomFactor);
    camlight('headlight');
    lighting gouraud;
    material dull;
    
    % 5-2. 왼쪽 hemisphere, 오른쪽 view
    subplot(2,2,3, 'Parent', fig_handle);
    patch('Vertices', lh_surface.vertices, 'Faces', lh_surface.faces, ...
          'FaceVertexCData', lh_activation, 'FaceColor', interp, 'EdgeColor', 'none');
    axis equal; axis off;
    view([90, 0]);   % 오른쪽 view
    camzoom(zoomFactor);
    camlight('headlight');
    lighting gouraud;
    material dull;
    
    % 5-3. 오른쪽 hemisphere, 왼쪽 view
    subplot(2,2,4, 'Parent', fig_handle);
    patch('Vertices', rh_surface.vertices, 'Faces', rh_surface.faces, ...
          'FaceVertexCData', rh_activation, 'FaceColor', interp, 'EdgeColor', 'none');
    axis equal; axis off;
    view([-90, 0]);  % 왼쪽 view
    camzoom(zoomFactor);
    camlight('headlight');
    lighting gouraud;
    material dull;
    
    % 5-4. 오른쪽 hemisphere, 오른쪽 view
    subplot(2,2,2, 'Parent', fig_handle);
    patch('Vertices', rh_surface.vertices, 'Faces', rh_surface.faces, ...
          'FaceVertexCData', rh_activation, 'FaceColor', interp, 'EdgeColor', 'none');
    axis equal; axis off;
    view([90, 0]);   % 오른쪽 view
    camzoom(zoomFactor);
    camlight('headlight');
    lighting gouraud;
    material dull;
    
    % 모든 subplot에 대해 색상 스케일 설정
    if isempty(min_scale) || isempty(max_scale)
        set(findall(fig_handle, 'Type', 'axes'), 'CLim', [plot_min_value, plot_max_value]);
    else
        set(findall(fig_handle, 'Type', 'axes'), 'CLim', [min_scale, max_scale]);
    end
    
    % 입력받은 colormap 적용
    colormap(color_map);
    
    % output_plot_path가 비어있지 않으면 해당 경로에 저장
    if ~isempty(output_plot_path)
        print(fig_handle, output_plot_path, '-djpeg', '-r300');
    end
end
