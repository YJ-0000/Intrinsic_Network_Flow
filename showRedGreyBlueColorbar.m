function showRedGreyBlueColorbar(minVal, maxVal, fontSize)
    % Custom Red → Grey → Blue colormap with 256 points
    n = 128; % number of points

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

    % Flip upside-down if needed
    color_map = flipud(cmap);

    % === Figure ===
    figure;
    colormap(color_map);

    % dummy image
    imagesc([0 1], [minVal maxVal], linspace(minVal, maxVal, size(color_map,1))');
    set(gca, 'Visible','off'); % 축 숨기기
    
    % 컬러바 생성
    cb = colorbar;
    cb.Ticks = [minVal maxVal];          
    cb.TickLabels = {num2str(minVal), num2str(maxVal)};
    cb.FontSize = fontSize;              

%     title('Custom Red → Grey → Blue Colorbar');
end
