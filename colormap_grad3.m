function cmap = colormap_grad3(n)

    if nargin < 1
        n = 256; % default colormap length
    end
    mid_point = n/2+1;

    % Key colors (RGB in [0,1])
    keyColors = [...
        1.00 0.00 0.00;   % red
        0.50 0.50 0.50;   % 
        0.10 0.10 1.00;   % blue
        % 0.50 0.50 0.50;   % 
        0.10 1.00 0.10;   % green
        0.50 0.50 0.50;   % 
        1.00 0.00 0.00];  % red

    % number of key colors
    m = size(keyColors,1);

    % interpolation positions
    x = linspace(0,1,m);
    xi = linspace(0,1,n);

    % interpolate each channel separately
    cmap = zeros(n,3);
    for k = 1:3
        cmap(:,k) = interp1(x, keyColors(:,k), xi, 'pchip');
    end
    % cmap(mid_point,:) = 0.5;
end
