function cmap = copper_s(n)
    % COPPER_S   Custom cyclic copper-style colormap with white ends
    %
    %   cmap = COPPER_S(n) returns an n-by-3 colormap matrix that
    %   fades white → copper shades → black → teal shades → white.
    %
    %   Example:
    %       imagesc(peaks(300));
    %       colormap(copper_s(256));
    %       colorbar;

    if nargin < 1
        n = 256; % default colormap length
    end

    % Key colors (RGB in [0,1])
    keyColors = [...
        0.95 0.90 0.70;   % 베이지
        0.80 0.45 0.35;   % 구리빛
        0.45 0.25 0.40;   % 보라톤
        0.40 0.40 0.40;   % 검정 (중앙)
        0.15 0.35 0.50;   % 청록빛
        0.80 0.95 0.85;   % 연한 민트
        0.95 0.90 0.70];

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
end
