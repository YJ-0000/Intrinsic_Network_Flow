function plot_cifti_on_cereb_flatmap(cifti, varargin)
% Plot a CIFTI dscalar onto the SUIT cerebellum flatmap (wb_command 불필요).
% - dscalar: 단일 프레임(3D)
%
% Name-Value (추가/변경된 것만 요약):
%   'Space'        : 'MNI' (default) or 'SUIT'
%   'ForceWarp'    : true/false (default: false)  % true면 suit_mni2suit로 항상 워핑 후 SUIT로 매핑
%
% 나머지 옵션: OutPng, Structure, Cropped, Dimension, CLim, Colormap, ShowColorbar, FigHandle, MapArgs

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
RB_color_map = flipud(cmap);

p = inputParser;
p.addParameter('OutPng','',@(s)ischar(s)||isstring(s));
p.addParameter('Structure','CEREBELLUM',@(s)ischar(s)||isstring(s));
p.addParameter('CLim',[],@(x)isnumeric(x)&&numel(x)==2 || isempty(x));
p.addParameter('Colormap',RB_color_map,@(x)true);
p.addParameter('ShowColorbar',true,@(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('FigHandle',[],@(h)isempty(h)||isgraphics(h));
p.addParameter('MapArgs',{'space','SUIT'},@(x)iscell(x));  % 기본값은 아래 Space 로직에서 덮어씀
p.addParameter('Space','MNI',@(s)any(strcmpi(s,{'MNI','SUIT'})));
p.addParameter('ForceWarp',false,@(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('CombineMode','max',@(s)any(strcmpi(s,{'max','sum','mean'})));
p.parse(varargin{:});
opt = p.Results;

assert(~isempty(which('suit_plotflatmap')), 'SUIT: suit_plotflatmap not found.');
assert(~isempty(which('suit_mni2suit')),    'SUIT: suit_mni2suit not found.');
assert(~isempty(which('suit_map2surf')),   'SUIT: suit_map2surf not found.');
assert(~isempty(which('spm_vol')),         'SPM12 not found.');

if isempty(opt.Colormap); opt.Colormap = RB_color_map; end

try
    [vol3d, sform, roiMask] = cifti_struct_dense_extract_volume_structure_data(cifti, opt.Structure);
catch
    nL = [opt.Structure,'_LEFT']; nR = [opt.Structure,'_RIGHT'];
    [vL, sL, rL] = cifti_struct_dense_extract_volume_structure_data(cifti, nL);
    [vR, sR, rR] = cifti_struct_dense_extract_volume_structure_data(cifti, nR);
    
    assert(all(all(sL==sR)));
    
    vol3d = vL + vR;
    sform = sL;
    roiMask = rL | rR;
end

% 2) 임시 NIfTI 저장
cwd = pwd;
tmpdir = fullfile(cwd,'work_dir'); 
if ~exist(tmpdir, 'dir')
    mkdir(tmpdir);
end
cd(tmpdir);
cerebNii = 'cereb_from_cifti.nii';
write_3d_nii(vol3d, sform, cerebNii);
suit_mni2suit(cerebNii,'Interp',1)
cd(cwd);

cerebNii_suit = fullfile(tmpdir,['Wm2s_',cerebNii]);
Data = suit_map2surf(cerebNii_suit, 'space','SUIT');

args = {};
if ~isempty(opt.Colormap), args = [args, {'cmap', opt.Colormap}]; end
if ~isempty(opt.CLim),     args = [args, {'cscale', opt.CLim}];  end
suit_plotflatmap(Data, args{:});
set(gca, 'xtick', [], 'ytick', [], 'ztick', [], 'xcolor', 'w', 'ycolor', 'w', 'zcolor', 'w', 'Visible', 0);

end

% --------- helpers ---------
function write_3d_nii(vol3d, sform1, outname)
    V = struct();
    V.fname = outname;
    V.dim   = size(vol3d); if numel(V.dim)==2, V.dim(3)=1; end
    V.dt    = [16 0]; % float32
    V.mat   = sform1;
    V.pinfo = [1;0;0];
    V.descrip = 'CIFTI->NIfTI (cerebellum)';
    Vn = spm_create_vol(V);
    spm_write_vol(Vn, single(vol3d));
end
