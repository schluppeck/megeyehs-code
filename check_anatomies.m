%% quick minimal code example for loading in MRI data
%  for megeyehs project
%
%  ds 2022-06-23

%% make sure fieldtrip code is added to path
%
% check that ft_read_mri is a function (if yes, assume ft is on path)



if exist('ft_read_mri') ~= 2
    disp('(uhoh) need to add fieldtrip to path')
    cwd = pwd();
    addpath(genpath(fullfile(cwd,'fieldtrip')));
    % check again
    assert( exist('ft_read_mri') == 2, 'still no good - whatsup?')
end

%% load a filename
%
% usign relative pathnames -- so this works on OneDrive folder which is
% mapped to different locations depmnding on machine... use your own
% convention to make this work with local data

% NB: the surfRelax folder is derived from freesurfer segmentation
% using a bunch of tools (mrTools) for fMRI analysis... denis uses this a
% lot, so if any qustions come up, check in...

participant = '15910';
datafolder = fullfile('..', participant, 'surfRelax');

M = ft_read_mri( fullfile(datafolder, sprintf('%s_mprage_pp.hdr', participant)));

%% and display - many ways to skin a cat 

robustRange = prctile(M.anatomy(:), [5 95]);
whichSlices = 1:5:size(M.anatomy,3);
imDims = size(M.anatomy);
montage(M.anatomy, 'DisplayRange', robustRange, 'Indices', whichSlices)


%% or a particular set of slices
figure()
midStack = round(imDims ./2);
s_ = slice(M.anatomy, midStack(2), midStack(1), midStack(3));

% make it look nice
shading('interp')
camproj('perspective')
caxis(robustRange)
colormap(gray())
% camlight()
axis('vis3d') % keep aspect ratios
axis('off')
rotate3d() % switch on interaction

%% can also try meshes... using Jonas Larsson's code
% which is wrapped into mrTools distro

f_skull_ = figure()
S = loadSurfOFF( fullfile(datafolder, sprintf('%s_left_WM.off', participant)));
curv = loadVFF( fullfile(datafolder, sprintf('%s__left_Curv.vff', participant)));

patch('vertices', S.vtcs, 'faces', S.tris,'facecolor', [0.8 0.2 0.2],'edgecolor','none');
colormap(gray);
light()
camlight()
daspect([1 1 1]);
axis('vis3d') % keep aspect ratios
axis('off')
rotate3d() % switch on interaction


%% adding a transparent skull surface to this

SKi = loadSurfOFF( fullfile(datafolder, sprintf('%s_outer_skin_surface.off', participant)));
SKo = loadSurfOFF( fullfile(datafolder, sprintf('%s_outer_skull_surface.off', participant)));

skin_patch_ = patch('vertices', SKi.vtcs, 'faces', SKi.tris,'facecolor', [0.1 0.3 0.3],'edgecolor','none');
skull_patch_ = patch('vertices', SKo.vtcs, 'faces', SKo.tris,'facecolor', [0.1 0.1 0.3],'edgecolor','none');
alpha(skull_patch_, 0.2)
alpha(skin_patch_, 0.4)
