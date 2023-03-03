%% ery_4a_masks_m6_s0_binary_mask_from_atlas.m
%
% This script creates a binary mask by combining regions from one or more
% atlases, and automatically saves the fmri_mask_image, and .nii versions
% of it in maskdir of your dataset
% 
% USAGE
%
% Script should be run from the root directory of the superdataset, e.g.
% /data/proj_discoverie
%
% DEPENDENCIES
%
% CANlab's CanlabCore and Neuroimaging_Pattern_Masks Github repos on your Matlab path
% if needed, clone from https://github.com/canlab
%
%__________________________________________________________________________
%
% authors: Aleksandra Budzinska, Lukas Van Oudenhove
% date:   KU Leuven, July, 2022
%
%__________________________________________________________________________
% @(#)% LaBGAScore_atlas_binary_mask_from_atlas.m         v1.2       
% last modified: 2022/07/26


% Define maskname and directory where mask will be written
%--------------------------------------------------------------------------

maskname = 'ery_4a_m6_mask_taste_cortex';
maskdir = '/data/proj_erythritol/proj_erythritol_4a/secondlevel/model_6m_long_12hmp_can_ar1/masks';

% NOTE: in default LaBGAS file organization this is a model-specific dir in
% rootdir/secondlevel/model_xxx, but you can obviously specify any dir here


% Load atlas of your choice
%--------------------------------------------------------------------------
canlab = load_atlas('canlab2018_2mm'); % you get a CANlab atlas object

% NOTE: you can not only load CANlab atlases using keywords as in this example,
% but any atlas you like by creating a new atlas object from a .nii atlas image, 
% e.g. aicha = atlas('/opt/KUL_apps/spm12/atlas/AICHA.nii');
% see help atlas for adding labels etc


% Select regions you want to include from the selected atlas 
%--------------------------------------------------------------------------
canlab_subset = select_atlas_subset(canlab, {'Ctx_AVI_L','Ctx_AVI_R','Ctx_AAIC_L','Ctx_AAIC_R','Ctx_MI_L','Ctx_MI_R','Ctx_FOP3_L','Ctx_FOP3_R','Ctx_FOP2_L','Ctx_FOP2_R'});
% NOTE: you can also use the numbers rather than the labels of regions to
% select as a vector
% see help atlas.select_atlas_subset


% Load and select the regions from another atlas (e.g. subcortical) like above
%--------------------------------------------------------------------------
%subcortical = load_atlas('cit168'); % you get a CANlab atlas object
%subcortical_subset = select_atlas_subset(subcortical, {'Put','Cau','NAC','BST_SLEA','GPe','GPi','SNc','RN','SNr','PBP','VTA','VeP','Hythal','STN'}); % still a CANlab atlas object


% Change data type of probability maps if needed
%--------------------------------------------------------------------------

% NOTE: Before merging the atlases, make sure that they have the same type of variables (look at the probability_maps. E.g. single/sparse double)
% This command converts probability_maps property variable type from single to double, and then to sparse double - may not be needed for all atlases

%subcortical_subset.probability_maps = double(subcortical_subset.probability_maps);
%subcortical_subset.probability_maps = sparse(subcortical_subset.probability_maps);


% Merge the atlases
%--------------------------------------------------------------------------

%combined_atlas = merge_atlases(subcortical_subset, canlab_subset); 

% NOTE: merge_atlases function automatically resamples the space of the second to the first, and adds consecutive labels


% Convert the atlas object to a CANlab fmri_mask_image object
%--------------------------------------------------------------------------

mask = fmri_mask_image(canlab_subset);

% NOTE: this automatically binarizes the mask!


% Save your binarized mask in maskdir
%--------------------------------------------------------------------------

write(mask, 'fname', fullfile(maskdir,[maskname,'.nii'])); % writes to .nii file

savefilenamedata = fullfile(maskdir,[maskname,'.mat']); % saves fmri_mask_image object as .mat file
save(savefilenamedata, 'mask', '-v7.3');
fprintf('\nSaved mask\n');

