%%% ery_4a_secondlevel_m6st_s7_e_run_rsa_st
%
%
% USAGE
%
% This script runs representational similarity analysis on an fmri_dat_st
% object created by the prep_3f script, calling Phil Kragel's
% rsa_regression() function under the hood
%
% For more info and an example, see
%
% https://github.com/labgas/CANlab_help_examples/tree/master/RSA_example
%
%
% OPTIONS FROM PREP_3f
%
% see first section of this script below and a2_set_default_options
%
%
% OPTIONS FOR CURRENT SCRIPT
%
% atlas = load_atlas('');            % keyword from load_atlas                       default 'canlab2018'
% regions = {};                      % vector of indices or cell array of labels for (contiguous!) regions you want to extract from atlas
%
%
%__________________________________________________________________________
%
% created by: Lukas Van Oudenhove
% date:   Leuven, January, 2023
%
%__________________________________________________________________________
% @(#)% e_run_rsa_st.m         v1.0
% last modified: 2023/01/25


%% GET AND SET OPTIONS
% -------------------------------------------------------------------------

% GET MODEL-SPECIFIC PATHS AND OPTIONS

ery_4a_secondlevel_m6st_s0_a_set_up_paths_always_run_first;

% NOTE: CHANGE THIS TO THE MODEL-SPECIFIC VERSION OF THIS SCRIPT
% NOTE: THIS WILL ALSO AUTOMATICALLY CALL A2_SET_DEFAULT_OPTIONS

% COPY MANDATORY OPTION FROM CORRESPONDING PREP_3f_ SCRIPT IF

results_suffix = ''; % suffix of your choice added to .mat file with saved results

% COPY OPTIONS FROM CORRESPONDING PREP_3f_ SCRIPT

% cons2exclude_dat_st = {}; % cell array of condition names to exclude, separated by commas (or blanks)
% behav_outcome_dat_st = 'rating'; % name of outcome variable in DAT.BEHAVIOR.behavioral_data_table_st
% subj_identifier_dat_st = 'participant_id'; % name of subject identifier variable in same table
% group_identifier_dat_st = 'group'; % name of group identifier variable in same table; leave commented out if you don't have groups

% SET CUSTOM OPTIONS FOR CURRENT SCRIPT

atlas_name = 'insula';                                  % keyword from load_atlas                       default 'canlab2018'
atlas = load_atlas(atlas_name);                     
regions = {atlas.labels{3}, atlas.labels{4}, atlas.labels{5}, atlas.labels{6}, atlas.labels{9}, atlas.labels{10}, atlas.labels{11}, atlas.labels{12}};               % vector of indices or cell array of labels for (contiguous!) regions you want to extract from atlas

% NOTE: the latter two option categories only need to be specified if you want to run a second version of your model 
% with different options than the defaults you set in your model-specific version of a2_set_default_options.m


%% LOAD FMRI_DATA_ST OBJECT AND OTHER NECESSARY VARIABLES IF NEEDED
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('LOADING DATA');
fprintf('\n\n');

if ~exist('DSGN','var') || ~exist('DAT','var')

    load(fullfile(resultsdir,'image_names_and_setup.mat'));

    if ~isfield(DAT,'BEHAVIOR')
        fprintf('\n');
        error('Behavioral data not yet added to DAT structure - run prep_1b script first');
    end

end

if ~exist('fmri_dat','var')
    
    if ~isempty(cons2exclude_dat_st)
        load(fullfile(resultsdir, ['single_trial_fmri_data_st_object_', behav_outcome_dat_st, '_exclude_cond_', char([cons2exclude_dat_st{:}]), '_', results_suffix, '.mat']));
        
    else
        load(fullfile(resultsdir, ['single_trial_fmri_data_st_object_', behav_outcome_dat_st, '_', results_suffix, '.mat']));

    end
end


%% DEFINE SUBJECT AND CONDITION IDENTIFIERS
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('GETTING SUBJECT AND CONDITION IDENTIFIERS');
fprintf('\n\n');

subject_id = fmri_dat.metadata_table.(subj_identifier_dat_st);
[uniq_subject_id, ~, subject_id] = unique(subject_id,'stable');
n_subj = size(uniq_subject_id,1);

cond_id = fmri_dat.metadata_table.trial_type;
[uniq_cond_id, ~, cond_id] = unique(cond_id,'stable');
n_cond = size(uniq_cond_id,1);


%% SCALE IMAGES ACCORDING TO OPTIONS
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('SCALING IMAGES IF REQUESTED IN OPTIONS');
fprintf('\n\n');

    switch myscaling_mvpa_reg_st

        case 'raw'

            fprintf('\nNo scaling of input images\n\n');

        case 'centerimages'

            fmri_dat = fmri_dat.rescale('centerimages'); 
            fprintf('\nCentering input images\n\n');

        case 'l2normimages'

            fmri_dat = fmri_dat.rescale('l2norm_images');
            fprintf('\nNormalizing input images by l2norm\n\n');

        case 'zscoreimages'

            fmri_dat = fmri_dat.rescale('zscoreimages');
            fprintf('\nZ-scoring input images\n\n');

        case 'zscorevoxels'

            fmri_dat = fmri_dat.rescale('zscorevoxels');
            fprintf('\nZ-scoring voxels across input images\n\n');

        otherwise 

            error('\ninvalid scaling option %s specified in myscaling_mvpa_reg_st variable defined in a2_set_default_options script, please correct\n\n', myscaling_mvpa_reg_st);

    end % switch scaling


%% CREATE MASK FROM ATLAS, AND APPLY TO FMRI_DATA_ST OBJECT
%--------------------------------------------------------------------------

fprintf('\n\n');
printhdr(['CREATING MASK FROM ATLAS ' atlas_name]);
fprintf('\n\n');

mask = select_atlas_subset(atlas,regions);

fmri_dat_mask = apply_mask(fmri_dat,mask);
fmri_dat_mask = trim_mask(fmri_dat);


%% RUN RSA REGRESSION AND VISUALIZE RESULTS
%--------------------------------------------------------------------------

fprintf('\n\n');
printhdr('CREATING REPRESENTATIONAL DISSIMILARITY MATRIX');
fprintf('\n\n');

rsa_stats = rsa_regression(fmri_dat_mask, cond_id, subject_id);


fprintf('\n\n');
printhdr('VISUALIZING RSA RESULTS');
fprintf('\n\n');

figure;
imagesc(rsa_stats.RDM); 
colorbar; 
xlabel 'Con Number'; 
title 'Brain RDM';

%show bootstrap distribution for generalization indices
figure;
histogram(rsa_stats.bs_gen_index(:,2));
title 'Condition';
ylabel 'Generalization Index';
orthviews(mask);