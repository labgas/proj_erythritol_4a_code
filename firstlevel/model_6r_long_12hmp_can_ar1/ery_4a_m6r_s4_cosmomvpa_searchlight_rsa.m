%% ery_4a_m6r_s4_cosmomvpa_searchlight_rsa.m
%
%
% *USAGE*
%
% This script uses functions of
% <https://www.cosmomvpa.org/index.html CoSMoMVPA> to
%
%
% # run first-level (i.e. within-subject) representational similarity
%   analysis comparing behavioral and neural dissimilarity matrices
%   for conditions of your choice in a mask of your choice, using
%   SPM.mat files as input and a leave-one-run-out crossvalidation
%   procedure
%
% # calculates and displays group-level results using permutation testing
%   and threshold-free cluster enhancement 
%
% # saves the results using standard naming and location
% 
% Run this script with Matlab's publish function to generate html report of results:
% publish('LaBGAScore_cosmomvpa_searchlight_rsa.m','outputDir',htmlsavedir)
%
% *OPTIONS*
%
% * results_suffix                  string to add to results files, useful if you want to run multiple analyses for example with different conditions
% * conds2include                   indices of conditions in DSGN.conditions (or SPM.Sess(x).U.name); assumes all conditions are present in all runs!
% * mask_name_rsa                   absolute path to mask, or name of mask file if already on Matlab path
% * atlas_name_rsa                  atlas object used for labeling regions
%                                       Two options:
%                                           1. atlas name from load_atlas.m for different atlas than
%                                               default canlab_2018 (used if you do not specify this option)
%                                           2. which('xxx.mat') name of .mat file in maskdir, generated by
%                                               https://github.com/labgas/LaBGAScore/blob/main/atlas_mask_tools/LaBGAScore_atlas_binary_mask_from_atlas.m
%                                               this file should contain a 'combined_atlas' var containing
%                                               an atlas object
%                                               USAGE: use this if you want to restrict the parcels to the
%                                                      ones included in maskname_rsa, or any atlas or a subset thereof
% * analysis_level                  'run', or 'subject' - the latter averages dissimilarity matrices over runs
%
% if analysis_level = 'subject'
% * phenofile_name                  name of file containing subject-level data for behavioral dissimilarity matrix as stored in BIDS/phenotype dir (in wide format, one row per subject, one var per condition)
% * phenofile_varnames              variable names containing behavioral data for different conditions
%
% if analysis_level = 'run'  
% * phenofile_st_name               name of file containing trial-level data for behavioral dissimilarity matrix as stored in BIDS/phenotype dir (in long format, one row per trial, one var containing behavioral data for all conditions)
% * phenofile_st_varname            single variable name containing behavioral data for different conditions
% * phenofile_st_subjid             single variable name containing subject identifiers
% * phenofile_st_runid              single variable name containing run identifiers
% * phenofile_st_condid             single variable name containing condition identifiers
% * do_null_data                    true/false - the former generates null data for every subject and uses it for the purpose of group level permutation stats, see docs cosmo_montecarlo_cluster_stats
%
% For more info on CoSMoMVPA, check out the
% <frontiersin.org/articles/10.3389/fninf.2016.00027/full accompanying paper>
%
% -------------------------------------------------------------------------
%
% Lukas Van Oudenhove
%
% date:   KU Leuven, June, 2023
%
% -------------------------------------------------------------------------
%
% LaBGAScore_cosmomvpa_searchlight_rsa.m         v1.1
%
% last modified: 2023/08/08
%
%
%% SET OPTIONS
% -------------------------------------------------------------------------

results_suffix = 'all_conds_centered';
conds2include = [1:4]; % indices of conditions in DSGN.conditions (or SPM.Sess(x).U.name); assumes all conditions are present in all runs!
mask_name_rsa = 'ery_4a_m6_mask_all_regions.nii'; % absolute path to mask, or name of mask file if already on Matlab path
atlas_name_rsa = which('ery_4a_m6_mask_all_regions.mat'); % atlas object used for labeling regions, see docs above for more info
analysis_level = 'run'; % alternative 'subject' - averages neural and behavioral dissimilarity matrices over runs prior to analysis; single trial not yet implemented

% if analysis level = 'subject'

phenofile_name = 'ratings_means.tsv'; % name of file containing data for average behavioral dissimilarity matrix as stored in BIDS/phenotype dir
phenofile_varnames = {'rating_sucrose','rating_erythritol','rating_sucralose','rating_water'}; % in same order as in DSGN.conditions!

% if analysis level = 'run'

phenofile_st_name = 'ratings_all.tsv'; % name of file containing data for single trial behavioral dissimilarity matrix as stored in BIDS/phenotype dir
phenofile_st_varname = 'rating';
phenofile_st_subjid = 'participant_id';
phenofile_st_runid = 'run_id';
phenofile_st_condid = 'trial_type';
do_null_data = false; % use null data in group-level permutation testing, see docs cosmo_montecarlo_cluster_stats


%% PREP WORK
% -------------------------------------------------------------------------

% set analysis name

analysis_name = 'rsa_searchlight';

% check whether LaBGAScore_prep_s0_define_directories has been run
% STUDY-SPECIFIC: replace LaBGAScore with study name in code below

if ~exist('rootdir','var') || ~exist('githubrootdir','var')
    fprintf('\n');
    warning('rootdir and/or githubrootdir variable not found in Matlab workspace, running LaBGAScore_prep_s0_define_directories before proceeding')
    fprintf('\n\n');
    ery_4a_prep_s0_define_directories;
    cd(rootdir);
else
    cd(rootdir);
end

% check whether LaBGAScore_firstlevel_s1_options_dsgn_struct.m has been run
% STUDY-SPECIFIC: replace LaBGAScore with study name and add model index in code below

if ~exist('DSGN','var')
    fprintf('\n');
    warning('DSGN variable not found in Matlab workspace, running LaBGAScore_firstlevel_s1_options_dsgn_struct.m before proceeding')
    fprintf('\n\n');
    ery_4a_firstlevel_m6r_s1_options_dsgn_struct;
end

% check whether LaBGAScore_firstlevel_s0_a_set_up_paths_always_run_first.m has been run
% STUDY-SPECIFIC: replace LaBGAScore with study name and add model index in code below

if ~exist('htmlsavedir','var')
    fprintf('\n');
    warning('DSGN variable not found in Matlab workspace, running LaBGAScore_s0_a_set_up_paths_always_run_first.m before proceeding')
    fprintf('\n\n');
    ery_4a_secondlevel_m6r_s0_a_set_up_paths_always_run_first;
end

% define phenotype dir & load behavioral data file

phenodir = fullfile(BIDSdir,'phenotype');

switch analysis_level
    
    case 'subject'
        
        phenofile = readtable(fullfile(phenodir, phenofile_name),'FileType', 'text', 'Delimiter', 'tab');
        phenofile_reduced = table2array(phenofile(:,phenofile_varnames));
        
    case 'run'
        
        phenofile_st = readtable(fullfile(phenodir, phenofile_st_name),'FileType', 'text', 'Delimiter', 'tab');
        phenofile_st.(phenofile_st_condid) = categorical(phenofile_st.(phenofile_st_condid));
        phenofile_st.(phenofile_st_subjid) = categorical(phenofile_st.(phenofile_st_subjid));
        
end

% get first level dir info

firstmodeldir = DSGN.modeldir;

firstlist = dir(fullfile(firstmodeldir,'sub-*'));
firstsubjs = cellstr(char(firstlist(:).name));
firstsubjdirs= cell(size(firstsubjs,1),1);

    for firstsub = 1:size(firstsubjs,1)
        firstsubjdirs{firstsub,1} = fullfile(firstlist(firstsub).folder,firstlist(firstsub).name);
    end
    
clear firstlist

% set firstlevel cosmo dir info & create dirs if needed

firstmodelcosmodir = fullfile(firstmodeldir,'cosmo');
if ~isfolder(firstmodelcosmodir)
    mkdir(firstmodelcosmodir);
end

firstmodelcosmodir_mask = fullfile(firstmodelcosmodir,'mask');
if ~isfolder(firstmodelcosmodir_mask)
    mkdir(firstmodelcosmodir_mask);
end

% set secondlevel TDT dir info & create dir if needed

secondmodelcosmodir = fullfile(resultsdir,'cosmo');
if ~isfolder(secondmodelcosmodir)
    mkdir(secondmodelcosmodir);
end

secondmodelcosmoanalysisdir = fullfile(secondmodelcosmodir, analysis_name);
if ~isfolder(secondmodelcosmoanalysisdir)
    mkdir(secondmodelcosmoanalysisdir);
end

% atlas prep work

if exist('atlas_name_rsa','var') && ~isempty(atlas_name_rsa) && exist(atlas_name_rsa, 'file')

    if contains(atlas_name_rsa,'.mat')
        
        load(atlas_name_rsa);
        [~, atlas_name_short] = fileparts(atlas_name_rsa);
        
        fprintf('\nLabeling regions using custom-made atlas %s\n\n', atlas_name_short);
        
        combined_atlas.atlas_name = atlas_name_short;
        combined_atlas.image_names = [atlas_name_short '.mat'];
        combined_atlas.fullpath = atlas_name_rsa;

    elseif ischar(atlas_name_rsa)

        combined_atlas = load_atlas(atlas_name_rsa);
        
        fprintf('\nLabeling regions using custom-made atlas %s\n\n', atlas_name_rsa);
        
        combined_atlas.atlas_name = atlas_name_rsa;
        combined_atlas.image_names = [atlas_name_rsa '.nii'];
        combined_atlas.fullpath = which([atlas_name_rsa '.nii']);

    else

         error('\ninvalid option "%s" defined in atlas_name_rsa variable, should be a keyword for load_atlas.m or a .mat file containing an atlas object, check docs"\n\n',atlas_name_rsa)

    end
    
else
    
    fprintf('\nLabeling regions using default canlab_2018 atlas\n\n');
    
end

% mask prep work

clear mask

mask_img = fmri_mask_image(which(mask_name_rsa),'noverbose');
target = fmri_data(fullfile(firstsubjdirs{1},'beta_0001.nii'),'noverbose'); % resample to space of functional images
mask2write = resample_space(mask_img,target);
mask2write.dat(mask2write.dat > 0) = 1; % binarize mask
mask_name_final = [mask_name_rsa(1:end-4) '_binary'];
mask = fullfile(firstmodelcosmodir_mask,[mask_name_final '.nii']);
write(mask2write,'fname',mask,'overwrite');

% pre-allocate loop vars

cosmo_input_datasets = cell(1,firstsub);
cosmo_behav_dsms = cell(1,firstsub);
cosmo_rsa_datasets = cell(1,firstsub);
cosmo_rsa_niifiles = cell(1,firstsub);
cosmo_random_input_datasets = cell(firstsub,100);
cosmo_null_datasets = cell(firstsub,100); % 100 null datasets per subjects recommended in cosmo_monte_carlo_cluster_stats docs
cosmo_stacked_null_datasets = cell(1,100);


%% LOOP COSMO CODE OVER SUBJECTS
% -------------------------------------------------------------------------

for sub = 1:firstsub
    
    fprintf('\n\n');
    printhdr(['SUBJECT #', num2str(sub)]);
    fprintf('\n\n');

    % CREATE COSMO DATASET, CLEAN, AND DISPLAY
    % ----------------------------------------
    ds = cosmo_fmri_dataset(fullfile(firstsubjdirs{sub},'SPM.mat'),'mask',mask);
    ds = cosmo_remove_useless_data(ds);
    cosmo_disp(ds);

    % LOAD CONDITION LABELS AND ADD TO DATASET
    % ----------------------------------------
    nr_runs = (size(ds.samples,1)/size(DSGN.conditions,2));
    ds.sa.labels = repmat(DSGN.conditions{1}',nr_runs,1); % assuming all conditions present in all runs!
    ds.sa.targets = repmat([1:size(DSGN.conditions,2)]',nr_runs,1);

    % SELECT CONDITIONS OF INTEREST AND AVERAGE OVER RUNS
    % ---------------------------------------------------
    ds_ofint = cosmo_slice(ds,(ds.sa.targets < (size(conds2include,2)+1)),1);
    
    if strcmpi(analysis_level,'run')
        cosmo_disp(ds_ofint);
    end
    
    if strcmpi(analysis_level,'subject')
        ds_ofint_runavg = cosmo_fx(ds_ofint, @(x)mean(x,1), 'targets', 1);
        cosmo_disp(ds_ofint_runavg);
    end

    % LOAD BEHAVIORAL RATINGS, CALCULATE, AND PLOT DISSIMILARITY MATRIX
    % -----------------------------------------------------------------
    
    switch analysis_level
        
        case 'subject'
    
            behav = phenofile_reduced(sub,:)';
            dsm_behav_vect = cosmo_pdist(behav); % row vector format
            dsm_behav_mat = cosmo_squareform(dsm_behav_vect); % matrix format
            
            figure;
            heatmap(ds_ofint_runavg.sa.labels,ds_ofint_runavg.sa.labels,dsm_behav_mat, 'Colormap', jet);
            figtitle = [firstsubjs{sub} ' dissimilarity matrix (Euclidean distance)'];
            set(gca,'Title',figtitle);
            set(gcf,'WindowState','maximized');
            drawnow, snapnow;
            
        case 'run'
    
            behav_sub_table = phenofile_st(phenofile_st.(phenofile_st_subjid) == firstsubjs{sub},:);

            behav_runs = [];
            behav_runx = [];

            for j = unique(behav_sub_table.(phenofile_st_runid))'
                behav_runx_table = behav_sub_table(behav_sub_table.(phenofile_st_runid) == j,:);
                    for k = 1:size(phenofile_varnames,2)
                        behav_runx(1,k) = mean(behav_runx_table.(phenofile_st_varname)(behav_runx_table.(phenofile_st_condid) == ds_ofint.sa.labels{k,1}));
                    end
                    behav_runs = [behav_runs,behav_runx];
            end

            dsm_behav_vect_runs = cosmo_pdist(behav_runs');
            dsm_behav_mat_runs = cosmo_squareform(dsm_behav_vect_runs');
            
            figure;
            heatmap(dsm_behav_mat_runs, 'Colormap', jet);
            figtitle = [firstsubjs{sub} ' dissimilarity matrix (Euclidean distance)'];
            set(gca,'Title',figtitle);
            set(gcf,'WindowState','maximized');
            drawnow, snapnow;
            
    end
        
    % SET UP COSMO MEASURE AND NEIGHBORHOOD
    % -------------------------------------
    measure = @cosmo_target_dsm_corr_measure;
    measure_args = struct();
    measure_args.center_data = true;
    
    switch analysis_level
        
        case 'subject'
            
            measure_args.target_dsm = dsm_behav_mat;
            nh = cosmo_spherical_neighborhood(ds_ofint_runavg,'count',100);
            
        case 'run'
            
            measure_args.target_dsm = dsm_behav_mat_runs;
            nh = cosmo_spherical_neighborhood(ds_ofint,'count',100);
            
    end

    % RUN SEARCHLIGHT RSA
    % -------------------
    
    switch analysis_level
        
        case 'subject'
            
            ds_rsa = cosmo_searchlight(ds_ofint_runavg,nh,measure,measure_args);
            
        case 'run'
            
            ds_rsa = cosmo_searchlight(ds_ofint,nh,measure,measure_args);
            
            if do_null_data
            
                parfor p = 1:100 % 100 null datasets per subjects recommended in cosmo_monte_carlo_cluster_stats docs

                    cosmo_random_input_datasets{sub,p} = ds_ofint;
                    cosmo_random_input_datasets{sub,p}.sa.targets = cosmo_randomize_targets(ds_ofint,'seed',12345);

                    cosmo_null_datasets{sub,p} = cosmo_searchlight(cosmo_random_input_datasets{sub,p},nh,measure,measure_args);

                end
                
            end
            
    end
    
    cosmo_disp(ds_rsa);

    % WRITE .NII FILE
    % ---------------
    firstsubjcosmodir = fullfile(firstsubjdirs{sub},'cosmo',analysis_name);
        if ~isfolder(firstsubjcosmodir)
            mkdir(firstsubjcosmodir);
        end
    niifile = fullfile(firstsubjcosmodir,[analysis_name '_' analysis_level '_' mask_name_final '_' results_suffix '.nii']);
    cosmo_map2fmri(ds_rsa, niifile);
    
    % PLOT MONTAGE
    % ------------
    
    corrmap = fmri_data(niifile);

    fprintf ('\nMONTAGE UNTHRESHOLDED SUBJECT LEVEL RESULTS OF %s ANALYSIS, MASK: %s, SUBJECT: %s', upper(analysis_name), upper(mask_name_final), firstsubjs{sub});

    figure;

    o2 = canlab_results_fmridisplay(region(corrmap), 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]}, 'cmaprange', [-1 1],'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');

    o2 = title_montage(o2, 5, [upper(analysis_name) ' UNTHRESHOLDED ' upper(mask_name_final) ' ' upper(firstsubjs{sub})]);
    
    o2 = legend(o2);

    figtitle = sprintf('%s_%s_unthresholded_montage_%s_%s', analysis_name, analysis_level, firstsubjs{sub}, mask_name_final);
    set(gcf, 'Tag', figtitle, 'WindowState','maximized');
    drawnow, snapnow;
    
    
    % STORE COSMO INPUT AND RESULTS DATASETS
    % --------------------------------------
    
    switch analysis_level
        
        case 'subject'
            
            cosmo_input_datasets{sub} = ds_ofint_runavg;
            cosmo_behav_dsms{sub} = dsm_behav_mat;
            
        case 'run'
            
            cosmo_input_datasets{sub} = ds_ofint;
            cosmo_behav_dsms{sub} = dsm_behav_mat_runs;
            
    end
    
    cosmo_rsa_datasets{sub} = ds_rsa;
    cosmo_rsa_niifiles{sub} = niifile;
    
    % CLEAR VARS
    % ----------
    clear ds* nr_runs behav* measure neighborhood firstsubjcosmodir niifile corrmap o2;

    
end % for loop over subjects


%% SAVE RESULTS
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('SAVING SUBJECT-LEVEL RSA RESULTS');
fprintf('\n\n');

savefilenamedata = fullfile(secondmodelcosmoanalysisdir, [analysis_name '_' analysis_level '_' mask_name_final '_' results_suffix '.mat']);
save(savefilenamedata, 'cosmo_input_datasets', 'cosmo_rsa_datasets','cosmo_rsa_niifiles','cosmo_behav_dsms', 'cosmo_random_input_datasets', 'cosmo_null_datasets', '-v7.3');

fprintf('\nSaved results for %s, analysis level %s\n', analysis_name, analysis_level);
fprintf('\nFilename: %s\n', savefilenamedata);


%% DO GROUP LEVEL STATS
% -------------------------------------------------------------------------

% CREATE STACKED COSMO DATASETS
% -----------------------------
ds_group = cosmo_stack(cosmo_rsa_datasets);
ds_group.sa.chunks = (1:size(cosmo_rsa_datasets,2))';
ds_group.sa.targets = ones(size(cosmo_rsa_datasets,2),1);

if do_null_data

    for p = 1:100
        cosmo_stacked_null_datasets{1,p} = cosmo_stack(cosmo_null_datasets(:,p));
    end

    for q = 1:100
        cosmo_stacked_null_datasets{1,q}.sa.chunks = (1:size(cosmo_rsa_datasets,2))';
        cosmo_stacked_null_datasets{1,q}.sa.targets = ones(size(cosmo_rsa_datasets,2),1);
    end
    
end

% CREATE NEIGHBORHOOD
% -------------------
nh_group = cosmo_cluster_neighborhood(ds_group);

% DEFINE OPTIONS
% --------------
opt_group = struct();
opt_group.niter = 10000; % 10k for publications
opt_group.h0_mean = 0;
opt_group.nproc = 1; % parallel throws warning, not sure whether this is a problem, hence we don't use parallel for now
opt_group.seed = 12345; % makes output deterministic

% RUN PERMUTATION TEST
% --------------------

if ~do_null_data
    
    z_ds_group = cosmo_montecarlo_cluster_stat(ds_group,nh_group,opt_group);
    
else
    
    z_ds_group = cosmo_montecarlo_cluster_stat(ds_group,nh_group,opt_group,'null',cosmo_stacked_null_datasets);
    
end

cosmo_disp(z_ds_group);

niifile_group = fullfile(secondmodelcosmoanalysisdir,[analysis_name '_' analysis_level '_' mask_name_final '_' results_suffix '.nii']);
cosmo_map2fmri(z_ds_group, niifile_group);

% PLOT MONTAGES
% -------------

corrmap_group = fmri_data(niifile_group);
corrmap_group_thresholded = threshold(corrmap_group,[-1.96 1.96],'raw-outside');

fprintf('\n\n');
printhdr('GROUP LEVEL RESULTS');
fprintf('\n\n');

% UNTHRESHOLDED

fprintf ('\nMONTAGE UNTHRESHOLDED CORRECTED GROUP LEVEL RESULTS OF %s ANALYSIS AT LEVEL %s, MASK: %s\n', upper(analysis_name), upper(analysis_level), upper(mask_name_final));

corrmap_group_region = region(corrmap_group,'noverbose');

figure;

o2 = canlab_results_fmridisplay(corrmap_group_region, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]}, 'cmaprange', [min(corrmap_group.dat) 0 0 max(corrmap_group.dat)],'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');

o2 = title_montage(o2, 5, [upper(analysis_name) ' UNTHRESHOLDED CORRECTED ' upper(analysis_level) ' LEVEL ' upper(mask_name_final)]);

o2 = legend(o2);

figtitle = sprintf('%s_%s_unthresholded_montage_%s', analysis_name, analysis_level, mask_name_final);
set(gcf, 'Tag', figtitle, 'WindowState','maximized');
drawnow, snapnow;

if ~isempty(corrmap_group_region)

    if exist('combined_atlas','var')
        [rpos, rneg, r_table_unthresh] = table(corrmap_group_region,'atlas_obj',combined_atlas);     % add labels from combined_atlas
        r_unthresh = [rpos rneg];                                                                    % re-concatenate labeled regions
    else
        [rpos, rneg, r_table_unthresh] = table(corrmap_group_region);                                % add labels from default canlab_2018 atlas 
        r_unthresh = [rpos rneg];                                                                    % re-concatenate labeled regions
    end

    % Montage of regions in table (plot and save)

    fprintf ('\nMONTAGE REGIONCENTERS UNTHRESHOLDED CORRECTED GROUP LEVEL RESULTS OF %s ANALYSIS AT LEVEL %s, MASK: %s\n', upper(analysis_name), upper(analysis_level), upper(mask_name_final));
    
    o3 = montage(r_unthresh, 'regioncenters', 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'cmaprange', [min(corrmap_group.dat) 0 0 max(corrmap_group.dat)]);

    % Activate, name, and save figure
    set(gcf, 'Tag', figtitle, 'WindowState','maximized');
    drawnow, snapnow;
    clear rpos rneg o3 figtitle
    
else
    
    r_unthresh = [];
    r_table_unthresh = [];

end % conditional montage plot if there are regions to show

clear o2

% THRESHOLDED

figure;

fprintf ('\nMONTAGE P < 0.05 CORRECTED GROUP LEVEL RESULTS OF %s ANALYSIS AT LEVEL %s, MASK: %s\n', upper(analysis_name), upper(analysis_level), upper(mask_name_final));

corrmap_group_thresh_region = region(corrmap_group_thresholded,'noverbose');

o2 = canlab_results_fmridisplay(corrmap_group_thresh_region,'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]}, 'cmaprange', [min(corrmap_group_thresholded.dat) -1.96 -1.96 -1.96], 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img'); 
% adapted since we only have negative values in this particular case, in case of positive values 'cmaprange', [min(corrmap_group_thresholded.dat) -1.96 1.96 max(corrmap_group_thresholded.dat)] 

o2 = title_montage(o2, 5, [upper(analysis_name) ' P < 0.05 CORRECTED ' upper(analysis_level) ' LEVEL ' upper(mask_name_final)]);

o2 = legend(o2);

figtitle = sprintf('%s_%s_corrected_montage_%s', analysis_name, analysis_level, mask_name_final);
set(gcf, 'Tag', figtitle, 'WindowState','maximized');
drawnow, snapnow;

if ~isempty(corrmap_group_thresh_region)

    if exist('combined_atlas','var')
        [rpos, rneg, r_table_thresh] = table(corrmap_group_thresh_region,'atlas_obj',combined_atlas);     % add labels from combined_atlas
        r_thresh = [rpos rneg];                                                                    % re-concatenate labeled regions
    else
        [rpos, rneg, r_table_thresh] = table(corrmap_group_thresh_region);                                % add labels from default canlab_2018 atlas 
        r_thresh = [rpos rneg];                                                                    % re-concatenate labeled regions
    end

    % Montage of regions in table (plot and save)

    fprintf ('\nMONTAGE REGIONCENTERS P < 0.05 CORRECTED GROUP LEVEL RESULTS OF %s ANALYSIS AT LEVEL %s, MASK: %s\n', upper(analysis_name), upper(analysis_level), upper(mask_name_final));
    
    o3 = montage(r_thresh, 'regioncenters', 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'cmaprange', [min(corrmap_group.dat) 0 0 max(corrmap_group.dat)]);

    % Activate, name, and save figure
    set(gcf, 'Tag', figtitle, 'WindowState','maximized');
    drawnow, snapnow;
    clear rpos rneg o3 figtitle
    
else

    r_thresh = [];
    r_table_thresh = [];

end % conditional montage plot if there are regions to show

clear o2


%% SAVE SECOND LEVEL RESULTS
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('SAVING CROSS-CLASSIFICATION ACCURACY GROUP LEVEL STATS');
fprintf('\n\n');

savefilenamedata = fullfile(secondmodelcosmoanalysisdir, ['group_level_stats_' analysis_name '_' analysis_level '_' mask_name_final '_' results_suffix '.mat']);
save(savefilenamedata, 'ds_group', 'cosmo_stacked_null_datasets', 'nh_group', 'opt_group', 'z_ds_group', 'r_unthresh', 'r_table_unthresh', 'r_thresh', 'r_table_thresh', '-v7.3');

fprintf('\nSaved results for %s\n', analysis_name);
fprintf('\nFilename: %s\n', savefilenamedata);