%% cosmovpa searchlight rsa on one subject

% load dataset, clean, and display
ds = cosmo_fmri_dataset('/data/proj_erythritol/proj_erythritol_4a/firstlevel/model_6r_long_12hmp_can_ar1/sub-003/SPM.mat');
ds = cosmo_remove_useless_data(ds);
cosmo_disp(ds)

% add load condition labels and add to dataset
load('DSGN.mat');
ds.sa.labels = repmat(DSGN.conditions{1}',6,1); % assuming all conditions present in all runs
ds.sa.targets = repmat([1:6]',6,1); % same, these two lines are still hardcoded but easy to adapt

% select conditions of interest only & average over runs
ds_ofint = cosmo_slice(ds,(ds.sa.targets < 5),1);
ds_ofint_runavg = cosmo_fx(ds_ofint, @(x)mean(x,1), 'targets', 1);

% load behavioral ratings averaged over runs & calculate dissimilarity matrix
behav = [6.3 10.9 7.3 0]'; % ratings in order of conditions/labels in dataset, hardcoded, can be easily read in by loading /data/proj_erythritol/proj_erythritol_4a/BIDS/phenotype/ratings_means.tsv
dsm_behav_vect = cosmo_pdist(behav); % row vector format
dsm_behav_mat = cosmo_squareform(dsm_behav_vect); % matrix format

% set up measure & neighborhood
measure = @cosmo_target_dsm_corr_measure;
measure_args = struct();
measure_args.target_dsm = dsm_behav_mat;
neighborhood = cosmo_spherical_neighborhood(ds_ofint_runavg,'count',100);

% run searchlight
ds_rsa = cosmo_searchlight(ds_ofint_runavg,neighborhood,measure,measure_args);

% write to .nii file
cosmo_map2fmri(ds_rsa,'ds_rsa.nii');