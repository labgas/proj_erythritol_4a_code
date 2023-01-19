%% ery_4a_secondlevel_m6_s7b_c2a_second_level_reg_rating.m
%
%
% USAGE
%
% This script displays second⁻level (i.e. across subjects) regression
% results generated by prep_3a_run_second_level_regression_and_save.m
%
% Run this script with Matlab's publish function to generate html report of results:
% publish('c2a_second_level_regression','outputDir',htmlsavedir)
% 
% See the documentation of prep_3a_run_second_level_regression_and_save.m
% for more info and options
%
% OPTIONS
% 
% NOTE: 
% defaults are specified in a2_set_default_options for any given model,
% but if you want to run the same model with different options, 
% you can make a copy of this script with a letter index (e.g. _s6a_) 
% and change the default option below
%
% save_figures_glm: true saves .svg files of all figures generated by c2a_second_level_regression.m (slow, takes up space) % added by @lukasvo76 May 2022
% q_threshold_glm: threshold for FDR-corrected display items
% p_threshold_glm: threshold for uncorrected display items
% k_threshold_glm: extent threshold for both corrected and uncorrected display items
% BF_threshold_glm: threshold for Bayes Factor maps, |BF| > 10 indicates strong evidence in favour of H1 (positive value) or H0 (negative value) - see help.statistic_image.estimateBayesFactor for details
%
% MANDATORY OPTIONS TO COPY FROM CORRESPONDING PREP_3a_ SCRIPT
%
% - mygroupfieldname: 'contrasts' or 'conditions'
% - results_suffix: name added to results file by prep_3a script in case of multiple versions of model, e.g. 'covariate_rating'
% - myscaling_glm: 'raw', 'scaled', or 'scaled_contrasts' (defined in a2_set_..., image scaling done in prep_2_... and prep_3_... data load)
% - maskname_glm: 
%       - default use of sparse gray matter mask
%       - model-specific maskdir defined in a_set_up_paths_always_run_first script
%       - if you do not want to mask, change to []
%       - if you want to use a custom mask, put it in maskdir and change name here
%       - applied before fdr correction in this script
%
%__________________________________________________________________________
%
% revamped by: Lukas Van Oudenhove
% date:   Dartmouth, May, 2022
%
%__________________________________________________________________________
% @(#)% c2a_second_level_regression.m         v4.0
% last modified: 2023/01/17


%% GET AND SET OPTIONS
% -------------------------------------------------------------------------

% GET MODEL-SPECIFIC PATHS AND OPTIONS

ery_4a_secondlevel_m6_s0_a_set_up_paths_always_run_first;

% NOTE: CHANGE THIS TO THE MODEL-SPECIFIC VERSION OF THIS SCRIPT
% NOTE: THIS WILL ALSO AUTOMATICALLY CALL A2_SET_DEFAULT_OPTIONS

% SET/COPY MANDATORY OPTIONS FROM CORRESPONDING PREP_3a_ SCRIPT

mygroupnamefield = 'contrasts'; 
results_suffix = 'test_rating'; % suffix of your choice added to .mat file with saved results
myscaling_glm = 'raw';
maskname_glm = which('gray_matter_mask_sparse.img'); % applied in this script before fdr thresholding

% SET CUSTOM OPTIONS

% NOTE: only specify if you want to run a second version of your model with different options
% than the defaults you set in your model-specific version of a2_set_default_options.m

% save_figures_glm = true/false;
% q_threshold_glm = .x;
% p_threshold_glm = .yyy;
% k_threshold_glm = zz;
% BF_threshold_glm = x;

% GET SCALING OPTION

switch myscaling_glm

    case 'raw'
        fprintf('\ncontrast calculated on raw (unscaled) condition images used in second-level GLM\n\n');
        scaling_string = 'no_scaling';

    case 'scaled'
        fprintf('\ncontrast calculated on z-scored condition images used in second-level GLM\n\n');
        scaling_string = 'scaling_z_score_conditions';

    case 'scaled_contrasts'
        fprintf('\nl2norm scaled contrast images used in second-level GLM\n\n');
        scaling_string = 'scaling_l2norm_contrasts';

    otherwise
        error('\ninvalid option "%s" defined in myscaling_glm variable in a2_set_default_options script, choose between "raw", "scaled", or "scaled_contrast" given option "%s" defined in mygroupnamefield variable\n', myscaling_glm, mygroupnamefield);

end


% GET MASKING OPTION

if exist('maskname_glm', 'var') && ~isempty(maskname_glm) && exist(maskname_glm, 'file')
    apply_mask_before_fdr = true;
    [~,maskname_short] = fileparts(maskname_glm);
    mask_string = sprintf('masked_with_%s', maskname_short);
    glmmask = fmri_mask_image(maskname_glm, 'noverbose'); 
else
    apply_mask_before_fdr = false;
    mask_string = sprintf('without_masking');
end  


%% LOAD RESULTS IF NEEDED
% -------------------------------------------------------------------------

% GLM
% -------------------------------------------------------------------------

if ~dorobfit_parcelwise
    resultsvarname = 'regression_stats_results';
    resultsstring = 'regression_stats_and_maps_';
    analysis_type = 'voxel-wise';
else
    resultsvarname = 'parcelwise_stats_results';
    resultsstring = 'parcelwise_stats_and_maps_';
    analysis_type = 'parcel-wise';
end

if ~exist(resultsvarname, 'var')
    
    fprintf('\n\n');
    printhdr('LOADING GLM DATA');
    fprintf('\n\n');

    savefilenamedata = fullfile(resultsdir, [resultsstring, mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']);

    if exist(savefilenamedata,'file')
        fprintf('\nLoading %s regression results and maps from %s\n\n', analysis_type, savefilenamedata);
        load(savefilenamedata, resultsvarname);
    else
        fprintf('\nNo saved results file %s. Skipping this analysis.', savefilenamedata);
        fprintf('\nRun prep_3a_run_second_level_regression_and_save.m to get %s regression results first.\n', analysis_type); 
        return
    end

else
    fprintf('\n%s %s found, displaying results\n\n', resultsvarname, analysis_type);

end

if ~dorobfit_parcelwise    
    results = regression_stats_results;

else
    results = parcelwise_stats_results;

end

% MVPA
% -------------------------------------------------------------------------

if domvpa_reg_cov
    mvpa_resultsvarname = 'mvpa_stats_results';
    mvpa_resultsstring = 'mvpa_stats_and_maps_';
    mvpa_analysis_type = algorithm_mvpa_reg_cov;
end

if ~exist('mvpa_resultsvarname','var')
    
    fprintf('\n\n');
    printhdr('LOADING MVPA DATA');
    fprintf('\n\n');
    
    savefilenamedata_mvpa = fullfile(resultsdir, ['mvpa_stats_and_maps_', mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']);

    if exist(savefilenamedata,'file')
        fprintf('\nLoading %s mvpa regression results and maps from %s\n\n', mvpa_analysis_type, savefilenamedata_mvpa);
        load(savefilenamedata_mvpa, mvpa_resultsvarname);
    else
        fprintf('\nNo saved results file %s. Skipping this analysis.', savefilenamedata_mvpa);
        fprintf('\nRun prep_3a_run_second_level_regression_and_save.m to get %s regression results first.\n', mvpa_analysis_type); 
        return
    end

else
    fprintf('\n%s %s found, displaying results\n\n', mvpa_resultsvarname, mvpa_analysis_type);

end
    
%% VISUALIZE GLM RESULTS FOR EACH CONTRAST
% -------------------------------------------------------------------------

for c = 1:size(results, 2) % number of contrasts or conditions

    analysisname = results{c}.analysis_name;
    names = results{c}.variable_names;
    
    names_string = names{1};
    
        for name = 2:size(names,2)
            names_string = [names_string,' ',names{name}];
        end
    
        if ~dorobfit_parcelwise
            t = results{c}.t;
        else
            t = results{c}.t_obj;
        end
        
    if doBayes
        
        BF = results{c}.BF;
        
    end
    
    if domvpa_reg_cov
        
        for covar = 1:size(mvpa_stats_results,2)
            mvpa_results{covar} = mvpa_stats_results{c,covar};
            mvpa_fmri_dats{covar} = mvpa_dats{c,covar};
        end
        
    end
    
    fprintf('\n\n');
    printhdr(['CONTRAST #', num2str(c), ': ', upper(analysisname)]);
    fprintf('\n\n');
    
    fprintf('\nREGRESSORS: %s\n\n', names_string);
    
        if isfield(results{c}, 'design_table')
            disp(results{c}.design_table);
            fprintf('\n\n');
        end
    
    num_effects = size(t.dat, 2); % number of regressors
    
    % BETWEEN-SUBJECT REGRESSORS & INTERCEPT: FDR corrected
    % ---------------------------------------------------------------------
    
    fprintf('\n\n');
    printhdr('FDR-CORRECTED GLM RESULTS');
    fprintf('\n\n');
    
    fprintf ('\nMONTAGE GLM RESULTS AT FDR q < %1.4f, k = %d, CONTRAST: %s, REGRESSORS: %s, MASK: %s, SCALING: %s\n\n', q_threshold_glm, k_threshold_glm, analysisname, names_string, mask_string, scaling_string);
    
    o2 = canlab_results_fmridisplay([], 'multirow', num_effects, 'outline', 'linewidth', 0.5, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]}, 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
    
        for j = 1:num_effects

            tj = get_wh_image(t, j);
                if apply_mask_before_fdr
                    tj = apply_mask(tj, glmmask);
                end
            tj = threshold(tj, q_threshold_glm, 'fdr', 'k', k_threshold_glm); 

            o2 = addblobs(o2, region(tj), 'wh_montages', (2*j)-1:2*j);
            o2 = title_montage(o2, 2*j, [analysisname ' FDR ' num2str(q_threshold_glm) ' ' names{j} ' ' mask_string ' ' scaling_string]);

        end % for loop over regressors in model
    
    figtitle = sprintf('%s_%s_%1.4f_FDR_montage_%s_%s_%s', analysisname, results_suffix, q_threshold_glm, names_string, mask_string, scaling_string);
    set(gcf, 'Tag', figtitle, 'WindowState','maximized');
    drawnow, snapnow;
    
        if save_figures_glm
            plugin_save_figure;
        end
        
    clear o2, clear figtitle, clear j, clear tj
    
    fprintf ('\nTABLES AND MONTAGE REGIONCENTERS GLM RESULTS AT FDR q < %1.4f, k = %d, CONTRAST: %s, REGRESSORS: %s, MASK: %s, SCALING: %s\n\n', q_threshold_glm, k_threshold_glm, analysisname, names_string, mask_string, scaling_string);
    
        for j = 1:num_effects

            fprintf ('\nTABLE GLM RESULTS AT FDR q < %1.4f, k = %d, CONTRAST: %s, REGRESSOR: %s, MASK: %s, SCALING: %s\n\n', q_threshold_glm, k_threshold_glm, analysisname, names{j}, mask_string, scaling_string);

            tj = get_wh_image(t, j);
                if apply_mask_before_fdr
                    tj = apply_mask(tj, glmmask);
                end
            tj = threshold(tj, q_threshold_glm, 'fdr', 'k', k_threshold_glm); 

            r = region(tj, 'noverbose');
            r(cat(1, r.numVox) < k_threshold_glm) = [];
            [rpos, rneg] = table(r);       % add labels
            r = [rpos rneg];               % re-concatenate labeled regions

            % Montage of regions in table (plot and save)
            if ~isempty(r)
                
                fprintf ('\nMONTAGE REGIONCENTERS GLM RESULTS AT FDR q < %1.4f, k = %d, CONTRAST: %s, REGRESSOR: %s, MASK: %s, SCALING: %s\n\n', q_threshold_glm, k_threshold_glm, analysisname, names{j}, mask_string, scaling_string);
                
                o3 = montage(r, 'colormap', 'regioncenters', 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});

                % Activate, name, and save figure
                figtitle = sprintf('%s_%s_%1.4f_FDR_regions_%s_%s_%s', analysisname, results_suffix, q_threshold_glm, names{j}, mask_string, scaling_string);
                set(gcf, 'Tag', figtitle, 'WindowState','maximized');
                drawnow, snapnow;
                    if save_figures_glm
                        plugin_save_figure;
                    end
                clear o3, clear figtitle, clear j, clear tj, clear r

            end % conditional montage plot if there are regions to show
            
        end % for loop over regressors in model

    
    % BETWEEN-SUBJECT REGRESSORS & INTERCEPT: uncorrected
    % ---------------------------------------------------------------------
    
    fprintf('\n\n');
    printhdr('UNCORRECTED GLM RESULTS');
    fprintf('\n\n');
    
    fprintf ('\nMONTAGE GLM RESULTS AT UNCORRECTED p < %1.4f, k = %d, CONTRAST: %s, REGRESSORS: %s, MASK: %s, SCALING: %s\n\n', p_threshold_glm, k_threshold_glm, analysisname, names_string, mask_string, scaling_string);
    
    o2 = canlab_results_fmridisplay([], 'multirow', num_effects, 'outline', 'linewidth', 0.5,'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]}, 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
    
        for j = 1:num_effects

            tj = get_wh_image(t, j);
                if apply_mask_before_fdr
                    tj = apply_mask(tj, glmmask);
                end
            tj = threshold(tj, p_threshold_glm, 'unc', 'k', k_threshold_glm); 

            o2 = addblobs(o2, region(tj), 'wh_montages', (2*j)-1:2*j);
            o2 = title_montage(o2, 2*j, [analysisname ' unc ' num2str(p_threshold_glm) ' ' names{j} ' ' mask_string ' ' scaling_string]);
        
        end % for loop over regressors in model
    
    figtitle = sprintf('%s_%s_%1.4f_unc_montage_%s_%s_%S', analysisname, results_suffix, p_threshold_glm, names_string, mask_string, scaling_string);
    set(gcf, 'Tag', figtitle, 'WindowState','maximized');
    drawnow, snapnow;
    
        if save_figures_glm
            plugin_save_figure;
        end
        
    clear o2, clear figtitle, clear j, clear tj
    
    fprintf ('\nTABLES AND MONTAGE REGIONCENTERS GLM RESULTS AT UNCORRECTED p < %1.4f, k = %d, CONTRAST: %s, REGRESSORS: %s, MASK: %s, SCALING: %s\n\n', p_threshold_glm, k_threshold_glm, analysisname, names_string, mask_string, scaling_string);
        
        for j = 1:num_effects

            fprintf ('\nTABLE GLM RESULTS AT UNCORRECTED p < %1.4f, k = %d, CONTRAST: %s, REGRESSOR: %s, MASK: %s, SCALING: %s\n\n', p_threshold_glm, k_threshold_glm, analysisname, names{j}, mask_string, scaling_string);

            tj = get_wh_image(t, j);
                if apply_mask_before_fdr
                    tj = apply_mask(tj, glmmask);
                end
            tj = threshold(tj, p_threshold_glm, 'unc', 'k', k_threshold_glm); 

            r = region(tj, 'noverbose');
            r(cat(1, r.numVox) < k_threshold_glm) = [];
            [rpos, rneg] = table(r);       % add labels
            r = [rpos rneg];               % re-concatenate labeled regions

            % Montage of regions in table (plot and save)
            if ~isempty(r)
                
                fprintf ('\nMONTAGE REGIONCENTERS GLM RESULTS AT UNCORRECTED p < %1.4f, k = %d, CONTRAST: %s, REGRESSOR: %s, MASK: %s, SCALING: %s\n\n', p_threshold_glm, k_threshold_glm, analysisname, names{j}, mask_string, scaling_string);
                
                o3 = montage(r, 'colormap', 'regioncenters', 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});

                % Activate, name, and save figure
                figtitle = sprintf('%s_%s_%1.4f_unc_regions_%s_%s_%s', analysisname, results_suffix, p_threshold_glm, names{j}, mask_string, scaling_string);
                set(gcf, 'Tag', figtitle, 'WindowState','maximized');
                drawnow, snapnow;
                    if save_figures_glm
                        plugin_save_figure;
                    end
                clear o3, clear figtitle, clear j, clear tj, clear r

            end % loop over regions in results
        
        end % loop over regressors
        
        
    % BETWEEN-SUBJECT REGRESSORS & INTERCEPT: Bayesian
    % ---------------------------------------------------------------------
    
    fprintf('\n\n');
    printhdr('BAYESIAN GLM RESULTS');
    fprintf('\n\n');
    
    fprintf ('\nMONTAGE BAYESIAN GLM RESULTS AT |BF| > %1.2f, CONTRAST: %s, REGRESSORS: %s, MASK: %s, SCALING: %s\n\n', BF_threshold_glm, analysisname, names_string, mask_string, scaling_string);
    
    o2 = canlab_results_fmridisplay([], 'multirow', num_effects, 'outline', 'linewidth', 0.5,'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]}, 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
    
        for j = 1:num_effects

            BFj = BF(1,j);
            BFj = threshold(BFj, [-2*(log(BF_threshold_glm)) 2*(log(BF_threshold_glm))], 'raw-outside'); 

            o2 = addblobs(o2, region(BFj), 'wh_montages', (2*j)-1:2*j);
            o2 = title_montage(o2, 2*j, [analysisname ' |BF| > ' num2str(BF_threshold_glm) ' ' names{j} ' ' mask_string ' ' scaling_string]);
        
        end % for loop over regressors in model
    
    figtitle = sprintf('%s_%s_%1.4f_unc_montage_%s_%s_%S', analysisname, results_suffix, BF_threshold_glm, names_string, mask_string, scaling_string);
    set(gcf, 'Tag', figtitle, 'WindowState','maximized');
    drawnow, snapnow;
    
        if save_figures_glm
            plugin_save_figure;
        end
        
    clear o2, clear figtitle, clear j, clear BFj
    
    fprintf ('\nTABLES AND MONTAGE REGIONCENTERS BAYESIAN GLM RESULTS AT |BF| > %1.2f, k = %d, CONTRAST: %s, REGRESSORS: %s, MASK: %s, SCALING: %s\n\n', BF_threshold_glm, k_threshold_glm, analysisname, names_string, mask_string, scaling_string);
        
        for j = 1:num_effects

            fprintf ('\nTABLE BAYESIAN GLM RESULTS AT |BF| > %1.2f, k = %d, CONTRAST: %s, REGRESSOR: %s, MASK: %s, SCALING: %s\n\n', BF_threshold_glm, k_threshold_glm, analysisname, names{j}, mask_string, scaling_string);

            BFj = BF(1,j);
            BFj = threshold(BFj, [-2*(log(BF_threshold_glm)) 2*(log(BF_threshold_glm))], 'raw-outside'); 

            r = region(BFj, 'noverbose');
            r(cat(1, r.numVox) < k_threshold_glm) = [];
            [rpos, rneg] = table(r);       % add labels
            r = [rpos rneg];               % re-concatenate labeled regions

            % Montage of regions in table (plot and save)
            if ~isempty(r)
                
                fprintf ('\nMONTAGE REGIONCENTERS BAYESIAN GLM RESULTS AT |BF| > %1.2f, k = %d, CONTRAST: %s, REGRESSOR: %s, MASK: %s, SCALING: %s\n\n', BF_threshold_glm, k_threshold_glm, analysisname, names{j}, mask_string, scaling_string);
                
                o3 = montage(r, 'colormap', 'regioncenters', 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});

                % Activate, name, and save figure
                figtitle = sprintf('%s_%s_%1.4f_unc_regions_%s_%s_%s', analysisname, results_suffix, p_threshold_glm, names{j}, mask_string, scaling_string);
                set(gcf, 'Tag', figtitle, 'WindowState','maximized');
                drawnow, snapnow;
                    if save_figures_glm
                        plugin_save_figure;
                    end
                clear o3, clear figtitle, clear j, clear tj, clear r

            end % if loop regions in results
        
        end % for loop over regressors
        
        
    % BETWEEN-SUBJECT REGRESSORS: MVPA
    % ---------------------------------------------------------------------
    
    if dobootstrap_mvpa_reg_cov
        
        mvpa_num_effects = size(mvpa_results,2);
        
        delete(gcp('nocreate'));
        c = parcluster('local'); % determine local number of cores, and initiate parallel pool with 80% of them
        nw = c.NumWorkers;
        parpool(round(0.8*nw));
        
        for j = 1:mvpa_num_effects
            
            fprintf('\n\n');
            printhdr('BOOTSTRAPPING MVPA WEIGHT MAP');
            fprintf('\n\n');

            t0_boot = tic;

            [~ , mvpa_bs_stats{j}] = predict(mvpa_fmri_dats{j}, 'algorithm_name', algorithm_mvpa_reg_cov,...
                'bootsamples', boot_n_mvpa_reg_cov, 'nfolds', 1, 'error_type', 'mse', ...
                parallelstr_mvpa_reg_cov, 'verbose', 0);

            t_end_boot = toc(t0_boot);
        
        end % for loop over regressors
        
        mvpa_bs_stats_results{c} = mvpa_bs_stats;
        
        fprintf('\n\n');
        printhdr('PLOTTING BOOTSTRAPPED MVPA WEIGHT MAPS');
        fprintf('\n\n');
        
        fprintf ('\nMONTAGE BOOTSTRAPPED %s RESULTS AT FDR q < %1.4f, k = %d, CONTRAST: %s, REGRESSORS: %s, MASK: %s, SCALING: %s\n\n', upper(algorithm_mvpa_reg_cov), q_threshold_mvpa_reg_cov, k_threshold_mvpa_reg_cov, analysisname, names_string, mask_string, scaling_string);
    
        o2 = canlab_results_fmridisplay([], 'multirow', mvpa_num_effects, 'outline', 'linewidth', 0.5, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]}, 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');

            for j = 1:mvpa_num_effects

                tj = mvpa_bs_stats{j}.weight_obj;
                    if apply_mask_before_fdr
                        tj = apply_mask(tj, glmmask);
                    end
                tj = threshold(tj, q_threshold_glm, 'fdr', 'k', k_threshold_glm); 

                o2 = addblobs(o2, region(tj), 'wh_montages', (2*j)-1:2*j);
                o2 = title_montage(o2, 2*j, [analysisname ' FDR ' num2str(q_threshold_glm) ' ' names{j} ' ' mask_string ' ' scaling_string]);

            end % for loop over regressors

        figtitle = sprintf('%s_%s_%s_%1.4f_FDR_montage_%s_%s_%s', analysisname, results_suffix, algorithm_mvpa_reg_cov, q_threshold_glm, names_string, mask_string, scaling_string);
        set(gcf, 'Tag', figtitle, 'WindowState','maximized');
        drawnow, snapnow;

            if save_figures_glm
                plugin_save_figure;
            end

        clear o2, clear figtitle, clear j, clear tj
        
    end % if loop bootstrap
   
end % for loop over contrasts/conditions


%% SAVE RESULTS
% -------------------------------------------------------------------------

if dobootstrap_mvpa_reg_cov
    
    fprintf('\n\n');
    printhdr('SAVING BOOTSTRAPPED MVPA RESULTS');
    fprintf('\n\n');
    
    savefilenamedata_mvpa = fullfile(resultsdir, ['mvpa_stats_and_maps_', mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']);
    save(savefilenamedata_mvpa, 'mvpa_bs_stats','-append');
    fprintf('\nSaved mvpa_stats_results for %s\n', mygroupnamefield);

    fprintf('\nFilename: %s\n', savefilenamedata_mvpa);

end