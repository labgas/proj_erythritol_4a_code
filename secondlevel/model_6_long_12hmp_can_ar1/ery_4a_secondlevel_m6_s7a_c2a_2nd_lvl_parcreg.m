%% ery_4a_secondlevel_m6_s7a_c2a_2nd_lvl_parcreg.m
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
% OPTIONS FOR CURRENT SCRIPT
% 
% NOTE: 
% defaults are specified in a2_set_default_options for any given model,
% but if you want to run the same model with different options, 
% you can make a copy of this script with a letter index (e.g. _s6a_) 
% and change the default option below
%
% GLM options
%
% save_figures_glm: true saves .svg files of all figures generated by c2a_second_level_regression.m (slow, takes up space) % added by @lukasvo76 May 2022
% q_threshold_glm: threshold for FDR-corrected display items
% p_threshold_glm: threshold for uncorrected display items
% k_threshold_glm: extent threshold for both corrected and uncorrected display items
% BF_threshold_glm: threshold for Bayes Factor maps, |BF| > 10 indicates strong evidence in favour of H1 (positive value) or H0 (negative value) - see help.statistic_image.estimateBayesFactor for details
%
% MVPA options
%
% dobootstrap_mvpa_reg_cov = true/false;
    % mvpa bootstrapping options
%     boot_n_mvpa_reg_cov = x;                                      
%     parallelstr_mvpa_reg_cov = 'parallel'/'noparallel';  
%     cons2boot = [y1,y2,...];
    % mvpa thresholding options
%     q_threshold_mvpa_reg_cov = .xx;                                  
%     k_threshold_mvpa_reg_cov = y;
%
% OPTIONS TO COPY FROM CORRESPONDING PREP_3a_ SCRIPT
%
% Mandatory options
%
% - mygroupfieldname: 'contrasts'/'conditions'
% - results_suffix: name added to results file by prep_3a script in case of multiple versions of model, e.g. 'covariate_rating'
%
% Options to copy if specified in prep_3a script
%
% covs2use = {'varname1, varname2,...'};
% group_id = {'varname'};
%
% Custom options from prep_3a script - see first section of this script below for overview
%
%__________________________________________________________________________
%
% revamped by: Lukas Van Oudenhove
% date:   Dartmouth, May, 2022
%
%__________________________________________________________________________
% @(#)% c2a_second_level_regression.m         v4.1
% last modified: 2023/01/20


%% GET AND SET OPTIONS
% -------------------------------------------------------------------------

% GET MODEL-SPECIFIC PATHS AND OPTIONS

ery_4a_secondlevel_m6_s0_a_set_up_paths_always_run_first;

% NOTE: CHANGE THIS TO THE MODEL-SPECIFIC VERSION OF THIS SCRIPT
% NOTE: THIS WILL ALSO AUTOMATICALLY CALL A2_SET_DEFAULT_OPTIONS

% COPY OPTIONS FROM CORRESPONDING PREP_3a_ SCRIPT

% Mandatory options

mygroupnamefield = 'contrasts'; 
results_suffix = 'no_cov'; % suffix of your choice added to .mat file with saved results

% Options to copy if specified in prep_3a script

% covs2use = {''};
% group_id = {'varname'};

% Custom options

% dorobust = true/false;
dorobfit_parcelwise = true;
%   csf_wm_covs = false;
%   remove_outliers = false;
%   atlasname_glm = 'atlas_name';
% maskname_glm = 'mask_name';
% myscaling_glm = 'raw'/'scaled'/'scaled_contrasts';
% design_matrix_type = 'onesample'/'group'/'custom';
% doBayes = true/false;
% domvpa_reg_cov = true/false;
%   algorithm_mvpa_reg_cov = 'cv_pcr';
%   holdout_set_method_mvpa_reg_cov = 'no_group'/'group';
%   nfolds_mvpa_reg_cov = x;
%   zscore_outcome_mvpa_reg_cov = true/false;

% SET CUSTOM OPTIONS FOR CURRENT SCRIPT

% NOTE: only specify if you want to run a second version of your model with different options
% than the defaults you set in your model-specific version of a2_set_default_options.m

% EXCEPTION: cons2boot, which is not in the a2 script, but can be used here
% if you do not want to bootstrap all contrasts

% GLM options

% save_figures_glm = true/false;
% q_threshold_glm = .x;
% p_threshold_glm = .yyy;
% k_threshold_glm = zz;
% BF_threshold_glm = x;

% MVPA options

% dobootstrap_mvpa_reg_cov = true/false;
    % mvpa bootstrapping options
%     boot_n_mvpa_reg_cov = x;                                      
%     parallelstr_mvpa_reg_cov = 'parallel'/'noparallel';   
%     cons2boot = [contrast_indices];
    % mvpa thresholding options
%     q_threshold_mvpa_reg_cov = .yy;                                  
%     k_threshold_mvpa_reg_cov = z; 

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

% NOTE: In case of parcelwise regression, masking is already done in prep_3a script as part of robfit_parcelwise, but we need the atlas object and define mask_string var

if ~dorobfit_parcelwise

    if exist('maskname_glm', 'var') && ~isempty(maskname_glm) && exist(maskname_glm, 'file')

        apply_mask_before_fdr = true;
        [~,maskname_short] = fileparts(maskname_glm);
        mask_string = sprintf('masked with %s', maskname_short);
        glmmask = fmri_mask_image(maskname_glm, 'noverbose'); 

    else

        apply_mask_before_fdr = false;
        mask_string = sprintf('without masking');

    end  
    
end

if exist('atlasname_glm','var') && ~isempty(atlasname_glm) && exist(atlasname_glm, 'file')

    if contains(atlasname_glm,'.mat')
        
        load(atlasname_glm);
        clear mask
        
        if dorobfit_parcelwise
            [~,maskname_short] = fileparts(atlasname_glm);
            mask_string = sprintf('masked with %s', maskname_short);
            
        end
        
        fprintf('\nLabeling regions using custom-made atlas %s\n\n', maskname_short);

    elseif ischar(atlasname_glm)

        combined_atlas = load_atlas(atlasname_glm);
        
        if dorobfit_parcelwise
            mask_string = sprintf('in atlas %s',atlasname_glm);
            
        end
        
        fprintf('\nLabeling regions using custom-made atlas %s\n\n', maskname_short);

    else

         error('\ninvalid option "%s" defined in atlasname_glm variable, should be a keyword for load_atlas.m or a .mat file containing an atlas object, check docs"\n\n',atlasname_glm)

    end

else

    if dorobfit_parcelwise
        mask_string = sprintf('without masking');
        
        fprintf('\nShowing parcelwise results without masking in 489 parcels of canlab_2018 atlas\n\n');
        
    end

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

if dobootstrap_mvpa_reg_cov
    
    mvpa_resultsvarname = 'mvpa_stats_results';
    mvpa_resultsstring = 'mvpa_stats_and_maps_';
    mvpa_analysis_type = algorithm_mvpa_reg_cov;

    if ~exist(mvpa_resultsvarname,'var')

        fprintf('\n\n');
        printhdr('LOADING MVPA DATA');
        fprintf('\n\n');

        savefilenamedata_mvpa = fullfile(resultsdir, ['mvpa_stats_and_maps_', mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']);

        if exist(savefilenamedata_mvpa,'file')
            fprintf('\nLoading %s mvpa regression results and maps from %s\n\n', mvpa_analysis_type, savefilenamedata_mvpa);
            load(savefilenamedata_mvpa);
        else
            fprintf('\nNo saved results file %s. Skipping this analysis.', savefilenamedata_mvpa);
            fprintf('\nRun prep_3a_run_second_level_regression_and_save.m to get %s regression results first.\n', mvpa_analysis_type); 
            return
        end

    else
        fprintf('\n%s %s found, displaying results\n\n', mvpa_resultsvarname, mvpa_analysis_type);

    end

end
    
%% VISUALIZE GLM RESULTS FOR EACH CONTRAST
% -------------------------------------------------------------------------
    
region_objs_fdr = cell(1,size(results,2));
region_tables_fdr = cell(1,size(results,2));

region_objs_unc = cell(1,size(results,2));
region_tables_unc = cell(1,size(results,2));

    if doBayes
        region_objs_Bayes = cell(1,size(results,2));
        region_tables_Bayes = cell(1,size(results,2));
    end

    
for c = 1:size(results, 2) % number of contrasts or conditions

    analysisname = results{c}.contrastname;
    names = results{c}.variable_names;
    
    names_string = names{1};
    
        for name = 2:size(names,2)
            names_string = [names_string,' ',names{name}];
        end
    
        if ~dorobfit_parcelwise
            
            t = results{c}.t; % NOTE: this statistic_object is unthresholded AND unmasked
            glmmask = resample_space(glmmask,t);
            
                if apply_mask_before_fdr
                    t = apply_mask(t, glmmask);
                end
            
        else
            
            t = results{c}.t_obj; % NOTE: this statistic_object is thresholded at FDR q < 0.05 AND masked!
            
        end
        
    if doBayes
        
        BF = results{c}.BF; % NOTE: unthresholded in both cases, masked in case of parcelwise only
        
        if ~dorobfit_parcelwise
            
            if apply_mask_before_fdr
                
                for j = 1:size(BF,2)
                    BF(j) = apply_mask(BF(j), glmmask);
                    
                end
            end
            
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
    
    o2 = canlab_results_fmridisplay([], 'multirow', num_effects, 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
    
        for j = 1:num_effects

            tj = get_wh_image(t, j);
            
                if ~dorobfit_parcelwise
                    
                    tj = threshold(tj, q_threshold_glm, 'fdr', 'k', k_threshold_glm); 
                    
                else
                    
                    if q_threshold_glm ~= .05 % already thresholded at q < 0.05
                        
                        tj = threshold(tj, q_threshold_glm, 'fdr', 'k', k_threshold_glm); 
                        
                    end
                        
                end

            o2 = addblobs(o2, region(tj), 'wh_montages', (2*j)-1:2*j, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
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
    
        region_fdr = cell(1,num_effects);
        table_fdr = cell(1,num_effects);
        
        for j = 1:num_effects

            fprintf ('\nTABLE GLM RESULTS AT FDR q < %1.4f, k = %d, CONTRAST: %s, REGRESSOR: %s, MASK: %s, SCALING: %s\n\n', q_threshold_glm, k_threshold_glm, analysisname, names{j}, mask_string, scaling_string);

            tj = get_wh_image(t, j);
            
                if ~dorobfit_parcelwise
                    
                    tj = threshold(tj, q_threshold_glm, 'fdr', 'k', k_threshold_glm); 
                    
                else
                    
                    if q_threshold_glm ~= .05
                        
                        tj = threshold(tj, q_threshold_glm, 'fdr', 'k', k_threshold_glm); 
                        
                    end
                        
                end

            r = region(tj, 'noverbose');
            r(cat(1, r.numVox) < k_threshold_glm) = [];
            
            if ~isempty(r)
            
                if exist('combined_atlas','var')
                    [rpos, rneg, r_table] = table(r,'atlas_obj',combined_atlas); % add labels from combined_atlas
                    r = [rpos rneg];                                    % re-concatenate labeled regions
                else
                    [rpos, rneg, r_table] = table(r);                            % add labels from default canlab_2018 atlas 
                    r = [rpos rneg];                                    % re-concatenate labeled regions
                end
            
            region_fdr{j} = r;
            table_fdr{j} = r_table;

            % Montage of regions in table (plot and save)
                
                fprintf ('\nMONTAGE REGIONCENTERS GLM RESULTS AT FDR q < %1.4f, k = %d, CONTRAST: %s, REGRESSOR: %s, MASK: %s, SCALING: %s\n\n', q_threshold_glm, k_threshold_glm, analysisname, names{j}, mask_string, scaling_string);
                
                o3 = montage(r, 'colormap', 'regioncenters', 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});

                % Activate, name, and save figure
                figtitle = sprintf('%s_%s_%1.4f_FDR_regions_%s_%s_%s', analysisname, results_suffix, q_threshold_glm, names{j}, mask_string, scaling_string);
                set(gcf, 'Tag', figtitle, 'WindowState','maximized');
                drawnow, snapnow;
                    if save_figures_glm
                        plugin_save_figure;
                    end
                clear o3, clear figtitle, clear j, clear tj, clear r, clear r_table

            end % conditional montage plot if there are regions to show
            
        end % for loop over regressors in model
        
        region_objs_fdr{c} = region_fdr;
        region_tables_fdr{c} = table_fdr;
            
   
    % BETWEEN-SUBJECT REGRESSORS & INTERCEPT: uncorrected
    % ---------------------------------------------------------------------
    
    fprintf('\n\n');
    printhdr('UNCORRECTED GLM RESULTS');
    fprintf('\n\n');
    
    fprintf ('\nMONTAGE GLM RESULTS AT UNCORRECTED p < %1.4f, k = %d, CONTRAST: %s, REGRESSORS: %s, MASK: %s, SCALING: %s\n\n', p_threshold_glm, k_threshold_glm, analysisname, names_string, mask_string, scaling_string);
    
    o2 = canlab_results_fmridisplay([], 'multirow', num_effects, 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
    
        for j = 1:num_effects

            tj = get_wh_image(t, j);
            
            tj = threshold(tj, p_threshold_glm, 'unc', 'k', k_threshold_glm); 

            o2 = addblobs(o2, region(tj), 'wh_montages', (2*j)-1:2*j, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
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
        
        region_unc = cell(1,num_effects);
        table_unc = cell(1,num_effects);
    
        for j = 1:num_effects

            fprintf ('\nTABLE GLM RESULTS AT UNCORRECTED p < %1.4f, k = %d, CONTRAST: %s, REGRESSOR: %s, MASK: %s, SCALING: %s\n\n', p_threshold_glm, k_threshold_glm, analysisname, names{j}, mask_string, scaling_string);

            tj = get_wh_image(t, j);

            tj = threshold(tj, p_threshold_glm, 'unc', 'k', k_threshold_glm); 

            r = region(tj, 'noverbose');
            r(cat(1, r.numVox) < k_threshold_glm) = [];
            
            if ~isempty(r)
            
                if exist('combined_atlas','var')
                    [rpos, rneg, r_table] = table(r,'atlas_obj',combined_atlas); % add labels from combined_atlas
                    r = [rpos rneg];                                    % re-concatenate labeled regions
                else
                    [rpos, rneg, r_table] = table(r);                            % add labels from default canlab_2018 atlas 
                    r = [rpos rneg];                                    % re-concatenate labeled regions
                end       
            
                region_unc{j} = r;
                table_unc{j} = r_table;

                % Montage of regions in table (plot and save)
                
                fprintf ('\nMONTAGE REGIONCENTERS GLM RESULTS AT UNCORRECTED p < %1.4f, k = %d, CONTRAST: %s, REGRESSOR: %s, MASK: %s, SCALING: %s\n\n', p_threshold_glm, k_threshold_glm, analysisname, names{j}, mask_string, scaling_string);
                
                o3 = montage(r, 'colormap', 'regioncenters', 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});

                % Activate, name, and save figure
                figtitle = sprintf('%s_%s_%1.4f_unc_regions_%s_%s_%s', analysisname, results_suffix, p_threshold_glm, names{j}, mask_string, scaling_string);
                set(gcf, 'Tag', figtitle, 'WindowState','maximized');
                drawnow, snapnow;
                    if save_figures_glm
                        plugin_save_figure;
                    end
                clear o3, clear figtitle, clear j, clear tj, clear r, clear r_table

            end % conditional montage plot if there are regions to show
        
        end % loop over regressors
        
    region_objs_unc{c} = region_unc;
    region_tables_unc{c} = table_unc;
        
        
    % BETWEEN-SUBJECT REGRESSORS & INTERCEPT: Bayesian
    % ---------------------------------------------------------------------
    
    if doBayes
    
        fprintf('\n\n');
        printhdr('BAYESIAN GLM RESULTS');
        fprintf('\n\n');

        fprintf ('\nMONTAGE BAYESIAN GLM RESULTS AT |BF| > %1.2f, CONTRAST: %s, REGRESSORS: %s, MASK: %s, SCALING: %s\n\n', BF_threshold_glm, analysisname, names_string, mask_string, scaling_string);

        o2 = canlab_results_fmridisplay([], 'multirow', num_effects, 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img'); 

            for j = 1:num_effects

                BFj = BF(1,j);
                
                BFj = threshold(BFj, [-2*(log(BF_threshold_glm)) 2*(log(BF_threshold_glm))], 'raw-outside'); 

                o2 = addblobs(o2, region(BFj), 'wh_montages', (2*j)-1:2*j, 'splitcolor',{[.25 0 0] [1 0 0] [0 0.25 0] [0 1 0]}); % red in favor of H0, green in favor of H1 for BF maps
                o2 = title_montage(o2, 2*j, [analysisname ' |BF| > ' num2str(BF_threshold_glm) ' ' names{j} ' ' mask_string ' ' scaling_string]);

            end % for loop over regressors in model

        figtitle = sprintf('%s_%s_%1.4f_unc_montage_%s_%s_%s', analysisname, results_suffix, BF_threshold_glm, names_string, mask_string, scaling_string);
        set(gcf,'Tag', figtitle, 'WindowState','maximized');
        drawnow, snapnow;

            if save_figures_glm
                plugin_save_figure;
            end

        clear o2, clear figtitle, clear j, clear BFj

        fprintf ('\nTABLES AND MONTAGE REGIONCENTERS BAYESIAN GLM RESULTS AT |BF| > %1.2f, k = %d, CONTRAST: %s, REGRESSORS: %s, MASK: %s, SCALING: %s\n\n', BF_threshold_glm, k_threshold_glm, analysisname, names_string, mask_string, scaling_string);

            region_Bayes = cell(1,num_effects);
            table_Bayes = cell(1,num_effects);
        
            for j = 1:num_effects

                fprintf ('\nTABLE BAYESIAN GLM RESULTS AT |BF| > %1.2f, k = %d, CONTRAST: %s, REGRESSOR: %s, MASK: %s, SCALING: %s\n\n', BF_threshold_glm, k_threshold_glm, analysisname, names{j}, mask_string, scaling_string);

                BFj = BF(1,j);
                BFj = threshold(BFj, [-2*(log(BF_threshold_glm)) 2*(log(BF_threshold_glm))], 'raw-outside'); 

                r = region(BFj, 'noverbose');  
                r(cat(1, r.numVox) < k_threshold_glm) = [];
                
                if ~isempty(r)
                
                    if exist('combined_atlas','var')
                        [rpos, rneg, r_table] = table(r,'atlas_obj',combined_atlas); % add labels from combined_atlas
                        r = [rpos rneg];                                    % re-concatenate labeled regions
                    else
                        [rpos, rneg, r_table] = table(r);                            % add labels from default canlab_2018 atlas 
                        r = [rpos rneg];                                    % re-concatenate labeled regions
                    end
                    
                    region_Bayes{j} = r;
                    table_Bayes{j} = r_table;

                    % Montage of regions in table (plot and save)

                    fprintf ('\nMONTAGE REGIONCENTERS BAYESIAN GLM RESULTS AT |BF| > %1.2f, k = %d, CONTRAST: %s, REGRESSOR: %s, MASK: %s, SCALING: %s\n\n', BF_threshold_glm, k_threshold_glm, analysisname, names{j}, mask_string, scaling_string);

                    o3 = montage(r, 'colormap', 'regioncenters', 'splitcolor',{[.25 0 0] [1 0 0] [0 0.25 0] [0 1 0]});

                    % Activate, name, and save figure
                    figtitle = sprintf('%s_%s_%1.4f_unc_regions_%s_%s_%s', analysisname, results_suffix, p_threshold_glm, names{j}, mask_string, scaling_string);
                    set(gcf, 'Tag', figtitle, 'WindowState','maximized');
                    drawnow, snapnow;
                        if save_figures_glm
                            plugin_save_figure;
                        end
                    clear o3, clear figtitle, clear j, clear tj, clear r, clear r_table

                end % conditional montage plot if there are regions to show

            end % for loop over regressors
            
        region_objs_Bayes{c} = region_Bayes;
        region_tables_Bayes{c} = table_Bayes;
            
    end % if loop doBayes
    
        
    % BETWEEN-SUBJECT REGRESSORS: MVPA
    % ---------------------------------------------------------------------
    
    if dobootstrap_mvpa_reg_cov
        
        mvpa_bs_stats_results = cell(1,size(results,2));
        
        for covar = 1:size(mvpa_stats_results,2)
            mvpa_results{covar} = mvpa_stats_results{c,covar};
            mvpa_fmri_dats{covar} = mvpa_dats{c,covar};
        end
        
        if isempty(cons2boot) || ismember(c,cons2boot)
        
            mvpa_num_effects = size(mvpa_results,2);

            mvpa_bs_stats = cell(1,mvpa_num_effects);

            if strcmp(parallelstr_mvpa_reg_cov,'parallel')

                delete(gcp('nocreate'));
                c = parcluster('local'); % determine local number of cores, and initiate parallel pool with 80% of them
                nw = c.NumWorkers;
                parpool(round(0.8*nw));

            end

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

            o2 = canlab_results_fmridisplay([], 'multirow', mvpa_num_effects, 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');

                for j = 1:mvpa_num_effects

                    tj = mvpa_bs_stats{j}.weight_obj;
                        if apply_mask_before_fdr
                            tj = apply_mask(tj, glmmask);
                        end
                    tj = threshold(tj, q_threshold_glm, 'fdr', 'k', k_threshold_glm); 

                    o2 = addblobs(o2, region(tj), 'wh_montages', (2*j)-1:2*j, 'splitcolor', {[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
                    o2 = title_montage(o2, 2*j, [analysisname ' FDR ' num2str(q_threshold_glm) ' ' names{j} ' ' mask_string ' ' scaling_string]);

                end % for loop over regressors

            figtitle = sprintf('%s_%s_%s_%1.4f_FDR_montage_%s_%s_%s', analysisname, results_suffix, algorithm_mvpa_reg_cov, q_threshold_glm, names_string, mask_string, scaling_string);
            set(gcf, 'Tag', figtitle, 'WindowState','maximized');
            drawnow, snapnow;

                if save_figures_glm
                    plugin_save_figure;
                end

            clear o2, clear figtitle, clear j, clear tj
        
        end % if loop exist cons2boot
        
    end % if loop bootstrap
   
end % for loop over contrasts/conditions


%% SAVE RESULTS
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('APPENDING REGION OBJECTS TO SAVED RESULTS');
fprintf('\n\n');
    
    if ~dorobfit_parcelwise
    
        savefilenamedata_region = fullfile(resultsdir, ['regression_stats_and_maps_', mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']);
        
        if ~doBayes
            save(savefilenamedata_region, 'region_objs_unc', 'region_objs_fdr', 'region_tables_unc', 'region_tables_fdr','-append');
            
        else
            save(savefilenamedata_region, 'region_objs_unc', 'region_objs_fdr', 'region_objs_Bayes', 'region_tables_unc', 'region_tables_fdr', 'region_tables_Bayes', '-append');
            
        end
        
    else
        
        savefilenamedata_region = fullfile(resultsdir, ['parcelwise_stats_and_maps_', mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']);
        
        if ~doBayes
            
                save(savefilenamedata_region, 'region_objs_unc', 'region_objs_fdr', 'region_objs_Bayes', '-append');
            
        else
            
                save(savefilenamedata_region, 'region_objs_unc', 'region_objs_Bayes', 'region_objs_fdr', 'region_tables_unc', 'region_tables_fdr', 'region_tables_Bayes','-append');
                  
        end
        
    end
        
        
    fprintf('\nAdded region objects for %s\n', mygroupnamefield);

    fprintf('\nFilename: %s\n', savefilenamedata_region);
    

if dobootstrap_mvpa_reg_cov
    
    fprintf('\n\n');
    printhdr('APPENDING BOOTSTRAPPED MVPA RESULTS');
    fprintf('\n\n');
    
    savefilenamedata_mvpa = fullfile(resultsdir, ['mvpa_stats_and_maps_', mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']);
    save(savefilenamedata_mvpa, 'mvpa_bs_stats_results','-append');
    fprintf('\nSaved mvpa_stats_results for %s\n', mygroupnamefield);

    fprintf('\nFilename: %s\n', savefilenamedata_mvpa);

end