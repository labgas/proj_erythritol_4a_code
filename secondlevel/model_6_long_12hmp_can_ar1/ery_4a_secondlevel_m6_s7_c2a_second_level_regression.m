%%% ery_4a_secondlevel_m6_s7_c2a_second_level_regression.m

% USAGE
%
% This script displays second⁻level (i.e. across subjects) regression
% results generated by prep_3a_run_second_level_regression_and_save.m
%
% Use Matlab's publish function to generate html report
% 
% See the documentation of prep_3a_run_second_level_regression_and_save.m
% for more info and options
%
%__________________________________________________________________________
%
% revamped by: Lukas Van Oudenhove
% date:   Dartmouth, May, 2022
%
%__________________________________________________________________________
% @(#)% c2a_second_level_regression.m         v2.0
% last modified: 2022/05/28


%% LOAD REGRESSION RESULTS IF NEEDED
%--------------------------------------------------------------------------

% options (from corresponding prep_3a script)

mygroupnamefield = 'contrasts';
results_suffix = 'OLS';

% check scaling

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

% check whether results are available, load if needed

if ~dorobfit_parcelwise
    resultsvarname = 'regression_stats_results';
    resultsstring = 'regression_stats_and_maps_';
    analysis_type = 'voxel-wise';
else
    resultsvarname = 'parcelwise_stats_results';
    resultsstring = 'parcelwise_stats_and_maps';
    analysis_type = 'parcel-wise';
end
    
if ~exist(resultsvarname, 'var')
    savefilenamedata = fullfile(resultsdir, [resultsstring, mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']);
        if exist(savefilenamedata,'file')
            fprintf('\nLoading %s regression results and maps from %s\n\n', analysis_type, savefilenamedata);
            load(savefilenamedata, resultsvarname);
        else
            fprintf('\nNo saved results file %s. Skipping this analysis.', savefilenamedata)
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


%% MASKING
%--------------------------------------------------------------------------

if exist(maskname_glm, 'file')
    apply_mask_before_fdr = true;
    [~,maskname_short] = fileparts(maskname_glm);
    mask_string = sprintf('within_%s', maskname_short);
    mask = fmri_data_st(maskname_glm, 'noverbose'); 
else
    apply_mask_before_fdr = false;
    mask_string = sprintf('without_masking');
end  


%% RUN MASS UNIVARIATE GLM
%--------------------------------------------------------------------------

for c = 1:size(results, 2) % number of contrasts or conditions

    analysisname = results{c}.analysis_name;
    names = results{c}.variable_names;
    
        if ~dorobfit_parcelwise
            t = results{c}.t;
        else
            t = results{c}.t_obj;
        end
    
    printhdr(analysisname)
    disp('Regressors: ')
    disp(names)
    
        if isfield(results{c}, 'design_table')
            disp(results{c}.design_table);
        end
    
    num_effects = size(t.dat, 2); % number of regressors
    
    % BETWEEN-SUBJECT REGRESSORS & INTERCEPT: FDR corrected
    % ---------------------------------------------------------------------
    o2 = canlab_results_fmridisplay([], 'multirow', num_effects, 'outline', 'linewidth', 0.5, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]}, 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
    
        for j = 1:num_effects

            fprintf ('\nShowing results at FDR q < %1.4f: %s\nEffect: %s, %s\n\n', q_threshold, analysisname, names{j}, mask_string);

            tj = get_wh_image(t, j);
                if apply_mask_before_fdr
                    tj = apply_mask(tj, mask);
                end
            tj = threshold(tj, q_threshold, 'fdr', 'k', k_threshold); 

            o2 = addblobs(o2, region(tj), 'wh_montages', (2*j)-1:2*j);
            o2 = title_montage(o2, 2*j, [analysisname ' FDR ' num2str(q_threshold) ' ' names{j} ' ' mask_string]);

        end % for loop over regressors in model
    
    figtitle = sprintf('%s_%s_%1.4f_FDR_montage_%s_%s', analysisname, results_suffix, q_threshold, scaling_string, mask_string);
    set(gcf, 'Tag', figtitle, 'WindowState','maximized');
    drawnow, snapnow;
    
        if save_figures
            plugin_save_figure;
        end
        
    clear o2, clear figtitle
    
        for j = 1:num_effects

            fprintf('\n\nTable of results for clusters %d contiguous voxels.\n', k_threshold);

            tj = get_wh_image(t, j);
                if apply_mask_before_fdr
                    tj = apply_mask(tj, mask);
                end
            tj = threshold(tj, q_threshold, 'fdr', 'k', k_threshold); 

            r = region(tj, 'noverbose');
            r(cat(1, r.numVox) < k_threshold) = [];
            [rpos, rneg] = table(r);       % add labels
            r = [rpos rneg];               % re-concatenate labeled regions

            % Montage of regions in table (plot and save)
            if ~isempty(r)
                o3 = montage(r, 'colormap', 'regioncenters', 'outline', 'linewidth', 0.5, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});

                % Activate, name, and save figure - then close
                figtitle = sprintf('%s_%s_%1.4f_FDR_regions_%s_%s_%s', analysisname, results_suffix, q_threshold, names{j}, scaling_string, mask_string);
                region_fig_han = activate_figures(o3);
                
                    if ~isempty(region_fig_han)
                        set(region_fig_han(1), 'Tag', figtitle, 'WindowState','maximized');
                        drawnow, snapnow;
                            if save_figures
                                plugin_save_figure;
                            end
                        close(region_fig_han(1)), clear o3, clear figtitle
                    else
                        fprintf('\n');
                        warning('Cannot find figure - Tag field was not set or figure was closed. Skipping save operation.');
                        fprintf('\n');
                    end

            end % conditional montage plot if there are regions to show
            
        end % for loop over regressors in model

    
    % BETWEEN-SUBJECT REGRESSORS & INTERCEPT: uncorrected
    % ---------------------------------------------------------------------    
    o2 = canlab_results_fmridisplay([], 'multirow', num_effects, 'outline', 'linewidth', 0.5, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]}, 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
    
        for j = 1:num_effects

            fprintf ('\nShowing results at uncorrected p < %1.4f: %s\nEffect: %s, %s\n\n', p_threshold, analysisname, names{j}, mask_string);

            tj = get_wh_image(t, j);
                if apply_mask_before_fdr
                    tj = apply_mask(tj, mask);
                end
            tj = threshold(tj, p_threshold, 'unc', 'k', k_threshold); 

            o2 = addblobs(o2, region(tj), 'wh_montages', (2*j)-1:2*j);
            o2 = title_montage(o2, 2*j, [analysisname ' unc ' num2str(p_threshold) ' ' names{j} ' ' mask_string]);
        
        end % for loop over regressors in model
    
    figtitle = sprintf('%s_%s_%1.4f_unc_montage_%s_%s', analysisname, results_suffix, p_threshold, scaling_string, mask_string);
    set(gcf, 'Tag', figtitle, 'WindowState','maximized');
    drawnow, snapnow;
    
        if save_figures
            plugin_save_figure;
        end
        
    clear o2, clear figtitle
        
        for j = 1:num_effects

            fprintf('\n\nTable of results for clusters >= %d contiguous voxels.\n', k_threshold);

            tj = get_wh_image(t, j);
                if apply_mask_before_fdr
                    tj = apply_mask(tj, mask);
                end
            tj = threshold(tj, p_threshold, 'unc', 'k', k_threshold); 

            r = region(tj, 'noverbose');
            r(cat(1, r.numVox) < k_threshold) = [];
            [rpos, rneg] = table(r);       % add labels
            r = [rpos rneg];               % re-concatenate labeled regions

            % Montage of regions in table (plot and save)
            if ~isempty(r)
                o3 = montage(r, 'colormap', 'regioncenters', 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});

                % Activate, name, and save figure - then close
                figtitle = sprintf('%s_%s_%1.4f_unc_regions_%s_%s_%s', analysisname, results_suffix, p_threshold, names{j}, scaling_string, mask_string);
                region_fig_han = activate_figures(o3);
                
                    if ~isempty(region_fig_han)
                        set(region_fig_han{1}, 'Tag', figtitle, 'WindowState','maximized');
                        drawnow, snapnow;
                            if save_figures
                                plugin_save_figure;
                            end
                        close(region_fig_han{1}), clear o3, clear figtitle
                    else
                        fprintf('\n');
                        warning('Cannot find figure - Tag field was not set or figure was closed. Skipping save operation.');
                        fprintf('\n');
                    end

            end % loop over regions in results
        
        end % loop over regressors
   
end % loop over contrasts/conditions
