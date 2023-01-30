%% ery_4a_secondlevel_m6_s11_c2_SVM_contrasts_masked.m
%
%
% USAGE
%
% This script displays second⁻level (i.e. across subjects) support vector
% machine (SVM) results generated by prep_3c_run_SVMs_on_contrasts_masked.m
%
% Run this script with Matlab's publish function to generate html report of results:
% publish('c2_SVM_contrasts_masked','outputDir',htmlsavedir)
% 
% See the documentation of prep_3c_run_SVMs_on_contrasts_masked.m
% for more info and options
%
% OPTIONS FOR CURRENT SCRIPT
%
% NOTE: 
% defaults are specified in a2_set_default_options for any given model,
% but if you want to run the same model with different options, 
% you can make a copy of this script with a letter index (e.g. _s6a_) 
% and change the default options below
%
% save_figures_svm: true saves .svg files of all figures generated by c2a_second_level_regression.m (slow, takes up space) % added by @lukasvo76 May 2022
% q_threshold_svm: threshold for FDR-corrected display items
% p_threshold_svm: threshold for uncorrected display items
% k_threshold_svm: extent threshold for both corrected and uncorrected display items
%
% OPTIONS TO COPY FROM CORRESPONDING PREP_3c_ SCRIPT
%
% Mandatory options
%
% - results_suffix: name added to results file by prep_3c script in case of multiple versions of model, e.g. 'covariate_rating'
%
% Options to copy if specified in prep_3c script
%
% - myscaling_svm: 'raw'/'subjectnorm'/'imagenorm'/'zscoreimages'/'zscorevoxels'
% - maskname_svm: which('')
% - dobootstrap_svm: true/false
%       cons2boot_svm: []; vector of indices for contrasts to bootstrap if you want to bootstrap a subset
% - dosearchlight_svm: true/false
%       cons2searchlight_svm: vector of indices for contrasts to perform a searchlight analysis on if you only want to do this in a subset
%
%__________________________________________________________________________
%
% revamped by: Lukas Van Oudenhove
% date:   KU Leuven, July, 2022
%
%__________________________________________________________________________
% @(#)% c2_SVM_contrasts_masked.m         v3.2
% last modified: 2022/09/02


%% GET AND SET OPTIONS
% -------------------------------------------------------------------------

% GET MODEL-SPECIFIC PATHS AND OPTIONS

ery_4a_secondlevel_m6_s0_a_set_up_paths_always_run_first;

% NOTE: CHANGE THIS TO THE MODEL-SPECIFIC VERSION OF THIS SCRIPT
% NOTE: THIS WILL ALSO AUTOMATICALLY CALL A2_SET_DEFAULT_OPTIONS

% COPY OPTIONS FROM CORRESPONDING PREP_3c_ SCRIPT

% Mandatory

results_suffix = ''; % suffix added to .mat file with saved results

% Options to copy if specified in prep_3a script

% myscaling_svm: 'raw'/'subjectnorm'/'imagenorm'/'zscoreimages'/'zscorevoxels'
% maskname_svm: which('maskname');
dobootstrap_svm = true;
      cons2boot_svm = [4:6];
% - dosearchlight_svm: true/false
%       cons2searchlight_svm: [x y z]

% SET CUSTOM OPTIONS FOR THIS SCRIPT

% NOTE: only specify if you want to run a second version of your model with different options
% than the defaults you set in your model-specific version of a2_set_default_options.m

% save_figures_svm = true/false;
% q_threshold_svm = .x;
% p_threshold_svm = .yyy;
% k_threshold_svm = zz;

% GET SCALING OPTION

switch myscaling_svm

    case 'raw'

        scaling_string = 'no_scaling';

    case 'subjectnorm'

        scaling_string = 'scaling_l2norm_subjects';
        fprintf('\nNormalized condition images for each subject by L2 norm of Condition 1 image before SVM\n\n')

    case 'imagenorm'

        scaling_string = 'scaling_l2norm_conditions';
        fprintf('\nNormalized each condition image for each subject by L2 norm before SVM\n\n')

    case 'zscoreimages'
        fprintf('\nZ-scored each condition image for each subject before SVM\n\n');
        scaling_string = 'scaling_z_score_conditions';

    case 'zscorevoxels'
        fprintf('\nZ-scored each condition image for each subject before SVM\n\n');
        scaling_string = 'scaling_z_score_conditions';

    otherwise
        error('\nincorrect scaling option %s specified in myscaling_svm option in a2_set_default_options.\nChoose between "raw", "subjectnorm", "imagenorm", "zscoreimages", or "zscorevoxels"\n\n', myscaling_svm);

end
    
% GET MASKING OPTION

if exist(maskname_svm, 'file')
    [~,maskname_short] = fileparts(maskname_svm);
    mask_string = sprintf('masked_with_%s', maskname_short);
else
    mask_string = sprintf('without_masking');
end 
    

%% LOAD NECESSARY VARIABLES IF NEEDED
%--------------------------------------------------------------------------

if ~exist('DSGN','var') || ~exist('DAT','var')
    
    load(fullfile(resultsdir,'image_names_and_setup.mat'));
    
end

if ~exist('DATA_OBJ','var')
    
    load(fullfile(resultsdir,'data_objects.mat'));
    
end


%% LOAD SVM RESULTS IF NEEDED
% -------------------------------------------------------------------------

resultsvarname = 'svm_stats_results';
resultsstring = 'svm_stats_results_contrasts_';

if ~exist(resultsvarname,'var')
    
    fprintf('\n\n');
    printhdr('LOADING DATA');
    fprintf('\n\n');
    
    if exist('maskname_short', 'var')
        savefilenamedata = fullfile(resultsdir, [resultsstring, scaling_string, '_', maskname_short, '_', results_suffix, '.mat']);
    else
        savefilenamedata = fullfile(resultsdir, [resultsstring, scaling_string, '_', results_suffix, '.mat']);
    end
    
    if exist(savefilenamedata, 'file')
        fprintf('\nLoading SVM results from %s\n\n', savefilenamedata);
        load(savefilenamedata);
    else
        fprintf('\nNo saved SVM results file %s. Skipping this analysis.\n\n', savefilenamedata)
        fprintf('\nRun prep_3c_run_SVMs_on_contrasts_masked.m to get svm results first.\n\n'); 
        return
    end
    
else
    fprintf('\n%s found, displaying results\n\n', resultsvarname);

end 

        
%% AVERAGE IMAGES FROM THE SAME PERSON WITHIN EACH SVM CLASS IF NEEDED
% -------------------------------------------------------------------------

% - This plugin function calculates (a) cross-validated distances from the
% SVM hyerplane, and (b) a cross-classification matrix across contrasts
% - It averages images within the same subject and condition (+ or -,
% on/off in the SVM analysis) into one image for testing purposes, so that
% it performs a subject-wise forced choice classification
% - It assumes that the image lists for each condition contain images for subjects 1:n
% in each condition, in the same order.

% The purpose of this plugin is to average over replicates of images with the
% same outcome collected on the same individuals.  We want one average
% image for outcome 1 and one average image for outcome -1 per person, and
% we want to test person-wise classification.
%
% Depending on how contrasts are specified, stats results structure from SVM training
% may involve multiple images coded
% with 1 or -1 per person, in which case the svm_stats_results stats
% structures will contain image-wise classification results, not
% person-wise averaging over images for each person.
% This function calculates an average pattern expression value (distance from hyperplane)
% for each person for images coded as 1 and those coded as -1.
%
% For example, if you are testing a [1 1 1 1 -1 -1 -1 -1] contrast, you
% will end up with 8 images per person in the SVM analysis, 4 for each
% condition. You want to test the average of the first 4 vs. the average of
% the last 4 when you caculate test accuracy.

switch ml_method_svm
    
    case 'predict'
        [dist_from_hyperplane, Y, svm_dist_pos_neg, svm_dist_pos_neg_matrix] = plugin_svm_contrasts_get_results_per_subject(DAT, svm_stats_results, DATA_OBJ);

    case 'oofmridataobj'
        [dist_from_hyperplane, Y, svm_dist_pos_neg, svm_dist_pos_neg_matrix] = plugin_oofmridataobj_svm_contrasts_get_results_per_subject(DAT, svm_stats_results, DATA_OBJ);
        
end


%% CHECK WHETHER WE HAVE PAIRED IMAGES AND ERROR OUT IF NOT
% -------------------------------------------------------------------------

kc = size(DAT.contrasts, 1);

ispaired = false(1, kc);

for i = 1:kc
    ispaired(i) = sum(Y{i} > 0) == sum(Y{i} < 0);
end

if ~all(ispaired)
    fprintf('\n');
    error('This script should only be run on paired, within-person contrasts. Check images and results. Skipping this analysis.');
end


%% DEFINE EFFECT SIZE FUNCTIONS AND ROC TYPE
% -------------------------------------------------------------------------

% Define paired and uppaired functions here for reference
% This script uses the paired option because it runs within-person
% contrasts

% ROC plot is different for paired samples and unpaired. Paired samples
% must be in specific order, 1:n for condition 1 and 1:n for condition 2.
% If samples are paired, this is set up by default in these scripts.
% But some contrasts entered by the user may be unbalanced, i.e., different
% numbers of images in each condition, unpaired. Other SVM scripts are set up
% to handle this condition explicitly and run the unpaired version.  

% Effect size, cross-validated, paired samples
dfun_paired = @(x, Y) mean(x(Y > 0) - x(Y < 0)) ./ std(x(Y > 0) - x(Y < 0));

% Effect size, cross-validated, unpaired sampled
dfun_unpaired = @(x, Y) (mean(x(Y > 0)) - mean(x(Y < 0))) ./ sqrt(var(x(Y > 0)) + var(x(Y < 0))); % check this. @lukasvo76 appears unused here hence we may not need it, but does not harm

rocpairstring = 'twochoice';  % 'twochoice' or 'unpaired'


%% CROSS-VALIDATED ACCURACY, ROC PLOTS, AND MONTAGES FOR EACH CONTRAST
% -------------------------------------------------------------------------

for c = 1:kc
    
    analysisname = DAT.contrastnames{c};
    
    fprintf('\n\n');
    printhdr(['CONTRAST #', num2str(c), ': ', upper(analysisname)]);
    fprintf('\n\n');
    
    if isempty(dist_from_hyperplane{c})
        
        warning('\nContrast %s not suitable for SVM. Skipping this contrast and continuing.\n\n', DAT.contrastnames{c});

        continue
        
    end
    
    % ROC PLOT
    % --------------------------------------------------------------------
    
    figtitle = sprintf('SVM ROC %s', upper(analysisname));
    create_figure(figtitle);
    set(gcf,'WindowState','maximized');
    
    ROC = roc_plot(dist_from_hyperplane{c}, logical(Y{c} > 0), 'color', DAT.contrastcolors{c}, rocpairstring);
    
    d_paired = dfun_paired(dist_from_hyperplane{c}, Y{c});
    fprintf('\nEffect size, cross-validated: Forced choice: d = %3.2f\n\n', d_paired);
    
    fprintf('\n\n');
    printhdr('ROC plot');
    fprintf('\n\n');
    drawnow, snapnow;
    
    if save_figures_svm
        plugin_save_figure
    end
    

    % PLOT THE THRESHOLDED SVM MAPS
    % --------------------------------------------------------------------    
    
    if dobootstrap_svm
        
        if isempty(cons2boot_svm) || ismember(c,cons2boot_svm)
        
            % get the stats results for this contrast, with weight map

            bs_stats = bootstrap_svm_stats{c};

            % FDR-corrected

            fprintf('\n\n');
            printhdr('FDR-corrected SVM results');
            fprintf('\n\n');

                % montage

                whmontage = 5; % montage to add title to

                fprintf ('\nMONTAGE SVM RESULTS AT FDR q < %1.4f, k = %d, CONTRAST: %s, %s, SCALING: %s\n\n', q_threshold_svm, k_threshold_svm, analysisname, mask_string, scaling_string);

                t = bs_stats.weight_obj;
                t = threshold(t, q_threshold_svm, 'fdr', 'k', k_threshold_svm); 
                r = region(t,'noverbose');

                o2 = montage(r, 'colormap', 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
                o2 = title_montage(o2, whmontage, [analysisname ' FDR ' num2str(q_threshold_svm) ' ' mask_string ' ' scaling_string]);

                figtitle = sprintf('%s_%s_%1.4f_FDR_montage_%s_%s', analysisname, results_suffix, q_threshold_svm, mask_string, scaling_string);
                set(gcf, 'Tag', figtitle, 'WindowState','maximized');
                drawnow, snapnow;
                    if save_figures_svm
                        plugin_save_figure;
                    end

                clear o2, clear figtitle

                % table and montage of regioncenters

                fprintf ('\nTABLE SVM RESULTS AT FDR q < %1.4f, k = %d, CONTRAST: %s, %s, SCALING: %s\n\n', q_threshold_svm, k_threshold_svm, analysisname, mask_string, scaling_string);

                r(cat(1, r.numVox) < k_threshold_svm) = [];
                [rpos, rneg] = table(r);       % add labels
                r = [rpos rneg];               % re-concatenate labeled regions

                if ~isempty(r)

                    fprintf ('\nMONTAGE REGIONCENTERS SVM RESULTS AT FDR q < %1.4f, k = %d, CONTRAST: %s, %s, SCALING: %s\n\n', q_threshold_svm, k_threshold_svm, analysisname, mask_string, scaling_string);

                    o3 = montage(r, 'colormap', 'regioncenters', 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});

                    % Activate, name, and save figure
                    figtitle = sprintf('%s_%s_%1.4f_FDR_regions_%s_%s', analysisname, results_suffix, q_threshold_svm, mask_string, scaling_string);
                    set(gcf, 'Tag', figtitle, 'WindowState','maximized');
                    drawnow, snapnow;
                        if save_figures_svm
                            plugin_save_figure;
                        end

                    clear o3, clear figtitle, clear t, clear r

                end % conditional montage plot if there are regions to show

            % uncorrected

            fprintf('\n\n');
            printhdr('uncorrected SVM results');
            fprintf('\n\n');

                % montage

                fprintf ('\nMONTAGE SVM RESULTS AT UNCORRECTED p < %1.4f, k = %d, CONTRAST: %s, %s, SCALING: %s\n\n', p_threshold_svm, k_threshold_svm, analysisname, mask_string, scaling_string);

                t = bs_stats.weight_obj;
                t = threshold(t, p_threshold_svm, 'unc', 'k', k_threshold_svm); 
                r = region(t,'noverbose');

                o2 = montage(r, 'colormap', 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
                o2 = title_montage(o2, whmontage, [analysisname ' unc ' num2str(p_threshold_svm) ' ' mask_string ' ' scaling_string]);

                figtitle = sprintf('%s_%s_%1.4f_unc_montage_%s_%s', analysisname, results_suffix, p_threshold_svm, mask_string, scaling_string);
                set(gcf, 'Tag', figtitle, 'WindowState','maximized');
                drawnow, snapnow;
                    if save_figures_svm
                        plugin_save_figure;
                    end

                clear o2, clear figtitle

                % table and montage of regioncenters

                fprintf ('\nTABLE SVM RESULTS AT UNCORRECTED < %1.4f, k = %d, CONTRAST: %s, %s, SCALING: %s\n\n', p_threshold_svm, k_threshold_svm, analysisname, mask_string, scaling_string);

                r(cat(1, r.numVox) < k_threshold_svm) = [];
                [rpos, rneg] = table(r);       % add labels
                r = [rpos rneg];               % re-concatenate labeled regions

                if ~isempty(r)

                    fprintf ('\nMONTAGE REGIONCENTERS SVM RESULTS AT UNCORRECTED p < %1.4f, k = %d, CONTRAST: %s, %s, SCALING: %s\n\n', q_threshold_svm, k_threshold_svm, analysisname, mask_string, scaling_string);

                    o3 = montage(r, 'colormap', 'regioncenters', 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});

                    % Activate, name, and save figure
                    figtitle = sprintf('%s_%s_%1.4f_unc_regions_%s_%s', analysisname, results_suffix, p_threshold_svm, mask_string, scaling_string);
                    set(gcf, 'Tag', figtitle, 'WindowState','maximized');
                    drawnow, snapnow;
                        if save_figures_svm
                            plugin_save_figure;
                        end

                    clear o3, clear figtitle, clear t, clear r

                end % conditional montage plot if there are regions to show
            
        end % if loop cons2boot
        
    end % if loop bootstrap
    
    
    % PLOT THE THRESHOLDED SEARCHLIGHT MAPS
    % --------------------------------------------------------------------    
    
    if dosearchlight_svm
        
        if isempty(cons2searchlight_svm) || ismember(c,cons2searchlight_svm)
        
            % get the stats results for this contrast, with accuracy map

            sl_stats = searchlight_svm_stats{c};

            % FDR-corrected
            
            fprintf('\n\n');
            printhdr('FDR-corrected searchlight results');
            fprintf('\n\n');

                % montage

                whmontage = 5; % montage to add title to

                fprintf ('\nMONTAGE SVM SEARCHLIGHT ACCURACY RESULTS AT FDR q < %1.4f, k = %d, CONTRAST: %s, %s, SCALING: %s\n\n', q_threshold_svm, k_threshold_svm, analysisname, mask_string, scaling_string);

                p = sl_stats.stat_img_obj;
                p = threshold(t, q_threshold_svm, 'fdr', 'k', k_threshold_svm); 
                r = region(p,'noverbose');

                o2 = montage(p, 'maxcolor', [1 1 0], 'mincolor', [1 0 0], 'cmaprange', [min(sl_stats.stat_img_obj.dat) max(sl_stats.stat_img_obj.dat)]);
                o2 = title_montage(o2, whmontage, [analysisname ' searchlight accuracy FDR ' num2str(q_threshold_svm) ' ' mask_string ' ' scaling_string]);

                figtitle = sprintf('%s_%s_%1.4f_FDR_searchlight_montage_%s_%s', analysisname, results_suffix, q_threshold_svm, mask_string, scaling_string);
                set(gcf, 'Tag', figtitle, 'WindowState','maximized');
                drawnow, snapnow;
                    if save_figures_svm
                        plugin_save_figure;
                    end

                clear o2, clear figtitle

                % table and montage of regioncenters

                fprintf ('\nTABLE SVM SEARCHLIGHT ACCURACY RESULTS AT FDR q < %1.4f, k = %d, CONTRAST: %s, %s, SCALING: %s\n\n', q_threshold_svm, k_threshold_svm, analysisname, mask_string, scaling_string);

                r(cat(1, r.numVox) < k_threshold_svm) = [];
                [rpos, rneg] = table(r);       % add labels
                r = [rpos rneg];               % re-concatenate labeled regions

                if ~isempty(r)

                    fprintf ('\nMONTAGE REGIONCENTERS SVM SEARCHLIGHT ACCURACY RESULTS AT FDR q < %1.4f, k = %d, CONTRAST: %s, %s, SCALING: %s\n\n', q_threshold_svm, k_threshold_svm, analysisname, mask_string, scaling_string);

                    o3 = montage(r, 'maxcolor', [1 1 0], 'mincolor', [1 0 0], 'regioncenters', 'cmaprange', [min(sl_stats.stat_img_obj.dat) max(sl_stats.stat_img_obj.dat)]);

                    % Activate, name, and save figure
                    figtitle = sprintf('%s_%s_%1.4f_FDR_searchlight_regions_%s_%s', analysisname, results_suffix, q_threshold_svm, mask_string, scaling_string);
                    set(gcf, 'Tag', figtitle, 'WindowState','maximized');
                    drawnow, snapnow;
                        if save_figures_svm
                            plugin_save_figure;
                        end

                    clear o3, clear figtitle, clear t, clear r

                end % conditional montage plot if there are regions to show

            % uncorrected

            fprintf('\n\n');
            printhdr('FDR-corrected searchlight results');
            fprintf('\n\n');

                % montage

                whmontage = 5; % montage to add title to

                fprintf ('\nMONTAGE SVM SEARCHLIGHT ACCURACY RESULTS AT UNCORRECTED p < %1.4f, k = %d, CONTRAST: %s, %s, SCALING: %s\n\n', p_threshold_svm, k_threshold_svm, analysisname, mask_string, scaling_string);

                p = sl_stats.stat_img_obj;
                p = threshold(t, p_threshold_svm, 'unc', 'k', k_threshold_svm); 
                r = region(p,'noverbose');

                o2 = montage(p, 'maxcolor', [1 1 0], 'mincolor', [1 0 0], 'cmaprange', [min(sl_stats.stat_img_obj.dat) max(sl_stats.stat_img_obj.dat)]);
                o2 = title_montage(o2, whmontage, [analysisname ' searchlight accuracy unc ' num2str(p_threshold_svm) ' ' mask_string ' ' scaling_string]);

                figtitle = sprintf('%s_%s_%1.4f_unc_searchlight_montage_%s_%s', analysisname, results_suffix, p_threshold_svm, mask_string, scaling_string);
                set(gcf, 'Tag', figtitle, 'WindowState','maximized');
                drawnow, snapnow;
                    if save_figures_svm
                        plugin_save_figure;
                    end

                clear o2, clear figtitle

                % table and montage of regioncenters

                fprintf ('\nTABLE SVM SEARCHLIGHT ACCURACY RESULTS AT UNCORRECTED p < %1.4f, k = %d, CONTRAST: %s, %s, SCALING: %s\n\n', p_threshold_svm, k_threshold_svm, analysisname, mask_string, scaling_string);

                r(cat(1, r.numVox) < k_threshold_svm) = [];
                [rpos, rneg] = table(r);       % add labels
                r = [rpos rneg];               % re-concatenate labeled regions

                if ~isempty(r)

                    fprintf ('\nMONTAGE REGIONCENTERS SVM SEARCHLIGHT ACCURACY RESULTS AT UNCORRECTED p < %1.4f, k = %d, CONTRAST: %s, %s, SCALING: %s\n\n', p_threshold_svm, k_threshold_svm, analysisname, mask_string, scaling_string);

                    o3 = montage(r, 'maxcolor', [1 1 0], 'mincolor', [1 0 0], 'regioncenters', 'cmaprange', [min(sl_stats.stat_img_obj.dat) max(sl_stats.stat_img_obj.dat)]);

                    % Activate, name, and save figure
                    figtitle = sprintf('%s_%s_%1.4f_unc_searchlight_regions_%s_%s', analysisname, results_suffix, p_threshold_svm, mask_string, scaling_string);
                    set(gcf, 'Tag', figtitle, 'WindowState','maximized');
                    drawnow, snapnow;
                        if save_figures_svm
                            plugin_save_figure;
                        end

                    clear o3, clear figtitle, clear t, clear r

                end % conditional montage plot if there are regions to show
            
        end % if loop cons2searchlight
        
    end % if loop searchlight
    
end  % loop over contrasts


%% CROSS-CLASSIFICATION MATRIX
% -------------------------------------------------------------------------

% uses svm_dist_pos_neg_matrix from plugin

fprintf('\n\n');
printhdr('SVM cross-classification matrix');
fprintf('\n\n');

% Get rid of empties (invalid contrasts)
isemptymtx = cellfun(@isempty, svm_dist_pos_neg_matrix);
wh_ok = ~all(isemptymtx);

%svm_dist_pos_neg_matrix(wh_ok, wh_ok)

diff_function = @(x) x(:, 1) - x(:, 2);         % should be positive for correct classification

iscorrect = @(x) sign(diff_function(x)) > 0;

cohens_d_function = @(x) mean(x) ./ std(x);

acc_function = @(corr_idx) 100 * sum(corr_idx) ./ length(corr_idx);

svm_dist_per_subject_and_condition = cellfun(diff_function, svm_dist_pos_neg_matrix(wh_ok, wh_ok), 'UniformOutput', false);

svm_cohens_d_train_transfer = cell2mat(cellfun(cohens_d_function, svm_dist_per_subject_and_condition, 'UniformOutput', false));

accuracy_by_subject_and_condition = cellfun(iscorrect, svm_dist_pos_neg_matrix(wh_ok, wh_ok), 'UniformOutput', false);

accuracy = cellfun(acc_function, accuracy_by_subject_and_condition, 'UniformOutput', false);
accuracy = cell2mat(accuracy);


%% FIGURE CROSS-VALIDATED DISTANCE FROM HYPERPLANE
% -------------------------------------------------------------------------

figtitle = sprintf('SVM cross-classification');
create_figure(figtitle);
set(gcf,'WindowState','maximized');

% pos = get(gcf, 'Position');
% pos(3) = pos(3) * 1.7;
% set(gcf, 'Position', pos)

fprintf('\n\n');
printhdr('SVM cross-validated distance from hyperplane');
fprintf('\n> 0 is correct classification\n');
fprintf('\n\n');

[ntrain, ntransfer] = size(svm_dist_per_subject_and_condition);
text_xval = [];
han = {};

trainnames = DAT.contrastnames(wh_ok);

for c = 1:ntrain  % for non-empty contrasts only
   
    dat = svm_dist_per_subject_and_condition(c, :);
    
    xvals = 1 + ntransfer * (c-1) : c * ntransfer;
    
    xvals = xvals + c - 1; % skip a space

    text_xval(c) = mean(xvals);
    mycolors = DAT.contrastcolors(wh_ok);
    
    trainname = trainnames{c};
    xtick_text{c} = sprintf('Train %s', trainname);
    
    sprintf('\nTrain on %s\n', trainname);
    
    han{c} = barplot_columns(dat, 'nofig', 'noviolin', 'colors', mycolors, 'x', xvals, 'names', trainnames);
    set(gca, 'XLim', [.5 xvals(end) + .5]);
    
end

xlabel(' ');
ylabel('Distance from hyperplane');

barhandles = cat(2, han{1}.bar_han{:});
legend(barhandles, trainnames)

set(gca, 'XTick', text_xval, 'XTickLabel', xtick_text, 'XTickLabelRotation', 45);

fprintf('\nAccuracy matrix - training (rows) by test contrasts (columns)\n');
print_matrix(accuracy, trainnames, trainnames);

if save_figures_svm
    plugin_save_figure;
end


%% FIGURE AND STATS ON EFFECT SIZES
% -------------------------------------------------------------------------

figtitle = sprintf('SVM cross-classification effect sizes');
create_figure(figtitle);
set(gcf,'WindowState','maximized');

% pos = get(gcf, 'Position');
% pos(3) = pos(3) * 1.7;
% set(gcf, 'Position', pos)

fprintf('\n\n');
printhdr(figtitle);
fprintf('\n\n');

imagesc(svm_cohens_d_train_transfer, [0 max(svm_cohens_d_train_transfer(:))])
set(gca, 'YDir', 'reverse', 'YTick', 1:ntrain,  'YTickLabel', xtick_text(1:ntrain), 'XTick', 1:ntrain, 'XTickLabel', trainnames, 'XTickLabelRotation', 45);
title(figtitle);
xlabel('Test condition');
ylabel('Training condition');
colorbar
cm = colormap_tor([1 1 1], [1 0 0]);
colormap(cm)

if save_figures_svm
    plugin_save_figure;
end

disp('Cohen''s d for training and transfer');
print_matrix(svm_cohens_d_train_transfer, trainnames, xtick_text(1:ntrain));

