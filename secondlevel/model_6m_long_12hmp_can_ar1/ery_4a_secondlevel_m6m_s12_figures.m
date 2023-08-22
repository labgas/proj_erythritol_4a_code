%% DEFINE PATHS AND LOAD VARIABLES
% -------------------------------------------------------------------------

% SET PATHS AND DEFINE DIRS

addpath(genpath('C:\Users\u0139539\OneDrive - KU Leuven\fmri course\proj_erythritol_4a_code')); % add proj-emosymp Github repo to your path
addpath(genpath('C:\Users\u0139539\OneDrive - KU Leuven\fmri course\RainCloudPlots')); % add RainCloudPlots Github repo to your path
addpath(genpath('C:\Users\u0139539\OneDrive - KU Leuven\fmri course\Robust_Statistical_Toolbox')); % add Robust Statistical Toolbox Github repo to your path
%addpath(genpath('C:\Users\lukas\Documents\MATLAB\cbrewer')); % add colorbrewer to your path for more color options

ery_4a_secondlevel_m6m_s0_a_set_up_paths_always_run_first

figspubdir = fullfile(resultsdir,'figures_publication');
if ~isfolder(figspubdir)
    mkdir(figspubdir);
end


% LOAD DATA AND DEFINE VARIABLES

% Behavioral

%load(fullfile(resultsdir,'image_names_and_setup.mat'));

% conditions
%idx = ~isnan(DAT.BEHAVIOR.behavioral_data_table.NA_neg);
% behdat = DAT.BEHAVIOR.behavioral_data_table(idx,:);
% behdat = sortrows(behdat,'patient','ascend');
% behdat.patient(behdat.patient == -1) = 2;
% group = [behdat.patient; behdat.patient; behdat.patient];
% condition = [ones(height(behdat),1); 2.*ones(height(behdat),1); 3.*ones(height(behdat),1)];
% NA = [behdat.NA_neg; behdat.NA_neu; behdat.NA_pos];
% symptoms = [behdat.symptoms_neg; behdat.symptoms_neu; behdat.symptoms_pos];
% outcomes = {NA,symptoms};
% outcome_names = {'negative affect','physical symptoms'};
% 
% for o = 1:size(outcomes,2)
%     D{o} = [outcomes{o},condition,group];
%     for i = 1:3
%         for j = 1:2
%         data{i,j} = D{o}(D{o}(:, 2) == i & D{o}(:, 3) ==j);
%         end
%     end
%     
%     data_all{o} = data;
% end
% 
% % contrasts
% for k = 1:max(unique(behdat.patient))
%     
%     NA_contrast{k} = [behdat.NA_neg_neu(behdat.patient==k),behdat.NA_neg_pos(behdat.patient==k)];
%     symptoms_contrast{k} = [behdat.symptoms_neg_neu(behdat.patient==k),behdat.symptoms_pos_neu(behdat.patient==k)];
%     
% end
% outcomes_contrast = {NA_contrast,symptoms_contrast};
% contrast_names = {'negative versus neutral','negative versus positive'};
% 
% % Signature responses
% 
% sigdat = DAT.SIG_contrasts.raw.dotproduct;
% behdat_full = DAT.BEHAVIOR.behavioral_data_table;
% behdat_full.patient(behdat_full.patient == -1) = 2;
% 
% for n = 1:max(unique(behdat_full.patient))
%     
%     NPS{n} = [sigdat.NPS.Negative_v_Neutral(behdat_full.patient==n),sigdat.NPS.Negative_v_Positive(behdat_full.patient==n)];
%     NPSpos{n} = [sigdat.NPSpos.Negative_v_Neutral(behdat_full.patient==n),sigdat.NPSpos.Negative_v_Positive(behdat_full.patient==n)];
%     NPSneg{n} = [sigdat.NPSneg.Negative_v_Neutral(behdat_full.patient==n),sigdat.NPSneg.Negative_v_Positive(behdat_full.patient==n)];
%     PINES{n} = [sigdat.PINES.Negative_v_Neutral(behdat_full.patient==n),sigdat.PINES.Negative_v_Positive(behdat_full.patient==n)];
%     SIIPS{n} = [sigdat.SIIPS.Negative_v_Neutral(behdat_full.patient==n),sigdat.SIIPS.Negative_v_Positive(behdat_full.patient==n)];
%     
% end
% 
% signatures = {NPS, NPSpos, NPSneg, PINES, SIIPS};
% signature_names = {'NPS', 'NPS positive', 'NPS negative', 'PINES', 'SIIPS'};
% contrast_names = {'negative versus neutral','negative versus positive'};
% 
% npsposregions = DAT.NPSsubregions.npspos_by_region_contrasts;
% npsposregions_names = DAT.NPSsubregions.posnames;
% npsposregions_names_full = {'vermis','R insula','R V1','R thalamus','L insula','R dorsal posterior insula','R S2 operculum','anterior midcingulate'};
% 
% for t = 1:max(unique(behdat_full.patient))
%     
%     for u = 1:size(npsposregions{1},2)
%        regions{u,t} = [npsposregions{1}(behdat_full.patient==t,u),npsposregions{2}(behdat_full.patient==t,u)];
%     end
%     
% end

% Robust parcelwise analysis

load(fullfile(resultsdir,'parcelwise_stats_and_maps_contrasts_no_scaling_cov_rating.mat'));
load(fullfile(resultsdir,'data_objects.mat'));
t_obj_suc_wat_delta_group = get_wh_image(parcelwise_stats_results{1,4}.t_obj,1); % loads statistic_image object for group effect on fourth contrast (sucrose-water) first variable (delta)
t_obj_suc_wat_intercept_group = get_wh_image(parcelwise_stats_results{1,4}.t_obj,2); % sucrose-water intercept
t_obj_ery_wat_delta_group = get_wh_image(parcelwise_stats_results{1,5}.t_obj,1); % ery-water delta
t_obj_ery_wat_intercept_group = get_wh_image(parcelwise_stats_results{1,5}.t_obj,2); %ery-water interecpt
t_obj_sl_wat_delta_group = get_wh_image(parcelwise_stats_results{1,6}.t_obj,1); % sucralose-water delta
t_obj_sl_wat_intercept_group = get_wh_image(parcelwise_stats_results{1,6}.t_obj,2); %sucralose-water intercept


% PDM analysis

% pdm_dir = fullfile(resultsdir,'mediation_analysis','brain_pdm');
% pdm_negneudir = fullfile(pdm_dir,'Negative_v_Neutral','somatic_symptoms');
% pdm_regiondir = fullfile(pdm_negneudir, 'region_objects_and_tables');
% 
% load(fullfile(pdm_regiondir,'region_objects_and_tables.mat'));
% load(fullfile(pdm_negneudir,'PDM_source_recon_Negative_v_Neutral_somatic_symptoms.mat'));

% Signature patterns for plotting

% gray_matter_mask = which('gray_matter_mask_sparse.img');
% nps_pat = apply_mask(fmri_data(which('weights_NSF_grouppred_cvpcr_FDR05_smoothed_fwhm05.img.gz')),gray_matter_mask);
% pines_pat = apply_mask(fmri_data(which('Rating_LASSO_PCR_boot5000_fdr05_k10_2.nii.gz')),gray_matter_mask);
% siips_pat = apply_mask(fmri_data(which('nonnoc_v11_4_subcluster_maps_fdr05_pattern_wttest.nii.gz')),gray_matter_mask);
% 
% pat_objs = {nps_pat,pines_pat,siips_pat};
% pat_obj_names = {'NPS','PINES','SIIPS'};


% DEFINE COLORS
% to be used in raincloudplots using the great cbrewer
% function - cbrewer() for overview of color schemes, and help cbrewer for
% info on function
try
    % get nice colours from colorbrewer
    % (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)
    [cb] = cbrewer('qual', 'Set1', 12, 'pchip');
catch
    % if you don't have colorbrewer, accept these far more boring colours
    cb = [0.5 0.8 0.9; 1 1 0.7; 0.7 0.8 0.9; 0.8 0.5 0.4; 0.5 0.7 0.8; 1 0.8 0.5; 0.7 1 0.4; 1 0.7 1; 0.6 0.6 0.6; 0.7 0.5 0.7; 0.8 0.9 0.8; 1 1 0.4];
end

cl(1, :) = cb(3, :);
cl(2, :) = cb(1, :);

fig_position = [200 200 600 400]; % coordinates for figures


%% BEHAVIORAL DATA
%--------------------------------------------------------------------------

% CONDITIONS - RM RAINCLOUD PLOT

% individual figures

% for o = 1:size(outcomes,2)
%     f  = figure('Position', fig_position,'WindowState','maximized');
%     h   = rm_raincloud(data_all{o}, cl, 0, 'rash');
%     ax{o} = gca;
%     ax{o}.FontSize = 14;
%     ax{o}.FontName = 'Cambria';
%     ax{o}.XAxis.LineWidth = 1;
%     ax{o}.YAxis.LineWidth = 1;
%     xlabel({strcat(outcome_names{o},' rating'),''},'FontSize',24,'FontWeight','bold');
%     ylabel({'','affective valence'},'FontSize',24,'FontWeight','bold');
%     yticklabels({'\fontsize{20} \bf positive','\fontsize{20} \bf neutral','\fontsize{20} \bf negative'});
%     legend([h.l(1,1) h.l(1,2)],{'FSS patients','healthy controls'},'Location','best','FontSize',24,'FontWeight','bold','Box','off');
%     for i = 1:3
%         for j = 1:2
%             h.s{i, j}.SizeData = 150;
%         end
%     end
%     
%     % save
%     print(f,fullfile(figspubdir,strcat('behav_',outcome_names{o},'.png')),'-dpng','-r600');
%     
%     f_all{o}= f;
%     h_all{o} = h;
%     
%     clear f h;
% end
% 
% % integrated two-panel figure for publication
% 
% fa  = figure('WindowState','maximized');
% subplot(1,2,1)
% ha   = rm_raincloud(data_all{1}, cl, 0, 'rash');
% ax1 = gca;
% ax1.FontSize = 12;
% ax1.FontName = 'Cambria';
% ax1.XAxis.LineWidth = 1;
% ax1.YAxis.LineWidth = 1;
% ax1.Position = [0.07 0.1100 0.42 0.8150];
% % title(['Figure 1' newline 'Repeated measures raincloud plot']);
% xlabel({strcat(outcome_names{1},' rating'),''},'FontSize',20,'FontWeight','bold');
% ylabel({'','affective valence'},'FontSize',20,'FontWeight','bold');
% yticklabels({'\fontsize{16} \bf positive','\fontsize{16} \bf neutral','\fontsize{16} \bf negative'});
% % legend([h.l(1,1) h.l(1,2)],{'FSS patients','healthy controls'},'Location','best','FontSize',20,'FontWeight','bold');
%     for i = 1:3
%         for j = 1:2
%             ha.s{i, j}.SizeData = 100;
%         end
%     end
% subplot(1,2,2)
% hb   = rm_raincloud(data_all{2}, cl, 0, 'rash');
% ax2 = gca;
% ax2.FontSize = 12;
% ax2.FontName = 'Cambria';
% ax2.XAxis.LineWidth = 1;
% ax2.YAxis.LineWidth = 1;
% ax2.Position = [0.5703 0.1100 0.42 0.8150];
% % title(['Figure 1' newline 'Repeated measures raincloud plot']);
% xlabel({strcat(outcome_names{2},' rating'),''},'FontSize',20,'FontWeight','bold');
% ylabel({'','affective valence'},'FontSize',20,'FontWeight','bold');
% yticklabels({'\fontsize{16} \bf positive','\fontsize{16} \bf neutral','\fontsize{16} \bf negative'});
% legend([hb.l(1,1) hb.l(1,2)],{'FSS patients','healthy controls'},'Location','best','FontSize',20,'FontWeight','bold','Box','off');
%     for i = 1:3
%         for j = 1:2
%             hb.s{i, j}.SizeData = 100;
%         end
%     end
%     
% % fa.Position = [-150 50 2250 971]; % can be used to change position of
% % entire figure
%     
% print(fa,fullfile(figspubdir,strcat('behav_',outcome_names{1},'_',outcome_names{2},'.png')),'-dpng','-r600');
% 
% 
% % CONTRASTS - CLASSIC RAINCLOUD PLOT
% 
% for o = 1:size(outcomes,2)
%     
%     for m = 1:size(NA_contrast,2)
%         
%         f2 = figure('Position', fig_position,'WindowState','maximized');
%         h1 = raincloud_plot(outcomes_contrast{o}{1}(:,m), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
%              'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
%              'box_col_match', 1, 'line_width', 3);
%         h2 = raincloud_plot(outcomes_contrast{o}{2}(:,m), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
%              'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 1, 'line_width', 3);
%         h1{1}.EdgeColor = 'none';
%         h2{1}.EdgeColor = 'none';
%         h1{2}.SizeData = 50;
%         h2{2}.SizeData = 50;
%         ax3{o,m} = gca;
%         ax3{o,m}.FontSize = 14;
%         ax3{o,m}.FontName = 'Cambria';
%         ax3{o,m}.YAxisLocation = 'origin';
%         ax3{o,m}.YTick = [];
%         ax3{o,m}.LineWidth = 0.25;
%         xlabel({'',[outcome_names{o},' rating']},'FontSize',24,'FontWeight','bold');
%         legend([h1{1} h2{1}], {'FSS patients', 'healthy controls'},'Location','best','FontSize',24,'FontWeight','bold','Box','off');
%         title(contrast_names{m},'FontSize',28,'FontWeight','bold');
%         ylim([(h2{3}.Position(2)+0.10.*h2{3}.Position(2)) (max([h1{1}.YData h2{1}.YData])+0.05.*max([h1{1}.YData h2{1}.YData]))]);
%         box off
%     
%     
%         % save
%         print(f2,fullfile(figspubdir,strcat('behav_contrasts_',outcome_names{o},'_',contrast_names{m},'.png')),'-dpng','-r600');
% 
%         f2_all{o,m}= f2;
%         h1_all{o,m} = h1;
%         h2_all{o,m} = h2;
% 
%         clear f2 h1 h2;
%         
%     end
%     
% end
% 
% 



%% ROBFIT PARCELWISE BRAIN FIGURES
%--------------------------------------------------------------------------


% NOTE: looping over elements of robfit_parcel_stats_results does not
% work for some weird reason, including using get_wh_image to get the first
% row of the dat

% EXAMPLE CODE FOR DIFFERENT CANLAB PLOTTING FUNCTIONS


% CANLAB_RESULTS_FMRIDISPLAY AND MONTAGE

% see also help for other fmridisplay methods, particularly addblobs and
% montage 

% suc-wat intercept (only negative values)
obj_suc_wat_interept_robust_parcelwise = canlab_results_fmridisplay(t_obj_suc_wat_intercept_group,'outline','linewidth',0.5,'montagetype','full hcp','splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
f7 = gcf;
f7.WindowState = 'maximized';
print(f7,fullfile(figspubdir,strcat('robfit_parcelwise_','sucrose-water_intercept','_montage.png')),'-dpng','-r600');
% 
    % RENDER_ON_SURFACE 
    figure;
    surface_handle_1 = addbrain('coronal_slabs_5'); % help addbrain for all options
    render_on_surface(t_obj_suc_wat_intercept_group, surface_handle_1,'colormap', 'winter');
    f8 = gcf;
    f8.WindowState = 'maximized';
    print(f8,fullfile(figspubdir,strcat('robfit_parcelwise_','sucrose-water_intercept','_coronal_slabs.png')),'-dpng','-r600');

%     surface_handle_2 = addbrain('hires'); % help addbrain for all options
%     render_on_surface(t_obj_suc_wat_intercept_group, surface_handle_2,'colormap','winter');


% ery-water intercept (only positive values)
obj_ery_wat_intercept_robust_parcelwise = canlab_results_fmridisplay(t_obj_ery_wat_intercept_group,'outline','linewidth',0.5,'montagetype','full hcp','splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
f7 = gcf;
f7.WindowState = 'maximized';
print(f7,fullfile(figspubdir,strcat('robfit_parcelwise_','erythritol-water_intercept','_montage.png')),'-dpng','-r600');

    % RENDER_ON_SURFACE 
    figure;
    surface_handle_1 = addbrain('coronal_slabs_5'); % help addbrain for all options
    render_on_surface(t_obj_ery_wat_intercept_group, surface_handle_1,'colormap', 'hot');
    f8 = gcf;
    f8.WindowState = 'maximized';
    print(f8,fullfile(figspubdir,strcat('robfit_parcelwise_','erythritol-water_intercept','_coronal_slabs.png')),'-dpng','-r600');

%     surface_handle_2 = addbrain('hires'); % help addbrain for all options
%     render_on_surface(t_obj_ery_wat_intercept_group, surface_handle_2,'full hcp','colormap','hot');

% sucralose-water intercept (has both positive and negative values)
obj_sl_wat_intercept_robust_parcelwise = canlab_results_fmridisplay(t_obj_sl_wat_intercept_group,'outline','linewidth',0.5,'montagetype','full hcp','splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
f7 = gcf;
f7.WindowState = 'maximized';
print(f7,fullfile(figspubdir,strcat('robfit_parcelwise_','sucralose-water_intercept','_montage.png')),'-dpng','-r600');

    % RENDER_ON_SURFACE 
    figure;
    surface_handle_1 = addbrain('coronal_slabs_5'); % help addbrain for all options
    render_on_surface(t_obj_sl_wat_intercept_group, surface_handle_1,'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
    f8 = gcf;
    f8.WindowState = 'maximized';
    print(f8,fullfile(figspubdir,strcat('robfit_parcelwise_','sucralose-water_intercept','_coronal_slabs.png')),'-dpng','-r600');
% 
% %     surface_handle_2 = addbrain('hires'); % help addbrain for all options
% %     render_on_surface(t_obj_sl_wat_intercept_group, surface_handle_2,'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
% 
% % suc-wat delta (positive and negative)
obj_suc_wat_delta_robust_parcelwise = canlab_results_fmridisplay(t_obj_suc_wat_delta_group,'outline','linewidth',0.5,'montagetype','full hcp','splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
f7 = gcf;
f7.WindowState = 'maximized';
print(f7,fullfile(figspubdir,strcat('robfit_parcelwise_','sucrose-water_delta_ratings','_montage.png')),'-dpng','-r600');

    % RENDER_ON_SURFACE 
    figure;
    surface_handle_1 = addbrain('coronal_slabs_5'); % help addbrain for all options
    render_on_surface(t_obj_suc_wat_delta_group, surface_handle_1,'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
    f8 = gcf;
    f8.WindowState = 'maximized';
    print(f8,fullfile(figspubdir,strcat('robfit_parcelwise_','sucrose-water_delta_ratings','_coronal_slabs.png')),'-dpng','-r600');

% %     surface_handle_2 = addbrain('hires'); % help addbrain for all options
% %     render_on_surface(t_obj_suc_wat_delta_group, surface_handle_2,'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
% 
% % ery-water delta (positive and negative)
obj_ery_wat_delta_robust_parcelwise = canlab_results_fmridisplay(t_obj_ery_wat_delta_group,'outline','linewidth',0.5,'montagetype','full hcp','splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
f7 = gcf;
f7.WindowState = 'maximized';
print(f7,fullfile(figspubdir,strcat('robfit_parcelwise_','erythritol-water_delta','_montage.png')),'-dpng','-r600');

    % RENDER_ON_SURFACE 
    figure;    
    surface_handle_1 = addbrain('coronal_slabs_5'); % help addbrain for all options
    render_on_surface(t_obj_ery_wat_delta_group, surface_handle_1,'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
    f8 = gcf;
    f8.WindowState = 'maximized';
    print(f8,fullfile(figspubdir,strcat('robfit_parcelwise_','erythritol-water_delta','_coronal_slabs.png')),'-dpng','-r600');

%     surface_handle_2 = addbrain('hires'); % help addbrain for all options
%     render_on_surface(t_obj_ery_wat_delta_group, surface_handle_2,'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
    
% sucralose-water delta (only positive)
obj_sl_wat_delta_robust_parcelwise = canlab_results_fmridisplay(t_obj_sl_wat_delta_group,'outline','linewidth',0.5,'montagetype','full hcp','splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
f7 = gcf;
f7.WindowState = 'maximized';
print(f7,fullfile(figspubdir,strcat('robfit_parcelwise_','sucralose-water_delta','_montage.png')),'-dpng','-r600');

    % RENDER_ON_SURFACE 
    figure;
    surface_handle_1 = addbrain('coronal_slabs_5'); % help addbrain for all options
    render_on_surface(t_obj_sl_wat_delta_group, surface_handle_1,'colormap', 'hot');
    f8 = gcf;
    f8.WindowState = 'maximized';
    print(f8,fullfile(figspubdir,strcat('robfit_parcelwise_','sucralose-water_delta','_coronal_slabs.png')),'-dpng','-r600');

%     surface_handle_2 = addbrain('hires'); % help addbrain for all options
%     render_on_surface(t_obj_sl_wat_delta_group, surface_handle_2,'colormap','hot');
    


% SURFACE

% can also be used for more flexibility in defining surfaces, cutaways, etc, 
% in combination with isosurf
% NOTE: see region.surface and image_vector.isosurface as well as
% canlab_canonical_brain_surface_cutaways for help
% 
% anat = fmri_data(which('keuken_2014_enhanced_for_underlay.img'), 'noverbose');
% p = isosurface(anat, 'thresh', 140, 'nosmooth', 'ylim', [-Inf -30]);
% p2 = isosurface(anat, 'thresh', 140, 'nosmooth', 'xlim', [-Inf 0], 'YLim', [-30 Inf]);
% alpha 1 ; lightRestoreSingle; view(135, 30); colormap gray;
% p3 = addbrain('limbic hires');
% set(p3, 'FaceAlpha', .6, 'FaceColor', [.5 .5 .5]);
% delete(p3(3)); p3(3) = [];
% lightRestoreSingle;
% surface_handles = [p p2 p3];
% 
% [surface_handles_1,pcl1,pcl2] = surface(t_obj_neg_neu_group,surface_handles);


% CLUSTER_SURF 

% does not seem to work for statistic_image objects, but check
% it out for atlas or region objects!
% NOTE: you can input any surface from the Github repo
% CanlabCore\canlab_canonical_brains\Canonical_brains_surfaces

% we use a region object for a single region here as example
%cluster_surf(region_13, 2, 'colors', {[1 1 0]}, 'surf_freesurf_inflated_Left.mat'); view(1,360);


% 
% %% PDM BRAIN FIGURES
% %--------------------------------------------------------------------------
% 
% % FULL MONTAGES
% 
% % single pdms
% for r = 1:size(reg_all_fdr,2)
%     o1 = canlab_results_fmridisplay(reg_all_fdr{1,r},'outline','linewidth',0.5,'montagetype','full hcp','splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img'); % o1 is an fmridisplay object - methods fmridisplay for help
% %     o1 = legend(o1);
% %     delete(o1.activation_maps{1,2}.legendhandle); % get rid of legend for contour activation map in object
%     fig1 = gcf;
%     fig1.WindowState = 'maximized';
%     montage_full_pdm{r} = o1;
%     print(fig1,fullfile(figspubdir,strcat('pdm',num2str(r),'_',contrast_names{1},'_montage_full.png')),'-dpng','-r600');
%     close gcf;
%     clear fig1;
% end % loop pdms
% 
% % source reconstruction maps for single pdms
% for s = 1:size(source_obj_j,2)
%    o1b = canlab_results_fmridisplay(source_obj_j{1,s},'outline','linewidth',0.5,'montagetype','full hcp','splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img'); % o1 is an fmridisplay object - methods fmridisplay for help
% %    o1b = canlab_results_fmridisplay(source_obj_j{1,s}.threshold([0 max(source_obj_j{1,s}.dat)],'raw-between'),'outline','linewidth',0.5,'montagetype','full hcp','splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img'); % o1 is an fmridisplay object - methods fmridisplay for help
% %     o1b = legend(o1b);
% %     delete(o1.activation_maps{1,2}.legendhandle); % get rid of legend for contour activation map in object
%     fig1b = gcf;
%     fig1b.WindowState = 'maximized';
%     montage_source_recon_pdm{s} = o1b;
%     print(fig1b,fullfile(figspubdir,strcat('pdm',num2str(s),'_source_recon_',contrast_names{1},'_montage_full.png')),'-dpng','-r600');
%     close gcf;
%     clear fig1b;
% end % loop pdms
% 
% 
% % BRAINSTEM RENDER FOR SINGLE PDMS
% 
% for r = 1:size(reg_all_fdr,2)
%     [cmap] = cbrewer('qual','Set1',24,'pchip');
%     [surfhan,~,~] = surface(reg_all_fdr{1,r},'brainstem_group');
%     for y = [1:10] % surfhan(1,11) is handle for the entire brainstem
%         surfhan(1,y).FaceColor = cmap(y,:);
%         surfhan(1,y).FaceAlpha = 0.10;
%         surfhan(1,y).FaceLighting = 'flat';
%     end
%     for y = 12:size(surfhan,2) % surfhan(1,11) is handle for the entire brainstem
%         surfhan(1,y).FaceColor = cmap(y,:);
%         surfhan(1,y).FaceAlpha = 0.10;
%         surfhan(1,y).FaceLighting = 'flat';
%     end
%     surfhan(1,11).FaceAlpha = 0.60;
%     surfhan(1,11).BackFaceLighting = 'reverselit';
%     surfhan(1,11).FaceLighting = 'gouraud';
%     surfhan(1,11).AmbientStrength = 0.5;
%     fig2 = gcf;
%     fig2.WindowState = 'maximized';
%     ax2 = fig2.Children(5);
%     ax2.XColor = 'none';
%     ax2.YColor = 'none';
%     ax2.ZColor = 'none';
%     ax2.View = [140 20];
%     brainstem_pdm{r} = fig2;
%     print(fig2,fullfile(figspubdir,strcat('pdm',num2str(r),'_',contrast_names{1},'_brainstem.png')),'-dpng','-r600');
%     close gcf;
%     clear fig2;
%     
% end % loop pdms
% 
% 
% % MULTIROW MONTAGE OF BOTH PDMS
% 
% o3 = canlab_results_fmridisplay([],'outline','linewidth',0.5,'montagetype','multirow',2,'overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
% for m = 1:size(o3.montage,2)
%     for a = 1:size(o3.montage{1,m}.axis_handles,2)
%         if m < (size(o3.montage,2)/2)+1
%             o3.montage{1,m}.axis_handles(a).Position(1,2) = o3.montage{1,m}.axis_handles(a).Position(1,2) - 0.10; % will depend on size of screen and figure, can be fixed more properly with screensize() function
%         else
%             o3.montage{1,m}.axis_handles(a).Position(1,2) = o3.montage{1,m}.axis_handles(a).Position(1,2) - 0.30; 
%         end
%     end
% end
% o3 = addblobs(o3,reg_all_fdr{1,1},'wh_montages',1:2,'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
% annotation('textbox',[.5 .475 .4 .5],'String','PDM #1','FontName','Cambria','FontSize',18,'FontWeight','bold','FitBoxToText','on','EdgeColor','none');
% o3 = addblobs(o3,reg_all_fdr{1,2},'wh_montages',3:4,'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
% annotation('textbox',[.5 0.025 .4 .5],'String','PDM #2','FontName','Cambria','FontSize',18,'FontWeight','bold','FitBoxToText','on','EdgeColor','none');
% o3 = legend(o3);
% delete(o3.activation_maps{1,2}.legendhandle); % get rid of legend for contour activation map in object
% fig3 = gcf;
% fig3.WindowState = 'maximized';
% print(fig3,fullfile(figspubdir,strcat('pdm_',contrast_names{1},'_montage_multirow.png')),'-dpng','-r600');
% 
% 
% % CUSTOM MONTAGES OF BOTH PDMS
% 
% % first example
% 
% o4 = canlab_results_fmridisplay(reg_all_fdr{1,1},'outline','splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
% for m = 1:size(o4.montage,2)
%     for a = 1:size(o4.montage{1,m}.axis_handles,2)
%         o4.montage{1,m}.axis_handles(a).Position(1,2) = o4.montage{1,m}.axis_handles(a).Position(1,2) + 0.20; % will depend on size of screen and figure, can be fixed more properly with screensize() function
%     end
% end
% annotation('textbox',[.5 .475 .4 .5],'String','PDM #1','FontName','Cambria','FontSize',18,'FontWeight','bold','FitBoxToText','on','EdgeColor','none');
% clear m a
% o4 = canlab_results_fmridisplay(reg_all_fdr{1,2},'outline','splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
% title('PDM #2');
% for m = 1:size(o4.montage,2)
%     for a = 1:size(o4.montage{1,m}.axis_handles,2)
%         o4.montage{1,m}.axis_handles(a).Position(1,2) = o4.montage{1,m}.axis_handles(a).Position(1,2) - 0.25; % will depend on size of screen and figure, can be fixed more properly with screensize() function
%     end
% end
% annotation('textbox',[.5 0.025 .4 .5],'String','PDM #2','FontName','Cambria','FontSize',18,'FontWeight','bold','FitBoxToText','on','EdgeColor','none');
% clear m a
% o4 = legend(o4);
% delete(o4.activation_maps{1,2}.legendhandle); % avoid duplicate legend with different limits
% o4.activation_maps{1,1}.legendhandle.Position(1,2) = o4.activation_maps{1,1}.legendhandle.Position(1,2) -0.05;
% brighten(.3);
% fig4 = gcf;
% fig4.WindowState = 'maximized';
% print(fig4,fullfile(figspubdir,strcat('pdm_',contrast_names{1},'_montage_custom.png')),'-dpng','-r600');
% 
% % second example
% 
% [fig5,~] = create_figure(strcat('pdm_',contrast_names{1},'_montage_custom'),2,10);
% fig5.WindowState = 'maximized';
% for c = 1:size(fig5.Children,1)
%     axh5{c} = fig5.Children(c);
% end
% o5 = fmridisplay('overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
% o5 = montage(o5,'coronal','slice_range', [-70 20], 'onerow', 'spacing', 10,'outline','existing_axes',[axh5{1,[12:20,1]}],'brighten',0.3); % weird order of axes in object, but no worries
% o5 = addblobs(o5,reg_all_fdr{1,1},'wh_montages',1,'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
% annotation('textbox',[0.45 .475 .4 .5],'String','PDM #1','FontName','Cambria','FontSize',18,'FontWeight','bold','FitBoxToText','on','EdgeColor','none');
% o5 = montage(o5,'coronal','slice_range', [-70 20], 'onerow', 'spacing', 10,'outline','existing_axes',[axh5{1,[2:11]}],'brighten',0.3); % weird order of axes in object, but no worries
% o5 = addblobs(o5,reg_all_fdr{1,2},'wh_montages',2,'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
% annotation('textbox',[0.45 0.055 .4 .5],'String','PDM #2','FontName','Cambria','FontSize',18,'FontWeight','bold','FitBoxToText','on','EdgeColor','none');
% enlarge_axes(gcf,1.75);
% for m = 1:size(o5.montage,2)
%     for a = 1:size(o5.montage{1,m}.axis_handles,2)
%         o5.montage{1,m}.axis_handles(a).Position(1,1) = o5.montage{1,m}.axis_handles(a).Position(1,1) -0.05;
%         if m < (size(o5.montage,2)/2)+1
%             o5.montage{1,m}.axis_handles(a).Position(1,2) = o5.montage{1,m}.axis_handles(a).Position(1,2) - 0.05; % will depend on size of screen and figure, can be fixed more properly with screensize() function
%         else
%             o5.montage{1,m}.axis_handles(a).Position(1,2) = o5.montage{1,m}.axis_handles(a).Position(1,2) + 0.025; 
%         end
%     end
% end
% o5 = legend(o5);
% delete(o5.activation_maps{1,2}.legendhandle); 
% print(fig5,fullfile(figspubdir,strcat('pdm_',contrast_names{1},'_montage_custom2.png')),'-dpng','-r600');
% 
%