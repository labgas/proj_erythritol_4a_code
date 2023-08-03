%% DEFINE PATHS AND LOAD VARIABLES
% -------------------------------------------------------------------------

% SET PATHS AND DEFINE DIRS

addpath(genpath('C:\Users\lukas\Documents\GitHub\proj-emosymp')); % add proj-emosymp Github repo to your path
addpath(genpath('C:\Users\lukas\Documents\GitHub\RainCloudPlots')); % add RainCloudPlots Github repo to your path
addpath(genpath('C:\Users\lukas\Documents\GitHub\Robust_Statistical_Toolbox')); % add Robust Statistical Toolbox Github repo to your path
addpath(genpath('C:\Users\lukas\Documents\MATLAB\cbrewer')); % add colorbrewer to your path for more color options

a_emosymp_m1_s1_set_up_paths_always_run_first

figspubdir = fullfile(resultsdir,'figures_publication');
if ~isfolder(figspubdir)
    mkdir(figspubdir);
end


% LOAD DATA AND DEFINE VARIABLES

% Behavioral

load(fullfile(resultsdir,'image_names_and_setup.mat'));

% conditions
idx = ~isnan(DAT.BEHAVIOR.behavioral_data_table.NA_neg);
behdat = DAT.BEHAVIOR.behavioral_data_table(idx,:);
behdat = sortrows(behdat,'patient','ascend');
behdat.patient(behdat.patient == -1) = 2;
group = [behdat.patient; behdat.patient; behdat.patient];
condition = [ones(height(behdat),1); 2.*ones(height(behdat),1); 3.*ones(height(behdat),1)];
NA = [behdat.NA_neg; behdat.NA_neu; behdat.NA_pos];
symptoms = [behdat.symptoms_neg; behdat.symptoms_neu; behdat.symptoms_pos];
outcomes = {NA,symptoms};
outcome_names = {'negative affect','physical symptoms'};

for o = 1:size(outcomes,2)
    D{o} = [outcomes{o},condition,group];
    for i = 1:3
        for j = 1:2
        data{i,j} = D{o}(D{o}(:, 2) == i & D{o}(:, 3) ==j);
        end
    end
    
    data_all{o} = data;
end

% contrasts
for k = 1:max(unique(behdat.patient))
    
    NA_contrast{k} = [behdat.NA_neg_neu(behdat.patient==k),behdat.NA_neg_pos(behdat.patient==k)];
    symptoms_contrast{k} = [behdat.symptoms_neg_neu(behdat.patient==k),behdat.symptoms_pos_neu(behdat.patient==k)];
    
end
outcomes_contrast = {NA_contrast,symptoms_contrast};
contrast_names = {'negative versus neutral','negative versus positive'};

% Signature responses

sigdat = DAT.SIG_contrasts.raw.dotproduct;
behdat_full = DAT.BEHAVIOR.behavioral_data_table;
behdat_full.patient(behdat_full.patient == -1) = 2;

for n = 1:max(unique(behdat_full.patient))
    
    NPS{n} = [sigdat.NPS.Negative_v_Neutral(behdat_full.patient==n),sigdat.NPS.Negative_v_Positive(behdat_full.patient==n)];
    NPSpos{n} = [sigdat.NPSpos.Negative_v_Neutral(behdat_full.patient==n),sigdat.NPSpos.Negative_v_Positive(behdat_full.patient==n)];
    NPSneg{n} = [sigdat.NPSneg.Negative_v_Neutral(behdat_full.patient==n),sigdat.NPSneg.Negative_v_Positive(behdat_full.patient==n)];
    PINES{n} = [sigdat.PINES.Negative_v_Neutral(behdat_full.patient==n),sigdat.PINES.Negative_v_Positive(behdat_full.patient==n)];
    SIIPS{n} = [sigdat.SIIPS.Negative_v_Neutral(behdat_full.patient==n),sigdat.SIIPS.Negative_v_Positive(behdat_full.patient==n)];
    
end

signatures = {NPS, NPSpos, NPSneg, PINES, SIIPS};
signature_names = {'NPS', 'NPS positive', 'NPS negative', 'PINES', 'SIIPS'};
contrast_names = {'negative versus neutral','negative versus positive'};

npsposregions = DAT.NPSsubregions.npspos_by_region_contrasts;
npsposregions_names = DAT.NPSsubregions.posnames;
npsposregions_names_full = {'vermis','R insula','R V1','R thalamus','L insula','R dorsal posterior insula','R S2 operculum','anterior midcingulate'};

for t = 1:max(unique(behdat_full.patient))
    
    for u = 1:size(npsposregions{1},2)
       regions{u,t} = [npsposregions{1}(behdat_full.patient==t,u),npsposregions{2}(behdat_full.patient==t,u)];
    end
    
end

% Robust parcelwise analysis

load(fullfile(resultsdir,'robfit_parcel_stats_and_maps_no_scaling.mat'));
load(fullfile(resultsdir,'data_objects.mat'));
t_obj_neg_neu_group = get_wh_image(robfit_parcel_stats_results{1,1}.t_obj,1); % loads statistic_image object for group effect on first contrast
region_13 = robfit_parcel_stats_results{1,1}.region_objects{1,1}(1,13); % loads region object for a single region for same effect on same contrast, to be used with cluster_surf below as an example
region_9 = robfit_parcel_stats_results{1,1}.region_objects{1,1}(1,9); % loads region object for a single region for same effect on same contrast, to be used to plot individual responses in this large somatosensory region
region_4 = robfit_parcel_stats_results{1,1}.region_objects{1,1}(1,4); % loads region object for a single region for same effect on same contrast, to be used to plot individual responses in smaller posterior insula region

% PDM analysis

pdm_dir = fullfile(resultsdir,'mediation_analysis','brain_pdm');
pdm_negneudir = fullfile(pdm_dir,'Negative_v_Neutral','somatic_symptoms');
pdm_regiondir = fullfile(pdm_negneudir, 'region_objects_and_tables');

load(fullfile(pdm_regiondir,'region_objects_and_tables.mat'));
load(fullfile(pdm_negneudir,'PDM_source_recon_Negative_v_Neutral_somatic_symptoms.mat'));

% Signature patterns for plotting

gray_matter_mask = which('gray_matter_mask_sparse.img');
nps_pat = apply_mask(fmri_data(which('weights_NSF_grouppred_cvpcr_FDR05_smoothed_fwhm05.img.gz')),gray_matter_mask);
pines_pat = apply_mask(fmri_data(which('Rating_LASSO_PCR_boot5000_fdr05_k10_2.nii.gz')),gray_matter_mask);
siips_pat = apply_mask(fmri_data(which('nonnoc_v11_4_subcluster_maps_fdr05_pattern_wttest.nii.gz')),gray_matter_mask);

pat_objs = {nps_pat,pines_pat,siips_pat};
pat_obj_names = {'NPS','PINES','SIIPS'};


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

for o = 1:size(outcomes,2)
    f  = figure('Position', fig_position,'WindowState','maximized');
    h   = rm_raincloud(data_all{o}, cl, 0, 'rash');
    ax{o} = gca;
    ax{o}.FontSize = 14;
    ax{o}.FontName = 'Cambria';
    ax{o}.XAxis.LineWidth = 1;
    ax{o}.YAxis.LineWidth = 1;
    xlabel({strcat(outcome_names{o},' rating'),''},'FontSize',24,'FontWeight','bold');
    ylabel({'','affective valence'},'FontSize',24,'FontWeight','bold');
    yticklabels({'\fontsize{20} \bf positive','\fontsize{20} \bf neutral','\fontsize{20} \bf negative'});
    legend([h.l(1,1) h.l(1,2)],{'FSS patients','healthy controls'},'Location','best','FontSize',24,'FontWeight','bold','Box','off');
    for i = 1:3
        for j = 1:2
            h.s{i, j}.SizeData = 150;
        end
    end
    
    % save
    print(f,fullfile(figspubdir,strcat('behav_',outcome_names{o},'.png')),'-dpng','-r600');
    
    f_all{o}= f;
    h_all{o} = h;
    
    clear f h;
end

% integrated two-panel figure for publication

fa  = figure('WindowState','maximized');
subplot(1,2,1)
ha   = rm_raincloud(data_all{1}, cl, 0, 'rash');
ax1 = gca;
ax1.FontSize = 12;
ax1.FontName = 'Cambria';
ax1.XAxis.LineWidth = 1;
ax1.YAxis.LineWidth = 1;
ax1.Position = [0.07 0.1100 0.42 0.8150];
% title(['Figure 1' newline 'Repeated measures raincloud plot']);
xlabel({strcat(outcome_names{1},' rating'),''},'FontSize',20,'FontWeight','bold');
ylabel({'','affective valence'},'FontSize',20,'FontWeight','bold');
yticklabels({'\fontsize{16} \bf positive','\fontsize{16} \bf neutral','\fontsize{16} \bf negative'});
% legend([h.l(1,1) h.l(1,2)],{'FSS patients','healthy controls'},'Location','best','FontSize',20,'FontWeight','bold');
    for i = 1:3
        for j = 1:2
            ha.s{i, j}.SizeData = 100;
        end
    end
subplot(1,2,2)
hb   = rm_raincloud(data_all{2}, cl, 0, 'rash');
ax2 = gca;
ax2.FontSize = 12;
ax2.FontName = 'Cambria';
ax2.XAxis.LineWidth = 1;
ax2.YAxis.LineWidth = 1;
ax2.Position = [0.5703 0.1100 0.42 0.8150];
% title(['Figure 1' newline 'Repeated measures raincloud plot']);
xlabel({strcat(outcome_names{2},' rating'),''},'FontSize',20,'FontWeight','bold');
ylabel({'','affective valence'},'FontSize',20,'FontWeight','bold');
yticklabels({'\fontsize{16} \bf positive','\fontsize{16} \bf neutral','\fontsize{16} \bf negative'});
legend([hb.l(1,1) hb.l(1,2)],{'FSS patients','healthy controls'},'Location','best','FontSize',20,'FontWeight','bold','Box','off');
    for i = 1:3
        for j = 1:2
            hb.s{i, j}.SizeData = 100;
        end
    end
    
% fa.Position = [-150 50 2250 971]; % can be used to change position of
% entire figure
    
print(fa,fullfile(figspubdir,strcat('behav_',outcome_names{1},'_',outcome_names{2},'.png')),'-dpng','-r600');


% CONTRASTS - CLASSIC RAINCLOUD PLOT

for o = 1:size(outcomes,2)
    
    for m = 1:size(NA_contrast,2)
        
        f2 = figure('Position', fig_position,'WindowState','maximized');
        h1 = raincloud_plot(outcomes_contrast{o}{1}(:,m), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
             'box_col_match', 1, 'line_width', 3);
        h2 = raincloud_plot(outcomes_contrast{o}{2}(:,m), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 1, 'line_width', 3);
        h1{1}.EdgeColor = 'none';
        h2{1}.EdgeColor = 'none';
        h1{2}.SizeData = 50;
        h2{2}.SizeData = 50;
        ax3{o,m} = gca;
        ax3{o,m}.FontSize = 14;
        ax3{o,m}.FontName = 'Cambria';
        ax3{o,m}.YAxisLocation = 'origin';
        ax3{o,m}.YTick = [];
        ax3{o,m}.LineWidth = 0.25;
        xlabel({'',[outcome_names{o},' rating']},'FontSize',24,'FontWeight','bold');
        legend([h1{1} h2{1}], {'FSS patients', 'healthy controls'},'Location','best','FontSize',24,'FontWeight','bold','Box','off');
        title(contrast_names{m},'FontSize',28,'FontWeight','bold');
        ylim([(h2{3}.Position(2)+0.10.*h2{3}.Position(2)) (max([h1{1}.YData h2{1}.YData])+0.05.*max([h1{1}.YData h2{1}.YData]))]);
        box off
    
    
        % save
        print(f2,fullfile(figspubdir,strcat('behav_contrasts_',outcome_names{o},'_',contrast_names{m},'.png')),'-dpng','-r600');

        f2_all{o,m}= f2;
        h1_all{o,m} = h1;
        h2_all{o,m} = h2;

        clear f2 h1 h2;
        
    end
    
end


%% SIGNATURE RESPONSES
%--------------------------------------------------------------------------

% INDIVIDUAL FIGURES

for s = 1:size(signatures,2)
    
    for p = 1:size(NPS,2)
        
        f3 = figure('Position', fig_position,'WindowState','maximized');
        h3 = raincloud_plot(signatures{s}{1}(:,p), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
             'box_col_match', 1, 'line_width', 3);
        h4 = raincloud_plot(signatures{s}{2}(:,p), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 1, 'line_width', 3);
        h3{1}.EdgeColor = 'none';
        h4{1}.EdgeColor = 'none';
        h3{2}.SizeData = 50;
        h4{2}.SizeData = 50;
        ax4{s,p} = gca;
        ax4{s,p}.FontSize = 14;
        ax4{s,p}.FontName = 'Cambria';
        ax4{s,p}.YAxisLocation = 'origin';
        ax4{s,p}.YTick = [];
        ax4{s,p}.LineWidth = 0.25;
        xlabel({'',[signature_names{s},' response']},'FontSize',24,'FontWeight','bold');
        legend([h3{1} h4{1}], {'FSS patients', 'healthy controls'},'Location','best','FontSize',24,'FontWeight','bold','Box','off');
        title(contrast_names{p},'FontSize',28,'FontWeight','bold');
        ylim([(h4{3}.Position(2)+0.10.*h4{3}.Position(2)) (max([h3{1}.YData h4{1}.YData])+0.05.*max([h3{1}.YData h4{1}.YData]))]);
        box off
    
    
        % save
        print(f3,fullfile(figspubdir,strcat(signature_names{s},'_',contrast_names{p},'.png')),'-dpng','-r600');

        f3_all{s,p}= f3;
        h3_all{s,p} = h3;
        h4_all{s,p} = h4;

        clear f3 h3 h4;
        
    end
    
end

clear p s;


% INTEGRATED PANEL FIGURES FOR PUBLICATION

% PINES, NPS, and SIIPS

f4  = figure('Position', fig_position,'WindowState','maximized');

for p = 1:size(NPS,2)
    subplot(3,2,p)
        h5 = raincloud_plot(signatures{4}{1}(:,p), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
                     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
                     'box_col_match', 1);
        h6 = raincloud_plot(signatures{4}{2}(:,p), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 1);
        h5{1}.EdgeColor = 'none';
        h6{1}.EdgeColor = 'none';
        h5{2}.SizeData = 10;
        h6{2}.SizeData = 10;
        ax4{p} = gca;
        ax4{p}.FontSize = 12;
        ax4{p}.FontName = 'Cambria';
        ax4{p}.YAxisLocation = 'origin';
        ax4{p}.YTick = [];
        ax4{p}.LineWidth = 0.25;
        xlabel({[signature_names{4},' response']},'FontSize',18,'FontWeight','bold');
        legend([h5{1} h6{1}], {'FSS patients', 'healthy controls'},'Location','best','FontSize',14,'FontWeight','bold','Box','off');
        title({(contrast_names{p}),''},'FontSize',20,'FontWeight','bold');
        ylim([(h6{3}.Position(2)+0.10.*h6{3}.Position(2)) (max([h5{1}.YData h6{1}.YData])+0.05.*max([h5{1}.YData h6{1}.YData]))]);
        box off
end

for q = 3:(size(NPS,2)+2)
    subplot(3,2,q)
        h7 = raincloud_plot(signatures{1}{1}(:,q-2), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
                     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
                     'box_col_match', 1);
        h8 = raincloud_plot(signatures{1}{2}(:,q-2), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 1);
        h7{1}.EdgeColor = 'none';
        h8{1}.EdgeColor = 'none';
        h7{2}.SizeData = 10;
        h8{2}.SizeData = 10;
        ax4{q} = gca;
        ax4{q}.FontSize = 12;
        ax4{q}.FontName = 'Cambria';
        ax4{q}.YAxisLocation = 'origin';
        ax4{q}.YTick = [];
        ax4{q}.LineWidth = 0.25;
        xlabel({[signature_names{1},' response']},'FontSize',18,'FontWeight','bold');
        legend([h7{1} h8{1}], {'FSS patients', 'healthy controls'},'Location','best','FontSize',14,'FontWeight','bold','Box','off');
%         title(contrast_names{q-2},'FontSize',20,'FontWeight','bold');
        ylim([(h8{3}.Position(2)+0.10.*h8{3}.Position(2)) (max([h7{1}.YData h8{1}.YData])+0.05.*max([h7{1}.YData h8{1}.YData]))]);
        box off
end

for r = 5:(size(NPS,2)+4)
    subplot(3,2,r)
        h7a = raincloud_plot(signatures{5}{1}(:,r-4), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
                     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
                     'box_col_match', 1);
        h8a = raincloud_plot(signatures{5}{2}(:,r-4), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 1);
        h7a{1}.EdgeColor = 'none';
        h8a{1}.EdgeColor = 'none';
        h7a{2}.SizeData = 10;
        h8a{2}.SizeData = 10;
        ax4{q} = gca;
        ax4{q}.FontSize = 12;
        ax4{q}.FontName = 'Cambria';
        ax4{q}.YAxisLocation = 'origin';
        ax4{q}.YTick = [];
        ax4{q}.LineWidth = 0.25;
        xlabel({[signature_names{5},' response']},'FontSize',18,'FontWeight','bold');
        legend([h7a{1} h8a{1}], {'FSS patients', 'healthy controls'},'Location','best','FontSize',14,'FontWeight','bold','Box','off');
%         title(contrast_names{q-2},'FontSize',20,'FontWeight','bold');
        ylim([(h8a{3}.Position(2)+0.10.*h8a{3}.Position(2)) (max([h7a{1}.YData h8a{1}.YData])+0.05.*max([h7a{1}.YData h8a{1}.YData]))]);
        box off
end
    
print(f4,fullfile(figspubdir,strcat(signature_names{4},'_',signature_names{1},'_',signature_names{5},'.png')),'-dpng','-r600');

clear p q r;

% NPS positive and negative

f5  = figure('Position', fig_position,'WindowState','maximized');

for p = 1:size(NPS,2)
    subplot(2,2,p)
        h9 = raincloud_plot(signatures{2}{1}(:,p), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
                     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
                     'box_col_match', 1);
        h10 = raincloud_plot(signatures{2}{2}(:,p), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 1);
        h9{1}.EdgeColor = 'none';
        h10{1}.EdgeColor = 'none';
        h9{2}.SizeData = 50;
        h10{2}.SizeData = 50;
        ax4{p} = gca;
        ax4{p}.FontSize = 12;
        ax4{p}.FontName = 'Cambria';
        ax4{p}.YAxisLocation = 'origin';
        ax4{p}.YTick = [];
        ax4{p}.LineWidth = 0.25;
        xlabel({[signature_names{2},' response']},'FontSize',18,'FontWeight','bold');
        legend([h9{1} h10{1}], {'FSS patients', 'healthy controls'},'Location','best','FontSize',14,'FontWeight','bold','Box','off');
        title({(contrast_names{p}),''},'FontSize',20,'FontWeight','bold');
        ylim([(h10{3}.Position(2)+0.10.*h10{3}.Position(2)) (max([h9{1}.YData h10{1}.YData])+0.05.*max([h9{1}.YData h10{1}.YData]))]);
        box off
end

for q = 3:(size(NPS,2)+2)
    subplot(2,2,q)
        h11 = raincloud_plot(signatures{3}{1}(:,q-2), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
                     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
                     'box_col_match', 1);
        h12 = raincloud_plot(signatures{3}{2}(:,q-2), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 1);
        h11{1}.EdgeColor = 'none';
        h12{1}.EdgeColor = 'none';
        h11{2}.SizeData = 50;
        h12{2}.SizeData = 50;
        ax4{q} = gca;
        ax4{q}.FontSize = 12;
        ax4{q}.FontName = 'Cambria';
        ax4{q}.YAxisLocation = 'origin';
        ax4{q}.YTick = [];
        ax4{q}.LineWidth = 0.25;
        xlabel({[signature_names{3},' response']},'FontSize',18,'FontWeight','bold');
        legend([h11{1} h12{1}], {'FSS patients', 'healthy controls'},'Location','best','FontSize',14,'FontWeight','bold','Box','off');
%         title(contrast_names{q-2},'FontSize',20,'FontWeight','bold');
        ylim([(h12{3}.Position(2)+0.10.*h12{3}.Position(2)) (max([h11{1}.YData h12{1}.YData])+0.05.*max([h11{1}.YData h12{1}.YData]))]);
        box off
end
    
print(f5,fullfile(figspubdir,strcat(signature_names{2},'_',signature_names{3},'.png')),'-dpng','-r600');

clear p q;


% NPS subregions

for c = 1:size(contrast_names,2)
    
    f6{c}  = figure('Position', fig_position,'WindowState','maximized');
    title(contrast_names{c},'FontSize',24,'FontWeight','bold');

    for p = 1:size(regions,1)
       subplot(3,3,p)
            h13 = raincloud_plot(regions{p,1}(:,c), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
                         'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
                         'box_col_match', 1, 'line_width', 1.5);
            h14 = raincloud_plot(regions{p,2}(:,c), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
                 'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 1, 'line_width', 1.5);
            h13{1}.EdgeColor = 'none';
            h14{1}.EdgeColor = 'none';
            h13{2}.SizeData = 20;
            h14{2}.SizeData = 20;
            ax5{p} = gca;
            ax5{p}.FontSize = 12;
            ax5{p}.FontName = 'Cambria';
            ax5{p}.YAxisLocation = 'origin';
            ax5{p}.YTick = [];
            ax5{p}.LineWidth = 0.25;
            xlabel({[npsposregions_names_full{p},' response']},'FontSize',14,'FontWeight','bold');
            legend([h13{1} h14{1}], {'FSS', 'controls'},'Location','best','FontSize',10,'FontWeight','bold','Box','off');
    %         title({(npsposregions_names{p}),''},'FontSize',20,'FontWeight','bold');
            ylim([(h14{3}.Position(2)+0.10.*h14{3}.Position(2)) (max([h13{1}.YData h14{1}.YData])+0.05.*max([h13{1}.YData h14{1}.YData]))]);
            box off
    end
    
    sgtitle({contrast_names{c},''},'FontSize',20,'FontWeight','bold','FontName','Cambria');
    
    print(f6{c},fullfile(figspubdir,strcat('npssubregions_',contrast_names{c},'.png')),'-dpng','-r600');
    
end


%% ROBFIT PARCELWISE BRAIN FIGURES
%--------------------------------------------------------------------------


% NOTE: looping over elements of robfit_parcel_stats_results does not
% work for some weird reason, including using get_wh_image to get the first
% row of the dat

% EXAMPLE CODE FOR DIFFERENT CANLAB PLOTTING FUNCTIONS


% CANLAB_RESULTS_FMRIDISPLAY AND MONTAGE

% see also help for other fmridisplay methods, particularly addblobs and
% montage

obj_robust_parcelwise = canlab_results_fmridisplay(t_obj_neg_neu_group,'outline','linewidth',0.5,'montagetype','full hcp','splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
f7 = gcf;
f7.WindowState = 'maximized';
print(f7,fullfile(figspubdir,strcat('robfit_parcelwise_',contrast_names{1},'_montage.png')),'-dpng','-r600');


% RENDER_ON_SURFACE

surface_handle_1 = addbrain('coronal_slabs_5'); % help addbrain for all options
render_on_surface(t_obj_neg_neu_group,surface_handle_1,'colormap','hot');
f8 = gcf;
f8.WindowState = 'maximized';
print(f8,fullfile(figspubdir,strcat('robfit_parcelwise_',contrast_names{1},'_coronal_slabs.png')),'-dpng','-r600');

surface_handle_2 = addbrain('hires'); % help addbrain for all options
render_on_surface(t_obj_neg_neu_group,surface_handle_2,'colormap','hot');

surface_handle_3 = addbrain('right_insula_slab'); % help addbrain for all options
render_on_surface(t_obj_neg_neu_group,surface_handle_3,'colormap','hot');
f9 = gcf;
f9.WindowState = 'maximized';
print(f9,fullfile(figspubdir,strcat('robfit_parcelwise_',contrast_names{1},'_right_insula_slab.png')),'-dpng','-r600');

surface_handle_4 = addbrain('right_cutaway'); % help addbrain for all options
render_on_surface(t_obj_neg_neu_group,surface_handle_4,'colormap','hot');


% SURFACE

% can also be used for more flexibility in defining surfaces, cutaways, etc, 
% in combination with isosurf
% NOTE: see region.surface and image_vector.isosurface as well as
% canlab_canonical_brain_surface_cutaways for help

anat = fmri_data(which('keuken_2014_enhanced_for_underlay.img'), 'noverbose');
p = isosurface(anat, 'thresh', 140, 'nosmooth', 'ylim', [-Inf -30]);
p2 = isosurface(anat, 'thresh', 140, 'nosmooth', 'xlim', [-Inf 0], 'YLim', [-30 Inf]);
alpha 1 ; lightRestoreSingle; view(135, 30); colormap gray;
p3 = addbrain('limbic hires');
set(p3, 'FaceAlpha', .6, 'FaceColor', [.5 .5 .5]);
delete(p3(3)); p3(3) = [];
lightRestoreSingle;
surface_handles = [p p2 p3];

[surface_handles_1,pcl1,pcl2] = surface(t_obj_neg_neu_group,surface_handles);


% CLUSTER_SURF 

% does not seem to work for statistic_image objects, but check
% it out for atlas or region objects!
% NOTE: you can input any surface from the Github repo
% CanlabCore\canlab_canonical_brains\Canonical_brains_surfaces

% we use a region object for a single region here as example
cluster_surf(region_13, 2, 'colors', {[1 1 0]}, 'surf_freesurf_inflated_Left.mat'); view(1,360);


% PLOT INDIVIDUAL RESPONSES IN LARGE SOMATOSENSORY CLUSTER

% extract data from region 9 for each condition
for cond = 1:size(DATA_OBJ,2)
    region_9_data{cond} = extract_data(region_9,DATA_OBJ{1,cond});
    region_9_dat{cond} = region_9_data{cond}.dat;
end

% prep data for plotting
fulldat = DAT.BEHAVIOR.behavioral_data_table;
fulldat = sortrows(fulldat,'patient','ascend');
fulldat.patient(fulldat.patient == -1) = 2;
fullgroup = [fulldat.patient; fulldat.patient; fulldat.patient];
fullcondition = [ones(height(fulldat),1); 2.*ones(height(fulldat),1); 3.*ones(height(fulldat),1)];
reg9 = [region_9_dat{1}; region_9_dat{2}; region_9_dat{2}];
fullD9 = [reg9,fullcondition,fullgroup];
somatomotor_data = table(fullD9(:,1),fullD9(:,2),fullD9(:,3),'VariableNames',{'somatomotor','valence','group'});
writetable(somatomotor_data,fullfile(figspubdir,'parcelwise_data.xlsx'),'Sheet',1);
    
for i = 1:2
    for j = 1:2
    reg9_data{i,j} = fullD9(fullD9(:, 2) == i & fullD9(:, 3) ==j);
    end
end

% plot
f9  = figure('Position', fig_position,'WindowState','maximized');
h9   = rm_raincloud(reg9_data, cl, 0, 'rash');
ax9 = gca;
ax9.FontSize = 14;
ax9.FontName = 'Cambria';
ax9.XAxis.LineWidth = 1;
ax9.XAxis.Limits = [-0.4 0.51];
ax9.YAxis.LineWidth = 1;
xlabel({'somatomotor cluster activity',''},'FontSize',24,'FontWeight','bold');
ylabel({'','affective valence'},'FontSize',24,'FontWeight','bold');
yticklabels({'\fontsize{20} \bf neutral','\fontsize{20} \bf negative'});
legend([h9.l(1,1) h9.l(1,2)],{'FSS patients','healthy controls'},'Location','best','FontSize',20,'FontWeight','bold','Box','off');
for i = 1:2
    for j = 1:2
        h9.s{i, j}.SizeData = 150;
    end
end

% save
print(f9,fullfile(figspubdir,strcat('parcelwise_region9_somatomotor.png')),'-dpng','-r600');


% PLOT INDIVIDUAL RESPONSES IN SMALL POSTERIOR INSULA CLUSTER

% extract data from region 4 for each condition
for cond = 1:size(DATA_OBJ,2)
    region_4_data{cond} = extract_data(region_4,DATA_OBJ{1,cond});
    region_4_dat{cond} = region_4_data{cond}.dat;
end

% prep data for plotting
reg4 = [region_4_dat{1}; region_4_dat{2}; region_4_dat{2}];
fullD4 = [reg4,fullcondition,fullgroup];
pINS_data = table(fullD4(:,1),fullD4(:,2),fullD4(:,3),'VariableNames',{'pINS','valence','group'});
writetable(pINS_data,fullfile(figspubdir,'parcelwise_data.xlsx'),'Sheet',2);
    
for i = 1:2
    for j = 1:2
    reg4_data{i,j} = fullD4(fullD4(:, 2) == i & fullD4(:, 3) ==j);
    end
end

% plot
f4  = figure('Position', fig_position,'WindowState','maximized');
h4   = rm_raincloud(reg4_data, cl, 0, 'rash');
ax4 = gca;
ax4.FontSize = 14;
ax4.FontName = 'Cambria';
ax4.XAxis.LineWidth = 1;
ax4.XAxis.Limits = [-0.45 0.6];
ax4.YAxis.LineWidth = 1;
xlabel({'posterior insula cluster activity',''},'FontSize',24,'FontWeight','bold');
ylabel({'','affective valence'},'FontSize',24,'FontWeight','bold');
yticklabels({'\fontsize{20} \bf neutral','\fontsize{20} \bf negative'});
legend([h4.l(1,1) h4.l(1,2)],{'FSS patients','healthy controls'},'Location','best','FontSize',20,'FontWeight','bold','Box','off');
for i = 1:2
    for j = 1:2
        h4.s{i, j}.SizeData = 150;
    end
end

% save
print(f4,fullfile(figspubdir,strcat('parcelwise_region4_pINS.png')),'-dpng','-r600');


%% PDM BRAIN FIGURES
%--------------------------------------------------------------------------

% FULL MONTAGES

% single pdms
for r = 1:size(reg_all_fdr,2)
    o1 = canlab_results_fmridisplay(reg_all_fdr{1,r},'outline','linewidth',0.5,'montagetype','full hcp','splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img'); % o1 is an fmridisplay object - methods fmridisplay for help
%     o1 = legend(o1);
%     delete(o1.activation_maps{1,2}.legendhandle); % get rid of legend for contour activation map in object
    fig1 = gcf;
    fig1.WindowState = 'maximized';
    montage_full_pdm{r} = o1;
    print(fig1,fullfile(figspubdir,strcat('pdm',num2str(r),'_',contrast_names{1},'_montage_full.png')),'-dpng','-r600');
    close gcf;
    clear fig1;
end % loop pdms

% source reconstruction maps for single pdms
for s = 1:size(source_obj_j,2)
   o1b = canlab_results_fmridisplay(source_obj_j{1,s},'outline','linewidth',0.5,'montagetype','full hcp','splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img'); % o1 is an fmridisplay object - methods fmridisplay for help
%    o1b = canlab_results_fmridisplay(source_obj_j{1,s}.threshold([0 max(source_obj_j{1,s}.dat)],'raw-between'),'outline','linewidth',0.5,'montagetype','full hcp','splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img'); % o1 is an fmridisplay object - methods fmridisplay for help
%     o1b = legend(o1b);
%     delete(o1.activation_maps{1,2}.legendhandle); % get rid of legend for contour activation map in object
    fig1b = gcf;
    fig1b.WindowState = 'maximized';
    montage_source_recon_pdm{s} = o1b;
    print(fig1b,fullfile(figspubdir,strcat('pdm',num2str(s),'_source_recon_',contrast_names{1},'_montage_full.png')),'-dpng','-r600');
    close gcf;
    clear fig1b;
end % loop pdms


% BRAINSTEM RENDER FOR SINGLE PDMS

for r = 1:size(reg_all_fdr,2)
    [cmap] = cbrewer('qual','Set1',24,'pchip');
    [surfhan,~,~] = surface(reg_all_fdr{1,r},'brainstem_group');
    for y = [1:10] % surfhan(1,11) is handle for the entire brainstem
        surfhan(1,y).FaceColor = cmap(y,:);
        surfhan(1,y).FaceAlpha = 0.10;
        surfhan(1,y).FaceLighting = 'flat';
    end
    for y = 12:size(surfhan,2) % surfhan(1,11) is handle for the entire brainstem
        surfhan(1,y).FaceColor = cmap(y,:);
        surfhan(1,y).FaceAlpha = 0.10;
        surfhan(1,y).FaceLighting = 'flat';
    end
    surfhan(1,11).FaceAlpha = 0.60;
    surfhan(1,11).BackFaceLighting = 'reverselit';
    surfhan(1,11).FaceLighting = 'gouraud';
    surfhan(1,11).AmbientStrength = 0.5;
    fig2 = gcf;
    fig2.WindowState = 'maximized';
    ax2 = fig2.Children(5);
    ax2.XColor = 'none';
    ax2.YColor = 'none';
    ax2.ZColor = 'none';
    ax2.View = [140 20];
    brainstem_pdm{r} = fig2;
    print(fig2,fullfile(figspubdir,strcat('pdm',num2str(r),'_',contrast_names{1},'_brainstem.png')),'-dpng','-r600');
    close gcf;
    clear fig2;
    
end % loop pdms


% MULTIROW MONTAGE OF BOTH PDMS

o3 = canlab_results_fmridisplay([],'outline','linewidth',0.5,'montagetype','multirow',2,'overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
for m = 1:size(o3.montage,2)
    for a = 1:size(o3.montage{1,m}.axis_handles,2)
        if m < (size(o3.montage,2)/2)+1
            o3.montage{1,m}.axis_handles(a).Position(1,2) = o3.montage{1,m}.axis_handles(a).Position(1,2) - 0.10; % will depend on size of screen and figure, can be fixed more properly with screensize() function
        else
            o3.montage{1,m}.axis_handles(a).Position(1,2) = o3.montage{1,m}.axis_handles(a).Position(1,2) - 0.30; 
        end
    end
end
o3 = addblobs(o3,reg_all_fdr{1,1},'wh_montages',1:2,'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
annotation('textbox',[.5 .475 .4 .5],'String','PDM #1','FontName','Cambria','FontSize',18,'FontWeight','bold','FitBoxToText','on','EdgeColor','none');
o3 = addblobs(o3,reg_all_fdr{1,2},'wh_montages',3:4,'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
annotation('textbox',[.5 0.025 .4 .5],'String','PDM #2','FontName','Cambria','FontSize',18,'FontWeight','bold','FitBoxToText','on','EdgeColor','none');
o3 = legend(o3);
delete(o3.activation_maps{1,2}.legendhandle); % get rid of legend for contour activation map in object
fig3 = gcf;
fig3.WindowState = 'maximized';
print(fig3,fullfile(figspubdir,strcat('pdm_',contrast_names{1},'_montage_multirow.png')),'-dpng','-r600');


% CUSTOM MONTAGES OF BOTH PDMS

% first example

o4 = canlab_results_fmridisplay(reg_all_fdr{1,1},'outline','splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
for m = 1:size(o4.montage,2)
    for a = 1:size(o4.montage{1,m}.axis_handles,2)
        o4.montage{1,m}.axis_handles(a).Position(1,2) = o4.montage{1,m}.axis_handles(a).Position(1,2) + 0.20; % will depend on size of screen and figure, can be fixed more properly with screensize() function
    end
end
annotation('textbox',[.5 .475 .4 .5],'String','PDM #1','FontName','Cambria','FontSize',18,'FontWeight','bold','FitBoxToText','on','EdgeColor','none');
clear m a
o4 = canlab_results_fmridisplay(reg_all_fdr{1,2},'outline','splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
title('PDM #2');
for m = 1:size(o4.montage,2)
    for a = 1:size(o4.montage{1,m}.axis_handles,2)
        o4.montage{1,m}.axis_handles(a).Position(1,2) = o4.montage{1,m}.axis_handles(a).Position(1,2) - 0.25; % will depend on size of screen and figure, can be fixed more properly with screensize() function
    end
end
annotation('textbox',[.5 0.025 .4 .5],'String','PDM #2','FontName','Cambria','FontSize',18,'FontWeight','bold','FitBoxToText','on','EdgeColor','none');
clear m a
o4 = legend(o4);
delete(o4.activation_maps{1,2}.legendhandle); % avoid duplicate legend with different limits
o4.activation_maps{1,1}.legendhandle.Position(1,2) = o4.activation_maps{1,1}.legendhandle.Position(1,2) -0.05;
brighten(.3);
fig4 = gcf;
fig4.WindowState = 'maximized';
print(fig4,fullfile(figspubdir,strcat('pdm_',contrast_names{1},'_montage_custom.png')),'-dpng','-r600');

% second example

[fig5,~] = create_figure(strcat('pdm_',contrast_names{1},'_montage_custom'),2,10);
fig5.WindowState = 'maximized';
for c = 1:size(fig5.Children,1)
    axh5{c} = fig5.Children(c);
end
o5 = fmridisplay('overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
o5 = montage(o5,'coronal','slice_range', [-70 20], 'onerow', 'spacing', 10,'outline','existing_axes',[axh5{1,[12:20,1]}],'brighten',0.3); % weird order of axes in object, but no worries
o5 = addblobs(o5,reg_all_fdr{1,1},'wh_montages',1,'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
annotation('textbox',[0.45 .475 .4 .5],'String','PDM #1','FontName','Cambria','FontSize',18,'FontWeight','bold','FitBoxToText','on','EdgeColor','none');
o5 = montage(o5,'coronal','slice_range', [-70 20], 'onerow', 'spacing', 10,'outline','existing_axes',[axh5{1,[2:11]}],'brighten',0.3); % weird order of axes in object, but no worries
o5 = addblobs(o5,reg_all_fdr{1,2},'wh_montages',2,'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
annotation('textbox',[0.45 0.055 .4 .5],'String','PDM #2','FontName','Cambria','FontSize',18,'FontWeight','bold','FitBoxToText','on','EdgeColor','none');
enlarge_axes(gcf,1.75);
for m = 1:size(o5.montage,2)
    for a = 1:size(o5.montage{1,m}.axis_handles,2)
        o5.montage{1,m}.axis_handles(a).Position(1,1) = o5.montage{1,m}.axis_handles(a).Position(1,1) -0.05;
        if m < (size(o5.montage,2)/2)+1
            o5.montage{1,m}.axis_handles(a).Position(1,2) = o5.montage{1,m}.axis_handles(a).Position(1,2) - 0.05; % will depend on size of screen and figure, can be fixed more properly with screensize() function
        else
            o5.montage{1,m}.axis_handles(a).Position(1,2) = o5.montage{1,m}.axis_handles(a).Position(1,2) + 0.025; 
        end
    end
end
o5 = legend(o5);
delete(o5.activation_maps{1,2}.legendhandle); 
print(fig5,fullfile(figspubdir,strcat('pdm_',contrast_names{1},'_montage_custom2.png')),'-dpng','-r600');


%% SIGNATURE FIGURES
%--------------------------------------------------------------------------

% MONTAGE OF NPS, PINES & SIIPS

for r = 1:size(pat_objs,2)
    o7 = canlab_results_fmridisplay(pat_objs{1,r},'outline','linewidth',0.5,'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img'); % o7 is an fmridisplay object - methods fmridisplay for help
    fig1 = gcf;
    fig1.WindowState = 'maximized';
    montage_sigs{r} = o7;
    print(fig1,fullfile(figspubdir,strcat(pat_obj_names{r},'_montage.png')),'-dpng','-r600');
    close gcf;
    clear fig1;
end % loop pdms
