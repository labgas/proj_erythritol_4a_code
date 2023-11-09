%% ery_4a_secondlevel_m6_s8_prep_4_apply_signatures
%
%
% USAGE
%
% This script 
% 1) calculates selected signature responses for conditions and
%       contrasts included in DAT, and saves them to new fields in DAT.
% 2) calculated responses per NPS subregion if NPS is selected, can/should
%       later be expanded to all signatures?
%
%
% OUTPUT
%
% These fields contain data tables whose columns are conditions or contrasts, 
% with variable names based on DAT.conditions or DAT.contrastnames, 
% but with spaces replaced with underscores:
% DAT.SIG_conditions.(myscaling_sigs).(similarity_metric_sigs).(keyword_sigs).signature_name
% DAT.SIG_contrasts.(myscaling_sigs).(similarity_metric_sigs).(keyword_sigs).signature_name
% 
%
% OPTIONS
%
% NOTE: 
% defaults are specified in a2_set_default_options, but can be changed below
% in case you want to add signature responses calculated with different
% options, e.g. scaling, and save your new version of the script with a
% letter index
%
% myscaling_sigs = 'raw'/'scaled';
% similarity_metric_sigs = 'dotproduct/'cosine_similarity','correlation';
% keyword_sigs = 'all'/any option from load_image_set;
% 
%
%__________________________________________________________________________
%
% revamped by: Lukas Van Oudenhove
% date:   Leuven, January, 2023
%
%__________________________________________________________________________
% @(#)% prep_4_apply_signatures_and_save.m         v2.0
% last modified: 2023/01/23


%% GET PATHS AND OPTIONS AND CHECK OPTIONS
% -------------------------------------------------------------------------

% GET MODEL-SPECIFIC PATHS AND OPTIONS

ery_4a_secondlevel_m6_s0_a_set_up_paths_always_run_first;

% NOTE: CHANGE THIS TO THE MODEL-SPECIFIC VERSION OF THIS SCRIPT
% NOTE: THIS WILL ALSO AUTOMATICALLY CALL A2_SET_DEFAULT_OPTIONS

% GET DEFAULT OPTIONS IF NOT SET IN A2_SET_DEFAULT_OPTIONS

options_needed = {'myscaling_sigs','similarity_metric_sigs','keyword_sigs'};  % Options we are looking for. Set in a2_set_default_options
options_exist = cellfun(@exist, options_needed); 

option_default_values = {'raw','dotproduct','all'};          % defaults if we cannot find info in a2_set_default_options at all 

if ~iscell(keyword_sigs)
    plugin_get_options_for_analysis_script
end

% SET CUSTOM OPTIONS

% NOTE: only specify if you want to run multiple versions of your model with different options
% than the defaults you set in your model-specific version of a2_set_default_options.m
% 
% use_scaled_images_sigs = true/false;
% similarity_metric_sigs = 'keyword';
% keyword_sigs = 'keyword';


%% LOAD NECESSARY VARIABLES IF NEEDED
% -------------------------------------------------------------------------

if ~exist('DSGN','var') || ~exist('DAT','var')
    
    load(fullfile(resultsdir,'image_names_and_setup.mat'));
    
end

if ~exist('DATA_OBJ','var') || ~exist('DATA_OBJsc','var')
    
    load(fullfile(resultsdir,'data_objects.mat'));
    load(fullfile(resultsdir,'data_objects_scaled.mat'));
    
end

if ~exist('DATA_OBJ_CON','var') || ~exist('DATA_OBJ_CONsc','var') || ~exist('DATA_OBJ_CONscc','var')
    
    load(fullfile(resultsdir,'contrast_data_objects.mat'));
    
end


%% CHECK REQUIRED DAT FIELDS
% -------------------------------------------------------------------------

% List required fields in DAT, in cell array
    
required_fields = {'contrastnames', 'contrasts' 'contrastcolors', 'conditions', 'colors'};

ok_to_run = plugin_check_required_fields(DAT, required_fields); % Checks and prints warnings
if ~ok_to_run
    return
end


%% SELECTED SIGNATURES
% -------------------------------------------------------------------------

fprintf('\n\n');
    if ~iscell(keyword_sigs)
        printhdr(['APPLYING SIGNATURE(S) ', upper(keyword_sigs), ' ON ', upper(myscaling_sigs), ' CONDITIONS AND CONTRASTS, SIMILARITY METRIC ', similarity_metric_sigs]);
    else
        [~,signame] = fileparts(char(keyword_sigs));
        printhdr(['APPLYING SIGNATURE(S) ', upper(signame), ' ON ', upper(myscaling_sigs), ' CONDITIONS AND CONTRASTS, SIMILARITY METRIC ', similarity_metric_sigs]);
    end
fprintf('\n\n');

switch myscaling_sigs
    
    case 'raw'
        
        data_object_conds = DATA_OBJ;
        data_object_conts = DATA_OBJ_CON;
        
    case 'scaled'
        
        data_object_conds = DATA_OBJsc;
        data_object_conts = DATA_OBJ_CONsc;
        
end

% CONDITIONS
% ----------

    if ~iscell(keyword_sigs)
        DAT.SIG_conditions.(myscaling_sigs).(similarity_metric_sigs).(keyword_sigs) = apply_all_signatures(data_object_conds, 'conditionnames', DAT.conditions, 'similarity_metric', similarity_metric_sigs, 'image_set', keyword_sigs);
    else
        DAT.SIG_conditions.(myscaling_sigs).(similarity_metric_sigs).(signame) = apply_all_signatures(data_object_conds, 'conditionnames', DAT.conditions, 'similarity_metric', similarity_metric_sigs, 'image_set', keyword_sigs);
    end

% CONTRASTS
% ---------
    if ~iscell(keyword_sigs)
        DAT.SIG_contrasts.(myscaling_sigs).(similarity_metric_sigs).(keyword_sigs) = apply_all_signatures(data_object_conts, 'conditionnames', DAT.contrastnames, 'similarity_metric', similarity_metric_sigs, 'image_set', keyword_sigs);
    else
        DAT.SIG_contrasts.(myscaling_sigs).(similarity_metric_sigs).(signame) = apply_all_signatures(data_object_conts, 'conditionnames', DAT.contrastnames, 'similarity_metric', similarity_metric_sigs, 'image_set', keyword_sigs);
    end

    
%% NPS SUBREGIONS
% -------------------------------------------------------------------------

if contains(keyword_sigs,'nps') || isequal(keyword_sigs,'all')

    nr_conds = size(DAT.conditions,2);

    % subregion names
    posnames = {'vermis' 'rIns' 'rV1' 'rThal' 'lIns' 'rdpIns' 'rS2_Op' 'dACC'};
    negnames = {'rLOC' 'lLOC' 'rpLOC' 'pgACC' 'lSTS' 'rIPL' 'PCC'};

    DAT.NPSsubregions.posnames = posnames;
    DAT.NPSsubregions.negnames = negnames;

    printhdr('Extracting NPS, adding to DAT')

    % CONDITIONS
    % ----------

    % NPS
    
    for i = 1:nr_conds

        switch similarity_metric_sigs
            
            case 'dotproduct'

                [DAT.npsresponse(i), ~, ~, DAT.NPSsubregions.npspos_by_region(i), DAT.NPSsubregions.npsneg_by_region(i)] = apply_nps(data_object_conds{i}, 'noverbose', 'notables');
                
            case 'cosine_similatiry'

                [DAT.npsresponse(i), ~, ~, DAT.NPSsubregions.npspos_by_region(i), DAT.NPSsubregions.npsneg_by_region(i)] = apply_nps(data_object_conds{i}, 'noverbose', 'notables', similarity_metric_sigs);

        end

    end

    % NPS subregions
    
    printhdr('Extracting NPS Subregions, adding to DAT.NPSsubregions')

    clear posdat negdat spos sneg xx
    
    for i = 1:nr_conds

        % Get averages
        DAT.NPSsubregions.posdat{i} = nanmean(DAT.NPSsubregions.npspos_by_region{i})'; % mean across subjects
        DAT.NPSsubregions.stepos{i} = ste(DAT.NPSsubregions.npspos_by_region{i})'; % ste

        DAT.NPSsubregions.negdat{i} = nanmean(DAT.NPSsubregions.npsneg_by_region{i})'; % mean across subjects
        DAT.NPSsubregions.steneg{i} = ste(DAT.NPSsubregions.npsneg_by_region{i})'; % ste

    end

    
    % CONTRASTS
    % ---------

    printhdr('Defining NPS contrasts, adding to DAT')

    nr_conts = size(DAT.contrasts, 1);

    DAT.npscontrasts = {};

    for c = 1:nr_conts
        
        mycontrast = DAT.contrasts(c, :);
        wh = find(mycontrast);

        DAT.npscontrasts{c} = cat(2, DAT.npsresponse{wh}) * mycontrast(wh)';

        % subregions
        DAT.NPSsubregions.npspos_by_region_contrasts{c} = zeros(size(DAT.NPSsubregions.npspos_by_region{wh(1)}));
        DAT.NPSsubregions.npsneg_by_region_contrasts{c} = zeros(size(DAT.NPSsubregions.npsneg_by_region{wh(1)}));

        for j = 1:length(wh)

            DAT.NPSsubregions.npspos_by_region_contrasts{c} = DAT.NPSsubregions.npspos_by_region_contrasts{c} + DAT.NPSsubregions.npspos_by_region{wh(j)} * mycontrast(wh(j));
            DAT.NPSsubregions.npsneg_by_region_contrasts{c} = DAT.NPSsubregions.npsneg_by_region_contrasts{c} + DAT.NPSsubregions.npsneg_by_region{wh(j)} * mycontrast(wh(j));

        end
    end

    for i = 1:nr_conts

        % Get averages
        DAT.NPSsubregions.posdat_contrasts{i} = nanmean(DAT.NPSsubregions.npspos_by_region_contrasts{i})'; % mean across subjects
        DAT.NPSsubregions.stepos_contrasts{i} = ste(DAT.NPSsubregions.npspos_by_region_contrasts{i})'; % ste

        DAT.NPSsubregions.negdat_contrasts{i} = nanmean(DAT.NPSsubregions.npsneg_by_region_contrasts{i})'; % mean across subjects
        DAT.NPSsubregions.steneg_contrasts{i} = ste(DAT.NPSsubregions.npsneg_by_region_contrasts{i})'; % ste

    end

end


%% SAVE RESULTS
% -------------------------------------------------------------------------

% NOTE: you may get permission errors here, in which case temporarily
% unannexing image_names_and_setup.mat using the command git annex unannex
% helps, followed by datalad saving the appended file

fprintf('\n\n');
printhdr('SAVING SIGNATURE RESPONSES TO DAT');
fprintf('\n\n');

savefilename = fullfile(resultsdir, 'image_names_and_setup.mat');
save(savefilename, 'DAT', '-append');
