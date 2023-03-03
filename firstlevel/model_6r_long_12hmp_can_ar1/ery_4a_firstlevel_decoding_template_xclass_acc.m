%% LaBGAScore_firstlevel_decoding_template_xclass_acc.m
%
% adapted from decoding_template.m

%% SET OPTIONS

results_suffix = 'all_conds';
conds2include = [1:4]; % indices of conditions in DSGN.conditions (or SPM.Sess(x).U.name); assumes all conditions are present in all runs!
mask_name = 'gray_matter_mask_sparse.img'; % absolute path to mask, or name of mask file if already on Matlab path

%% PREP WORK

% set analysis name

analysis_name = 'xclass_acc';

% check whether LaBGAScore_firstlevel_s1_options_dsgn_struct has been run

if ~exist('DSGN','var')
    warning('\nDSGN structure variable not found in Matlab workspace, running ery_4a_firstlevel_m6r_s1_options_dsgn_struct before proceeding')
    ery_4a_firstlevel_m6r_s1_options_dsgn_struct;
    cd(rootdir);
else
    cd(rootdir);
end

% get first level dir info

firstmodeldir = DSGN.modeldir;

firstlist = dir(fullfile(firstmodeldir,'sub-*'));
firstsubjs = cellstr(char(firstlist(:).name));

    for firstsub = 1:size(firstsubjs,1)
        firstsubjdirs{firstsub,1} = fullfile(firstlist(firstsub).folder,firstlist(firstsub).name);
    end
    
clear firstlist

% set TDT dir info

firstmodelTDTdir = fullfile(firstmodeldir,'TDT');
if ~isfolder(firstmodelTDTdir)
    mkdir(firstmodelTDTdir);
end

firstmodelTDTdir_mask = fullfile(firstmodelTDTdir,'mask');
if ~isfolder(firstmodelTDTdir_mask)
    mkdir(firstmodelTDTdir_mask);
end

% pre-allocate loop vars

spmdotmats = cell(1,firstsub);
results_combined = cell(1,firstsub);
group_results = cell(1,size(conds2include,2));

for cond = 1:size(group_results,2)
    group_results{cond} = [];
end

% confusion_matrices = zeros(size(conds2include,2),size(conds2include,2),firstsub);
confusion_matrices = [];

%% LOOP TEMPLATE CODE OVER SUBJECTS

for sub = 1:firstsub
    
    spmdotmats{sub} = load(fullfile(firstsubjdirs{sub},'SPM.mat'));

% This script is a template that can be used for a decoding analysis on 
% brain image data. It is for people who have betas available from an 
% SPM.mat and want to automatically extract the relevant images used for
% classification, as well as corresponding labels and decoding chunk numbers
% (e.g. run numbers). If you don't have this available, then use
% decoding_template_nobetas.m

% Make sure the decoding toolbox and your favorite software (SPM or AFNI)
% are on the Matlab path (e.g. addpath('/home/decoding_toolbox') )
% TDT
% addpath('$ADD FULL PATH TO TDT TOOLBOX AS STRING OR MAKE THIS LINE A COMMENT IF IT IS ALREADY$')
% assert(~isempty(which('decoding_defaults.m', 'function')), 'TDT not found in path, please add')
% SPM/AFNI
% addpath('$ADD FULL PATH TO SPM/AFNI (if you need them) AS STRING OR MAKE THIS LINE A COMMENT IF IT IS ALREADY$')
% assert((~isempty(which('spm.m', 'function')) || ~isempty(which('BrikInfo.m', 'function'))) , 'Neither SPM nor AFNI found in path, please add (or remove this assert if you really dont need to read brain images)')

% Set defaults
cfg = decoding_defaults;

% Set the analysis that should be performed (default is 'searchlight')
cfg.analysis = 'wholebrain'; % standard alternatives: 'wholebrain', 'ROI' (pass ROIs in cfg.files.mask, see below)

% Set the filename of your brain mask (or your ROI masks as cell matrix) 
% for searchlight or wholebrain e.g. 'c:\exp\glm\model_button\mask.img' OR 
% for ROI e.g. {'c:\exp\roi\roimaskleft.img', 'c:\exp\roi\roimaskright.img'}
% You can also use a mask file with multiple masks inside that are
% separated by different integer values (a "multi-mask")
if sub == 1
    mask_img = fmri_mask_image(which(mask_name),'noverbose');
    target = fmri_data(fullfile(firstsubjdirs{sub},'beta_0001.nii'),'noverbose');
    mask2write = resample_space(mask_img,target);
    if contains(mask_name,'.nii')
        mask = fullfile(firstmodelTDTdir_mask,mask_name);
    else
        mask = fullfile(firstmodelTDTdir_mask,[mask_name(1:end-4) '.nii']);
    end
    write(mask2write,'fname',mask,'overwrite');
end

cfg.files.mask = mask;

% Set the output directory where data will be saved, e.g. 'c:\exp\results\buttonpress'
cfg.results.dir = fullfile(firstsubjdirs{sub},'TDT',analysis_name,[cfg.analysis '_' mask_name '_' results_suffix]);

if ~isfolder(cfg.results.dir)
    mkdir(cfg.results.dir);
end

% Set the filepath where your SPM.mat and all related betas are, e.g. 'c:\exp\glm\model_button'
beta_loc = firstsubjdirs{sub};

% Set the label names to the regressor names which you want to use for 
% decoding, e.g. 'button left' and 'button right'
% don't remember the names? -> run display_regressor_names(beta_loc)
% infos on '*' (wildcard) or regexp -> help decoding_describe_data
for label = conds2include
     labelnames{label} = char(spmdotmats{1,sub}.SPM.Sess(1).U(label).name); % labels in SPM.Sess(1).U.name, which weirdly have a blank space added to names in DSGN.conditions
end

%% Set additional parameters
% Set additional parameters manually if you want (see decoding.m or
% decoding_defaults.m). Below some example parameters that you might want 
% to use a searchlight with radius 12 mm that is spherical:

% cfg.searchlight.unit = 'mm';
% cfg.searchlight.radius = 12; % if you use this, delete the other searchlight radius row at the top!
% cfg.searchlight.spherical = 1;
% cfg.verbose = 2; % you want all information to be printed on screen
% cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; 

% Enable scaling min0max1 (otherwise libsvm can get VERY slow)
% if you dont need model parameters, and if you use libsvm, use:
cfg.scale.method = 'min0max1';
cfg.scale.estimation = 'all'; % scaling across all data is equivalent to no scaling (i.e. will yield the same results), it only changes the data range which allows libsvm to compute faster

% if you like to change the decoding software (default: libsvm):
% cfg.decoding.software = 'liblinear'; % for more, see decoding_toolbox\decoding_software\. 
% Note: cfg.decoding.software and cfg.software are easy to confuse.
% cfg.decoding.software contains the decoding software (standard: libsvm)
% cfg.software contains the data reading software (standard: SPM/AFNI)

% Some other cool stuff
% Check out 
%   combine_designs(cfg, cfg2)
% if you like to combine multiple designs in one cfg.

%% Decide whether you want to see the searchlight/ROI/... during decoding
cfg.plot_selected_voxels = 500; % 0: no plotting, 1: every step, 2: every second step, 100: every hundredth step...

%% Add additional output measures if you like
% See help decoding_transform_results for possible measures

cfg.results.output = {'confusion_matrix'}; % 'accuracy_minus_chance' by default

% You can also use all methods that start with "transres_", e.g. use
%   cfg.results.output = {'SVM_pattern'};
% will use the function transres_SVM_pattern.m to get the pattern from 
% linear svm weights (see Haufe et al, 2015, Neuroimage)

%% Nothing needs to be changed below for a standard leave-one-run out cross
%% validation analysis.

% The following function extracts all beta names and corresponding run
% numbers from the SPM.mat
regressor_names = design_from_spm(beta_loc);

% Extract all information for the cfg.files structure (labels will be [1 -1] if not changed above)
cfg = decoding_describe_data(cfg,labelnames,conds2include ,regressor_names,beta_loc);

% This creates the leave-one-run-out cross validation design:
cfg.design = make_design_cv(cfg); 

% Run decoding
results = decoding(cfg);
results_combined{sub} = results;
confusion_matrices = cat(3, confusion_matrices, results.confusion_matrix.output{1});

%% Make figure
figure;
heatmap(results.confusion_matrix.output{1}, 'Colormap', jet);
figtitle = firstsubjs{sub};
set(gca,'Title',figtitle);
drawnow, snapnow;

%% Extract accuracy for each condition from confusion matrix
for cond = 1:size(conds2include,2)
    group_results{cond} = [group_results{cond};results.confusion_matrix.output{1}(cond,cond)]; 
end

end


%% CALCULATE AVERAGE AND DO GROUP LEVEL STATS

mean_matrix = mean(confusion_matrices,3);
std_matrix = std(confusion_matrices,[],3);
figure;
heatmap(mean_matrix, 'Colormap', jet);
figure;
heatmap(std_matrix, 'Colormap', jet);
group_results = cell2mat(group_results);
group_results_tbl = array2table(group_results);
varnames = group_results_tbl.Properties.VariableNames;
for var = 1:size(varnames,2)
    group_results_tbl.Properties.VariableNames{var} = ['c' num2str(var)];
    group_results_tbl.Properties.VariableDescriptions{var} = labelnames{var};
end
conds = table([1:size(labelnames,2)]','VariableNames',{'Conditions'});
rm = fitrm(group_results_tbl,[group_results_tbl.Properties.VariableNames{1} '-' group_results_tbl.Properties.VariableNames{end} ' ~ 1'],'WithinDesign',conds);
ranovatbl = ranova(rm);
margmeanstbl = margmean(rm,'Conditions');
posthoctbl = multcompare(rm,'Conditions'); % default Tukey-Kramer corrected
figure;
plot = plot(rm);
figure;
margmeansplot = plotprofile(rm,'Conditions');