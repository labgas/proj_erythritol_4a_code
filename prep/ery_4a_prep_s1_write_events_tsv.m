%% ery_4a_prep_s1_write_events_tsv
%
% This script reads logfiles, extract the onsets and durations for
% different conditions, and write events.tsv files to the BIDS dir for
% each subject, for two different models
% 
% DEPENDENCIES
% None
% 
% INPUTS
% Presentation .log files in sourcedata dir for each subject
%
% OUTPUT
% events.tsv files for each run in BIDS dir for each subject
%
%__________________________________________________________________________
%
% author: Lukas Van Oudenhove
% date:   December, 2021
%
%__________________________________________________________________________
% @(#)% ery_4a_prep_s1_write_events_tsv.m         v1.0        
% last modified: 2021/12/23


%% DEFINE DIRECTORIES, SUBJECTS, RUNS, CONDITIONS, AND IMPORT OPTIONS
%--------------------------------------------------------------------------

ery_4a_prep_s0_define_directories;

subjs2write = {};

runnames = {'run-1','run-2','run-3','run-4','run-5','run-6'};
logfilenames = {'*_run1.log','*_run2.log','*_run3.log','*_run4.log','*_run5.log','*_run6.log'};
taskname = 'sweettaste_';
sweet_labels = {'sucrose delivery';'erythritol delivery';'sucralose delivery';'control delivery'};
swallow_rinse_labels = {'sucrose_swallowing';'erythritol_swallowing';'sucralose_swallowing';'control_swallowing'};
rating_labels = {'Sucrose','Erythritol','Sucralose','Control'};
fixation_labels = {'fixation_cross','Sucrose fixation cross','Erythritol fixation cross','Sucralose fixation cross','Control fixation cross'};
events_interest = {'sucrose','erythritol','sucralose','control'};
events_nuisance = {'swallow_rinse','rating'};

varNames = {'Trial','Event Type','Code','Time','TTime','Uncertainty','Duration','Uncertainty','ReqTime','ReqDur','Stim Type','Pair Index'};
selectedVarNames = [1:5 7 9];
varTypes = {'double','categorical','categorical','double','double','double','double','double','double','char','char','double'};
delimiter = '\t';
dataStartLine = 5;
extraColRule = 'ignore';

opts = delimitedTextImportOptions('VariableNames',varNames,...
                                'SelectedVariableNames',selectedVarNames,...
                                'VariableTypes',varTypes,...
                                'Delimiter',delimiter,...
                                'DataLines', dataStartLine,...
                                'ExtraColumnsRule',extraColRule); 


%% LOOP OVER SUBJECTS TO READ LOGFILES, CREATE TABLE WITH ONSETS AND DURATIONS, AND SAVE AS EVENTS.TSV TO BIDSSUBJDIRS
%-------------------------------------------------------------------------------------------------------------------
if ~isempty(subjs2write)
    [C,ia,~] = intersect(sourcesubjs,subjs2write);
    
    if ~isequal(C',subjs2write)
        error('subject defined in subjs2write not present in sourcedir, please check before proceeding');
    
    else
        
        for sub = ia'
        
        % DEFINE SUBJECT LEVEL DIRS
        subjsourcedir = sourcesubjdirs{sub};
        subjBIDSdir = fullfile(BIDSsubjdirs{sub},'func');

            % LOOP OVER RUNS
            for run = 1:size(logfilenames,2)
                
                logfilename = dir(fullfile(subjsourcedir,'logfiles',logfilenames{run}));
                logfilename = char(logfilename(:).name);
                logfilepath = fullfile(subjsourcedir,'logfiles',logfilename);
                
                if ~isfile(logfilepath)
                    warning(strcat(logfilepath,' does not exist, please check'));
                    continue
                
                elseif size(logfilepath,1) > 1
                    error(strcat('more than one logfile with run index_',run,' for subject_',sourcesubjs{sub}, ' please check'))
                
                else
                    log = readtable(logfilepath,opts);
                    log = log(~isnan(log.Trial),:);
                    time_zero = log.Time(log.Trial == 0 & log.EventType == 'Pulse'); % time for onsets and durations is counted from the first scanner pulse onwards
                        
                        if size(time_zero,1) > 1
                            error('ambiguity about time zero in %s%s, please check logfile',subjs{sub},logfilenames{run});
                        end
                        
                    log.TimeZero = log.Time - time_zero;
                    log.onset = log.TimeZero ./ 10000; % convert to seconds 
                    log(log.EventType == 'Pulse',:) = [];
                    log.trial_type = cell(height(log),1);
                        
                        for k = 1:height(log)
                            if ismember(log.Code(k),swallow_rinse_labels)
                                log.trial_type{k} = events_nuisance{1};
                            
                            elseif ismember(log.Code(k),rating_labels)
                                log.trial_type{k} = events_nuisance{2};
                            
                            elseif ismember(log.Code(k),sweet_labels)
                                idx = (log.Code(k) == sweet_labels);
                                log.trial_type{k} = events_interest{idx'};
                            
                            elseif ismember(log.Code(k),fixation_labels)
                                log.trial_type{k} = 'fixation';
                            
                            else log.trial_type{k} = '';
                            
                            end
                        end
                    
                    log.trial_type = categorical(log.trial_type);
                    log = log(~isundefined(log.trial_type),:);        
                    log = removevars(log,{'Trial','EventType','Code','Time','TTime','Duration','ReqTime','TimeZero'}); % get rid of junk variables from logfile we don't need
                    log.duration = zeros(height(log),1);
                        
                        for m = 1:height(log)
                            if ~isequal(log.trial_type(m),'fixation')
                                log.duration(m) = log.onset(m+1) - log.onset(m);
                            
                            else log.duration(m) = NaN;
                            
                            end
                        end
                        
                    log = log(~isnan(log.duration),:);
                        
                    filename = strcat(subjBIDSdir,'/',sourcesubjs{sub},'_task-',taskname,runnames{run},'_model_1_events.tsv');
                    writetable(log,filename,'Filetype','text','Delimiter','\t');
                    clear logfile log time_zero filename

                end % if loop checking whether logfile exists

            end % for loop runs
            
        end % for loop subjects
    
    end % if loop checking subjs2write present in sourcesubjs        
    
else
    
    if ~isequal(sourcesubjs, BIDSsubjs)
        error('subjects in sourcedir and BIDSdir do not match, please check before proceeding');
        
    else
    
        for sub = 1:size(sourcesubjs,1)

        % DEFINE SUBJECT LEVEL DIRS
        subjsourcedir = sourcesubjdirs{sub};
        subjBIDSdir = fullfile(BIDSsubjdirs{sub},'func');

            % LOOP OVER RUNS
            for run = 1:size(logfilenames,2)
                
                logfilename = dir(fullfile(subjsourcedir,'logfiles',logfilenames{run}));
                logfilename = char(logfilename(:).name);
                logfilepath = fullfile(subjsourcedir,'logfiles',logfilename);
                
                if ~isfile(logfilepath)
                    warning(strcat(logfilepath,' does not exist, please check'));
                    continue
                
                elseif size(logfilepath,1) > 1
                    error(strcat('more than one logfile with run index_',run,' for subject_',sourcesubjs{sub}, ' please check'))
                
                else
                    log = readtable(logfilepath,opts);
                    log = log(~isnan(log.Trial),:);
                    time_zero = log.Time(log.Trial == 0 & log.EventType == 'Pulse'); % time for onsets and durations is counted from the first scanner pulse onwards
                        
                        if size(time_zero,1) > 1
                            error('ambiguity about time zero in %s%s, please check logfile',subjs{sub},logfilenames{run});
                        end
                        
                    log.TimeZero = log.Time - time_zero;
                    log.onset = log.TimeZero ./ 10000; % convert to seconds 
                    log(log.EventType == 'Pulse',:) = [];
                    log.trial_type = cell(height(log),1);
                        
                        for k = 1:height(log)
                            
                            if ismember(log.Code(k),swallow_rinse_labels)
                                log.trial_type{k} = events_nuisance{1};
                            
                            elseif ismember(log.Code(k),rating_labels)
                                log.trial_type{k} = events_nuisance{2};
                            
                            elseif ismember(log.Code(k),sweet_labels)
                                idx = (log.Code(k) == sweet_labels);
                                log.trial_type{k} = events_interest{idx'};
                            
                            elseif ismember(log.Code(k),fixation_labels)
                                log.trial_type{k} = 'fixation';
                            
                            else log.trial_type{k} = '';
                            
                            end
                        end
                    
                    log.trial_type = categorical(log.trial_type);
                    log = log(~isundefined(log.trial_type),:);        
                    log = removevars(log,{'Trial','EventType','Code','Time','TTime','Duration','ReqTime','TimeZero'}); % get rid of junk variables from logfile we don't need
                    log.duration = zeros(height(log),1);
                        
                        for m = 1:height(log)
                            
                            if ~isequal(log.trial_type(m),'fixation')
                                log.duration(m) = log.onset(m+1) - log.onset(m);
                            
                            else log.duration(m) = NaN;
                            
                            end
                        end
                        
                    log = log(~isnan(log.duration),:);
                        
                    filename = strcat(subjBIDSdir,'/',sourcesubjs{sub},'_task-',taskname,runnames{run},'_model_1_events.tsv');
                    writetable(log,filename,'Filetype','text','Delimiter','\t');
                    clear logfile log time_zero filename

                end % if loop checking whether logfile exists

            end % for loop runs

        end % for loop subjects
    
    end % if loop checking sourcesubjs == BIDSsubjs
    
end % if loop checking writing option