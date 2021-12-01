%% WAD_prep_s2_write_events_tsv
%
% This script reads logfiles, extract the onsets and durations for
% different conditions, and write events.tsv files in the rawdata dir for
% each subject, for two different models
% 
% DEPENDENCIES
% None
% 
% INPUTS
% Presentation .log files in sourcedata dir for each subject
%
% OUTPUT
% events.tsv files for each run in rawdata dir for each subject
%
%__________________________________________________________________________
%
% author: Lukas Van Oudenhove and Iris Coppieters
% date:   April, 2021
%
%__________________________________________________________________________
% @(#)% WAD_prep_s2_write_events_tsv.m         v1.0        
% last modified: 2021/04/06


%% DEFINE DIRECTORIES, SUBJECTS, RUNS, AND CONDITIONS
%--------------------------------------------------------------------------

rootdir = 'C:\Users\lukas\Dropbox (Dartmouth College)\4.Iris_WAD_MRI';
sourcedir = fullfile(rootdir,'sourcedata');
rawdir = fullfile(rootdir,'rawdata');

tasknames = categorical({'picA'; 'picB'});
logfilenames = {'*partA*.log';'*partB*.log'};
event_cats = categorical({'high_fear';'moderate_fear';'neutral_fear';'imagine'});
high_fear_pics = {'pic_fear1','pic_fear2','pic_fear3','pic_fear4','pic_fear5'};
mod_fear_pics = {'pic_fear6','pic_fear7','pic_fear8','pic_fear9','pic_fear10'};
neutral_pics = {'pic_neutral1','pic_neutral2','pic_neutral3','pic_neutral4','pic_neutral5'};
imagine_cats = categorical({'imagine_high';'imagine_moderate';'imagine_neutral'});

subjs = dir(fullfile(sourcedir,'sub-*'));
idx = [subjs.isdir]';
subjs = {subjs(idx).name}';

subjs_raw = dir(fullfile(rawdir,'sub-*'));
idx_raw = [subjs_raw.isdir]';
subjs_raw = {subjs_raw(idx_raw).name}';

if ~isequal(subjs,subjs_raw)
    error('subject directories in sourcedata and rawdata do not match; check and correct before proceeding');
end

clear idx idx_raw subjs_raw % script hygiene


%% LOOP OVER SUBJECTS TO READ LOGFILES, CREATE TABLE WITH ONSETS AND DURATIONS, AND SAVE AS EVENTS.TSV TO SUBJRAWDIR
%-------------------------------------------------------------------------------------------------------------------

for i = 1:size(subjs,1)
    
    % DEFINE SUBJECT LEVEL DIRS
    subjsourcedir = fullfile(sourcedir,'\',subjs{i});
    subjrawdir = fullfile(rawdir,'\',subjs{i},'\ses-1\func');
    
    % LOOP OVER RUNS
    for j = 1:size(tasknames,1)
        % Model 1
        % NOTE: models the picture and imagine periods together as one long event
        logfile = ls(fullfile(subjsourcedir,logfilenames{j}));
        if ~isfile(fullfile(subjsourcedir,logfile))
            warning(strcat(logfilenames{j},' does not exist for_',subjs{i},', please check'));
            continue
        else
            log = readtable(fullfile(subjsourcedir,logfile),'FileType','text','Delimiter','tab');
            log.EventType = categorical(log.EventType);
            log.Code = categorical(log.Code);
            time_zero = log.Time(log.Trial == 1 & log.EventType == 'Pulse'); % time for onsets and durations is counted from the first scanner pulse onwards
                if size(time_zero,1) > 1
                    error('ambiguity about time zero in %s%s, please check logfile',subjs{i},logfilenames{j});
                end
            log.TimeZero = log.Time - time_zero;
            log.onset = log.TimeZero ./ 10000; % convert to seconds (which equals TRs in our particular case). This is important for spm input. 
            log(log.EventType == 'Pulse',:) = [];
            log2 = log;
                for k = 1:height(log)
                    if ismember(log.Code(k),high_fear_pics)
                        log.trial_type(k) = event_cats(1);
                    elseif ismember(log.Code(k),mod_fear_pics)
                        log.trial_type(k) = event_cats(2);
                    elseif ismember(log.Code(k),neutral_pics)
                        log.trial_type(k) = event_cats(3);
                    end
                end
            log = removevars(log,{'Subject','Trial','EventType','Code','Time','TTime','Uncertainty','Duration','Uncertainty_1','ReqTime','ReqDur','StimType','PairIndex','TimeZero'}); % get rid of junk variables from logfile we don't need
            log.duration = zeros(height(log),1);
            log.duration(isundefined(log.trial_type)) = NaN;
                for m = 1:height(log)
                    if ~isundefined(log.trial_type(m))
                        log.duration(m) = log.onset(m+2) - log.onset(m); % we want to include the duration of the imagine period in this model
                    end
                end
            log(isundefined(log.trial_type),:) = [];
            filename = strcat(subjrawdir,'\',subjs{i},'_ses-1_task-',char(tasknames(j)),'_model_1_events.tsv');
            writetable(log,filename,'Filetype','text','Delimiter','\t');
            clear logfile log time_zero filename

            % Model 2
            % NOTE: models picture and imagine periods as separate events
                for k = 1:height(log2)
                    if ismember(log2.Code(k),high_fear_pics)
                        log2.trial_type(k) = event_cats(1);
                    elseif ismember(log2.Code(k),mod_fear_pics)
                        log2.trial_type(k) = event_cats(2);
                    elseif ismember(log2.Code(k),neutral_pics)
                        log2.trial_type(k) = event_cats(3);
                    elseif isequal(log2.Code(k),event_cats(4)) % we take advantage of the fact that imagine always follows a pic
                        if ismember(log2.Code(k-1),high_fear_pics)
                            log2.trial_type(k) = imagine_cats(1);
                        elseif ismember(log2.Code(k-1),mod_fear_pics)
                            log2.trial_type(k) = imagine_cats(2);
                        elseif ismember(log2.Code(k-1),neutral_pics)
                            log2.trial_type(k) = imagine_cats(3);
                        else
                            error('imagine row does not follow a picture row in %s%s, please check logfile',subjs{i},logfilenames{j});
                        end
                    end
                end
            log2 = removevars(log2,{'Subject','Trial','EventType','Code','Time','TTime','Uncertainty','Duration','Uncertainty_1','ReqTime','ReqDur','StimType','PairIndex','TimeZero'}); % get rid of junk variables from logfile we don't need
            log2.duration = zeros(height(log2),1);
            log2.duration(isundefined(log2.trial_type)) = NaN;
                for m = 1:height(log2)
                    if ~isundefined(log2.trial_type(m))
                        log2.duration(m) = log2.onset(m+1) - log2.onset(m);
                    end
                end
            log2(isundefined(log2.trial_type),:) = [];
            filename2 = strcat(subjrawdir,'\',subjs{i},'_ses-1_task-',char(tasknames(j)),'_model_2_events.tsv');
            writetable(log2,filename2,'Filetype','text','Delimiter','\t');
            clear log2 filename2
        end % if loop checking whether logfile exists
    end % for loop runs
end % for loop subjects