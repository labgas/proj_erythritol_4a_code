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
% last modified: 2021/12/01


%% DEFINE DIRECTORIES, SUBJECTS, RUNS, AND CONDITIONS
%--------------------------------------------------------------------------

ery_4a_prep_s0_define_directories;

subjs2write = {};

runnames = {'run-1','run-2','run-3','run-4','run-5','run-6'};
logfilenames = {'*_run1.log','*_run2.log','*_run3.log','*_run4.log','*_run5.log','*_run6.log'};
taskname = 'sweettaste_';
event_cats = categorical({'sucrose delivery';'erythritol delivery';'sucralose delivery';'neutral solution delivery';'rating'});
rating_labels = {'Sucrose','Erythritol','Sucralose','Control'};
% mod_fear_pics = {'pic_fear6','pic_fear7','pic_fear8','pic_fear9','pic_fear10'};
% neutral_pics = {'pic_neutral1','pic_neutral2','pic_neutral3','pic_neutral4','pic_neutral5'};
% imagine_cats = categorical({'imagine_high';'imagine_moderate';'imagine_neutral'});


%% LOOP OVER SUBJECTS TO READ LOGFILES, CREATE TABLE WITH ONSETS AND DURATIONS, AND SAVE AS EVENTS.TSV TO BIDSSUBJDIRS
%-------------------------------------------------------------------------------------------------------------------
if ~isempty(subjs2write)
    [C,ia,~] = intersect(sourcesubjs,subjs2write);
    
    if ~isequal(C',subjs2write)
        error('subject defined in subjs2write not present in sourcedir, please check before proceeding');
    else
        for sub = ia'
        end    
    end 
    
else
    
    if ~isequal(sourcesubjs, BIDSsubjs)
        error('subjects in sourcedir and BIDSdir do not match, please check before proceeding');
        
    else
    
        for sub = 1:size(sourcesubjs,1)

        % DEFINE SUBJECT LEVEL DIRS
        subjsourcedir = sourcesubjdirs{sub};
    %     subjlogfiles = dir(fullfile(subjsourcedir,'/logfiles/*.log'));
    %     subjlogfiles = cellstr(char(subjlogfiles(:).name));
        subjBIDSdir = BIDSsubjdirs{sub};

            % LOOP OVER RUNS
            for run = 1:size(logfilenames,1)
                % Model 1
                % NOTE: models the picture and imagine periods together as one long event
                logfilename = dir(fullfile(subjsourcedir,'logfiles',logfilenames{run}));
                logfilename = char(logfilename(:).name);
                logfilepath = fullfile(subjsourcedir,'logfiles',logfilename);
                
                if ~isfile(logfilepath)
                    warning(strcat(logfilepath,' does not exist, please check'));
                    continue
                
                elseif size(logfile,1) > 1
                    error(strcat('more than one logfile with run index_',run,' for subject_',sourcesubjs{sub}, ' please check'))
                
                else
                    log = readtable(logfilepath,'FileType','text','Delimiter','tab');
                    log = log(~isnan(log.Trial),:);
                    log.EventType = categorical(log.EventType);
                    log.Code = categorical(log.Code);
                    time_zero = log.Time(log.Trial == 0 & log.EventType == 'Pulse'); % time for onsets and durations is counted from the first scanner pulse onwards
                        
                        if size(time_zero,1) > 1
                            error('ambiguity about time zero in %s%s, please check logfile',subjs{sub},logfilenames{run});
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
                    filename = strcat(subjBIDSdir,'/',sourcesubjs{sub},'_task-',taskname,runnames{run},'_model_1_events.tsv');
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
                                    error('imagine row does not follow a picture row in %s%s, please check logfile',subjs{sub},logfilenames{run});
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
                    filename2 = strcat(subjrawdir,'\',subjs{sub},'_ses-1_task-',char(tasknames(run)),'_model_2_events.tsv');
                    writetable(log2,filename2,'Filetype','text','Delimiter','\t');
                    clear log2 filename2
                end % if loop checking whether logfile exists

            end % for loop runs

        end % for loop subjects
    
    end % if loop checking sourcesubjs == BIDSsubjs
    
end % if loop checking writing option