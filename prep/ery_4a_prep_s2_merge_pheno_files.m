%% ery_4a_prep_s2_merge_pheno_files
%
% This quick and dirty script merges the phenotype.tsv files with online
% and offline ratings into one single ratings_all.tsv file, and saves it
%
% NOTE: script needs to be run from BIDS/phenotype dir, NOT from root dir
%
%
%__________________________________________________________________________
%
% author: Lukas Van Oudenhove
% date:   April, 2022
%
%__________________________________________________________________________
% @(#)% ery_4a_prep_s2_merge_pheno_files         v1.0        
% last modified: 2022/04/25


%% READ OFFLINE FILE

varNames = {'trial_type','participant_id','hunger','intensity','liking'};
varTypes = {'categorical','categorical','double','double','double'};
delimiter = '\t';
filetype = 'text';
datastartline = 2;

opts = delimitedTextImportOptions('VariableNames',varNames,...
    'VariableTypes',varTypes,...
    'Delimiter',delimiter,...
    'DataLines',datastartline);

offline = readtable('ratings_offline.tsv',opts);


%% READ ONLINE FILE

online = readtable('ratings_online.tsv','Delimiter',delimiter,'FileType',filetype);
online.participant_id = categorical(online.participant_id);
online.trial_type = categorical(online.trial_type);


%% MERGE AND WRITE

all = join(online,offline);

writetable(all,'ratings_all.tsv','FileType',filetype,'Delimiter',delimiter);