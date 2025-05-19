% code for z-scoring LFP and EEG Data. Additionally you can also
% re-reference your data to your desires
clc; clear; close all;
[scriptPath,scriptName] = fileparts(which(mfilename));
addpath(fullfile(scriptPath,'functions'))

ftpath = 'Z:\Ivan\scriptsForCastor\fieldtrip-master';
addpath(ftpath);
% 2 run fieldtrip startup (no output = good!)
ft_defaults;

% re-feference?
reference=false;
ref= {'all'}; % here you can put your desired channel for referencing in form of a string inside a cell. 'all' means all available
               % channels will be used to create a common average reference
%PARAMETERS
%----------
%FILES, data and protocol (checks which protocol fits to which data)
Files.data = import_fileList(... opens file browser if not exist
    ... 'Z:\1 TIDIS Lab\Ida\Script\filelists\Opto Stimulation\*.txt');
    'Z:\1 TIDIS Lab\Ida\Script\Audio Analysis_PulsePal_3\filelists\*.txt');
Files.prot = fullfile(... script automatically use appropriate protocol
    'Z:\1 TIDIS Lab\Ida\Script\Audio Analysis_PulsePal_3\protocols',...
    {...'protocol_laser_3600s_80dB_108x3beeps.mat',... we now have two
    ...'protocol_laser_3600s_80dB_108x3beeps_2.mat',...
    ...'protocol_laser_60s_80dB_2x3beeps.mat',...
    ...'FC_protocol_45s_14beeps.mat',...
    ...'FC_protocol_1h_342beeps.mat',... 
    ...'protocol_laser_600s_80dB_27x2beeps.mat'...
    ...'protocol_laser_120s_80dB_5x2beeps.mat'...
    ...'protocol_morning.mat',...
    ...'protocol_2min_80dB_15beeps.mat'...
    ...'protocol_laser_3600s_80dB_108x3beeps_2.mat'
    'protocol_laser_600s_80dB_27x2beeps.mat'});

%% MAIN PROGRAM
%---------------
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)));
%number of ...
noFIL = numel(Files.data);
%noSTA=size(opt.Stages,1);  %stages

%COMMON PATH (to shorten output)
strs1 = regexp(fileparts(Files.data{1}),filesep,'split');
n = numel(strs1);
str = [newline,char(0)]; %string that sholud not appear in any path name
for fil = 2:noFIL
    strs2 = [regexp(fileparts(Files.data{fil}),filesep,'split'),str];
    n   = min([n,find(~ismember(lower(strs2),lower(strs1)),1,'first')-1]);
end
mPath = fullfile(strs1{1:n}); clear strs1 strs2
fprintf('COMMON PATH: %s\n',mPath);

%CHECK DATA FILES
%files that not exist
ind = cellfun(@(x)exist(x,'file')~=2,Files.data);
if any(ind)
    fprintf('[\bFiles NOT found (removed)]\b\n')
    fprintf('  - %s\n',Files.data{ind})
    Files.data(ind) = [];
end
%files duplicate
ind = true(size(Files.data));
[~,tmp] = unique(lower(Files.data));
ind(tmp) = false;
if any(ind)
    fprintf('[\bFiles duplicate (kept only once)]\b\n')
    fprintf('  - %s\n',Files.data{ind})
    Files.data(ind) = [];
end
%corresponding hypnogram files needed
ind = cellfun(@(x)exist(x,'file')~=2,...
    cellfun(@(x)fullfile(fileparts(x),'Hypnogram.mat'),...
    Files.data,'uniformoutput',false));
if any(ind)
    fprintf('[\bFiles removed (missing hypnogram)]\b\n')
    fprintf('  - %s\n',Files.data{ind})
    Files.data(ind) = [];
end
%files remaining
recs={};
for i= 1 : size(Files.data,1)
    recs{i}= Files.data{i,1}(1:end-14);
end
recs=unique(recs);
noFIL = numel(recs);
if noFIL==0
    fprintf(2,'NO FILE FOUND FOR ANALYSIS!\n')
    return
else
    fprintf('FILES N = %i\n',noFIL)
end
%files must have equal sampling rates
load(Files.data{1},'SampRate');
fs = SampRate;
for fil = 2:noFIL
    load(Files.data{fil},'SampRate');
    if fs~=SampRate
        error('Sampling rate must be the same in all files (%i vs %i)',...
            fs,SampRate)
    end
end

%% LOAD DATA
fprintf('\nLOAD DATA:\n'); 
DATA = cell(noFIL,1); %structure per data file
nnFIL  = numel(num2str(noFIL));
indent = blanks(2*nnFIL+1); %indent for fprintf

for fil=1:noFIL
    
    cd([recs{fil}])
    
    allchans=[];
    for q=1: size(Files.data,1)
        allchans(q)= strcmp(recs{fil}, Files.data{q}(1:end-14));
    end
    allchans= find(logical(allchans));
    
    chanNames={};
    for q=1:length(allchans)
        chanNames{q}= Files.data{allchans(q)}(end-12:end);
    end
    chanNames{end+1}= 'EEG1.mat';
    chanNames{end+1}= 'EEG2.mat';
    
    trial=[];
    dat=[];
    
    for i= 1:length(chanNames)
        
        load([chanNames{i}],'resampled_data_mV', 'SampRate', 'stimTimes')
        if size(resampled_data_mV, 1) == 1
            resampled_data_mV= resampled_data_mV';
        end
        trial = vertcat (trial, resampled_data_mV');
        dat.label{i,1} = chanNames{i};
        
    end
    
    dat.trial= {trial};
    dat.time  = {[1:numel(resampled_data_mV)]/SampRate};
    dat = ft_checkdata(dat);
    
    cfg = [];
    cfg.standardize = 'yes';
    dat_rawz = ft_preprocessing(cfg,dat);
    
    if reference==1
        cfg = [];
        cfg.reref= 'yes';
        cfg.refchannel = ref;
        dat_rawzr = ft_preprocessing(cfg,dat_rawz);
        dat_rawz= dat_rawzr;
        
        for i= 1:size (dat_rawzr.label,1)
            resampled_data_mV= dat_rawzr.trial{1,1}(i,:);
            save (['z-scored-referenced-' dat_rawzr.label{i}], ...
                'resampled_data_mV', 'SampRate', 'stimTimes')
        end
    else
        
        for i= 1:size (dat_rawz.label,1)
            resampled_data_mV= dat_rawz.trial{1,1}(i,:);
            save (['z-scored-' dat_rawz.label{i}], ...
                'resampled_data_mV', 'SampRate', 'stimTimes')
        end
        
    end
    
end