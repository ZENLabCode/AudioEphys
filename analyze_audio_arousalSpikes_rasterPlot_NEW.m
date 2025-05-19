% analyze_audio_arousalSpikes
%------------------------------------------------------------------------
% NEW 2022
% - Read several spike files, average all data
% - Add 2nd Figure with Rasterplot
% !!! if it's requestng the rhd fie just go to one folder with recording
%tmp = dir(['Z:\1 TIDIS Lab\Ida\Recording\Tetrode Recordings\',...
%     'FC\**\',...
%     '*HP_SUA.txt']);
% Averaging data & spike ratio
% 1: Averages data across channels of a spike file
%    (spikes are equal for all tetrodes, no averaging needed)
% 2: Averages data/ratio of selected time region
%    (across beeps per transition type & sleep stage)
% 3: Averages data/ratio across all files
% '*spike_SUA.txt']);
%     '*Au1_SUA.txt']);
%     '*0.txt']);

clc; clear; close all;
[scriptPath,scriptName] = fileparts(which(mfilename));
addpath(fullfile(scriptPath,'functions'))
addpath(genpath('Z:\OptoLab_v4.1\function')); %f_spikeRate_ISI.m


%PARAMETERS
%----------
%FILES
%spike data
tmp = dir(['Z:\1 TIDIS Lab\Ida\Recording\',...
    'Tetrode Recordings\**\',...
    '*CMT.txt']);
Files.spikes = selectionList(fullfile({tmp.folder},{tmp.name}));
%protocols (auto find appropriate one, error if no one fits)
% PS: script selects appropriate protocol from given ones by comparing
%     recorded stimTimes from data channels with protocol beep times
Files.protocol = fullfile(scriptPath,'protocols',{... possible protocols
    ...'protocol_morning.mat','protocol_afternoon.mat',...
    ...'protocol_laser_3600s_80dB_108x3beeps.mat',...
    ...'FC_protocol_1h_342beeps.mat',...
    ...'protocol_morning.mat',...
    'protocol_morning.mat'});


%ANALYSIS OPTIONS
%stimTimes label (stimTimes' or 'stimTimes2' or ...)
ana.stimTimes = 'stimTimes'; %to find appropriate protocol
%beep selection, power and frequency (pooled together, no single analysis)
ana.powers      = 80; %if empty, all available ones
ana.frequencies = 5000; %if empty, all available ones
%stages to analyze
% - transitions from selected stage to Wake (Sleep-Wake)
%   or staying in sleep (Sleep-Sleep, e.g. NREM-NREM or NREM-REM)
% - transitions to wake if switch to wake before next beep, or if transtion
%   to wake within delay limit (delayLIM)
ana.Stages = {... {number, label}, all stages! may have multiple numbers
    1,'Wake';...
    2,'NREM';...
    3,'REM';...
    };
ana.delayLIM = 10; %[s]
ana.Trans = {... {stage from, {stages to}, {not stages}, label, plot color}
    ... PS: transitions to must have all possible stages!!!
    ...     so transition to wake within 10s have {'NREM','REM','Wake'}
    'NREM' ,{'NREM','REM'} ,{'Wake'} ,'Sleep-Sleep' ,'b';...
    'NREM' ,{'Wake'}       ,{}       ,'Sleep-Wake'  ,'r';...
    'REM'  ,{'NREM','REM'} ,{'Wake'} ,'Sleep-Sleep' ,'b';...
    'REM'  ,{'Wake'}       ,{}       ,'Sleep-Wake'  ,'r';...
    };
%margin to plot
ana.margin = [0.7,0.6]; %[s], before beep start & after beep end
%sleep-wake transition: wake before next beep AND wake delay <= delayLIM
ana.delayLIM = 10; %[s]
ana.labels.spikeData={'Channel','Unit','Timestamp','PC 1','PC 2','PC 3'};
%channel data filter (optional, uses passband_fourier.m)
%  - lower is 0 if NaN
%  - upper is fs/2 if NaN
ana.bandpass = [... [lower,upper; lower,upper; ...]
    0,45;... remove 50 Hz,
    55,NaN;...
    ];
%select spikes units
% - if empty, takes all available units
% - 1:10 takes units 1 to 10 (if available)
%   PS: - unit zero are normally unsorted spikes
%         For having all units but 0, just select number high enough
%         Normally spike sorting have less than 10. However, you should
%         know that from spike sorting.
ana.units = 1:100;

%RATIO OPTIONS
% Using function f_spikeRate_ISI.m
%  --> uses movmean across winLen, results per [s]
%  --> output in time steps of winLen/2
ana.ratio.winLen = 0.1; %[s]

%TABLE
% Counts number of spikes
%  - counts for each file, averaged across specific transition
%    (so might be floating numbers)
% Tcnt is matrix with time start and time end in [s]
%  - time start/end is relative to beep start (at t==0)
%    E.g: Tcnt = [-0.2,0] counts spikes 0.2 s before beep start
%    PS : - uses  spikeTimes >= Tcnt(1) & spikeTimes < Tcnt(2)
%           (for to have spikes real before or from inclusive beep start)
%         - must be within margin (beep duration is)
%  - counts seoarately for each row of Tcnt --> separate xls sheet
% Opens file browser to ask where to save
% Does nothing if Tcnt is empty
Tcnt = [... [s], 
    -0.2 , 0  ;...
     0   , 0.2;...
     ];

%PLOT/TABLE OR NOT (true or false)
plo.spikes = true; %spike data (all files,channels & selected units)
plo.data   = true; %channel data & spike rate
plo.raster = true; %raster plot / bar

%PLOT PROPERTIES (must not be empty!)
%figure / axes
props.figure = {'position',[440,330,880,450]};
props.axes   = {'visible','on'}; %y-axes are beeps
%fig_createAxes
props.fig_createAxes = {[65,65,12],[50,10,80],'pixel'}; %basic
props.axisRatioY     = [5,1]; %signal, spindle ratio
%plot
props.plot = {'linewidth',1};
props.fill = {'facecolor',[1,1,1]*.6,'edgecolor','none'};
props.bar  = {'BarWidth',1,'facealpha',.5}; %PS: replaced eith plot
props.timeFac = 1;   %time factor to [s]
props.timeUni = 's'; %time unit based on factor
%text
props.title  = {'fontweight','normal'};
props.xlabel = {'visible','on'};
props.ylabel = {'visible','on'};
props.text   = {'visible','on'};
props.legend = {'visible','on'};
props.superTitle = {'visible','on'};

%MAIN PROGRAM
%--------------
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)));

%INIT
if isempty(Files.spikes)
    fprintf(2,'No File Selected!\n')
    return
end
%check filenames
tmp = fieldnames(Files);
for k = 1:numel(tmp)
    field = tmp{k};
    %must be cell
    if ~iscell(Files.(field))
        Files.(field) = {Files.(field)};
    end
    %exclude inexistent files
    ind = cellfun(@(x)exist(x,'file')~=2,Files.(field));
    if any(ind)
        if all(ind)
            fprintf(2,'Files.%s, no file exist!\n',field)
            return
        else
            fprintf('[\bFiles.%s, inexistent files removed:]\b\n',field)
            fprintf(' - %s\n',Files.(field){ind})
            Files.(field)(ind) = [];
        end
    end
end
%margin (relative to time)
ana.margin(1) = -abs(ana.margin(1));
ana.margin(2) = +abs(ana.margin(2));

%FIND HYPNOGRAM & CHANNEL FILES
noFIL = numel(Files.spikes);
Files.hypnogram = cell(noFIL,1);
Files.channels  = cell(noFIL,1);
for fil = 1:noFIL
    [rPath,rFile] = fileparts(Files.spikes{fil});
    %hypnogram file
    file = fullfile(rPath,'Hypnogram.mat');
    if exist(file,'file')==2
        Files.hypnogram{fil} = file;
    end
    %channel files
    tmp = regexp(rFile,'_','split');
    ind = contains(tmp,'amp-');
    if sum(ind)==0
        Files.channels{fil} = 'Cannot extract channel names from filename';
        continue
    end
    tmp = tmp{ind};
    ind = regexp(tmp,'[0-9]*');
    num = cellfun(@str2double,regexp(tmp,'[0-9]*','match'),...
        'uniformoutput',false);
    channels = cellfun(@(x)sprintf('%s%03i.mat',tmp(1:ind(1)-1),x),num,...
        'uniformoutput',false);
    files = cellfun(@(x)fullfile(rPath,x),channels,...
        'uniformoutput',false);
    ind = cellfun(@(x)exist(x,'file')~=2,files);
    if any(ind)
        Files.channels{fil} = sprintf('Channels not found:%s',...
            strjoin(channels(ind),', '));
        continue
    end
    Files.channels{fil} = files;
end
%remove 'no hypnogram'
ind = cellfun(@isempty,Files.hypnogram);
if any(ind)
    fprintf('[\bFiles removed (missing Hypnogram.mat):]\b\n')
    fprintf(' - %s\n',Files.spikes{ind})
    Files.spikes(ind)    = [];
    Files.hypnogram(ind) = [];
    Files.channels(ind)  = [];
end
%remove 'channel issue'
ind = cellfun(@ischar,Files.channels);
if any(ind)
    fprintf('[\bFiles removed (channel issue):]\b\n')
    tmp = unique(Files.channels(ind));
    for k = 1:numel(tmp)
        str = tmp{k};
        ind =  cellfun(@(x)ischar(x)&&strcmpi(x,str),Files.channels);
        fprintf('  %s\n',str)
        fprintf('   - %s\n',Files.spikes{ind})
        Files.spikes(ind)    = [];
        Files.hypnogram(ind) = [];
        Files.channels(ind)  = [];
    end
end

%NUMBER OF ...
noFIL = numel(Files.spikes);
noPTC = numel(Files.protocol);
noSTA = size(ana.Stages,1);
noTRA = size(ana.Trans,1);

%LOAD PROTOCOLS
ind = cellfun(@(x)exist(x,'file')~=2,Files.protocol);
if any(ind)
    fprintf('[\Protocol-files removed (do not exist)]\b\n')
    fprintf(' - %s\n',Files.protocol{ind});
    Files.protocol(ind) = [];
end
PP = f_protocols_load(Files.protocol);

%CHECK PROTOCOLS
%must have equal signals (different timing allowed)
for ptc = 2:noPTC
    if ~isequaln(PP(1).beep,PP(ptc).beep) || ...
            ~isequaln(PP(1).laser,PP(ptc).laser)
        error(sprintf(['All protocols must have same beeps/laser\n',...
            '(for averaging, different beep timing is OK)']))
    end
    if numel(PP(1).stimSTA)~=numel(PP(ptc).stimSTA)
        error('All protocols must have the same amount of beeps')
    end
end
%check variables
tmp = unique([PP(1).beep.Signals{:,2}]); %powers
if isempty(ana.powers)
    ana.powers = tmp;
else
    ind = ~ismember(ana.powers,tmp);
    if any(ind)
        str = sprintf('%g, ',ana.powers(ind));
        fprintf('[\bPower %s dB not available (removed)]\b\n',...
            str(1:end-2))
        ana.powers(ind) = [];
    end
end
tmp = unique([PP(1).beep.Signals{:,1}]); %frequencies
if isempty(ana.frequencies)
    ana.frequencies = tmp;
else
    ind = ~ismember(ana.frequencies,tmp);
    if any(ind)
        str = sprintf('%g, ',ana.frequencies(ind));
        fprintf('[\bFrequency %s Hz not available (removed)]\b\n',...
            str(1:end-2))
        ana.frequencies(ind) = [];
    end
end

%% READ/PREPARE DATA (and plot spikes)
fig = 0; %init
fprintf('\nREAD DATA\n'); tic
nnFIL  = numel(num2str(noFIL));
indent = blanks(2*nnFIL+2); %indent for fprintf
clear DATA
for fil = 1:noFIL
    fileSPK  = Files.spikes{fil};
    fileHYP  = Files.hypnogram{fil};
    filesCHA = Files.channels{fil};
    rPath = fileparts(fileSPK);
    fprintf('%*i/%i: %s\n',nnFIL,fil,noFIL,rPath)
    %channels
    noCHA = numel(filesCHA);
    channels = cell(1,noCHA);
    for k = 1:numel(channels)
        [~,channels{k}] = fileparts(filesCHA{k});
    end
    
    %% LOAD SPIKE DATA
    file = fileSPK;
    [~,rFile,rExt] = fileparts(file);
    fprintf('%s Load Spikes:\n',indent);
    fprintf('%s   - %s\n',indent,[rFile,rExt]);
    try
        spikes = f_readSpikes(file,ana.labels.spikeData);
    catch
        fprintf(2,'%s Not a regular spike-file\n',indent)
        continue
    end
    %remove unwanted units
    if ~isempty(ana.units)
        ind = ismember(spikes.unit,ana.units);
        tmp = fieldnames(spikes);
        for k = 1:numel(tmp)
            spikes.(tmp{k}) = spikes.(tmp{k})(ind,:);
        end
    end
    
    %% PLOT
    if plo.spikes
        tmp = fullfile(rPath,ls(fullfile(rPath,'*.rhd')));
        tmp = read_Intan_RHD2000_file_zenlab(tmp,false); %for sampling rate
        x = (1:size(spikes.spikes,2))/...
            tmp.frequency_parameters.amplifier_sample_rate * 10^3;
        units = unique(spikes.unit);
        for k = 1:numel(units)
            unit = units(k);
            fig = fig+1;
            hf(fig) = figure; hold on
            ind = spikes.unit==unit;
            m = mean(spikes.spikes(ind,:),1);   %mean
            s = std(spikes.spikes(ind,:),[],1); %STD
            plot(x,m,'k','linewidth',2);
            plot(x,m-s,':k','linewidth',2)
            plot(x,m+s,':k','linewidth',2)
            %settings & text
            yl = [min(m-s) max(m+s)];
            yl = yl+[-.1 .1]*diff(yl);
            yt = linspace(0,x(end),5);
            set(gca,'xlim',[0,x(end)],'ylim',yl,'box','on',...
                'xtick',yt,'xticklabel',[],'xgrid','on',...
                'ytick',[]);
            ylabel('Mean \pm STD')
            title({strrep(rFile,'_','\_'),sprintf('Spikes Unit %i',unit)})
            %line
            xx = [-1 0]+x(end);
            yy = [1,1]*min(yl)+0.1*diff(yl);
            line(xx,yy,'color','k','linewidth',2)
            text(sum(xx)/2,yy(1),sprintf('%g ms',diff(xx)),...
                'horizontalalignment','center','verticalalignment','top')
            %channel text
            for q = 1:4
                text(sum(yt(q:q+1))/2,yl(1),sprintf('Channel %i',q),...
                    'horizontalalignment','center',...
                    'verticalalignment','top')
            end
        end
    end
    
    %% LOAD CHANNEL DATA
    fprintf('%s Load Channel Data:\n',indent);
    for cha = 1:noCHA
        file = filesCHA{cha};
        [~,rFile,rExt] = fileparts(file);
        fprintf('%s   - %s\n',indent,[rFile,rExt]);
        tmp = load(file);
        %init/check
        if cha==1
            fs        = tmp.SampRate;
            stimTimes = tmp.(ana.stimTimes);
            noSAM     = numel(tmp.resampled_data_mV);
            Data      = NaN(noSAM,noCHA);
        elseif ~isequaln(fs,tmp.SampRate) || ...
                ~isequaln(stimTimes,tmp.(ana.stimTimes)) || ...
                noSAM~= numel(tmp.resampled_data_mV)
            error('uups') %must be all the same
        end
        %append
        Data(:,cha) = tmp.resampled_data_mV(:);
    end
    %filter
    if ~isempty(ana.bandpass)
        Data = passband_fourier(Data,ana.bandpass,fs);
        fprintf('%s   Data filtered\n',indent);
    end
    %mean
    if noCHA>1
        Data = mean(Data,2);
        fprintf('%s   Data averaged across channels\n',indent);
    end
    
    %% LOAD HYPNOGRAM
    file = fileHYP;
    [~,rFile,rExt] = fileparts(file);
    fprintf('%s Load Hypnogram:\n',indent);
    fprintf('%s   - %s\n',indent,[rFile,rExt]);
    tmp = load(file);
    hypnogram = tmp.Hypnogram(:);
    if isfield(tmp,'fs')
        fac = fs/tmp.fs;
    elseif isfield(tmp,'SampRate')
        fac = fs/tmp.SampRate;
    else
        fac = floor(noSAM/numel(hypnogram));
    end
    if fac<1
        error('uups')
    elseif fac>1
        hypnogram = repmat(hypnogram',fac,1);
        hypnogram = hypnogram(:);
        fprintf('%s   Upsampled by factor %i\n',fac)
    end
    hypnogram(numel(hypnogram)+1:noSAM) = NaN;
    
    %% GET APPROPRIATE PROTOCOL
    warning off
    P = f_protocol_check(PP,stimTimes,fs);
    warning on
    if isempty(P)
        fprintf('%s [\bProtocol does not fit to data]\b\n',indent)
        %continue
    end
    file = P.protocol.file;
    [~,rFile,rExt] = fileparts(P.protocol.file);
    fprintf('%s Protocol:\n',indent)
    fprintf('%s   - %s\n',indent,[rFile,rExt]);
    fprintf('%s     repetitions, N = %i\n',indent,...
        numel(unique(P.repetition)))
    %selected frequencies & power
    ind = isnan(P.repetition) | ... remove NaNs
        ~ismember(P.power,ana.powers) | ...
        ~ismember(P.frequency,ana.frequencies);
    tmp = fieldnames(P);
    for k = 1:numel(tmp)
        if size(P.(tmp{k}),1) == numel(ind)
            P.(tmp{k})(ind,:) = [];
        end
    end
    %print out
    fprintf('%s     beeps excluded, N = %i\n',indent,sum(ind));
    fprintf('%s     beeps analyzed, N = %i\n',indent,sum(~ind))
    
    %% INIT (new, to average later the needed)
    if fil==1
        beepDUR = P.protocol.beep.duration;
        margin  = round((ana.margin+[0,beepDUR])*fs); %index
        noSTI   = numel(P.stimSTA);
        ind0    = margin(1):margin(2); %to add at stim starts
    elseif beepDUR~=P.protocol.beep.duration || noSTI~=numel(P.stimSTA)
        error('uups') %should be checked before file loop
    end
    
    %% INIT
    if ~exist('DATA','VAR')
        N = noFIL*noSTI;
        DATA.Data  = NaN(N,numel(ind0));
        DATA.File  = NaN(N,1);
        DATA.Trans = NaN(N,1);
        DATA.SpikeTime = cell(N,1);
        %spike vector
        rProps = {ana.ratio.winLen*1000,round((diff(margin)+1)/fs*1000)}; %[ms]
        [sRate,tRate,sVec] = f_spikeRate_ISI([],rProps{:}); %init, no spikes
        DATA.SpikeInd = NaN(N,numel(sVec));
        DATA.Ratio  = NaN(N,numel(tRate));
        %time vectors
        DATA.tData  = ind0/fs; %[s]
        DATA.tRate  = tRate/1000 + margin(1)/fs; %[s]
        DATA.tSpike = linspace(margin(1),margin(2),numel(sVec));
    end
    
    %BEEP LOOP
    for sti = 1:noSTI
        stimSTA = P.stimSTA(sti);
        stimEND = stimSTA+ana.delayLIM*fs-1; %to get transition state
        if sti==1
            stimEND = min([stimEND,P.stimSTA(sti+1)-1]);
        else
            stimEND = min([stimEND,find(~isnan(hypnogram),1,'last')]);
        end
        
        %get row
        row = find(isnan(DATA.File),1,'first');
        
        %CHANNEL DATA
        indDAT = ind0 + stimSTA;
        DATA.File(row,1) = fil;
        DATA.Data(row,:) = Data(indDAT); %care of index!!!
        
        %SPIKE TIME/RATE
        ind = ... spike index (timestap in [s])
            spikes.timestamp >= indDAT(1)/fs & ...
            spikes.timestamp <= indDAT(end)/fs;
        %time (zero at beep start [s])
        DATA.SpikeTime{row,1} = spikes.timestamp(ind)-stimSTA/fs; 
        %rate (zero at margin(1), [ms])
        tmp = ceil((spikes.timestamp(ind)-indDAT(1)/fs)*1000);
        tmp(tmp==0) = []; %because rounding
        [sRate,tRate,sVec] = f_spikeRate_ISI(tmp,rProps{:});
        if ~isequal(DATA.tRate,tRate/1000+margin(1)/fs) %check
            error('uups')
        end
        DATA.SpikeInd(row,:)  = sVec;
        DATA.Ratio(row,:)     = sRate;
        
        %TRANSITION
        stages = hypnogram(stimSTA:stimEND);
        ind1   = ismember([ana.Stages{:,1}],stages(1));
        ind2   = ismember([ana.Stages{:,1}],stages(2:end));
        stage1  = ana.Stages{ind1,2}; %stage beep starts
        stages2 = ana.Stages(ind2,2); %from beep start till ctrl. dur
        ind = find(ismember(ana.Trans(:,1),stage1) & ...
            cellfun(@(x) any(ismember(stages2,x)),ana.Trans(:,2)) & ...
            cellfun(@(x)~any(ismember(stages2,x)),ana.Trans(:,3)));
        if ~isempty(ind)
            if numel(ind)~=1
                error('uups')
            end
            DATA.Trans(row,1) = ind;
        end
    end %beep loop
end %file loop
if ~exist('DATA','var')
    fprintf(2,'NO DATA LOADED (check protocols)\n')
    toc
    return
end
%remove inexistent data
if any(ind)
    fields = fieldnames(DATA);
    for k = 1:numel(fields)
        field = fields{k};
        if ~strcmp(field(1),'t') %time vector
            try
                DATA.(field)(ind,:) = [];
            catch
            end
        end
    end
end

%% MEAN DATA ACROSS FILES & TRANSITIONS (CONCAT SPIKES)
DATA.mData  = NaN(noTRA,size(DATA.Data,2));
DATA.sData  = NaN(noTRA,size(DATA.Data,2));
DATA.mRatio = NaN(noTRA,size(DATA.Ratio,2));
DATA.cSpikeTime = cell(noTRA,1);
DATA.cSpikeInd  = cell(noTRA,1);
%DATA.cntStim    = NaN(noTRA,1);
clear dat
for tra = 1:noTRA
    indTRA = DATA.Trans==tra;
    %1st mean across files
    dat.mData  = NaN(noFIL,size(DATA.Data,2));
    dat.mRatio = NaN(noFIL,size(DATA.Ratio,2));
    %dat.cSpikeInd  = repmat(1:size());
    dat.cNoSPK = 0;
    for fil = 1:noFIL
        ind = indTRA & DATA.File==fil;
        if any(ind)
            dat.mData  = nanmean(DATA.Data(ind,:),1);
            dat.mRatio = nanmean(DATA.Ratio(ind,:),1);
        end
    end
    %add mean/concat
    DATA.mData(tra,:)  = nanmean(dat.mData,1);
    DATA.sData(tra,:)  = nanstd(dat.mData,[],1);
    DATA.mRatio(tra,:) = nanmean(dat.mRatio,1);
    %DATA.cntStim(tra) = max(sum(~isnan(DATA.Data(indTRA,:)),1));
    %concat
    tmp = cellfun(@(x)x(:),DATA.SpikeTime(indTRA),...
        'UniformOutput',false);
    DATA.cSpikeTime{tra} = vertcat(tmp);
    DATA.cSpikeInd{tra} = DATA.SpikeInd(indTRA,:);
end


%% FIGURE CHANNEL DATA & SPIKE RATES
% clc; close all; fig = 0; clear hf
if plo.data
    fig = fig+1;
    hf(fig) = figure(props.figure{:}); drawnow
    transStages = unique(ana.Trans(:,1));
    noLAB = numel(transStages);
    movegui(hf(fig),'center')
    HA = NaN(2,noLAB);
    HA(1,:) = fig_createAxes(hf(fig),[1,noLAB],props.fig_createAxes{:});
    for k = 1:noLAB
        set(HA(1,k),'tag',transStages{k})
        pos0 = get(HA(1,k),'position'); uni = get(HA(1,k),'unit');
        dy = props.fig_createAxes{2}(2);
        y = (pos0(4)-dy)/sum(props.axisRatioY).*props.axisRatioY;
        pos1 = [pos0(1) sum(pos0([2,4]))-y(1) pos0(3) y(1)];
        pos2 = [pos0(1:3) y(2)];
        set(HA(1,k),'position',pos1)
        HA(2,k) = axes('unit',uni,'position',pos2);
    end
    set(HA,'unit','normalized')
    
    %INIT
    tmp = [min(DATA.mData(:)-DATA.sData(:)),...
        max(DATA.mData(:)+DATA.sData(:))];    
    yl     = tmp+diff(tmp).*[-1,1]*0.11;
    LegSTR = cell(1,noLAB);
    HP     = cell(1,noLAB);
    
    %PLOT DATA
    %beep
    Z = zeros(1,4);
    X = [0,0,beepDUR,beepDUR];
    Y = [min(yl),max(yl),max(yl),min(yl)];
    for k = 1:noLAB
        set(hf(fig),'CurrentAxes',HA(1,k))
        fill(X,Y,Z,props.fill{:}); hold on
    end
    %data
    tDat = DATA.tData*props.timeFac;
    tRat = DATA.tRate*props.timeFac;
    for tra = 1:noTRA
        [stage,~,~,lab,col] = ana.Trans{tra,:};
        indAXI = find(strcmpi(transStages,stage));
        set(hf(fig),'CurrentAxes',HA(1,indAXI))
        %data
        m = DATA.mData(tra,:);
        s = DATA.sData(tra,:);
        r = DATA.mRatio(tra,:);
        legStr = sprintf('%s, N = %i',lab,numel(DATA.cSpikeTime{tra}));
        
        %PLOT DATA
        if any(s~=0) %when noFIL>1
            ind = numel(tDat):-1:1;
            dat = [m-s,m(ind)+s(ind)];
            h = fill([tDat,tDat(ind)],dat,zeros(size(dat)),...
                'edgecolor','none','facecolor',col,'facealpha',0.5);
            uistack(h,'bottom')
        end
        hp = plot(tDat,m,'color',col,props.plot{:});
        HP(indAXI) = {[HP{indAXI},hp]};
        if isempty(LegSTR{indAXI})
            LegSTR{indAXI} = legStr;
        else
            LegSTR{indAXI} = [LegSTR(indAXI),legStr];
        end
        
        %PLOT RATIO
        set(hf(fig),'CurrentAxes',HA(2,indAXI))
        plot(tRat,r,'color',col,props.plot{:}); hold on
        ylabel({'Spike','Ratio'},props.ylabel{:});
        xlabel(sprintf('Time [%s]',props.timeUni),props.xlabel{:})
    end
    
    %TEXT
    %legend/title
    for k = 1:noLAB
        set(hf(fig),'CurrentAxes',HA(1,k))
        legend(HP{k},LegSTR{k})
        title(transStages{k},props.title{:})
    end
    %super title
    str=cell(1,3);
    if noFIL==1
        [rPath,rFile,rExt]=fileparts(Files.spikes{1});
        str{1}=regexprep(rPath,{'\','_'},{'\\\\','\\_'});
        str{2}=regexprep([rFile,rExt],{'\','_'},{'\\\\','\\_'});
    else
        str{1}=sprintf('Files N = %i',noFIL);
    end
    str{3}=sprintf('%sdB,  %sHz,  max wake delay %g s',...
        sprintf('%g ',ana.powers),sprintf('%g ',ana.frequencies),...
        ana.delayLIM);
    h = fig_superLabel(HA(1,:),'title',str);
    set(h,props.superTitle{:})
    
    %SETTINGS
    linkaxes(HA,'x')
    set(HA,'xlim',margin/fs*props.timeFac,'unit','normalized',props.axes{:})
    set(HA(1,:),'ylim',yl,'xticklabel',[])
    set(HA(2,:),'ylim',max(DATA.mRatio(:))*[0,1.1])
end

%% FIGURE
% clc; close all; fig = 0; clear hf
% clc; fig = 0; clear hf
%plot options
%  - 1 bar-plot new calculated mean firing rate
%  - 2 bar-plot ration calculated above (quite simular)
%  - 3 is 1 & line plot ratio
%  - 4 is 1 & bar-plot ratio
plotOpt = 2;
maxR = 0; %max riring rate
if plo.raster
    transStages = unique(ana.Trans(:,1));
    for k = 1:numel(transStages)
        stage = transStages{k};
        indTRA = find(ismember(ana.Trans(:,1),stage));
        
        %figure & axes
        fig = fig+1;
        hf(fig) = figure(props.figure{:}); drawnow
        movegui(hf(fig),'center')
        noCOL = numel(indTRA);
        HA = fig_createAxes(hf(fig),[2,noCOL],props.fig_createAxes{:});
        set(HA,'unit','normalized')
        %init
        xl = [ana.margin(1),ana.margin(2)+beepDUR]*props.timeFac;
        tRat = DATA.tRate*props.timeFac;
        
        
        for axi = 1:noCOL
            trans = indTRA(axi);
            label = ana.Trans{trans,4};
            %spike/ratio data           
            spikes = cellfun(@(x)(x')*props.timeFac,...
                DATA.cSpikeTime{trans},'uniformoutput',false);
            N = numel(spikes); %beeps this stage & type
            rates = DATA.mRatio(trans,:);
            SpikeInd = DATA.cSpikeInd{trans};
            
            %PLOT RASTER PLOT
            set(hf(fig),'CurrentAxes',HA(1,axi))
            ind = ~cellfun(@isempty,spikes);
            if any(ind)
                plotSpikeRaster(spikes,'PlotType','vertline',...
                    'RelSpikeStartTime',0.0);
            else
            end
            %text
            set(gca,'visible','off')
            title({label,sprintf('N = %i',N)},'visible','on')
            
            %FIRING RATE
            %count spikes
            dt = mean(diff(tRat));
            tt = linspace(xl(1),xl(2),numel(tRat)+1);
            n   = numel(tt)-1;
            CNT = NaN(N,n);
            for spk = 1:N
                s = spikes{spk};
                for k = 1:n
                    CNT(spk,k) = sum(s>tt(k) & s<=tt(k+1));
                end
            end
            %mean across transitions
            FR = mean(CNT,1)/dt;
            tFR  = tt(1:end-1)+diff(tt)/2;
            if ismember(plotOpt,[1,3,4])
                maxR = max([maxR;FR(:);rates(:)]);
            else
                maxR = max([maxR;rates(:)]);
            end
            
            %PLOT RATIO
            set(hf(fig),'CurrentAxes',HA(2,axi))
            switch plotOpt
                case 1
                    bar(tFR,FR,'BarWidth',1);
                    ylabel('Mean Firing Rate [Hz]')
                case 2
                    bar(tRat,rates,'BarWidth',1)
                    ylabel('Spike Ratio')
                case 3
                    bar(tFR,FR,'BarWidth',1); hold on
                    plot(tRat,rates,'r')
                case 4
                    bar(tFR,FR,'BarWidth',1); hold on
                    bar(tRat,rates,'BarWidth',1,...
                        'FaceAlpha',0.5);
            end
            xlabel(sprintf('Time [%s]',props.timeUni))
        end
        %super title
        if noFIL==1
            [~,rFile] = fileparts(Files.spikes{1});
            sgtitle({strrep(rFile,'_','\_'),stage},'fontweight','bold');
        else
            sgtitle({sprintf('Average of %i Files',noFIL),stage},...
                'fontweight','bold');
        end
        %settings
        linkaxes(HA,'x')
        set(HA,'xlim',xl)
        if maxR==0
            maxR = 1;
        end
        set(HA(2,:),'ylim',[0,maxR*1.2],'box','off')
    end
end

%% TABLE
clc %%%

if ~isempty(Tcnt)
    %init
    labels = ['File',cell(1,noTRA),'Time Range [s]','No Transitions'];
    offC    = 1; %offset column to transition data
    Table  = cell(0,numel(labels));
    
    %count loop
    for lim = 1:size(Tcnt,1)
        tLim = Tcnt(lim,:);
        Tab  = cell(noFIL+1,numel(labels)); %+ one empty row
        
        %file loop
        for fil = 1:noFIL
            str = cell(1,noTRA); %transition count (per transition)
            
            %transition loop
            for tra = 1:noTRA
                col = tra+offC; %column
                %transition label
                [stage,~,~,lab,~] = ana.Trans{tra,:};
                if isempty(labels{col})
                    labels{col} = sprintf('%s %s',stage,lab);
                elseif ~strcmpi(labels{col},sprintf('%s %s',stage,lab))
                    error('uups')
                end
                
                %count spikes
                ind = DATA.Trans==tra & DATA.File==fil;
                str{tra} = num2str(sum(ind));
                tmp = cellfun(@(x)x(x>=tLim(1)&x<tLim(2)),...
                    DATA.SpikeTime(ind),'uniformoutput',false);
                %append
                Tab{fil,col} = mean(cellfun(@numel,tmp));
            end %transition loop
            
            %append
            [~,Tab{fil,strcmpi(labels,'File')}] = ...
                fileparts(Files.spikes{fil});
            Tab{fil,strcmpi(labels,'no Transitions')} = ...
                strjoin(str,', ');
            Tab{fil,strcmpi(labels,'Time Range [s]')} = ...
                sprintf('%g to %g',tLim);
        end %file loop
        %append
        Table = [Table;labels;Tab];
    end %count loop
    
    %SAVE
    %default savepath (common path)
    rPaths = cellfun(@fileparts,Files.spikes,'uniformoutput',false);
    sPath = rPaths{1};
    for fil = 2:numel(rPaths)
        tmp2 = regexp(rPaths{fil},filesep,'split');
        tmp1 = regexp(sPath,filesep,'split');
        tmp1(numel(tmp2)+1:end) = [];
        ind = 1;
        while ind<numel(tmp1) && strcmpi(tmp1{ind},tmp2{ind})
            ind = ind+1;
        end    
        sPath = strjoin(tmp1(1:ind-1),filesep);
    end    
    %savename
    [sFile,sPath,indx] = uiputfile('*.xls;*.xlsx','Save xls-file as',...
        fullfile(sPath,'spikeCount.xlsx'));
    if ischar(sFile)
        sname = fullfile(sPath,sFile);
        %proper saving xls-file
        if exist(sname,'file')==2
            delete(sname) 
        end
        %save
        xlswrite(sname,Table)
        xls_cellFit(sname)
        fprintf('SAVED: %s\n',sname)
    end
end