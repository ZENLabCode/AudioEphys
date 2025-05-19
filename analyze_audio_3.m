% analyze_audio_??
%------------------------------------------------------------------------
% Plots mean of recording signals at audio beeps (plus/minus margin) per
% beep frequency, beep power and scored sleep stage.
%
% There are options for averaging data:
%  1: mean of pooled data; mean across all signals at audio beeps of all
%     individuals
%     Note: individuals having more beeps in a special sleep stage (per
%           frequency and power) contributes more to the mean
%           (higher weighted)
%  2: mean of individual means
%     Note: all individuals are equally weighted, unless no beep occured
%           in a specific sleep stage, frequency, power
%
% PS: - use audioBeeps_createProtocol.m for creating beep protocols
%       (beep duration bust be equal) and audioBeeps_runProtocol.m for
%       recording with PulsePal.
%     - protocol files are needed, as stimulus of recording only containes
%       beep start and end but not frequency and power.
%     - Difference analyze_audio_2 to analyze_audio. Latter only averages
%       equal channels (equal channel names). In this can also average
%       data from different channelnames. However, one has to define the
%       channel names for each read path. Also save code is excluded
%       (needs to be done manually now).
%
%
% Thomas Rusterholz, Mai 2019
%-------------------------------------------------------------------------      

% Data Ida:
%   EEG1 - 28
%   EEG2 - 31
%   EMG  - 3, 6
%   CMT  - 20-23
%   AU1  - 14-15, 17-18
%   MG   - 8-10
% PS: Use CTX AU or THA MG channels

clc; clear; close all;
[scriptPath,scriptName] = fileparts(which(mfilename));
addpath(fullfile(scriptPath,'functions'))


%PARAMETERS
%----------
%FILES, data and protocol (checks which protocol fits to which data)
Files.data = import_fileList(... opens file browser if not exist
    ... 'Z:\1 TIDIS Lab\Ida\Script\filelists\Opto Stimulation\*.txt');
    'Z:\1 TIDIS Lab\Ida\Script\Audio Analysis_PulsePal_3\filelists\*.txt');
Files.prot = fullfile(... script automatically use appropriate protocol
    'Z:\1 TIDIS Lab\Ida\Script\Audio Analysis_PulsePal_3\protocols',{...
    ...'protocol_laser_3600s_80dB_108x3beeps.mat',... we now have two
    ...'protocol_laser_3600s_80dB_108x3beeps_2.mat',...
    ...'protocol_laser_600s_80dB_27x2beeps.mat',...
    ...'FC_protocol_45s_14beeps.mat',...
    ...'protocol_laser_120s_80dB_5x2beeps.mat',...
    ...'FC_protocol_1h_342beeps.mat',... 
    ...'protocol_morning.mat',...
    ...'protocol_laser_3600s_80dB_108x3beeps_2.mat'
    ...'protocol_45s_80dB.mat'
    'protocol_morning.mat'
    });

% TO TEST:
% Files.data(5:end)=[];


%PROPERTIES / OPTIONS
%figure to plot
opt.signal.plot   = true;
opt.signal.margin = [0.7,0.5]; %margin [s], before and after beep
opt.minmax.plot   = true;
opt.minmax.fun    = @min; %@min or @max for plot
opt.minmax.margin = [0,0]; %not larger than opt.signal.margin
% table options (50 for min, 100 for max is fine)
opt.table.dur = [50 100];  %[ms] range from beep
opt.table.fun = opt.minmax.fun; %@min or @max for .xls 
%average type: 1 = mean across all, 1 = mean across individual mean
tmp = {... {type, text}, text in title unless only one rPaths
    1,'Mean Across Pooled Data';...
    2,'Mean Across Individual Mean';...
    };
[opt.average.type,opt.average.str] = tmp{2,:}; %SELECT ONE ROW!
%stim variable ('stimTimes', 'stimTimes2', ...)
opt.stimulus.str = 'stimTimes';
%stages definition
opt.Stages = {... {number in hypnogram, label}
    1,'Wake';...
    2,'NREM';...
    3,'REM';...
    };
%signal filter (lower is 0 if NaN, upper is fs/2 if NaN);
opt.bandpass = [... [low,up; low,up; ...]
    0,45;... remove 50 Hz,
    55,NaN;...
    ];

%PLOT PROPERTIES
%---------------
%GENERAL, ALL FIGURES & PLOTS
%colors (numeric row for each beep power; auto if not enough)
props.all.color.power  = [0 0 1; 0.9290 0.6940 0.1250; 1 0 0]; %colors
props.all.time.factor  = 1000; %time factor to unit seconds
props.all.time.unit    = 'ms';   %time unit based on time factor
props.all.ylimPerStage = false;

%SIGNAL PLOT
%figure axes
props.signal.figure = {'position',[440,330,880,600]};
props.signal.fig_createAxes = {[60,30,10],[50,30,63],'pixel'};
props.signal.axes = {'visible','on','ylim',[-1,1]*1,'ytick',-50:25:50};
%plottings
props.signal.plot = {'visible','on'}; %average data
props.signal.fill = {'edgecolor','none','facecolor',[1,1,1]*.7}; %beep area
%text
props.signal.xlabel = {'visible','on'};
props.signal.ylabel = {'visible','on'};
props.signal.title  = {'fontweight','normal'};
props.signal.superTitle = {'fontweight','bold'};
props.signal.text   = {'visible','on'}; %all text(...)

%MINMAX PLOT
props.minmax.figure={'position',[440,330,880,600]};
props.minmax.fig_createAxes={[60,30,10],[50,30,63],'pixel'};
props.minmax.axes={'visible','on',...
    'xtick',0:25:100}; %'ylim',[-1,1]*300
%plottings
props.minmax.bar={'visible','on','BarWidth',2,'edgecolor','k'};
%text
props.minmax.xlabel={'visible','on'};
props.minmax.ylabel={'visible','on'};
props.minmax.title={'fontweight','normal'};
props.minmax.superTitle={'fontweight','bold'};
props.minmax.text={'visible','on'}; %all text(...)

%% MAIN PROGRAM
%---------------
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)));
%number of ...
noFIL = numel(Files.data);
noSTA=size(opt.Stages,1);  %stages

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
noFIL = numel(Files.data);
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

%LOAD PROTOCOLS
Protocols = f_protocols_load(Files.prot);
noPTC = numel(Files.prot);
for ptc = 2:noPTC %must have equal signals (different timing allowed)
    if ~isequaln(Protocols(1).beep,Protocols(ptc).beep) || ...
            ~isequaln(Protocols(1).laser,Protocols(ptc).laser)
        error(sprintf(['All Protocols must have same beeps/laser\n',...
            '(different beep timing is OK)']))
    end
end
%variable definitions from protocol
beeps = Protocols(1).beep;
laser = Protocols(1).laser;
freqs  = unique([beeps.Signals{:,1}]);
powers = unique([beeps.Signals{:,2}]);
%number of ...
noBIN  = beeps.duration*fs; %beep bins
noPOW = numel(powers);
noFRQ = numel(freqs);

%CHANNEL NAMES (only for title)
channels = cellfun(@(x) x(...
    find(x==filesep,1,'last')+1:find(x=='.',1,'last')-1),...
    Files.data,'uniformoutput',false);

%% LOAD DATA
fprintf('\nLOAD DATA:\n'); tic
DATA = cell(noFIL,1); %structure per data file
nnFIL  = numel(num2str(noFIL));
indent = blanks(2*nnFIL+1); %indent for fprintf
for fil=1:noFIL
    clear data Hypnogram P %clean up
    fileDAT = Files.data{fil};
    fileHYP = fullfile(fileparts(fileDAT),'Hypnogram.mat');
    fprintf('%*i/%i File: %s\n',nnFIL,fil,noFIL,...
        strrep(fileDAT,mPath,'...'));

    %DATA
    data = load(fileDAT);
    %stimFile= fileDAT(end-12:end);
    %stimFile= [mPath '\' stimFile];
    % load(stimFile, 'stimTimes');
    % data.stimTimes=stimTimes;
    %filter data
    if ~isempty(opt.bandpass)
        signal = passband_fourier(data.resampled_data_mV,...
            opt.bandpass,fs);
        fprintf('%s Data filtered\n',indent)
    else
        signal = data.resampled_data_mV;
    end
    noSAM = numel(signal); %number of samples
    %stimuli
    stimTimes = round(data.(opt.stimulus.str)); %rounded index
    stimSTA = stimTimes(1:2:end);
    stimEND = stimTimes(2:2:end);
    
    %HYPNOGRAM (and upsampling)
    load(fileHYP,'Hypnogram');
    fac = round(noSAM/numel(Hypnogram)); %upsample factor
    if fac>1
        fprintf('%s Hypnogram upsampled by factor %i\n',indent,fac)
        Hypnogram = repmat(Hypnogram(:)',fac,1);
    elseif fac<1
        error('uups')
    end
    Hypnogram = Hypnogram(:);
    Hypnogram(end+1:noSAM)=NaN;
    stages = Hypnogram(round((stimSTA+stimEND)/2));
    if any(isnan(stages))
        warning(['Hypnogram not fully scored\n',...
            '%i stimuli occurs in stage NaN'],sum(isnan(stages)))
    end
    
    %PROTOCOL (find correct protocol)
    P = f_protocol_check(Protocols,stimTimes,fs);
    [~,tmp] = fileparts(P.protocol.file);
    fprintf('%s Protocol: %s.mat:\n',indent,tmp)
    %remove NaN's
    ind = isnan(P.repetition);
    if any(ind)
        tmp = fieldnames(P);
        for k = 1:numel(tmp)
            if isequal(size(P.(tmp{k})),size(ind))
                P.(tmp{k})(ind,:) = [];
            end
        end
    end
    %append
    P.noBEP = numel(P.stimSTA); %valid beeps
    P.noREP = numel(unique(P.repetition));
    fprintf('%s   %i beeps, repeated %ix\n',indent,P.noBEP,P.noREP)
    
    %APPEND SIGNAL / BEEP STAGES
    P.stages  = Hypnogram(round((P.stimSTA+P.stimEND)/2));
    P.data = signal;
    if any(isnan(P.stages))
        warning(['Hypnogram not fully scored\n',...
            '%i stimuli occurs in stage NaN'],sum(isnan(stages)))
    end
    DATA(fil)={P};
end %file loop
toc

%% BEEP DATA
fprintf('\nGET BEEP DATA & AVERAGE\n'); tic
margin    = abs(opt.signal.margin)*fs; %in bins
tBIN      = -margin(1):margin(2)+noBIN; %in bins (bin 0 in addition)
indMIMA   = find(tBIN>=0&tBIN<noBIN);
t         = (tBIN)/fs*props.all.time.factor;
DataBeeps = cell(noFIL,noSTA,noFRQ,noPOW);
for fil = 1:noFIL
    P = DATA{fil};
    %power, freq and stages loop
    for pow = 1:noPOW
        indPOW = P.power==powers(pow);
        for frq = 1:noFRQ
            indFRQ = P.frequency==freqs(frq);
            for sta = 1:noSTA
                stage = opt.Stages{sta,1};
                ind = indPOW & indFRQ & P.stages==stage;
                datSTA = P.stimSTA(ind)-margin(1);
                datEND = P.stimEND(ind)+margin(2);
                M = NaN(numel(datSTA),numel(tBIN));
                for k = 1:numel(datSTA)
                    M(k,:) = P.data(datSTA(k):datEND(k));
                end
                DataBeeps{fil,sta,frq,pow} = M;
            end
        end
    end
end

%% MEAN DATA (two different means) & MIN TABLE
tmp=size(DataBeeps);
mDATA1=cell(tmp(2:end)); %mean across all files
mDATA2=cell(tmp(2:end)); %mean across individual means
mERR1=cell(tmp(2:end)); %mean across all files
mERR2=cell(tmp(2:end)); %mean across individual means
%index to first N ms after beep start
indPOI = (opt.table.dur(1)/1000*fs+1:opt.table.dur(2)/1000*fs)+margin(1);
tmp = {'Stage','Frequency','Power','Peak Delay [ms]','Peak Amp','Comment','File'};
table = cell(noFIL*noSTA*noFRQ*noPOW+1,numel(tmp));
table(1,:) = tmp;
cnt = 1;
%loops
for sta=1:noSTA
    stage = opt.Stages{sta,2};
    for frq=1:noFRQ
        freq = freqs(frq);
        for pow=1:noPOW
            power = powers(pow);
            mTOT2=[];
            mTOT1=[];
            %recording loop
            for fil=1:noFIL
                file = Files.data{fil};
                %data
                data  = DataBeeps{fil,sta,frq,pow};
                if isempty(data)
                    NaN(size(data));
                else
                    mData = nanmean(data,1); %mean
                end
                mTOT1=[mTOT1;data];
                mTOT2=[mTOT2;mData];
                %table
                cnt = cnt+1;
                [m,ind] = opt.table.fun(mData(indPOI));
                tmp = (ind+indPOI(1)-1)/fs*1000-margin(1);
                table(cnt,:) = {stage,freq,power,tmp,m,'',file};
            end
            mDATA1{sta,frq,pow}=nanmean(mTOT1,1);
            mDATA2{sta,frq,pow}=nanmean(mTOT2,1);
            mERR1{sta,frq,pow}=nanstd(mTOT1,[],1)/sqrt(size(mTOT1,1));
            mERR2{sta,frq,pow}=nanstd(mTOT2,[],1)/sqrt(size(mTOT2,1));
        end
    end
end
table(2,6) = {sprintf(['time point of %s value from beep start in ',...
    '[%g %g] after beeps start'],func2str(opt.table.fun),opt.table.dur)};
%select mean data
switch opt.average.type
    case 1
        mDATA=mDATA1;
        mERR=mERR1;
        str='Average across pooled data of all recordings';
    case 2
        mDATA=mDATA2;
        mERR=mERR2;
        str='Average across average data of individual recordings';
    otherwise
        error('opt.average.type = %i is not supported!',opt.average.type)
end
%minmax
fprintf('  %s\n',str)
toc

% FIGURES
%% PLOT SIGNALS
fprintf('\nPLOT\n'); tic
%plot colors per power
if size(props.all.color.power,1)<noPOW
    tmp=lines(noPOW);
    props.all.color.power=[props.all.color.power;...
        NaN(noPOW-size(props.all.color.power,1),3)];
    for k=find(isnan(props.all.color.power(:,1)),1,'first'):noPOW
        props.all.color.power(k,:)=tmp(find(~ismember(tmp,...
            props.all.color.power,'rows'),1,'first'),:);
    end
end
HF=[];

%PLOT SIGNAL
if opt.signal.plot
    fprintf('  - Signals\n')
    %FIGURE
    hf=figure(props.signal.figure{:}); drawnow
    HF(end+1)=hf;
    fig_size(hf);
    ha=fig_createAxes(hf,[noFRQ,noSTA],props.signal.fig_createAxes{:});
    MiMa=NaN(noSTA,noFRQ,noPOW,2);
    %DATA LOOPS /PLOT
    for sta=1:noSTA
        stage=opt.Stages{sta,2}; %stage label
        for frq=1:noFRQ
            frequency=freqs(frq);
            set(hf,'CurrentAxes',ha(frq,sta))
            %init
            hp=NaN(noPOW,1);
            legSTR=cell(noPOW,1);
            for pow=1:noPOW
                power=powers(pow);
                legSTR{pow}=sprintf('%g dB',power);
                mDAT=mDATA{sta,frq,pow};
                err=mERR{sta,frq,pow};
                MiMa(sta,frq,pow,:)=[min(mDAT),max(mDAT)];
                %plot
                col=props.all.color.power(pow,:);
                hp(pow)=plot(t,mDAT,'color',col,props.signal.plot{:});
                hold on
                %error
                x=[t,t(end:-1:1)];
                y=[mDAT-err,mDAT(end:-1:1)+err(end:-1:1)];
                h=fill(x,y,zeros(size(x)),'edgecolor','none',...
                    'facecolor',col,'facealpha',0.2);
                uistack(h,'bottom')
                
            end
            %text
            if frq==1
                title(stage,props.signal.title{:})
            end
            if frq==noFRQ
                xlabel(sprintf('Time [%s]',...
                    props.all.time.unit),props.signal.xlabel{:})
            end
            if sta==1
                ylabel(sprintf('%g Hz',frequency),...
                    props.signal.ylabel{:})
            end
            drawnow
        end
    end
    %SETTINGS / TEXT
    %linkaxes and specific settings
    if props.all.ylimPerStage
        for k=1:noSTA
            linkaxes(ha(:,k),'xy')
            tmp=MiMa(k,:,:,:);
            set(ha(:,k),'ylim',[-1,1]*ceil(max(abs(tmp(:)))))
        end
    else
        linkaxes(ha,'xy')
        set(ha,'ylim',[-1,1]*ceil(max(abs(MiMa(:)))))
        set(ha(:,2:end),'yticklabel',[])
    end
    %text & settings
    if noFIL==1
        str=strrep(Files.data{:},mPath,'...');
        h=fig_superLabel(ha,'title',{'Indivudual Data',regexprep(str,...
            {'\\','\_'},{'\\\\','\\_'})});
    else
        str=sprintf(', %s',channels{:});
        h=fig_superLabel(ha,'title',{opt.average.str,str(3:end)});
    end
    set(h,props.signal.superTitle{:})
    set(ha(1:end-1,:),'xticklabel',[])
    set(ha,'xlim',[floor(tBIN(1)),ceil(tBIN(end))],...
        'unit','normalized',...
        props.signal.axes{:})
    %beep areas
    yl=get(ha(1),'ylim'); %actual axis from superLabel
    X=[0,0,beeps.duration,beeps.duration]*props.all.time.factor;
    Y=[yl(1),yl(2),yl(2),yl(1)];
    Z=zeros(1,4);
    for axi=1:numel(ha)
        set(hf,'CurrentAxes',ha(axi));
        h=fill(X,Y,Z,props.signal.fill{:});
        uistack(h,'bottom')
        text(tBIN(margin(1)+noBIN/2),yl(1),'beep',...
            'horizontalalignment','center',...
            'verticalalignment','bottom',props.signal.text{:})
    end
    %legend
    legend(hp,legSTR)
end %plot signal

%% PLOT MIN OR MAX
if opt.minmax.plot
    funSTR=func2str(opt.minmax.fun);
    fprintf('  - %s data',funSTR)
    [M,IND]=cellfun(@(x)opt.minmax.fun(x(indMIMA)),mDATA);
    tM=t(indMIMA);
    xl=t([indMIMA(1),indMIMA(end)+1]);
    xl=xl+[-1,1]*0.1*diff(xl);
    
    %FIGURE
    hf=figure(props.minmax.figure{:}); drawnow
    HF(end+1)=hf;
    fig_size(hf);
    ha=fig_createAxes(hf,[noFRQ,noSTA],props.minmax.fig_createAxes{:});
    %DATA LOOPS /PLOT
    for sta=1:noSTA
        stage=opt.Stages{sta,2}; %stage label
        for frq=1:noFRQ
            frequency=freqs(frq);
            set(hf,'CurrentAxes',ha(frq,sta))
            %init
            [~,inds]=sort(squeeze(M(sta,frq,:)),'descend'); %sort index
            hp=NaN(noPOW,1);
            legSTR=cell(noPOW,1);
            for pow=1:noPOW
                ind=inds(pow);
                col=props.all.color.power(ind,:);
                power=powers(ind);
                legSTR{ind}=sprintf('%g dB',power);
                hp(ind)=bar(tM(IND(sta,frq,ind)),M(sta,frq,ind),...
                    'facecolor',col,'edgecolor',col,props.minmax.bar{:});
                hold on
            end
            %text
            if frq==1
                title(stage,props.minmax.title{:})
            end
            if frq==noFRQ
                xlabel(sprintf('Time [%s]',props.all.time.unit),props.minmax.xlabel{:})
            end
            if sta==1
                ylabel(sprintf('%g Hz',frequency),props.minmax.ylabel{:})
            end
            drawnow
        end
    end
    
    %SETTINGS / TEXT
    %linkaxes and specific settings
    if props.all.ylimPerStage
        for k=1:noSTA
            linkaxes(ha(:,k),'xy')
            tmp=M(k,:,:);
            if strcmpi(funSTR,'min')
                yl=[min(tmp(:)),max([0;tmp(:)])];
            else
                yl=[min([0;tmp(:)]),max(tmp(:))];
            end
            dy=0.02;
            if yl(1)~=0
                yl(1)=yl(1)+sign(yl(1))*dy*diff(yl);
            end
            if yl(2)~=0
                yl(2)=yl(2)+sign(yl(2))*dy*diff(yl);
            end
            set(ha(:,k),'ylim',[floor(yl(1)),ceil(yl(2))])
        end
    else
        linkaxes(ha,'xy')
        if strcmpi(funSTR,'min')
            yl=[min(M(:)),max([0;max(M(:))])];
        else
            yl=[min([0;min(M(:))]),max(M(:))];
        end
        dy=0.05;
        if yl(1)~=0
            yl(1)=yl(1)+sign(yl(1))*dy*diff(yl);
        end
        if yl(2)~=0
            yl(2)=yl(2)+sign(yl(2))*dy*diff(yl);
        end
        set(ha(:,1),'ylim',[floor(yl(1)),ceil(yl(2))])
        set(ha(:,2:end),'yticklabel',[])
    end
    set(ha(1:end-1,:),'xticklabel',[])
    set(ha,'xlim',xl,'xgrid','on','ygrid','on',...
        'xtick',[0,noBIN/fs*props.all.time.factor],...
        'unit','normalized',props.minmax.axes{:})
    %text
    str=sprintf('%s - PEAK',upper(funSTR));
    if noFIL==1
        h=fig_superLabel(ha,'title',{'Indivudual Data',str});
    else
        h=fig_superLabel(ha,'title',{opt.average.str,str});
    end
    set(h,props.minmax.superTitle{:})
    %legend
    legend(hp,legSTR)
    fprintf('\n')
end %plot minmax
toc

%% SAVE table
fprintf('\n')
[sFile,sPath] = uiputfile('*.xls;*.xlsx','Save min peak delay',...
    ['Z:\1 TIDIS Lab\Ida\Recording\\Analysis\Max Min Peaks\Z-scored',...
    'Max Min Peaks\Z-scored\PeakDelay.xls']);
if ischar(sFile)
    sname = fullfile(sPath,sFile);
    if exist(sname,'file')==2
        delete(sname);
    end
    xlswrite(sname,table)
    try
        xls_cellFit(sname)
    catch
    end
    
    fprintf('Min peak delays saved as:\n')
    fprintf(' %s\n',sname)
    
else
    fprintf('[\bMin Peak Delays NOT Saved!]\b\n')
end

return %SAVING EXLUDED (see analyze_audio.m, it's in there)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE
if saveFIG
    fprintf('Save figures\n')
    for cha=1:numel(channels)
        hf=HF(cha); figure(hf);
        fprintf(' Figure %i:\n',get(hf,'Number'))
        %save name
        [sPath,sFile,sExt]=fileparts(saveNAM(channels{cha},noFIL));
        if isempty(sPath)
            sPath=rPaths{1};
            ind=cellfun(@(x) strncmpi(sPath,x,numel(sPath)),rPaths);
            while ~all(ind)
                sPath=fileparts(sPath);
                ind=cellfun(@(x) strncmpi(sPath,x,numel(sPath)),rPaths);
            end
        end
        if isempty(sExt)
            sExt='.fig'; %for file-selection, real defined in saveFUNs
        end
        [sFile,sPath]=uiputfile('*.*','Save Figure As (no extension)',...
            fullfile(sPath,sFile));
        if ~ischar(sFile)
            fprintf('\b [\bNOT saved!]\b\n')
            continue
        end
        [~,sFile,~]=fileparts(sFile);
        sname=fullfile(sPath,sFile);
        fprintf('\b %s.*\n',sname)
        %save
        for k=1:numel(saveFUNs)
            fun=saveFUNs{k};
            fun(hf,sname);
        end
    end
else
    fprintf('[\bFigures NOT saved]\b\n')
end