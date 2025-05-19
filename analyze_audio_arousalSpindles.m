% analyze_audio_arousalSpindles
%------------------------------------------------------------------------
% ...
%tabel & figure
% - per delay

clc; clear; close all;
[scriptPath,scriptName] = fileparts(which(mfilename));
addpath(fullfile(scriptPath,'functions'))

%PARAMETERS
%----------
%FILES DATA
% fileList = ['Z:\1 TIDIS Lab\Ida\Script\Audio Analysis_PulsePal_3\',...
%     'filelists\Audio\Baseline\CMT_A10.txt'];
fileList = ['Z:\1 TIDIS Lab\Ida\Script\Audio Analysis_PulsePal_3\',...
    'filelists\Audio\CMT_z-scored.txt'];
[Files.data,fileList] = import_fileList(fileList);
% Files.data(4:end) = []


%FILES HYPNOGRAMS/SPINDLES
tmp = cellfun(@fileparts,Files.data,'uniformoutput',false); %read paths
Files.hypno = fullfile(tmp,'Hypnogram.mat');
EEG = 'EEG1';
Files.spind = fullfile(tmp,sprintf('spindle_%s.mat',EEG));
%Files.spind = []; %test, no spindles if empty

%FILES PROTOCOL
%  - cell; finds the one that fits to data stimTimes
%    But it's more save when you know the correct one
Files.proto = fullfile(scriptPath,'protocols',{... possible protocols
    ...'protocol_morning.mat';'protocol_afternoon.mat';...
    ...'protocol_laser_3600s_80dB_108x3beeps.mat';...
    ...'protocol_laser_3600s_80dB_108x3beeps_2.mat';...
    ...'FC_protocol_1h_342beeps.mat';...
    ...'FC_protocol_45s_14beeps.mat';...
    ...'protocol_morning.mat';...
    ...'protocol_laser_3600s_80dB_108x3beeps_2.mat'    
    ...'protocol_laser_600s_80dB_27x2beeps.mat'
    ...'protocol_laser_3600s_80dB_108x3beeps.mat'
    'protocol_morning.mat';...
    });

%FILE SAVENAME, count transitions (askes again)
[~,tmp] = fileparts(fileList);
snameCNT = fullfile(['Z:\1 TIDIS Lab\Ida\',...
    'Analysis\Arousal'],sprintf('Arousals_%s.xlsx',tmp));

%ANALYSIS OPTIONS
ana.stimTimes   = 'stimTimes'; %'stimTimes' or 'stimTimes2' or ...
ana.powers      = 80; %all if empty (pooled together, no single analysis)
ana.frequencies = [5000]; %all if empty (pooled together)
%transitions
ana.Transitions = {...{stage from, stages to, label, color, stages perc}
    ... stages perc: transitions contributing to percentage calc
    ... one axis per transition start state
    2 ,1 ,'NREM-Wake' ,'r',[1,2];...
    2 ,2 ,'NREM-NREM' ,'b',[1,2];...
%     3 ,1 ,'REM-Wake' ,'r',[1,3];...
%     3 ,3 ,'REM-REM'  ,'b',[1,3];...
    ...2 ,1:2 ,'NREM-Sleep'  ,'g',[1,3];...
    };
%transitions count by delay time window
ana.windows = [10]; %[s]
%plot margin, plotting ± from beep start/end
ana.margin = [-0.7,0.6]; %[s]
%signal filter (PS: uses fs/2 for NaNs in 2nd column)
ana.bandpass = [... [lower,upper; lower,upper; ...]
    0,45;... remove components around 50 Hz
    55,NaN;...
    ];

%AVERAGE DATA
%  - 1 for mean across all available beeps (pooled data)
%    --> individuals with more beeps (per conditions) contributes more
%        (higher weighted individuals)
%  - 2 for mean across individual mean
%    --> all individuals contributes equally (unless individuals that have
%        no beep for a specific condition)
%  NOTE: the same for standard error. E.g. if 2 and only one individual
%        --> standard error will be zero, even it's an average across
%            several beeps.
tmp = {... {type, label}, text in title unless only one rPaths
    1,'Mean Across Pooled Data';...
    2,'Mean Across Individual Mean';...
    };
[ana.average.type,ana.average.label] = tmp{2,:}; %SELECT ONE ROW!

%RATIO OPTIONS
% Calculates ratio per second (using function f_ratio.m)
% If x is the data vector for calculating the ratio, sampled with fs:
%  - winLen defines the window length in [s]
%    Calculating for each window: sum(x(window))/winLen
%    --> Ratio is sum of events per second for that window
%        (i.h., if you need the ratio per minute, multiply ratio with 60)
%  - stepLen defines moving window step in [s]
%    It's somehow the resolution ratio
%    Note: - if stepLen < sample duration of x, uses stepLen = 1/fs.
%            So if you set stepLen = 0, uses highest possible resolution.
%          - if stepLen is empty, NaN or > winLen, uses stepLen = winLen
% PS: - in old scripts, stepLen was fixed to winLen!
ana.ratio.type = 'count'; %count or ratio [1/s]
ana.ratio.winLen = 0.1; %[s]
ana.ratio.stepLen = 0.1; %[s] (for ratio, count uses winLen)
ana.ratio.fac = 60; %factor to ratio in [1/s]
ana.ratio.unit = '[min^{-1}]'; %corresponding unit of fac (for plot labels)

%PLOT PROPERTIES (look&feel)
%figure / axes
props.figure = {'visible','on'};
Ylim = [-1,1]; %auto if empty, for ylim
props.axes   = {'xtick',-3:3}; %y-axes are beeps without ylim
%fig_createAxes
props.fig_createAxes = {[75,40,12],[50,30,80],'pixel'}; %basic
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
%number of ...
noFIL = numel(Files.data);
%noSTA = size(ana.Stages,1);
noTRA = size(ana.Transitions,1);
noWIN = numel(ana.windows);
%sort transition stages
ana.Transitions(:,[1,2,end]) = cellfun(@sort,...
    ana.Transitions(:,[1,2,end]),'uniformoutput',false);


%CHECK FILES
%inexisting files
withSPN = ~isempty(Files.spind);
if withSPN
    ind = ...
        cellfun(@(x)exist(x,'file')~=2,Files.data)  | ...
        cellfun(@(x)exist(x,'file')~=2,Files.spind) | ...
        cellfun(@(x)exist(x,'file')~=2,Files.hypno);
else
    ind = ...
        cellfun(@(x)exist(x,'file')~=2,Files.data)  | ...
        cellfun(@(x)exist(x,'file')~=2,Files.hypno);
end
if any(ind)
    fprintf('[\bData excluded because of missing files (%i of %i)]\b\n',...
        sum(ind),noFIL)
    fprintf('  - %s\n',Files.data{ind})
    Files.data(ind)  = [];
    Files.hypno(ind) = [];
    if withSPN
        Files.spind(ind) = [];
    end
end
%duplicate files (data files only)
ind = true(size(Files.data));
[~,tmp] = unique(lower(Files.data));
ind(tmp) = false;
if any(ind)
    fprintf('[\bFiles duplicate (kept only once)]\b\n')
    fprintf('  - %s\n',Files.data{ind})
    Files.data(ind)  = [];
    Files.hypno(ind) = [];
    if withSPN
        Files.spind(ind) = [];
    end
end
%remaining files
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
Protocols = f_protocols_load(Files.proto);
noPTC = numel(Files.proto);
for ptc = 2:noPTC %must have equal signals (different timing allowed)
    if ~isequaln(Protocols(1).beep.duration,Protocols(ptc).beep.duration) || ...
            ...~isequaln(Protocols(1).beep.Signals,Protocols(ptc).beep.Signals) || ...
            ~isequaln(Protocols(1).laser,Protocols(ptc).laser)
        error(sprintf(['All Protocols must have same beeps/laser\n',...
            '(different beep timing is OK)']))
    end
end
%check variables
tmp = unique([Protocols(1).beep.Signals{:,2}]); %powers
if isempty(ana.powers) || any(isnan(ana.powers))
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
tmp = unique([Protocols(1).beep.Signals{:,1}]); %frequencies
if isempty(ana.frequencies) || any(isnan(ana.frequencies))
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
%check frequencies & power
for ptc = 1:noPTC
     p = Protocols(ptc);
     tmp = cell2mat(p.beep.Signals);
     if ~all(ismember(ana.frequencies,tmp)&ismember(ana.powers,tmp))
         error('uups')
     end
end

%% READ DATA
fprintf('\nREAD DATA\n'); tic
nnFIL = numel(num2str(noFIL));
indent = blanks(2*nnFIL+2); %indent for fprintf
for fil=1:noFIL
    fileDAT = Files.data{fil};
    fileHYP = Files.hypno{fil};
    [rPath.dat,rFile.dat,rExt.dat] = fileparts(fileDAT);
    [rPath.hyp,rFile.hyp,rExt.hyp] = fileparts(fileHYP);
    if ~strcmpi(rPath.dat,rPath.hyp)
        error('uups')
    end
    fprintf('%*i/%i: %s\n',nnFIL,fil,noFIL,rPath.dat)
    %spindles
    if withSPN
        fileSPN = Files.spind{fil};
        [rPath.spn,rFile.spn,rExt.spn]=fileparts(fileSPN);
        if ~strcmpi(rPath.dat,rPath.spn)
            error('uups')
        end
    end
    
    %DATA (recording)
    fprintf('%s File Data      : %s\n',indent,[rFile.dat,rExt.dat])
    data      = load(fileDAT);
    signal    = data.resampled_data_mV(:);
    stimTimes = round(data.(ana.stimTimes)); %rounded index
    noSAM     = numel(signal); %number of samples
    %filter data
    if ~isempty(ana.bandpass)
        tmp = ana.bandpass;
        %for print out (basically passband_fourier would made it self)
        tmp(isnan(tmp(:,1)),1) = 0; %made
        tmp(isnan(tmp(:,2)),2) = fs/2;
        %filter
        signal = passband_fourier(signal,tmp,fs);
        %print out
        str = sprintf(' ,%g-%g',tmp');
        fprintf('%s   Bandpass %s Hz\n',indent,str(3:end))
    end
    
    %HYPNOGRAM (and upsampling)
    fprintf('%s File Hypnogram : %s\n',indent,[rFile.hyp,rExt.hyp])
    clear Hypnogram
    load(fileHYP,'Hypnogram');
    %upsample
    fac = round(noSAM/numel(Hypnogram));
    if fac>1
        fprintf('%s   Upsampled by factor %i\n',indent,fac)
        Hypnogram = repmat(Hypnogram(:)',fac,1);
    elseif fac<1
        error('uups')
    end
    Hypnogram = Hypnogram(:);
    Hypnogram(end+1:noSAM) = NaN;
    %test scoring
    tmp = Hypnogram(stimTimes(1:2:end)); %at stim start
    if any(isnan(tmp))
        warning(['Hypnogram not fully scored\n',...
            '%i stimuli occurs in stage NaN'],sum(isnan(tmp)))
    end
    
    %PROTOCOL (find correct protocol)
    P = f_protocol_check(Protocols,stimTimes,fs);
    [~,rFile.pro,rExt.pro] = fileparts(P.protocol.file);
    fprintf('%s File Protocol  : %s\n',indent,[rFile.pro,rExt.pro])
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
    noSTI = numel(P.stimSTA);
    %print out
    fprintf('%s   Total beeps, N = %i\n',indent,numel(ind))
    fprintf('%s   Excl. beeps, N = %i\n',indent,sum(ind));
    fprintf('%s   Remaining  , N = %i\n',indent,noSTI);
    fprintf('%s   Protocol repetitions, N = %i\n',indent,...
        numel(unique(P.repetition)))
    
    %SPINDLES (center)
    if withSPN
        fprintf('%s File Spindle   : %s\n',indent,[rFile.spn,rExt.spn])
        tmp = load(fileSPN);
        spindles = false(size(signal));
        %index spindle midpoints, in sample rate of signal
        ind = ([tmp.spindle.end]+[tmp.spindle.str])/2/tmp.SampRate*fs;
        spindles(round(ind)) = true;
    end
    
    %INIT
    fieldsV = {'stimSTA','stimEND',... fields vectors
        'power','frequency','stage','nextStim','nextStage',...
        'nextIndex','nextDelay'};
    fieldsM = {'Data','Spindles'}; %fields matrix
    if fil==1
        margin = abs(ana.margin);
        timeIND = -round(margin(1)*fs):...
            round((margin(2)+P.protocol.beep.duration)*fs);
        t = timeIND/fs;
        for k = 1:numel(fieldsV)
            DATA.(fieldsV{k}) = NaN(noSTI,noFIL);
        end
        for k = 1:numel(fieldsM)
            DATA.(fieldsM{k}) = NaN(noSTI,numel(timeIND),noFIL);
        end
    elseif size(DATA.Data,1)<noSTI %increase data size
        for k = 1:numel(fieldsV)
            DATA.(fieldsV{k})(end+1:noSTI,:) = NaN;
        end
        for k = 1:numel(fieldsM)
            DATA.(fieldsM{k})(end+1:noSTI,:,:) = NaN;
        end
    end
    
    %APPEND
    %from protocol
    ind = 1:noSTI; %if there are less stims as in previous ones
    DATA.stimSTA(ind,fil)   = P.stimSTA;
    DATA.stimEND(ind,fil)   = P.stimEND;
    DATA.power(ind,fil)     = P.power;
    DATA.frequency(ind,fil) = P.frequency;
    DATA.stage(ind,fil)     = Hypnogram(P.stimSTA);
    DATA.nextStim(ind,fil) = [P.stimSTA(2:end);noSAM];
    for sti = 1:noSTI
        %stim index & stage
        ind1 = DATA.stimSTA(sti,fil);
        indE = DATA.nextStim(sti,fil);
        stage1 = Hypnogram(ind1);
        if ~isequaln(stage1,DATA.stage(sti,fil))
            error('uups')
        end
        %next ind/stage
        ind2 = find(Hypnogram(ind1+1:indE)~=stage1,1,'first')+ind1;
        if isempty(ind2)
            ind2 = ind1;
        end
        %append
        DATA.nextIndex(sti,fil) = ind2;
        DATA.nextStage(sti,fil) = Hypnogram(ind2);
        DATA.nextDelay(sti,fil) = (ind2-ind1)/fs;
        
        %SIGNALS
        ind = timeIND+ind1;
        DATA.Data(sti,:,fil) = signal(ind);
        if withSPN
            DATA.Spindles(sti,:,fil)= spindles(ind);
        end
    end
end
toc

%% AVERAGE/COUNT DATA
fprintf('\nAVERAGE/GROUP DATA\n'); tic
mDATA(1:noTRA,1:noWIN) = {NaN(0,numel(timeIND))};
mSPIN = mDATA;
mCNT  = zeros(noTRA,noWIN); %count transitions
mPER  = zeros(noTRA,noWIN); %count for percentage calculation
MCNT  = zeros(noTRA,noWIN,noFIL); %count transitions
MPER  = zeros(noTRA,noWIN,noFIL); %count for percentage calculation
for fil = 1:noFIL
    %index selected power & frequency
    indSEL = ... PS: NaN's will be excluded
        ismember(DATA.power(:,fil),ana.powers) & ...
        ismember(DATA.frequency(:,fil) ,ana.frequencies);
    %transition loop
    for tra = 1:noTRA
        [stages1,stages2,label,col,stagesP] = ana.Transitions{tra,:};
        indTRA = ismember(DATA.stage(:,fil),stages1); %start
        %window loop
        for win = 1:noWIN
            indWIN = ... add start stagees2 (excl with indTRA)
                (DATA.nextDelay(:,fil)<=ana.windows(win) & ...
                ismember(DATA.nextStage(:,fil),stages2)) | ...
                (DATA.nextDelay(:,fil)>ana.windows(win) & ...
                ismember(DATA.stage(:,fil),stages2));
            %   NOTE: - 1st part, transition to stages2 within time
            %           window (for stages1 ~= stages2).
            %         - 2nd part, if no transition within time window but
            %           stages1 is equal to stages2. If stages1 ~= stages2
            %           indTRA will remove this again.
            indPER = ...
                (DATA.nextDelay(:,fil)<=ana.windows(win) & ...
                ismember(DATA.nextStage(:,fil),stagesP)) | ...
                (DATA.nextDelay(:,fil)>ana.windows(win) & ...
                ismember(DATA.stage(:,fil),stagesP));
            ind  = indSEL & indTRA & indWIN;
            indP = indSEL & indTRA & indPER;
            mCNT(tra,win) = mCNT(tra,win)+sum(ind);
            mPER(tra,win) = mPER(tra,win)+sum(indP);
            MCNT(tra,win,fil) = sum(ind);
            MPER(tra,win,fil) = sum(indP);
            
            %spindles
            if withSPN
                mSPIN{tra,win} = [mSPIN{tra,win};...
                    DATA.Spindles(ind,:,fil)];
            end
            %data
            switch ana.average.type
                case 1 %pool data
                    mDATA{tra,win} = [mDATA{tra,win};...
                        DATA.Data(ind,:,fil)];
                case 2 %average individual data
                    mDATA{tra,win} = [mDATA{tra,win};...
                        mean(DATA.Data(ind,:,fil),1)];
                otherwise
                    error('ana.average.type = %i is not implemented!',...
                        ana.average.type)
            end
        end %window loop
    end  %transition loop
end
%grande average of signale
for k=1:numel(mDATA)
    %data
    if k==1
        mDAT = cell(size(mDATA));
        mERR = cell(size(mDATA));
    end
    mDAT{k} = nanmean(mDATA{k},1);
    mERR{k} = nanstd(mDATA{k},[],1)/sqrt(size(mDATA{k},1));
    %spindle
    if withSPN
        if k==1
            mSPN = cell(size(mSPIN));
            cSPN = cell(size(mCNT));
        end
        tmp = mSPIN{k};
        mSPN{k} = any(tmp,1);
        cSPN{k} = sum(tmp,1);
    end
end
toc

%% PLOT
fprintf('\nPLOT DATA\n'); tic
if true %for testing, exclude plotting
    % clc; close all
    %axes for transisions
    axesNUM = zeros(noTRA,1);
    for tra = 1:noTRA
        num = ana.Transitions{tra,1}; %sort if more than one stage
        ind = find(cellfun(@(x)isequal(sort(x),sort(num)),...
            ana.Transitions(1:tra-1,1)),1,'first');
        if isempty(ind)
            axesNUM(tra) = max(axesNUM)+1;
        else
            axesNUM(tra) = axesNUM(ind);
        end
    end
    noAXI = numel(unique(axesNUM));
    %figure position
    dx = props.fig_createAxes{1}(:);
    posF = [440,330,[1,noAXI-1,1]*dx+380*noAXI,440];
    %init
%     ind = find(cellfun(@(x)ischar(x)&&strcmpi(x,'ylim'),props.axes),...
%         1,'last');
    if ~isempty(Ylim)
        yl = Ylim;
    else
        mima = [min([mDAT{:}]-[mERR{:}]) , max([mDAT{:}]+[mERR{:}])];
        yl = mima+diff(mima)*[-0.1,0.1];
    end
    x = t*props.timeFac;
    X = [0,0,repmat(P.protocol.beep.duration,1,2)]*props.timeFac;
    Y = [min(yl),max(yl),max(yl),min(yl)];
    Z = zeros(1,4)-1;
    maxRAT = 1;
    
    %plots
    for win = 1:noWIN
        hf = figure('position',posF,props.figure{:}); drawnow
        movegui(hf,'center')
        ha = fig_createAxes(hf,[1,noAXI],props.fig_createAxes{:});
        if withSPN %add spindle axes
            ha(2,:) = NaN;
            for k = 1:noAXI
                pos0 = get(ha(1,k),'position'); uni=get(ha(1,k),'unit');
                dy = pos0(4)/sum(props.axisRatioY)*props.axisRatioY;
                %set ha1
                pos = pos0;
                pos([2,4])=[pos(2)+sum(dy(2:end)),dy(1)];
                
                set(ha(1,k),'position',pos)
                %additional axis
                pos([2,4])=[pos(2)-dy(2),dy(2)];
                ha(2,k) = axes('unit',uni,'position',pos);
            end
        end
        
        %PLOT BEEPS
        hb = NaN(noAXI,1);
        for k = 1:noAXI
            set(hf,'CurrentAxes',ha(1,k));
            hb(k) = fill(X,Y,Z,props.fill{:}); hold on
        end
        
        %PLOT DATA/SPINDLES
        hp = NaN(noTRA,2); %mean signal and spindles
        LegSTR = cell(noTRA,2);
        for tra = 1:noTRA
            [str,col] = ana.Transitions{tra,3:4};
            LegSTR(tra,1) = {sprintf('%s, N = %i',str,mCNT(tra,win))};
            set(hf,'CurrentAxes',ha(1,axesNUM(tra)))
            
            %DATA
            y = mDAT{tra,win};
            s = mERR{tra,win};
            hp(tra,1) = plot(x,y,'color',col,props.plot{:});
            %error
            xs = [t,t(end:-1:1)];
            ys = [y-s,y(end:-1:1)+s(end:-1:1)];
            h = fill(xs,ys,zeros(size(xs)),'edgecolor','none',...
                'facecolor',col,'facealpha',0.2);
            uistack(h,'bottom')
            %spindles
            if withSPN
                tmp = mSPN{tra,win};
                y = NaN(size(tmp));
                y(tmp==true) = yl(1)+0.01*diff(yl)*tra;
                hp(tra,2) = plot(x,y,'.','color',col,props.plot{:});
                LegSTR{tra,2} = sprintf('spindles N = %i',...
                    sum(cSPN{tra,win}));
            end
            %text
            title(sprintf('Transitions within %g s',ana.windows(win)),...
                props.title{:})
            ylabel('',props.ylabel{:});
            
            %SPINLDE RATIO
            if withSPN
                set(hf,'CurrentAxes',ha(2,axesNUM(tra)))
                y = cSPN{tra,win};
                %ratio
                switch lower(ana.ratio.type)
                    case 'count'
                        ny = fs*ana.ratio.winLen;
                        nx = floor(numel(y)/ny);
                        ratio = sum(reshape(y(1:nx*ny),[ny,nx]),1);
                        tRat  = mean(reshape(x(1:nx*ny),[ny,nx]),1);
                    case 'ratio'
                        [ratio,tRat,~,stepLen] = f_eventRatio(~isnan(y),...
                            fs,ana.ratio.winLen,ana.ratio.stepLen);
                        ratio = ratio*ana.ratio.fac; %desired unit
                        tRat = tRat-margin(1); %correct beep at zero
                    otherwise
                        error('Event ratio type ''%s'' not implemented',...
                            ana.ratio.type)
                end
                plot(tRat,ratio,'color',col,props.plot{:}); hold on
                maxRAT = max([ratio(:);maxRAT]);
                %text
                if tra==1
                    xlabel(sprintf('Time [%s]',props.timeUni),...
                        props.xlabel{:})
                    if axesNUM(tra)==1
                        switch lower(ana.ratio.type)
                            case 'count'        
                                ylabel({EEG,'# spindles',...
                                    sprintf('%g-s bins',ana.ratio.winLen)});
                            case 'ratio'
                                ylabel({sprintf('Ratio %s',ana.ratio.unit),...
                                    sprintf('win: %g-s',ana.ratio.winLen),...
                                    sprintf('steps: %g-s',stepLen)})
                        end
                    end
                end
            end
        end %transition loop
        
        %LEGEND
        for k = 1:noAXI
            ind = axesNUM==k;
            set(hf,'CurrentAxes',ha(1,k))
            if withSPN
                a = hp(ind,:); b = LegSTR(ind,:);
                legend(a(:),b(:))
            else
                legend(hp(ind,1),LegSTR(ind,1))
            end
        end
        
        %SETTINGS
        for k = 1:noAXI
            uistack(hb(k),'bottom')
        end
        linkaxes(ha,'x')
        if withSPN
            set(ha(2,:),'ylim',[0,1.2*maxRAT])
        end
        set(ha(1,:),'ylim',yl,'xtick',[])
        set(ha,'unit','normalized','xlim',[x(1),x(end)],props.axes{:})
        
        %SUPER TITLE
        str=cell(1,3);
        if noFIL==1
            [rPath,rFile,rExt]=fileparts(Files.data{1});
            str{1}=regexprep(rPath,{'\','_'},{'\\\\','\\_'});
            str{2}=regexprep([rFile,rExt],{'\','_'},{'\\\\','\\_'});
        else
            str{1}=sprintf('Files N = %i',noFIL);
            str{2}=ana.average.label;
        end
        str{3}=sprintf('%sdB,  %sHz',...
            sprintf('%g ',ana.powers),sprintf('%g ',ana.frequencies));
        h=fig_superLabel(ha(1,:),'title',str);
        set(h,props.superTitle{:})
    end
    toc
else
    fprintf(' [\bexcluded !!!]\b\n')
end

%% TABLE
fprintf('\nTRANSITION COUNT\n')
labelTab = ['Window [s]',num2cell(ana.windows)];
labelFIL = ['File',repmat({'Average'},1,noWIN)];
Tab = cell(noTRA,numel(labelTab),2);
Tab(:,1,1) = ana.Transitions(:,3);
Tab(:,1,2) = cellfun(@(x)sprintf('%s [%%]',x),ana.Transitions(:,3),...
    'uniformoutput',false);
Tab(:,2:end,1) = num2cell(mCNT);
Tab(:,2:end,2) = num2cell(mCNT./mPER*100);
for fil = 1:noFIL
    ind = size(Tab,2)+(1:size(MCNT,2));
    Tab(:,ind,:) = {NaN};
    Tab(:,ind,1) = num2cell(MCNT(:,:,fil));
    Tab(:,ind,2) = num2cell(MCNT(:,:,fil)./MPER(:,:,fil)*100);
    labelFIL(ind) = Files.data(fil);
    labelTab(ind) = num2cell(ana.windows);
end
table = [labelFIL;labelTab;Tab(:,:,1);Tab(:,:,2)]';
% disp(table)
% fprintf('\nTRANSITION COUNT\n')
% labelTab = ['Window [s]',num2cell(ana.windows)];
% Tab = cell(noTRA,numel(labelTab),2);
% Tab(:,1,1) = ana.Transitions(:,3);
% Tab(:,1,2) = cellfun(@(x)sprintf('%s [%%]',x),ana.Transitions(:,3),...
%     'uniformoutput',false);
% Tab(:,2:end,1) = num2cell(mCNT);
% Tab(:,2:end,2) = num2cell(mCNT./mPER*100);
% table = [labelTab;Tab(:,:,1);Tab(:,:,2)]';
% disp(table)


%%


%EXPORT EVENT COUNT
if exist('snameCNT','var') %to get rid for testing
    [sFile,sPath] = uiputfile('*.xls;*xlsx','Save Transition Count',...
        snameCNT);
    if ischar(sFile)
        sname = fullfile(sPath,sFile);
        if exist(sname,'file')==2
            delete(sname)
        end
        xlswrite(sname,table)
        try %when function available
            xls_cellFit(sname)
        catch
        end
    end
end

