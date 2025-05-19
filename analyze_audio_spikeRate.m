% analyze_audio_spikeRate
%------------------------------------------------------------------------
% Plots mean spike rates at beeps (bar plot), per sleep stage and
% frequency.

%
% Thomas Rusterholz, Mai 2019
%-------------------------------------------------------------------------
clc; clear; %close all;
scriptPath=fileparts(which(mfilename));
addpath(fullfile(scriptPath,'functions'))
addpath('Z:\OptoLab_v4.1\function\spike'); %f_spikeRate_ISI.m

%PARAMETERS
%----------
%FILES
Files.spikes    = ... opens file browser to select one or more
    ['Z:\1 TIDIS Lab\Ida\Recording\Tetrode Recordings\',...
    'Opto Stimulation\Baseline No Audio Opto\Morning\CMT18_portB\old laser too high\*.txt'];
% Files.protocol  = fullfile(scriptPath,'protocols',{... protocols (cell)
%     'protocol_laser_3600s_80dB_108x3beeps.mat'});
Files.protocol  = {['Z:\1 TIDIS Lab\Ida\Script\',...
'Audio Analysis_PulsePal_3\protocols\',...
 ...'FC_protocol_45s_14beeps.mat',...
 ...'FC_protocol_1h_342beeps.mat',...
 ...'protocol_laser_3600s_80dB_108x3beeps.mat',...
 'protocol_laser_600s_80dB_27x2beeps.mat']};
Files.hypnogram = 'Hypnogram.mat'; %adds read path of spikes

%ANALYSIS OPTIONS
ana.power=80; %beep power [dB], one only
ana.Stages={... {number in hypnogram, label}
    1,'Wake';...    
    2,'NREM';...
    3,'REM';...
    };
%labels protocol (hypnogram sub-path) and stimTimes (data files channels)
ana.labels.protocol='Morning'; %appropriate for Files.protocol
ana.labels.stimTimes='stimTimes'; %'stimTimes' or 'stimTimes2' or ...
%column labels in spike file (correct order!, no label for waveforms!)
ana.labels.spikeData={'Channel','Unit','Timestamp','PC 1','PC 2','PC 3'};
%plot margin: pre beep & post beep [s]
ana.margin=[0.5,0.4];
%RATIO (bar plot or line plot using f_spikeRate_ISI)
ana.rate.select = 'bar'; %'bar' or 'line' plot
ana.rate.barWidth = 0.03; %bar width for bar-plot
ana.rate.winLen = 50; %[ms], sliding window length for line plot
%  PS:  margin(1) AND margin(2) + beep duration must be dividable by
%       winLen. Else an error might occure


%PLOT PROPERTIES
%figure / axes
props.figure={'position',[440,330,880,600]};
props.fig_createAxes={[60,30,10],[50,30,80],'pixel'}; %own function
props.axes={'xtick',-700:100:700}; %y-axes are beeps
%plot
props.fill={'edgecolor','none','facecolor',[1,1,1]*.7}; %beep area
props.bar={'barwidth',1}; %spike rate bar-plot
props.line={'linewidth',2}; %spike rate line-plot
props.time.factor=1000; %time factor to unit seconds
props.time.unit='ms';   %time unit based on time factor
%text
props.xlabel={};
props.ylabel={};
props.title={'fontweight','normal'};
props.superTitle={'fontweight','bold'};
props.text={};


%MAIN PROGRAM
%--------------
str=sprintf('%s',upper(mfilename));
fprintf('%s\n%s\n',str,repmat('-',size(str)));
%props must not be empty
tmp = fieldnames(props);
for k = 1:numel(tmp)
    if isempty(props.(tmp{k}))
        props.(tmp{k}) = {'visible','on'};
    end
end


%LOAD PROTOCOL FILES
ind=cellfun(@(x)exist(x,'file')==2,Files.protocol);
if all(~ind)
    error('No protocol file found')
end
if any(~ind)
    fprintf('[\bProtocols excluded (file NOT found):\n%s]\b',...
        sprintf(' - %s\n',Files.protocol{~ind}))
    Files.protocol(~ind)=[];
end
for k = 1:numel(Files.protocol)
    pp = f_protocols_load(Files.protocol{k});
    PP(k) = pp;
end

%GET SPIKE FILES
[rFiles,rPath]=uigetfile('.txt','Select Spike Files',Files.spikes,...
    'Multiselect','on');
if ~ischar(rPath)
    fprintf('[\bNo File Selected!]\b\n')
    return
elseif ischar(rFiles)
    rFiles={rFiles};
end

%number of ...
noSTA=size(ana.Stages,1);
noFIL=numel(rFiles);

%LOAD DATA (ID, for each spike-file)
fprintf('LOAD DATA\n')
cntFIL=0; %init, count files with valid data
nnFIL = numel(num2str(noFIL));
indent = blanks(2*nnFIL+2);
DATA = cell(2,noFIL);
for fil=1:noFIL
    %FILES
    rFile=rFiles{fil};
    ID.fileSPK  = fullfile(rPath,rFile);
    ID.fileHYP  = fullfile(rPath,Files.hypnogram);
    ID.filePTC  = ''; %init, find later
    ID.str.protocol = ana.labels.protocol; %label    
    %info from read file string
    tmp=regexp(rFile,'_','split');
    [ID.str.day, ID.str.unknown, ID.str.port, ID.str.location] = tmp{1:4};
    ID.channels = cellfun(@(x)sprintf('%s-%03s',tmp{5}(1:5),x),...
        regexp(tmp{5},'[0-9]*','match'),'UniformOutput',false);
    
%     %FIND HYPNOGRAM FILE
%     tmp=regexprep(Files.hypnogram,strREP.str,eval(strREP.rep));   
%     file=dir(tmp);
%     switch numel(file)
%         case 0
%             error('Hypnogram file NOT found for:\n%s\n]\b',tmp)
%         case 1
%             fileHYP=fullfile(file.folder,file.name);
%         otherwise
%             for k=1:numel(file)
%                 list{k}=fullfile(file(k).folder,file(k).name);
%             end
%             ind=listdlg('PromptString',...
%                 {sprintf('Several files found for ''%s''',rFile),...
%                 'Select one'},...
%                 'SelectionMode','single',...
%                 'ListSize', [600 80*numel(file)],...
%                 'Name','Select a Hypnogram',...
%                 'ListString',list);
%             if isempty(ind)
%                 error('You must select a Hypnogram file for %s\n',...
%                     rFile)
%             end
%             fileHYP=list{ind};
%     end
%     ID.fileHYP=fileHYP;
    
    %FIND PROTOCOL
    tmp=load(fullfile(fileparts(ID.fileHYP),[ID.channels{1},'.mat']));
    noSAM=numel(tmp.resampled_data_mV);
    ID.fs=tmp.SampRate;
    tmp = f_protocol_check(PP,tmp.(ana.labels.stimTimes),ID.fs);
    ind = isnan(tmp.frequency); %remove NaNs
    if any(ind)
        tmp.stimSTA(ind) = [];
        tmp.stimEND(ind) = [];
        tmp.frequency(ind) = [];
        tmp.power(ind) = [];
        tmp.repetition(ind) = [];
    end
    ID.prot = tmp;
    
    if isempty(ID.prot)
         fprintf('  [\bOMITTED! No appropriate protocol found]\b\n\n')
        continue %excluded
    end
    ID.filePTC=ID.prot.protocol.file;
    
    %PRINT OUT
    fprintf('%*i/%i: Files\n',nnFIL,fil,noFIL)
    fprintf('%s protocol  : %s\n',indent,ID.filePTC)
    fprintf('%s hypnogram : %s\n',indent,ID.fileHYP)
    fprintf('%s spikes    : %s\n',indent,ID.fileSPK)

    %LOAD HYPNOGRAM
    load(ID.fileHYP,'Hypnogram')
    %upsamling or not
    fac=round(noSAM/numel(Hypnogram));
    if fac>1
        fprintf('  [\bHypnogram - upsampled by factor %i!]\b\n',...
            fac)
        Hypnogram=repmat(Hypnogram',fac,1);
        Hypnogram=Hypnogram(:);
    end
    %fill with NaNs
    if numel(Hypnogram)<noSAM
        Hypnogram(numel(Hypnogram)+1:noSAM)=NaN;
    elseif numel(Hypnogram)>noSAM
        error('uups')
    end
    ID.hypnogram=Hypnogram;
    
    %LOAD SPIKES
    try
        ID.spikes=f_readSpikes(ID.fileSPK,ana.labels.spikeData);
    catch
        fprintf('%s [\bNot a regular spike-file]\b\n',indent)
        continue
    end
    
    %apply data
    cntFIL=cntFIL+1;
    DATA(cntFIL)={ID};
end
DATA(cellfun(@isempty,DATA)) = [];
noFIL=numel(DATA);
if noFIL==0
   fprintf(2,'No Data Loaded!\n') 
end

%% FILE LOOP
fprintf('\nPLOTS, FILES:\n')
nnFIL=numel(num2str(noFIL));
for fil=1:noFIL
    ID=DATA{fil};
    [rPath,fileSPK]=fileparts(ID.fileSPK); 
    pathDAT=fileparts(ID.fileHYP);
    fprintf('- %s\n',fileSPK);
    
    %SPECIAL DATA TO WORK WIDTH (often used)
    Hypnogram=(ID.hypnogram(:));
    fs=ID.fs;
    beepDUR=ID.prot.protocol.beep.duration;
    %init time
    tLim=[-abs(ana.margin(1)),abs(ana.margin(2))+beepDUR]; %time limit
    switch ana.rate.select
        case 'bar'
            tRat=tLim(1):ana.rate.barWidth:tLim(2); %time boarders spike rate
            tBar=tRat(1:end-1)+ana.rate.barWidth/2; %time bar midpoints
        case 'line'
            tLin = linspace(tLim(1),tLim(2),...
                diff(tLim)*ID.fs/ana.rate.winLen*2+1);
        otherwise
            error('Selected ratio ''%s'' not defined',ana.rate.select)
    end
    %spike data
    units=unique(ID.spikes.unit);
    units(units==0)=[]; %artifacted spikes
    noUNI=numel(units);
    %protocol data
    frequencies=unique(ID.prot.frequency);
    noFRQ=numel(frequencies);
  
    %LOOPS
    %unit loop
    for uni=1:noUNI
        unit=units(uni);
        clear spikes
        ind=ID.spikes.unit==unit;
        spikes.time=ID.spikes.timestamp(ind);
        spikes.index=round(spikes.time*fs);
        spikes.hypnogram=Hypnogram(spikes.index);
        if strcmp(ana.rate.select,'line')
            [sRate,t_sRate,sVec] = f_spikeRate_ISI(...
                spikes.time*1000,ana.rate.winLen);
        end
        
        %FIGURE
        hf=figure(props.figure{:}); drawnow
        fig_size(hf);
        ha=fig_createAxes(hf,[noFRQ,noSTA],props.fig_createAxes{:});
        set(ha,'unit','normalized')
        %init
        maxRate=-inf; %max spike rate
        ht=[]; %text handle
       
        
        
        %freq loop
        for frq=1:noFRQ
            frequency=frequencies(frq);
            ind=ID.prot.frequency==frequency & ...
                ID.prot.power==ana.power; %one power only!
            clear stims
            stims.index=ID.prot.stimSTA(ind); %start index
            stims.time=stims.index/fs;
            stims.hypnogram=Hypnogram(stims.index);
            noBEP=sum(ind);
            if noBEP~=ID.prot.protocol.prot.beepsPerSignal
                error('uups')
            end
            
            %stage loop
            for sta=1:noSTA
                [stageNUM,stage]=ana.Stages{sta,:};
                indHYP = find(stims.hypnogram==stageNUM);
                
                %SPIKE RATE / PLOT I
                set(hf,'CurrentAxes',ha(frq,sta))
                N=numel(indHYP);
                switch ana.rate.select
                    case 'bar'
                        SpikeRate=zeros(N,numel(tBar));
                        for k=1:N
                            tSTA=stims.time(indHYP(k));
                            for q=1:numel(tBar)
                                indSPK = spikes.hypnogram==stageNUM & ...
                                    spikes.time >= tSTA+tRat(q) & ...
                                    spikes.time <= tSTA+tRat(q+1);
                                SpikeRate(k,q)=sum(indSPK)/...
                                    ana.rate.barWidth; %[Hz]
                            end
                        end
                        spikeRate=mean(SpikeRate,1);
                        maxRate=max([maxRate;spikeRate(:)]);
                        %plot
                        bar(tBar*props.time.factor,spikeRate,...
                            props.bar{:});
                    case 'line'
                        SpikeRate=NaN(N,numel(tLin));
                        for k=1:N
                            tSTA = stims.time(indHYP(k));
                            tmp = (stims.time(indHYP(k))+tLim)*1000; %[ms]
                            ind1 = round(tmp(1)/ana.rate.winLen*2)+1;
                            ind2 = round(tmp(2)/ana.rate.winLen*2)+1;
                            ind = ind1:ind2;
                            if any(ind<1)
                                ind(ind<1) = [];
                                SpikeRate(k,end-numel(ind)+1:numel(ind))...
                                    = sRate(ind);
                            elseif any(ind)>numel(sRate)
                                ind(ind>numel(sRate)) = [];
                                SpikeRate(k,1:numel(ind)) = sRate(ind);
                            else
                                SpikeRate(k,:) = sRate(ind);
                            end
                        end
                        spikeRate=nanmean(SpikeRate,1);
                        maxRate=max([maxRate;spikeRate(:)]);
                        %plot
                        plot(tLin*props.time.factor,spikeRate,...
                            props.line{:});
                end
                
                %PLOT II
                %settings
                hold on
                set(gca,'xlim',tLim*props.time.factor)
                %text
                if frq==1
                    title(stage,props.title{:})
                end
                if frq==noFRQ
                    xlabel(sprintf('Time [%s]',props.time.unit),...
                        props.xlabel{:})
                end
                if sta==1
                    ylabel(sprintf('%g Hz',frequency),props.ylabel{:})
                end
                ht(end+1)=text(beepDUR*props.time.factor/2,max(ylim),...
                    sprintf('N = %i',N),...
                    'horizontalalignment','center',...
                    'verticalalignment','top',props.text{:});
                drawnow
            end %stage loop
        end %freq loop
        h=fig_superLabel(ha,'title',{...
            sprintf('Mean Spike Rate [Hz]'),...
            strrep(fileSPK,'_','\_'),...
            sprintf('Unit %i, %g dB',unit,ana.power)});
        set(h,props.superTitle{:})
        linkaxes(ha,'xy')
        m=1.2*maxRate;
        for k=1:numel(ha)
            set(hf,'CurrentAxes',ha(k))
            %beep area
            h=fill([0,0,beepDUR,beepDUR]*props.time.factor,...
                [0,m,m,0],zeros(1,4),props.fill{:});
            uistack(h,'bottom')
            %correct text position
            pos=get(ht(k),'position');
            pos(2)=m;
            set(ht(k),'position',pos)
            %axes settings
            set(ha(k),'ylim',[0,m],props.axes{:})
        end
    end %unit loop
end %file loop