% analyze_audio_spike
%------------------------------------------------------------------------
% Plots spikes before, while after beeps per frequency and sleep stage.
% One figure per selected spike file and spike unit (all but unit zero).
%
% So far, only one beep power selectable and only one protocol selectable.
% Script will check whether selected protocol is appropriate (by comparing
% defined beep timing with recorded stimuli)
%
% PS: hypnogram will not be upsampled! For actual data, it's already saved
%     upsampled. If not, one has to change the script!

%Files.spikes = ... opens file browser to select one or more
%     ['Z:\1 TIDIS Lab\Ida\Recording\Tetrode Recordings\',...
%     'Audio Stimulation\1d\*.txt'];
% Thomas Rusterholz, Mai 2019
%-------------------------------------------------------------------------
clc; clear; close all;
scriptPath = fileparts(which(mfilename));
addpath(fullfile(scriptPath,'functions'))

%PARAMETERS
%----------
%FILES
Files.spikes = ... opens file browser to select one or more
    ['Z:\1 TIDIS Lab\Ida\Recording\Tetrode Recordings\',...
    'Opto Stimulation\1d\Afternoon\CMT19_portC\*SUA.txt'];
Files.protocol = fullfile(scriptPath,'protocols',{... protocols (cell)
 ...'FC_protocol_45s_14beeps.mat',...
 ...'FC_protocol_1h_342beeps.mat',...
 ...'protocol_laser_3600s_80dB_108x3beeps.mat',...
 ...'protocol_laser_600s_80dB_27x2beeps.mat'...
 ...'protocol_morning.mat'...
 ...'protocol_laser_3600s_80dB_108x3beeps.mat'...
 ...'protocol_laser_3600s_80dB_108x3beeps_2.mat'...
 'protocol_laser_600s_80dB_27x2beeps.mat'});
Files.hypnogram = 'Hypnogram.mat'; %adds read path of spikes
%string replacments for Files.hypnogram, splits spike filename at '_'
% E.g  4d_A_portB_Au1_amp-B-0-1-2-3_Sorted.mat
%      DAY_UNKNOWN_PORT_LOCATION_CHANNELS*
% strREP.str={'ID-DAYTIME','ID-PORT','ID-DAY'}; %to replace
% strREP.rep='{ID.str.daytime,ID.str.port,ID.str.day}'; %replace with this, ID.*


%ANALYSIS OPTIONS
ana.power  = 80; %beep power [dB], one only
ana.Stages = {... {number in hypnogram, label}
    1,'Wake';...
    2,'NREM';...
    3,'REM';...
    };
%labels protocol (hypnogram sub-path) and stimTimes (data files channels)
ana.labels.protocol  = 'Morning'; %appropriate for Files.protocol
ana.labels.stimTimes = 'stimTimes'; %'stimTimes' or 'stimTimes2' or ...
%column labels in spike file (correct order!, no label for waveforms!)
ana.labels.spikeData = {'Channel','Unit','Timestamp','PC 1','PC 2','PC 3'};
%accuracy protocol to stimTimes, sum(|prot-stim|)/sum(prot)
ana.accuracy = 10^-3; %PS: zero is an utopic value :-)
%plot margin: pre beep & post beep
ana.margin = [0.7,0.7]; %[s]


%PLOT PROPERTIES (must not be empty!)
%figure / axes
props.figure = {'position',[440,330,880,600]};
props.fig_createAxes = {[60,30,10],[50,30,80],'pixel'}; %own function
props.axes = {'xtick',-700:100:700}; %y-axes are beeps, 'ytick',0:10:40
%plot
props.fill = {'edgecolor','none','facecolor',[1,1,1]*.7}; %beep area
props.line = {'color','k','linewidth',1}; %spike lines
props.time.factor = 1000; %time factor to unit seconds
props.time.unit   = 'ms';   %time unit based on time factor
%text
props.xlabel = {'visible','on'};
props.ylabel = {'visible','on'};
props.title  = {'fontweight','normal'};
props.superTitle = {'fontweight','bold'};
props.text   = {'visible','on'};


%MAIN PROGRAM
%--------------
str = sprintf('%s',upper(mfilename));
fprintf('%s\n%s\n',str,repmat('-',size(str)));

%LOAD PROTOCOL FILES
ind = cellfun(@(x)exist(x,'file')==2,Files.protocol);
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
[rFiles,rPath] = uigetfile('.txt','Select Spike Files',Files.spikes,...
    'Multiselect','on');
if ~ischar(rPath)
    fprintf('[\bNo File Selected!]\b\n')
    return
elseif ischar(rFiles)
    rFiles = {rFiles};
end

%number of ...
noSTA = size(ana.Stages,1);
noFIL = numel(rFiles);


%LOAD DATA (ID, for each spike-file)
fprintf('LOAD DATA\n')
cntFIL = 0; %init, count files with valid data
for fil = 1:noFIL
    %FILES
    rFile = rFiles{fil};
    ID.fileSPK  = fullfile(rPath,rFile);
    tmp = dir(fullfile(rPath,Files.hypnogram));
    ID.fileHYP = fullfile(tmp.folder,tmp.name);
    ID.filePTC = ''; %init, find later
    ID.str.protocol = ana.labels.protocol; %label
    fprintf(' File: %s\n',ID.fileSPK)
    %info from read file string
    tmp = regexp(rFile,'_','split');
    [ID.str.day, ID.str.unknown, ID.str.port, ID.str.location] = tmp{1:4};
    ID.channels = cellfun(@(x)sprintf('%s-%03s',tmp{5}(1:5),x),...
        regexp(tmp{5},'[0-9]*','match'),'UniformOutput',false);
    switch lower(ID.str.unknown)
        case 'a'
            ID.str.daytime = 'Afternoon';
        case 'm'
            ID.str.daytime = 'Morning';
        otherwise
            error('Fatal error')
    end
     
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
    tmp   = load(fullfile(fileparts(ID.fileHYP),[ID.channels{1},'.mat']));
    noSAM = numel(tmp.resampled_data_mV);
    ID.fs = tmp.SampRate;
    ID.prot = f_protocol_check(PP,tmp.(ana.labels.stimTimes),ID.fs);
    if isempty(ID.prot)
         fprintf('  [\bOMITTED! No appropriate protocol found]\b\n\n')
        continue %excluded
    end
    ID.filePTC = ID.prot.protocol.file;
    
    %LOAD HYPNOGRAM
    load(ID.fileHYP,'Hypnogram')
    %upsamling or not
    fac = round(noSAM/numel(Hypnogram));
    if fac>1
        fprintf('  [\bHypnogram - upsampled by factor %i!]\b\n',...
            fac)
        Hypnogram = repmat(Hypnogram',fac,1);
        Hypnogram = Hypnogram(:);
    end
    %fill with NaNs
    if numel(Hypnogram)<noSAM
        Hypnogram(numel(Hypnogram)+1:noSAM) = NaN;
    elseif numel(Hypnogram)>noSAM
        error('uups')
    end
    ID.hypnogram = Hypnogram;
    
    %LOAD SPIKES
    ID.spikes = f_readSpikes(ID.fileSPK,ana.labels.spikeData);
    
    %apply data
    cntFIL = cntFIL+1;
    DATA(cntFIL) = ID;
end
fprintf('\n')

%FILE LOOP
noFIL = numel(DATA);
nnFIL = numel(num2str(noFIL));
for fil = 1:noFIL
    ID = DATA(fil);
    [rPath,fileSPK] = fileparts(ID.fileSPK); 
    pathDAT = fileparts(ID.fileHYP);
    fprintf('PLOTS, FILES:\n')
    fprintf(' spikes    : %s\n',fileSPK);
    fprintf(' hypnogram : %s\n',ID.fileHYP);
    fprintf(' protocol  : %s\n',ID.filePTC);
    
    %SPECIAL DATA TO WORK WIDTH (often used)
    Hypnogram = (ID.hypnogram(:));
    fs = ID.fs;
    beepDUR = ID.prot.protocol.beep.duration;
    t = (-ana.margin(1)*fs:ana.margin(2)*fs + beepDUR*fs)/fs; %plot time
    %spike data
    units = unique(ID.spikes.unit);
    units(units==0) = []; %artifacted spikes
    noUNI = numel(units);
    %protocol data
    frequencies = unique(ID.prot.frequency);
    noFRQ = numel(frequencies);
  
    %LOOPS
    %unit loop
    for uni = 1:noUNI
        unit = units(uni);
        clear spikes
        ind = ID.spikes.unit==unit;
        spikes.time  = ID.spikes.timestamp(ind);
        spikes.index = round(spikes.time*fs);
        spikes.hypnogram = Hypnogram(spikes.index);
        
        %FIGURE
        hf = figure(props.figure{:}); drawnow
        fig_size(hf);
        ha = fig_createAxes(hf,[noFRQ,noSTA],props.fig_createAxes{:});
        set(ha,'unit','normalized')
        
        %freq loop
        for frq = 1:noFRQ
            frequency = frequencies(frq);
            ind = find(...
                ID.prot.frequency==frequency & ...
                ID.prot.power==ana.power); %one power only!
            clear stims
            stims.index = ID.prot.stimSTA(ind); %start index
            stims.time  = stims.index/fs;
            stims.hypnogram = Hypnogram(stims.index);
            noBEP = numel(stims.index);
            
            %stage loop
            for sta = 1:noSTA
                [stageNUM,stage] = ana.Stages{sta,:};
                indHYP = find(stims.hypnogram==stageNUM);
                
                %PLOT DATA
                set(hf,'CurrentAxes',ha(frq,sta))
                %beep area                
                fill([0,0,beepDUR,beepDUR]*props.time.factor,...
                    [0,noBEP+.5,noBEP+.5,0],zeros(1,4),...
                    props.fill{:})
                hold on
                %spikes
                N = numel(indHYP);
                for k = 1:N
                    tSTA   = stims.time(indHYP(k));
                    indSPK = spikes.hypnogram==stageNUM & ...
                        spikes.time >= tSTA+t(1) & ...
                        spikes.time <= tSTA+t(end);
                    if sum(indSPK)==0
                        continue
                    end
                    tSPK = (spikes.time(indSPK)-tSTA)*props.time.factor;
                    x = repmat((tSPK(:))',2,1);
                    y = repmat([-.5;.5]+k,1,numel(tSPK));
                    line(x,y,props.line{:})
                    drawnow
                end
                %settings
                set(gca,'xlim',[t(1),t(end)]*props.time.factor,...
                    'ylim',[0,noBEP+.5],props.axes{:})
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
                text(beepDUR*props.time.factor/2,max(ylim),...
                    sprintf('N = %i',N),...
                    'horizontalalignment','center',...
                    'verticalalignment','top',props.text{:})
                drawnow
            end %stage loop
        end %freq loop
        h = fig_superLabel(ha,'title',{'Raster Plot',...
            strrep(fileSPK,'_','\_'),...
            sprintf('Unit %i, %g dB',unit,ana.power)});
        set(h,props.superTitle{:})
                
    end %unit loop
end %file loop