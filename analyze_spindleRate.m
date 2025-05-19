% analyze_spindleRate
%------------------------------------------------------------------------
% Calculates the spindle rate in NREM sleep
% Stage based on spindle center


clc; clear; close all;
[scriptPath,scriptName] = fileparts(which(mfilename));
addpath(fullfile(scriptPath,'functions'))

%PARAMETERS
%----------
%PATHS AND FILES
%read paths (get paths from file-list)
[tmp,fileList] = import_fileList(... reads files!
    ['Z:\1 TIDIS Lab\Ida\Script\Audio Analysis_PulsePal_3\filelists\',...
    '\*.txt']);
rPaths = cellfun(@fileparts,tmp,'uniformoutput',false); %paths only
rPaths = natsort(unique(rPaths)); %unique & natsort
%spindle files (without path)
% Ex.'spindle_EEG1.mat','spindle_EEG2.mat','spindle_amp-B-024.mat'
Files.spind = {'spindle_amp-A-028.mat','spindle_amp-B-028.mat','spindle_amp-C-028.mat','spindle_amp-D-028.mat','spindle_amp-A-000.mat','spindle_amp-B-000.mat','spindle_amp-C-000.mat','spindle_amp-D-000.mat','spindle_amp-A-004.mat','spindle_amp-B-004.mat','spindle_amp-C-004.mat','spindle_amp-D-004.mat'}; %one or more
Files.hypno = 'Hypnogram.mat';
fsHYP = 1000; %sampling rate of saved hypnogram

%TIME UNIT (and corresponding factor to seconds)
%   for tUNI = 's', tFAC = 1
%   for tUNI = 'min, tFAC = 1/60
%   ...
tUNI = 'min';
tFAC = 1/60;

%NREM stage number in hypnogram (may also be a vector)
stageNREM = 2;

%SAVE NAME (excel-file, asks again for re-confirm)C:\Users\i0325207\Desktop\spindleRate.xlsx
saveName = 'Z:\1 TIDIS Lab\Ida\Analysis\Sleep Analysis\Spindles\Opto Row Files\*.xlsx';
saveName = 'spindleRate_EEG.xlsx';
% 'spindleRate_Tetrode.xlsx' or 'spindleRate_EEG.xlsx'

%MAIN PROGRAM
%--------------
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)));

%CHECK
if ischar(rPaths) %must be cell
    rPaths = {rPaths};
end
if ischar(Files.spind) %must be cell
    rPaths = {rPaths};
    Files.spind = {Files.spind};
end
if iscell(Files.hypno) %must be char
    Files.hypno = Files.hypno{1};
end

%NUMBER OF ...
noPAT = numel(rPaths);
noFIL = numel(Files.spind);

%INIT TABLE
labels = ['Path',Files.spind(:)','Unit'];
tab = cell(noPAT,numel(labels));
tab{1,end} = sprintf('[1/%s]',tUNI);

%PATH LOOP
nnPAT = numel(num2str(noPAT));
indent = blanks(2*nnPAT+2);
for pat = 1:noPAT
    rPath = rPaths{pat};
    tab(pat,1) = {rPath};
    fprintf('%*i/%i: %s\n',nnPAT,pat,noPAT,rPath)
    
    %LOAD HYPNOGRAM
    rFile = Files.hypno;
    fname = fullfile(rPath,rFile);
    fprintf('%s Read %s\n',indent,rFile);
    if exist(fname,'file')~=2
        fprintf('\b - [\bNOT found]\b\n')
        continue
    end
    tmp = load(fname);
    Hypnogram = tmp.Hypnogram(:);
    %duration
    durTOT  = numel(Hypnogram)/fsHYP*tFAC;
    durNREM = sum(ismember(Hypnogram,stageNREM))/fsHYP*tFAC;
    fprintf('%s   duration total [%s] : %g\n',indent,tUNI,durTOT)
    fprintf('%s   duration NREM [%s]  : %g\n',indent,tUNI,durNREM)
    
    %FILE LOOP
    for fil = 1:noFIL
        rFile = Files.spind{fil};
        fprintf('%s Read %s\n',indent,rFile)
        
        %LOAD DATA
        fname = fullfile(rPath,rFile);
        if exist(fname,'file')~=2
            fprintf('\b - [\bNOT found]\b\n')
            continue
        end
        clear SampRate spindle
        load(fname,'SampRate','spindle');
        
        %SPINDLES
        spnIND = ([spindle.str]+[spindle.end])/2; %spindle center index
        hypnogram = Hypnogram;
        if fsHYP>SampRate %index for hypnogram samples
            spnIND = spnIND/SampRate*fsHYP; %index in fsHYP
        elseif SampRate>fsHYP %upsample hypnogram
            hypnogram = repmat(hypnogram',SampRate/fsHYP,1);
            hypnogram = hypnogram(:);
        end
        spnIND = round(spnIND);
        nTOT = numel(spnIND);
        %remove spindles
        spnIND(spnIND>numel(hypnogram)) = []; %beyond scoring
        spnIND(~ismember(hypnogram(spnIND),stageNREM)) = [];
        noSPN = numel(spnIND);
        spnRate = noSPN/durNREM;
        %append & print out
        tab(pat,1+fil) = {spnRate};
        fprintf('%s   total  removed  in_NREM \n',indent)
        fprintf('%s   %5i  %7i  %7i \n',indent,nTOT,nTOT-noSPN,noSPN)
        fprintf('%s   spindle rate [1/%s] : %f\n',indent,tUNI,spnRate)
        
%         % test plot
%         close all
%         figure
%         t = (1:numel(hypnogram))/fsHYP*tFAC;
%         plot(t,hypnogram); hold on;
%         plot(t(spnIND),hypnogram(spnIND),'rx')
%         set(gca,'xlim',[0,t(end)],'ylim',[0.5,3.5])
%         xlabel(sprintf('[%s]',tUNI))
%         return
        
    end
end

%SAVE
[sFile,sPath] = uiputfile('*.xls;*.xlsx','Save Tabel As',saveName);
if ~ischar(sFile)
    fprintf(2,'Data NOT saved!\n')
    return
end
sname = fullfile(sPath,sFile);
xlswrite(sname,[labels;tab]);
try %if function is available
    xls_cellFit(sname)
catch
end
fprintf('Saved: %s\n',sname)