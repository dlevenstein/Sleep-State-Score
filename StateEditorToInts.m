function StateIntervals = StateEditorToInts(basepath,basename)
% Takes output from stateeditor and puts it back into StateIntervals.  Also
% puts metadata into apppropriate places as well, including saving any
% manual thresholds
% Brendon Watson July 2016

if ~exist('basepath','var')
    [~,basename] = fileparts(cd);
    basepath = cd;
end

load(fullfile(basepath,[basename '-states.mat']));
I = IDXtoINT(states);

NREMints = I{3};
WAKEints = cat(1,I{1},I{2});%recombine MA and WAKE so they can be re-separated by length... to remain consistent
WAKEints = sortrows(WAKEints,1);
REMints = I{5};

StateIntervals = StatesToFinalScoring(NREMints,WAKEints,REMints);

% gsipath = fullfile(basepath,[basename '_GoodSleepInterval.mat']);
% if exist(gsipath,'file')
%     load(gsipath)
% end


SWLFPpath = fullfile(basepath,[basename '_SWLFP.mat']);
if exist(SWLFPpath,'file')
    load(SWLFPpath,'SWchannum')
    StateIntervals.metadata.SWchannum = SWchannum;
else
    disp('No _SWLFP.mat') %will be useful for reminding to change code after file type conversion
end


THLFPpath = fullfile(basepath,[basename '_ThetaLFP.mat']);
if exist(THLFPpath,'file')
    load(THLFPpath,'THchannum')
    StateIntervals.metadata.THchannum = THchannum;
else
    disp('No _ThetaLFP.mat') 
end
save(fullfile(basepath,[basename '_SleepScoreManualReviewed.mat']),'StateIntervals');


% save manualy-set thresholds (and associated histograms) created in
% StateEditor... if they exist
SSFSEpath = fullfile(basepath,[basename '_SleepScore_FromStateEditor.mat']);
if exist(SSFSEpath,'file')
    load(SSFSEpath,'histsandthreshs')
    ManualHistsAndThreshs = histsandthreshs;
    save(fullfile(basepath,[basename '_StateScoreMetrics.mat']),'ManualHistsAndThreshs','-append');
end

