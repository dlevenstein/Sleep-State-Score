function SleepScoreMaster(datasetfolder,recordingname,varargin)
%SleepScoreMaster(datasetfolder,recordingname)
%This is the master function for sleep state scoring.
%
%INPUT (optional)   If no inputs included, select folder containing .eeg
%                   and .xml file in prompt.
%                   
%   datasetfolder   Top level folder in which the dataset resides. 
%                   For example:
%                   '/Users/dlevenstein/Dropbox/Research/Datasets/BWData/'
%                   -if not included, prompt comes up to
%                   navigate to the folder holding your recording.
%   recordingname   (optional)
%                   Name of the recording, this will be the name of the
%                   folder in which the .lfp file and other files reside.
%                   For example, the .lfp file should be:
%                   'datasetfolder/recordingname/recordingname.lfp'
%                   ... it is also assumed that this serves as the basename
%                   for the files for instance data will be at
%                   /datasetfolder/recordingname/recordingname.lfp
%   'savedir'       Default: datasetfolder
%   'overwrite'     Default: false
%   'savebool'      Default: true
%   'spindledelta'  Default: false (spindle detector not ready yet...)
%   'scoretime'     Default: [0 Inf]
%   'badchannels'   file datasetfolder/recordingname/'bad_channels.txt'
%                   that lists channels will omit certain channels from EMG
%                   detection and LFP selection
%
%OUTPUT
%   StateIntervals  structure containing start/end times (seconds) of
%                   NREM, REM, WAKE states and episodes. states is the 
%                   "raw" state scoring. episodes are joined episodes of 
%                   extended (40s) time in a given states, allowing for 
%                   brief interruptions. also contains NREM packets, 
%                   unitary epochs of NREM as described in Watson et al 2016.
%                   saved in a .mat file:
%                   recordingname_SleepScore.mat 
%   
%
% DLevenstein and BWatson 2015/16
%% Recording Selection
%if recname is 'select' or something
%use uigetfile to pick and get list of filenames
%if recname is 'all', get all recordings in a folder and
%then run SleepScoreMaster on each of the filenames'
%if no input arguements... select uigetfile

%Select from no input
if ~exist('datasetfolder','var')
    DIRECTORYNAME = uigetdir('',...
        'Which recording(s) would you like to state score?');
    if isequal(DIRECTORYNAME,0);return;end  
    [datasetfolder,recordingname] = fileparts(DIRECTORYNAME); 
end
  
%Select from dataset folder
switch recordingname
    case 'select'
        foldercontents = dir(datasetfolder);
        possiblerecordingnames = {foldercontents([foldercontents.isdir]==1).name};
        [s,v] = listdlg('PromptString','Which recording(s) would you like to state score?',...
                        'ListString',possiblerecordingnames);
        recordingname = possiblerecordingnames(s);
end

%If multiple recordings, loop
numrecs = length(recordingname);
if numrecs > 1 & iscell(recordingname)
    display(['Multiple Recordings (',num2str(numrecs),')'])
    for rr = 1:numrecs
        SleepScoreMaster(datasetfolder,recordingname{rr},varargin{:})
        close all
    end
    return
elseif numrecs == 1 & iscell(recordingname)
        recordingname = recordingname{1};
end

display(['Scoring Recording: ',recordingname]);

%% Deal with input options from varargin
%none yet, but will do this with inputParser when we have input options

%possible variable input options
%'timewin' - only state score a subset of the recording
%'HPCsites' - site indices for HPC probes - will only check these for theta
%           if applicable
%'figloc' - secondardy folder to save figures to
%'spikegroups' - if not in the .xml file
%'SWChannel', 'ThetaChannel' - can enter manually instead of determining
%                               algorithmically

%% inputParse for Optional Inputs and Defaults
p = inputParser;

defaultOverwrite = false;    %Pick new and Overwrite existing ThLFP, SWLFP?
defaultSavebool = true;    %Save Stuff (EMG, LFP)
defaultSpindledelta = false; %Detect spindles/delta?

defaultSavedir = datasetfolder;

defaultScoretime = [0 Inf];

addParameter(p,'overwrite',defaultOverwrite,@islogical)
addParameter(p,'savebool',defaultSavebool,@islogical)
addParameter(p,'spindledelta',defaultSpindledelta,@islogical)
addParameter(p,'savedir',defaultSavedir)
addParameter(p,'scoretime',defaultScoretime)


parse(p,varargin{:})
%Clean up this junk...
overwrite = p.Results.overwrite; 
savebool = p.Results.savebool;
spindledelta = p.Results.spindledelta;
savedir = p.Results.savedir;
scoretime = p.Results.scoretime;

%% Database File Management 
savefolder = fullfile(savedir,recordingname);

if ~exist(savefolder,'dir')
    mkdir(savefolder)
end

%Figure locations
figloc = [fullfile(savefolder,'StateScoreFigures'),'/'];
if ~exist(figloc,'dir')
    mkdir(figloc)
end

%Filenames for EMG, thLFP, and swLFP .mat files in the database.
EMGpath = fullfile(savefolder,[recordingname '_EMGCorr.mat']);
thetalfppath = fullfile(savefolder,[recordingname,'_ThetaLFP.mat']);
swlfppath = fullfile(savefolder,[recordingname,'_SWLFP.mat']);
%Filenames for State and Event .mat files.
sleepstatepath = fullfile(savefolder,[recordingname,'_SleepScore.mat']);
sleepeventpath = fullfile(savefolder,[recordingname,'_SleepEvents.mat']);
spindlestatspath = fullfile(savefolder,[recordingname,'_SpindleStats.mat']);
%Filenames for StateCluster Metrics (broadband/theta)
scoremetricspath = fullfile(savefolder,[recordingname,'_StateScoreMetrics.mat']);

%Filename for .lfp file
if exist (fullfile(datasetfolder,recordingname,[recordingname,'.lfp']),'file')
    rawlfppath = fullfile(datasetfolder,recordingname,[recordingname,'.lfp']);
elseif exist (fullfile(datasetfolder,recordingname,[recordingname,'.eeg']),'file')
    rawlfppath = fullfile(datasetfolder,recordingname,[recordingname,'.eeg']);
elseif ~overwrite
    display('No .lfp file... but using saved files so maybe it''s ok!')
else
    display('No .lfp file')
    return
end

%% CALCULATE EMG FROM HIGH-FREQUENCY COHERENCE
% Load/Calculate EMG based on cross-shank correlations 
% (high frequency correlation signal = high EMG).  
% Schomburg E.W. Neuron 84, 470?485. 2014)
% Do this before lfp load because finding proper lfps will depend on this.
% If EMG is already calculated and in it's own .mat, then load, otherwise
% calculate this

if ~exist(EMGpath,'file') || overwrite;
    display('Calculating EMG')
    EMGCorr = EMGCorrForSleepscore(rawlfppath,scoretime);%BW modify this to have different dependencies, currently assumes presence of: 
    % eeg filename - ok
    % .xml filename - ok
    %     Save ..._EMGCorr file
    if savebool
        save(EMGpath,'EMGCorr')
    end

else
    display('EMG aleady calculated: Loading...')
    load(EMGpath,'EMGCorr')
end
EMG = EMGCorr(:,2);
clear EMGCorr


%% DETERMINE BEST SLOW WAVE AND THETA CHANNELS

sf_LFP = 1250;
sf_EMG = 2;

if ~exist(thetalfppath,'file') && ~exist(swlfppath,'file') || overwrite; % if no lfp file already, load lfp and make lfp file?

    display('Picking SW and TH Channels')
    [SWchannum,THchannum,swLFP,thLFP,t_LFP] = PickSWTHChannel(datasetfolder,recordingname,figloc,scoretime);
    
    if savebool
        %Transfer this into scoremetricspath, predownsampled to what it
        %needs to be for ClusterStates.
        save(swlfppath,'swLFP','SWchannum','t_LFP','sf_LFP');
        save(thetalfppath,'thLFP','THchannum','t_LFP','sf_LFP');
    end
else
    display('SW and TH Channels Already Extracted, Loading...')
    load(swlfppath,'swLFP','SWchannum')
    load(thetalfppath,'thLFP','THchannum')
end



%% CLUSTER STATES BASED ON SLOW WAVE, THETA, EMG

display('Clustering States Based on EMG, SW, and TH LFP channels')
[stateintervals,~,~,~,~,broadbandSlowWave,thratio,EMG,t_clus,badtimes,reclength] = ClusterStates(swLFP,thLFP,EMG,sf_LFP,sf_EMG,figloc,recordingname);

    if savebool
        %Should save (downsampled to what's used in clusterstates...)
        %sw/thLFP in scoremetricspath here!
        save(scoremetricspath,...
            'broadbandSlowWave','thratio','EMG','t_clus',...
            'SWchannum','THchannum','badtimes','reclength');
    end

%% JOIN STATES INTO EPISODES

NREMints = stateintervals{2};
REMints = stateintervals{3};
WAKEints = stateintervals{1};

StateIntervals = StatesToFinalScoring(NREMints,WAKEints,REMints);
StateIntervals.metadata.SWchannum = SWchannum;
StateIntervals.metadata.THchannum = THchannum;

save(sleepstatepath,'StateIntervals');


%% Find Slow Waves and Spindle Times
if spindledelta

    [ pSpindleInts,cycletimemap,deltapeaks,SpindleStats ] = FindSpindlesAndSWs(datasetfolder,recordingname,figloc,StateIntervals);

    SleepEvents.Spindles = pSpindleInts;
    SleepEvents.DeltaPeaks = deltapeaks;
    SleepEvents.SpindleCycleTime = cycletimemap;


    save(sleepeventpath,'SleepEvents');
    save(spindlestatspath,'SpindleStats')
end
end

