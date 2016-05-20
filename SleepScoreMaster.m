function [stateintervals,episodeintervals] = SleepScoreMaster(datasetfolder,recordingname,varargin)
%[stateintervals,episodeintervals] = SleepScoreMaster(datasetfolder,recordingname)
%This is the master function for sleep state scoring.
%
%INPUT
%   datasetfolder   Top level folder in which the dataset resides. 
%                   For example
%                   '/Users/dlevenstein/Dropbox/Research/Datasets/BWData/'
%   recordingname   Name of the recording, this will be the name of the
%                   folder in which the .lfp file and other files reside.
%                   For example, the .lfp file should be:
%                   'datasetfolder/recordingname/recordingname.lfp'
%                   ... it is also assumed that this serves as the basename
%                   for the files for instance data will be at
%                   /datasetfolder/recordingname/recordingname.lfp
%
%OUTPUT
%   stateintervals  {Nstates} Cell array of state start and end times, each
%                   member of the array will be [Nints x 2]
%   episodeintervals    similar for episodes
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
    [datasetfolder,recordingname] = fileparts(DIRECTORYNAME); 
    if isequal(DIRECTORYNAME,0) 
       disp('User pressed cancel...')
       return
    end

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
        SleepScoreMaster(datasetfolder,recordingname{rr},varargin)
        close all
    end
elseif numrecs == 1 & iscell(recordingname)
        recordingname = recordingname{1};
end
	
    

%%
display(['Scoring Recording: ',recordingname]);
%
%% DEV
% Load the necessary files as needed for development
% datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/DTData/';
%datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/BWData/~updated/Recordings (1)/';
%recordingname = 'c3po_160202';
% recordingname = 'DT2_rPPC_rCCG_362um_218um_20160209_160209_183610';
sessionfolder = fullfile(datasetfolder,recordingname);

sf_LFP = 1250;
sf_EMG = 2;

savebool = 1;

figloc = [fullfile(sessionfolder,'StateScoreFigures'),'/'];
%figloc = fullfile(sessionfolder,'StateScoreFigures');

if ~exist(figloc,'dir')
    mkdir(figloc)
end
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
%'savefiles'    - save the EMG,LFP files to .mats?

%Overwrite pickes new and writes over existing ThLFP ans SWLFP
overwrite = 1;


%% Database File Management 

%Filenames for EMG, thLFP, and swLFP .mat files in the database.
EMGpath = fullfile(datasetfolder,recordingname,[recordingname '_EMGCorr.mat']);
thetalfppath = fullfile(datasetfolder,recordingname,[recordingname,'_ThetaLFP.mat']);
swlfppath = fullfile(datasetfolder,recordingname,[recordingname,'_SWLFP.mat']);

sleepstatepath = fullfile(datasetfolder,recordingname,[recordingname,'_SleepScore.mat']);

if exist (fullfile(datasetfolder,recordingname,[recordingname,'.lfp']),'file')
    rawlfppath = fullfile(datasetfolder,recordingname,[recordingname,'.lfp']);
elseif exist (fullfile(datasetfolder,recordingname,[recordingname,'.eeg']),'file')
    rawlfppath = fullfile(datasetfolder,recordingname,[recordingname,'.eeg']);
else 
    display('No .lfp file')
end

%% CALCULATE EMG FROM HIGH-FREQUENCY COHERENCE
% Load/Calculate EMG based on cross-shank correlations 
% (high frequency correlation signal = high EMG).  
% Schomburg E.W. Neuron 84, 470?485. 2014)
% Do this before lfp load because finding proper lfps will depend on this.
% If EMG is already calculated and in it's own .mat, then load, otherwise
% calculate this

if ~exist(EMGpath,'file')
    display('Calculating EMG')
    EMGCorr = EMGCorrForSleepscore(rawlfppath);%BW modify this to have different dependencies, currently assumes presence of: 
    % eeg filename - ok
    % .xml filename - ok
    %     Save ..._EMGCorr file

else
    display('EMG aleady calculated: Loading')
    load(EMGpath,'EMGCorr')
end
EMG = EMGCorr(:,2);
clear EMGCorr


%% DETERMINE BEST SLOW WAVE AND THETA CHANNELS
%Possibility - use multiple channels for concensus, find one good SWChannel
%from each shank? - could even use all SW for clustering - will improve
%both time resolution and reliability

%Possbile cases:
%   -SW,TH channel picked and have .mat
%   -SW,TH channel picked but need to be loaded from .eeg -> save a .mat
%   -SW,TH channel not picked

%To Do: get LoadLFP (or other) to load already downsampled LFP....
%Then return channel and load for .mat and Clustering



if ~exist(thetalfppath,'file') && ~exist(swlfppath,'file') || overwrite; % if no lfp file already, load lfp and make lfp file?
%     if exist (fullfile(datasetfolder,recordingname,[recordingname,'.lfp']),'file')
%         rawlfppath = fullfile(datasetfolder,recordingname,[recordingname,'.lfp']);
%     elseif exist (fullfile(datasetfolder,recordingname,[recordingname,'.eeg']),'file')
%         rawlfppath = fullfile(datasetfolder,recordingname,[recordingname,'.eeg']);
%     end

%     Par = LoadPar(fullfile(datasetfolder,recordingname,[recordingname,'.xml']));
    display('Picking SW and TH Channels')
    [SWchannum,THchannum,swLFP,thLFP] = PickSWTHChannel(datasetfolder,recordingname,figloc);
    
    %open lfp, 
%         lfp = readmulti(eegloc, nChannels, xcorr_chs) * bmd.voltsperunit*1000; %read and convert to mV    
%         or lfp = LoadBinary...
    %[SWChannel] = FindSWChannel();
    %[ThetaChannel] = FindThetaChannel();
    %swLFP = lfp(SWChannel);
    %thLFP = lfp(ThetaChannel);
    %store _ThetaLFP.mat
    %store _SWLFP.mat
    if savebool
    % save...
        save(swlfppath,'swLFP');
        save(thetalfppath,'thLFP');
    end
else
    display('SW and TH Channels Already Extracted, Loading...')
    load(swlfppath,'swLFP')
    load(thetalfppath,'thLFP')
end



%% CLUSTER STATES BASED ON SLOW WAVE, THETA, EMG

display('Clustering States Based on EMG, SW, and TH LFP channels')
[stateintervals,states] = ClusterStates(swLFP,thLFP,EMG,sf_LFP,sf_EMG,figloc,[]);


%% JOIN STATES INTO EPISODES

SWSints = stateintervals{2};
REMints = stateintervals{3};
WAKEints = stateintervals{1};

minPACKdur = 20;
SWSlengths = SWSints(:,2)-SWSints(:,1);
packetintervals = SWSints(SWSlengths>=minPACKdur,:);

minintdur = 40;
minSWSdur = 20;
[episodeintervals{2}] = IDStateEpisode(SWSints,minintdur,minSWSdur);

minintdur = 40;
minWAKEdur = 20;
[episodeintervals{1}] = IDStateEpisode(WAKEints,minintdur,minWAKEdur);

minintdur = 40;
minREMdur = 20;
[episodeintervals{3}] = IDStateEpisode(REMints,minintdur,minREMdur);

% episodeidx = INTtoIDX(episodeintervals,ceil(t_LFP(end)));
% episodeintervals=IDXtoINT(episodeidx);
% episodeintervals = episodeintervals(2:4);

% BW: I added some other code, have to look at this too

%% Save
StateIntervals.SWSstate = SWSints;
StateIntervals.REMstate = REMints;
StateIntervals.WAKEstate = WAKEints;
StateIntervals.SWSepisode = episodeintervals{2};
StateIntervals.REMepisode = episodeintervals{3};
StateIntervals.WAKEeposode = episodeintervals{1};
StateIntervals.SWSpacket = packetintervals;

save(sleepstatepath,'stateintervals','episodeintervals','StateIntervals');
end

