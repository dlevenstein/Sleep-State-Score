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
%% DEV
% Load the necessary files as needed for development
datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/BWData/';
recordingname = 'Dino_061814_mPFC';
sessionfolder = fullfile(datasetfolder,recordingname);

sf_LFP = 1250;
sf_EMG = 2;

figloc = fullfile('StateScoreFigures',recordingname);
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

%% Database File Management 

%Filenames for EMG, thLFP, and swLFP .mat files in the database.
EMGpath = fullfile(datasetfolder,recordingname,[recordingname '_EMGCorr.mat']);
thetalfppath = fullfile(datasetfolder,recordingname,[recordingname,'_ThetaLFP.mat']);
swlfppath = fullfile(datasetfolder,recordingname,[recordingname,'_SWLFP.mat']);


%% CALCULATE EMG FROM HIGH-FREQUENCY COHERENCE
% Load/Calculate EMG based on cross-shank correlations 
% (high frequency correlation signal = high EMG).  
% Schomburg E.W. Neuron 84, 470?485. 2014)
% Do this before lfp load because finding proper lfps will depend on this.
% If EMG is already calculated and in it's own .mat, then load, otherwise
% calculate this

if ~exist(EMGpath,'file')

    EMGCorr = EMGCorrForSleepscore(sessionfolder,recordingname);%BW modify this to have different dependencies, currently assumes presence of: 
    % eeg filename - ok
    % .xml filename - ok
    %     Save ..._EMGCorr file

else
    load(EMGpath,'EMGCorr')
end
EMG = EMGCorr(:,2);
clear EMGCorr


%% DETERMINE BEST SLOW WAVE AND THETA CHANNELS
%Possibility - use multiple channels for concensus, find one good SWChannel
%from each shank? - could even use all SW for clustering - will improve
%both time resolution and reliability

if ~exist(thetalfppath,'file')% if no lfp file already, load lfp and make lfp file?
    if exist (fullfile(datasetfolder,recordingname,[recordingname,'.lfp']),'file')
        rawlfppath = fullfile(datasetfolder,recordingname,[recordingname,'.lfp']);
    elseif exist (fullfile(datasetfolder,recordingname,[recordingname,'.eeg']),'file')
        rawlfppath = fullfile(datasetfolder,recordingname,[recordingname,'.eeg']);
    end

    Par = LoadPar(fullfile(datasetfolder,recordingname,[recordingname,'.xml']));
    %open lfp, 
%         lfp = readmulti(eegloc, nChannels, xcorr_chs) * bmd.voltsperunit*1000; %read and convert to mV    
%         or lfp = LoadBinary...
    %[SWChannel] = FindSWChannel();
    %[ThetaChannel] = FindThetaChannel();
    %swLFP = lfp(SWChannel);
    %thLFP = lfp(ThetaChannel);
    %store _ThetaLFP.mat
    %store _SWLFP.mat
else
    load(swlfppath,'swLFP')
    load(thetalfppath,'thLFP')
end


%% CLUSTER STATES BASED ON SLOW WAVE, THETA, EMG

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

end

