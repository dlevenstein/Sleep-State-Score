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

%Currently, LFP are held in their own .mats... will make this to extract
%directly from .lfp files
load([datasetfolder,recordingname,'/',recordingname,'_LFP.mat'])
swLFP = LFP;
load([datasetfolder,recordingname,'/',recordingname,'_ThetaLFP.mat'])
thLFP = LFP;
clear LFP

%EMG is already calculated and in it's own .mat, need to add the ability to
%extract this.
load([datasetfolder,recordingname,'/',recordingname,'_EMGCorr.mat'])
EMG = EMGCorr(:,2);
clear EMGCorr

sf_LFP = 1250;
sf_EMG = 2;

figloc = ['StateScoreFigures/',recordingname];
%% Deal with input options from varargin
%none yet, but will do this with inputParser when we have input options

%possible variable input options
%'timewin' - only state score a subset of the recording
%'HPCsites' - site indices for HPC probes
%'figloc' - secondardy folder to save figures to
%'spikegroups' - if not in the .xml file

%% CALCULATE EMG FROM HIGH-FREQUENCY COHERENCE
%BW: you have all the code for this

%Select channels to use for EMG.... for now just use one from each spike
%group? or whatever you did previously


%[EMG] = GetEMGCorr(EMGchannels);


%% DETERMINE BEST SLOW WAVE AND THETA CHANNELS
%Possibility - use multiple channels for concensus, find one good SWChannel
%from each shank? - could even use all SW for clustering - will improve
%both time resolution and reliability
%[SWChannel] = FindSWChannel();
%[ThetaChannel] = FindThetaChannel();


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

end

