function [ stateintervals,episodeintervals ] = SleepScoreMaster(datasetfolder,recordingname,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%
%
% DLevenstein and BWatson 2015/16
%% DEV
datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/BWData/';
recordingname = 'Dino_061814_mPFC';

load([datasetfolder,recordingname,'/',recordingname,'_LFP.mat'])
swLFP = LFP;
load([datasetfolder,recordingname,'/',recordingname,'_ThetaLFP.mat'])
thLFP = LFP;
clear LFP

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

%% CALCULATE EMG FROM HIGH-FREQUENCY COHERENCE
%BW: you have all the code for this

%[EMG] = GetEMGCorr();


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

