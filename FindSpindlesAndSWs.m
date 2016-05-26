function [ pSpindleInts,cycletimemap,deltapeaks ] = FindSpindlesAndSWs(datasetfolder,recname,figfolder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%




%% LOAD THE SESSION IN FMA
%recname = 'c3po_160202';
%recname = '20140526_277um';
%datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/BWData/~updated/Recordings (1)/';
%recname = '20140526_277um';
%datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/DTData/';
%recname = 'DT2_rPPC_rCCG_362um_218um_20160209_160209_183610';
%figfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/TransitionsAndSpindles/AnalysisFigsMats/SpindleID/';
% probegroups = {1:5,7:12}; %Group 1: Cingulate.  %Group2: Parietal DT2, spike group 6,13 are top of shank
% %probegroups = {7:12,1:6}; %Group 1: Cingulate.  %Group2: HPC c3po
% probegroups = {1:4,5:6}; %Group 1: Cingulate.  %Group2: HPC jenn1
%numprobes = length(probegroups);



xmlfilename = [datasetfolder,'/',recname,'/',recname,'.xml'];
if exist (fullfile(datasetfolder,recname,[recname,'.lfp']),'file')
    rawlfppath = fullfile(datasetfolder,recname,[recname,'.lfp']);
elseif exist (fullfile(datasetfolder,recname,[recname,'.eeg']),'file')
    rawlfppath = fullfile(datasetfolder,recname,[recname,'.eeg']);
else 
    display('No .lfp file')
end


SetCurrentSession(xmlfilename)
global DATA;

spikegroups = DATA.spikeGroups.groups;

numsites = DATA.nChannels;
numgroups = length(spikegroups);


% Par = LoadPar(xmlfilename);
% Fs = Par.lfpSampleRate; % Hz, LFP sampling rate
% nChannels = Par.nChannels;
% spikegroups = {Par.SpkGrps(:).Channels};
numgroups = length(spikegroups);

%% Load the channel map to get the cortical channels
CTXlabels = {'mPFC','ACC','MotorCtx','OFC'};

spikegroupanatomyfilename = fullfile(datasetfolder,recname,[recname,'_SpikeGroupAnatomy.csv']);
spkgroupanatomy=readtable(spikegroupanatomyfilename);

ctxgroups = ismember(spkgroupanatomy.AnatomicalSite,CTXlabels);
ctxgroups = spkgroupanatomy.SpikeGroup(ctxgroups);
ctxchannels = [spikegroups{ctxgroups}];

%%

load(fullfile(datasetfolder,recname,[recname,'_SleepScore.mat']));

NREMint = StateIntervals.NREMstate;


%% Step 2b: Identidy Spindle Intervals by amplitude envelope peak.
chanavg = true;
%broadspband = [7 18]; %gardner
broadspband = [8 18]; 
%broadspband = [10 15]; %peyrache
[pSpindleInts_env] = SpindleEnvelopeInts(ctxchannels,NREMint,broadspband,chanavg,figfolder,recname);


%% Step 2c: Wavelet Method
%random site from each shank.
%ctxprobechannel = cellfun(@(X) datasample(X,1),spikegroups(probegroups{1}));
%[pSpindleInts_wav] = SpindleWaveletInts(ctxprobechannel,NREMint,figloc,recname);


%% Adjust Start End Time to peaks, return spindle cycle normalized time
frange = [8 18];
[pSpindleInts,cycletimemap] = SetSpindleIntTime(pSpindleInts_env,frange,ctxchannels);


%% Delta Peak Times
[deltapeaks] = DeltaPeakTimes(ctxchannels,NREMint,figfolder,recname);


%% Characteriation of pSpindle LFP and peaks (add peaks characterization from below)
SpindleIntLFP(pSpindleInts,cycletimemap,deltapeaks,ctxchannels,NREMint,figfolder,recname)



end

