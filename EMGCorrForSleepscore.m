function EMGCorr = EMGCorrForSleepscore(basenamepath,scoretime,specialchannels,specialshanks)
% Based on Erik Schomburg's work and code.  Grabs channels and calculates
% their correlations in the 300-600Hz band over sliding windows of 0.5sec.
% Channels are automatically selected and are a combination first channels
% on each shank (not using directly adjacent shanks, they must be at least 
% 2 shanks away.  
% First channels on "specialshanks" are always used
% "specialchannels" are also always used
% Requires .eeg/lfp and .xml.  Assumes each spikegroup in the .xml
% represents a "shank"
% 
% Mean pairwise correlations are calculated for each time point.
% 
% Adapted 2015 Brendon Watson

%% Parameters
savebool = 0;

%% input handling
if ~exist('specialchannels','var')
    specialchannels = [];
end
if ~exist('specialshanks','var')
    specialshanks = [];
end

%% get basics about eeg/lfp file
if strcmp(basenamepath(end-3:end),'.lfp') || strcmp(basenamepath(end-3:end),'.eeg')
    eegloc = basenamepath;
    xmlloc = [basenamepath(1:end-4),'.xml'];
    saveloc = [basenamepath(1:end-4),'_EMGCorr.mat'];
else
    if ~isempty(dir('*.eeg'))
        eegloc = [basenamepath '.eeg'];
    elseif ~isempty(dir('*.lfp'))
        eegloc = [basenamepath '.lfp'];
    else
        return
    end
    xmlloc = [basenamepath,'.xml'];
    saveloc = [basenamepath,'_EMGCorr.mat'];
end

Par = LoadPar_SleepScore(xmlloc);

Fs = Par.lfpSampleRate; % Hz, LFP sampling rate
nChannels = Par.nChannels;

if isfield(Par,'SpkGrps')
    SpkGrps = Par.SpkGrps;
elseif isfield(Par,'AnatGrps')
    SpkGrps = Par.AnatGrps;
    display('No SpikeGroups, Using AnatomyGroups')
else
    display('No SpikeGroups...')
end
    

xcorr_halfwindow_s = 0.5;%specified in s
% downsampleFs = 125;
% downsampleFactor = round(Fs/downsampleFs);
binScootS = 0.5;
binScootSamps = Fs*binScootS;
corrChunkSz = 20;%for batch-processed correlations

%% Pick shanks to analyze
% get spike groups,
% pick every other one... unless specialshanks, in which case pick non-adjacent
spkgrpstouse = 1:2:length(SpkGrps);
spkgrpstouse = spkgrpstouse(:);

if ~isempty(specialshanks) 
    spkgrpstouse = unique(cat(1,spkgrpstouse,specialshanks));
end
if ~isempty(specialchannels)
    for a = 1:length(specialchannels)
        for b = 1:length(SpkGrps)
            if ismember(specialchannels(a),SpkGrps(b).Channels)
                spkgrpstouse = cat(1,spkgrpstouse,b);
                continue
            end
        end
    end
end

xcorr_chs = [];
for a = 1:length(spkgrpstouse)
    [lia,lib]=ismember(specialchannels,SpkGrps(a).Channels+1);
    if any(lia)
        xcorr_chs(end+1) = SpkGrps(a).Channels(lib)+1;
    else
        xcorr_chs(end+1) = SpkGrps(a).Channels(1)+1;%correct for base 0/1 difference
    end
end

xcorr_chs = unique(xcorr_chs);

%% Read and filter channel
% read channel
eeg = readmulti(eegloc, nChannels, xcorr_chs); %read and convert to mV    
% Filter first in high frequency band to remove low-freq physiologically
% correlated LFPs (e.g., theta, delta, SPWs, etc.)

eeg = LoadBinary_Down(eegloc,'frequency',Fs,...
    'nchannels',nChannels,'channels',xcorr_chs,...
    'start',scoretime(1),'duration',diff(scoretime));


xcorr_freqband = [275 300 600 625]; % Hz
eeg = filtsig_in(eeg, Fs, xcorr_freqband);

%% xcorr 'strength' is the summed correlation coefficients between channel
% pairs for a sliding window of 25 ms
xcorr_window_samps = round(xcorr_halfwindow_s*Fs);
xcorr_window_inds = -xcorr_window_samps:xcorr_window_samps;%+- that number of ms in samples

% xcorr 'strength' is the averaged correlation coefficients between channel
% old version... single correlation calculated at once
% % pairs for a sliding window of 25 ms
% xcorr_window_samps = round(xcorr_window_s*Fs);
% xcorr_window_inds = -xcorr_window_samps:xcorr_window_samps;%+- that number of ms in samples
% 
% numbins = (size(eeg,1) - xcorr_window_samps*2)/binScootSamps;
% xcorrStrength = zeros(numbins, 1);
% binind = 0;
% for i=(1+xcorr_window_inds(end)):binScootSamps:(size(eeg,1)-xcorr_window_inds(end))
%     binind = binind+1;
%     for j=1:(length(xcorr_chs)-1)
%         for k=(j+1):length(xcorr_chs)
%             s1 = eeg(i + xcorr_window_inds, j);
%             s2 = eeg(i + xcorr_window_inds, k);
%             tmp = corrcoef(s1,s2);
%             xcorrStrength(binind) = xcorrStrength(binind) + tmp(1,2);
%         end
%     end
% end
% xcorrStrength = xcorrStrength/(length(xcorr_chs)*(length(xcorr_chs)-1)/2);


% new version... batches of correlation calculated at once
timestamps = (1+xcorr_window_inds(end)):binScootSamps:(size(eeg,1)-xcorr_window_inds(end));
numbins = length(timestamps);
EMGCorr = zeros(numbins, 1);
% tic
for j=1:(length(xcorr_chs)-1)
    for k=(j+1):length(xcorr_chs)
        c1 = [];
        c2 = [];
        binind = 0;
        binindstart = 1;
        for i = timestamps
            binind = binind+1;
            s1 = eeg(i + xcorr_window_inds, j);
            s2 = eeg(i + xcorr_window_inds, k);
            c1 = cat(2,c1,s1);
            c2 = cat(2,c2,s2);
            if size(c1,2) == corrChunkSz || i == timestamps(end)
                binindend = binind;
                tmp = corr(c1,c2);
                tmp = diag(tmp);
                EMGCorr(binindstart:binindend) = EMGCorr(binindstart:binindend) + tmp;
                c1 = [];
                c2 = [];
                binindstart = binind+1;
            end
        end
    end
end
% toc

EMGCorr = EMGCorr/(length(xcorr_chs)*(length(xcorr_chs)-1)/2);

EMGCorr = cat(2,timestamps'/Fs,EMGCorr);
ChannelsCompared = xcorr_chs;
% EMGCorrData = v2struct(EMGCorr,ChannelsCompared,AnatShankSite);
EMGCorrData = v2struct(EMGCorr,ChannelsCompared);
if savebool
    % save...
    save(saveloc,'EMGCorrData','EMGCorr');
end

function [filt_sig, Filt] = filtsig_in(sig, Fs, filtband_or_Filt)
% [filt_sig, Filt] = filtsig(sig, dt_ms, filtband_or_Filt)
%
% Created by: Erik Schomburg, 2011

if isnumeric(filtband_or_Filt)
    h  = fdesign.bandpass(filtband_or_Filt(1), filtband_or_Filt(2), filtband_or_Filt(3), filtband_or_Filt(4), ...
        60, 1, 60, Fs);
    Filt = design(h, 'butter', 'MatchExactly', 'passband');
else
    Filt = filtband_or_Filt;
end

if ~isempty(sig)
    if iscell(sig)
        filt_sig = cell(size(sig));
        for i=1:length(sig(:))
            filt_sig{i} = filter(Filt, sig{i});
            filt_sig{i} = filter(Filt, filt_sig{i}(end:-1:1));
            filt_sig{i} = filt_sig{i}(end:-1:1);
        end
    elseif ((size(sig,1) > 1) && (size(sig,2) > 1))
        filt_sig = zeros(size(sig));
        for i=1:size(filt_sig,2)
            filt_sig(:,i) = filter(Filt, sig(:,i));
            filt_sig(:,i) = filter(Filt, filt_sig(end:-1:1,i));
            filt_sig(:,i) = filt_sig(end:-1:1,i);
        end
    else
        filt_sig = filter(Filt, sig);
        filt_sig = filter(Filt, filt_sig(end:-1:1));
        filt_sig = filt_sig(end:-1:1);
    end
else
    filt_sig = [];
end

