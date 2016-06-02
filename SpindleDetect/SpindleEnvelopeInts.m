function [ pSpindleInts ] = SpindleEnvelopeInts(  ctxchannels,NREMint,broadspband,chanavg,figloc,recname )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%% Load LFP

% Load downsampled LFP
downsamplefactor = 5;
allLFP = GetLFP_Down(ctxchannels+1,'intervals',NREMint,'downsample',downsamplefactor);
sf_LFP = 1250./downsamplefactor;

%% mean over LFP
if chanavg
    allLFP(:,2) = mean(allLFP(:,2:end),2);
    allLFP(:,3:end) = [];
end

%% filter in spindle band, half rectify
%broadspband = [10 17]; %Looks pretty good...
%broadspband = [7 16];
%frange = [6 18];
LFPspindle = FilterLFP(allLFP,'passband',broadspband,'nyquist',sf_LFP./2);
LFPspindle(:,2:end) = zscore(LFPspindle(:,2:end));
allLFP(:,2:end) = zscore(allLFP(:,2:end));

%% Amplitude
[phase,amplitude] = Phase(LFPspindle);

%% find peaks and intervals over threshold on any sites
peakthresh = 3; %STD
intthresh = 1.5;
minpeakdist = 0.1;

minintdur = 0.4;
maxintdur = 4;
minintseparation = 0.1;

%Shape for findpeaks
numchans = length(allLFP(1,2:end));

peakchan = [];
peakt = [];
peakheight = [];
intt = [];
for cc = 1:numchans
    [cpeak,cpeakt] = findpeaks(amplitude(:,cc+1),amplitude(:,1),...
        'MinPeakHeight',peakthresh,'MinPeakDistance',minpeakdist);
    cchan = cc*ones(size(cpeakt));
    peakchan = [peakchan; cchan];
    peakt = [peakt; cpeakt];
    peakheight = [peakheight; cpeak];
    
    
    [periods] = Threshold(amplitude(:,[1 cc+1]),'>',intthresh,'min',minintdur);
    intt = [intt ; periods];
end

[ megergints ] = MergeSeparatedInts( intt,minintseparation );

%% Find Intervals that have a peak within them

 [~,~,pSpindleInts] = RestrictInts(peakt,megergints);
 pSpindleInts = megergints(pSpindleInts,:);
 
 spdur = pSpindleInts(:,2)-pSpindleInts(:,1);
 pSpindleInts(spdur>maxintdur,:) = [];
 
 numsp = length(pSpindleInts(:,1));
 
 %% Return normalized spindle cycle time...
%%

figure
for ss = 1:2
spnum = randi(numsp);
    subplot(4,2,ss)
    plot(amplitude(:,1),amplitude(:,2),'k')
    hold on
    plot(LFPspindle(:,1),LFPspindle(:,2),'b')
    plot(allLFP(:,1),allLFP(:,2),'k')
    plot(peakt,peakheight,'r.')
    plot(pSpindleInts(spnum,:),3*[1 1],'r')
    xlim(pSpindleInts(spnum,:)+[-1 1])
end
    subplot(4,1,2)
    plot(amplitude(:,1),amplitude(:,2),'k')
    hold on
    plot(LFPspindle(:,1),LFPspindle(:,2),'b')
    plot(allLFP(:,1),allLFP(:,2),'k')
    plot(peakt,peakheight,'r.')
    plot(pSpindleInts(spnum,:),3*[1 1],'r')
    %xlim(pSpindleInts(spnum,:)+[-1 1])
end

