function [ pSpindleInts ] = SpindleWaveletInts( ctxchannels,NREMint,figloc,recname )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Load LFP

% Load downsampled LFP
downsamplefactor = 5;
allLFP = GetLFP_Down(ctxchannels,'intervals',NREMint,'downsample',downsamplefactor);
sf_LFP = 1250./downsamplefactor;

%% mean over LFP
chanavg = true;
if chanavg
     allLFP(:,2:end) = zscore(allLFP(:,2:end));
     allLFP(:,end+1) = mean(allLFP(:,2:end),2);
%     allLFP(:,3:end) = [];
end

%% fparms
frange = [7 18];
nfreqs = 50;
ncyc = 8;

%peakthresh = 4;
peakthresh = 3.5;
%peakthresh = 1.5;
minpeakdist = 0.1;
intthresh = 2.5;
%intthresh = 2.5;

minintdur = 0.4;
maxintdur = 4;
maxintseparation = 0.1;

peakchan = [];
peakt = [];
peakheight = [];
intt = [];

numchans = length(allLFP(1,2:end));
%%
for cc = 1:numchans
 %   cc = 1;
    [freqs,t,spec] = WaveSpec(allLFP(:,cc+1),frange,nfreqs,ncyc,1/sf_LFP,'lin');
    spec = zscore(abs(spec)'.^2);
    %spec = zscore(abs(spec)');

    %%    
    [maxval,maxfreqidx] = max(spec,[],2);
    maxfreq = freqs(maxfreqidx);
    
  % maxval(maxfreqidx==1 | maxfreqidx==nfreqs) = 0;
   % maxval = zscore(maxval);
    
    [cpeak,cpeakt] = findpeaks(maxval,allLFP(:,1),...
        'MinPeakHeight',peakthresh,'MinPeakDistance',minpeakdist);
    cchan = cc*ones(size(cpeakt));
    peakchan = [peakchan; cchan];
    peakt = [peakt; cpeakt];
    peakheight = [peakheight; cpeak];
    
    
    [periods] = Threshold([allLFP(:,1),maxval],'>',intthresh,'min',minintdur./2);
    
     [~,~,withpeakidx] = RestrictInts(cpeakt,periods);
     periods = periods(withpeakidx,:);

    
    intt = [intt ; periods];
    
    

%     thresh = 2.5;
%     overthresh = maxval>=peakthresh;
%     maxfreq(~overthresh) = nan;

end

%%

[ pSpindleInts ] = MergeSeparatedInts( intt,maxintseparation );

%  [~,~,pSpindleInts] = RestrictInts(peakt,megergints);
%  pSpindleInts = megergints(pSpindleInts,:);


 
 spdur = pSpindleInts(:,2)-pSpindleInts(:,1);
 pSpindleInts(spdur>maxintdur | spdur<minintdur ,:) = [];
 
 [~,intIDX] = RestrictInts(allLFP(:,1),pSpindleInts);

 maxfreq(~intIDX) = nan;
 
 numsp = length(pSpindleInts(:,1));
%%
figure
for ss = 1:4
    spnum = randi(numsp);
    trange = pSpindleInts(spnum,:)+[-1 1];
    specints = allLFP(:,1) >= trange(1) & allLFP(:,1) <= trange(2);
subplot(2,2,ss)
    imagesc(trange,(freqs),spec(specints,:)')
    hold on
    plot(allLFP(:,1),(maxfreq),'k')
    plot(allLFP(:,1),allLFP(:,end)+4,'k','Linewidth',1)
    ylim([0 20])
    caxis([0 3])
    axis xy
   % LogScale_ss('y',2)
    plot(pSpindleInts(spnum,:),3*[1 1],'r')
    xlim(pSpindleInts(spnum,:)+[-1 1])
end
%%
    figure
subplot(2,2,3)
    hist(maxval,15)
subplot(2,2,4)
    hist(maxfreq,15)
%%
    
    
    
% for ss = 1:2
% spnum = randi(numsp);
%     subplot(4,2,ss)
%     plot(amplitude(:,1),amplitude(:,2),'k')
%     hold on
%     plot(LFPspindle(:,1),LFPspindle(:,2),'b')
%     plot(allLFP(:,1),allLFP(:,2),'k')
%     plot(peakt,peakheight,'r.')
%     plot(pSpindleInts(spnum,:),3*[1 1],'r')
%     xlim(pSpindleInts(spnum,:)+[-1 1])
% end
% 
% 


end

