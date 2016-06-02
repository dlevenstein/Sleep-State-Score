function [ peakboundedints,cycletimemap ] = SetSpindleIntTime(pSpindleInts,frange,ctxchannels);
%[ pSpindleInts,cycletimemap ] = SetSpindleIntTime(pSpindleInts,frange,ctxchannels,NREMint,figloc,recname);
%changes spindle interval times to start/end at the previous/next peak in
%the spindle-filtered signal, and returns cycle-normalized time
%% Load LFP
intsupport = 1;
loadints = bsxfun(@(X,Y) X+Y,pSpindleInts,intsupport*[-1 1]);
[ loadints ] = MergeSeparatedInts( loadints,0 );
%%
% Load downsampled LFP
downsamplefactor = 1;
allLFP = GetLFP_Down(ctxchannels,'intervals',loadints,'downsample',downsamplefactor);
sf_LFP = 1250./downsamplefactor;

%% mean over LFP
chanavg = true;
if chanavg
     allLFP(:,2:end) = zscore(allLFP(:,2:end));
     allLFP(:,2) = mean(allLFP(:,2:end),2);
     allLFP(:,3:end) = [];
end

%% Find Spindle Peaks

LFPspindle = FilterLFP(allLFP,'passband',frange,'nyquist',sf_LFP./2,'order',4);

minpeakdist = 0.9/frange(2);
peakthresh = 0;
    [~,Sppeakt] = findpeaks(LFPspindle(:,2),LFPspindle(:,1),...
        'MinPeakHeight',peakthresh,'MinPeakDistance',minpeakdist);
    
%% Find peaks immediately before/after spindle intervals
peakboundedints = pSpindleInts;
peakboundedints(:,1) = interp1(Sppeakt,Sppeakt,pSpindleInts(:,1),'previous');
peakboundedints(:,2) = interp1(Sppeakt,Sppeakt,pSpindleInts(:,2),'next');


%% real time to indextime
t_idx = 1:length(allLFP(:,1));
SpIntIdx = interp1(allLFP(:,1),t_idx,peakboundedints,'nearest');

%%
numspindles = length(peakboundedints(:,1));
pad = 20; %samples

for ss = 1:numspindles
    idxwithpad = SpIntIdx(ss,1)-pad:SpIntIdx(ss,2)+pad;
    samples = LFPspindle(idxwithpad,:);
    [phases,~,unwrapped] = Phase(samples);
    unwrapped(1:pad,:) = [];
    unwrapped(end-pad:end,:) = [];
    if unwrapped(1,2) > pi;
        unwrapped(:,2) = unwrapped(:,2)-2*pi;
    end
    
    normtime = (1:length(unwrapped(:,1)))./length(unwrapped(:,1));
    
    cycletimemap{ss} = [unwrapped normtime'];
end


