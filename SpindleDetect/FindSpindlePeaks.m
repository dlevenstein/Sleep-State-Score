function [ joinedpeaks,joinedpeaksfreq,joinedpeakssign,joinedpeaksheight ] = FindSpindlePeaks( ctxchannels,NREMint,frange,thresh,chanavg,figloc,recname )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
%% Load LFP

% Load downsampled LFP
downsamplefactor = 1;
allLFP = GetLFP_Down(ctxchannels,'intervals',NREMint,'downsample',downsamplefactor);
sf_LFP = 1250./downsamplefactor;

%% mean over LFP
if chanavg
    allLFP(:,2) = mean(allLFP(:,2:end),2);
    allLFP(:,3:end) = [];
end

%% filter in spindle band, half rectify
%broadspband = [10 17]; %Looks pretty good...
%broadspband = [7 16];
%broadspband = [8 20];
LFPspindle = FilterLFP(allLFP,'passband',frange,'nyquist',sf_LFP./2,'order',4);
LFPspindle(:,2:end) = zscore(LFPspindle(:,2:end));
%LFPspindle(LFPspindle<0) = -LFPspindle(LFPspindle<0);
%% find peaks over threshold
%thresh = 2.5; %STD
%thresh = 2; %STD
%minpeakdist = 0.01;
minpeakdist = 0.9/frange(2);

%Shape for findpeaks
numchans = length(ctxchannels);
numchans = length(allLFP(1,2:end));
%t4peaks = repmat(LFPspindle(:,1),[numchans,1]);
%LFP4peaks = LFPspindle(:,2:end);%LFP4peaks=LFP4peaks(:);

peakchan = [];
peakt = [];
peakheight = [];
peaksign = [];
peakfreq = [];
for cc = 1:numchans
    [cpeak,cpeakt] = findpeaks(LFPspindle(:,cc+1),LFPspindle(:,1),...
        'MinPeakHeight',thresh,'MinPeakDistance',minpeakdist);
    cchan = cc*ones(size(cpeakt));
    psign = 1*ones(size(cpeakt));
    peakchan = [peakchan; cchan];
    peaksign = [peaksign; psign];
    peakt = [peakt; cpeakt];
    peakheight = [peakheight; cpeak];
    peakfreq = [peakfreq; 1./(cpeakt(2:end)-cpeakt(1:end-1)); nan];
    
    [cpeak,cpeakt] = findpeaks(-LFPspindle(:,cc+1),LFPspindle(:,1),...
        'MinPeakHeight',thresh,'MinPeakDistance',minpeakdist);
    cchan = cc*ones(size(cpeakt));
    psign = -1*ones(size(cpeakt));
    peakchan = [peakchan; cchan];
    peaksign = [peaksign; psign];
    peakt = [peakt; cpeakt];
    peakheight = [peakheight; cpeak];
    peakfreq = [peakfreq; 1./(cpeakt(2:end)-cpeakt(1:end-1)); nan];
end


%% sorted peaks
[peakt_sorted,sortpeaks] = sort(peakt);
peakchan_sorted = peakchan(sortpeaks);
peakheight_sorted = peakheight(sortpeaks);
peaksign_sorted = peaksign(sortpeaks);
peakfreq_sorted = peakfreq(sortpeaks);
%IPIs

%% Joining Peaks
minpeakdist = 0.02;

[joinedpeaks,joindex] = unique(peakt_sorted);
joinedpeakssign = peaksign_sorted(joindex);
joinedpeaksheight = peakheight_sorted(joindex);
joinedpeaksfreq = peakfreq_sorted(joindex);
pp = 1;
while pp < length(joinedpeaks)
    if joinedpeaks(pp+1)-joinedpeaks(pp)<=minpeakdist
        joinedpeaks(pp+1) = [];
        joinedpeakssign(pp+1) = [];
        joinedpeaksfreq(pp+1) = [];
        joinedpeaksheight(pp) = max(joinedpeaksheight([pp pp+1]));
        joinedpeaksheight(pp+1) = [];
    else
        pp = pp+1;
    end
end

peaks_pos = joinedpeaks(joinedpeakssign==1);
peaks_neg = joinedpeaks(joinedpeakssign==-1);
IPIs_pos = diff(peaks_pos);
IPIs_neg = diff(peaks_neg);
IPIs_all = diff(joinedpeaks);
IPIs_same = [IPIs_pos; IPIs_neg];
IPIs_all(IPIs_all<0)=[];

numpeaks = length(joinedpeaks);
%%
%% InterPeak Temporal Statistics - summary figure
if exist('figloc','var')

negcolor = [0 0 0.8];
poscolor = [0.6 0 0];

figure
    subplot(4,2,[3 5])
        plot(log10(IPIs_pos(1:end-1)),log10(IPIs_pos(2:end)),'.','color',poscolor,'MarkerSize',0.1)
        hold on
        plot(log10(IPIs_neg(1:end-1)),log10(IPIs_neg(2:end)),'.','color',negcolor,'MarkerSize',0.1)
       LogScale('xy',10)
       text(0.1,1,['# Peaks: ',num2str(numpeaks)])
       xlabel('Inter-Peak Interval n (s)');
       ylabel('Inter-Peak Interval n+1 (s)');
        xlim([-1.3 1.2]);ylim([-1.3 1.2])
     
	subplot(4,2,1)
        [peakhist,histbins] = hist(log10(IPIs_pos(IPIs_pos>0.02 & IPIs_pos<10)),100);
        [troughhist] = hist(log10(IPIs_neg(IPIs_neg>0.02 & IPIs_neg<10)),histbins);
        b = bar(histbins,[troughhist' peakhist'],'stacked','linestyle','none');
        b(1).FaceColor = poscolor;
        b(2).FaceColor = negcolor;
            xlim([-1.3 1.2]);
            LogScale('x',10)
            legend('Positive Peaks','Negative Peaks')
        
      fmin = 6;
      fmax = 18;
    subplot(4,2,[4 6])
        %     imagesc(freqbins(freqbins>fmin&freqbins<fmax),freqbins(freqbins>fmin&freqbins<fmax),...
        %         log10(N_freq(freqbins>fmin&freqbins<fmax,freqbins>fmin&freqbins<fmax)))
        plot(1./(IPIs_pos(1:end-1)),1./(IPIs_pos(2:end)),'.','color',poscolor,'MarkerSize',8)
        hold on
        plot(1./(IPIs_neg(1:end-1)),1./(IPIs_neg(2:end)),'.','color',negcolor,'MarkerSize',8)
        axis xy
         xlim([fmin fmax])
       %LogScale('xy',10)
       xlabel('1/Inter-Peak Interval n (Hz)');
       %colorbar
      % caxis([0 100])
       ylabel('1/Inter-Peak Interval n+1 (Hz)');
       ylim([fmin fmax])

    subplot(4,2,2)
        [peakhist,histbins] = hist(1./(IPIs_pos(IPIs_pos>0.05 & IPIs_pos<0.16)),25);
        [troughhist] = hist(1./(IPIs_neg(IPIs_neg>0.05 & IPIs_neg<0.16)),histbins);
        b = bar(histbins,[troughhist' peakhist'],'stacked','linestyle','none');
        b(1).FaceColor = poscolor;
        b(2).FaceColor = negcolor;
%         [peakhist,histbins] = hist(1./(IPIs_same(IPIs_same>0.05 & IPIs_same<0.16)),25);
%         b = bar(histbins,peakhist','stacked','linestyle','none');
%         b(1).FaceColor = poscolor;
%         b(2).FaceColor = negcolor;
        xlim([fmin fmax])
    
    
    saveas(gcf,[figloc,recname,'_peaktimingmasterfig'],'jpeg')    
end
end

