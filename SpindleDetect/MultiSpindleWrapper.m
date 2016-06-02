function [ output_args ] = MultiSpindleWrapper(recname,recfolder,corticalchannels)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%
%
%% LOAD THE SESSION IN FMA
%recname = 'c3po_160202';
%recname = '20140526_277um';
datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/BWData/~updated/Recordings (1)/';
recname = '20140526_277um';
%datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/DTData/';
%recname = 'DT2_rPPC_rCCG_362um_218um_20160209_160209_183610';
figloc = '/Users/dlevenstein/Dropbox/Research/Current Projects/TransitionsAndSpindles/AnalysisFigsMats/SpindleID/';
probegroups = {1:5,7:12}; %Group 1: Cingulate.  %Group2: Parietal DT2, spike group 6,13 are top of shank
%probegroups = {7:12,1:6}; %Group 1: Cingulate.  %Group2: HPC c3po
probegroups = {1:4,5:6}; %Group 1: Cingulate.  %Group2: HPC jenn1
numprobes = length(probegroups);



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
%%
load([datasetfolder,recname,'/',recname,'_SleepScore.mat'])

NREMint = StateIntervals.SWSpacket;
%NREMint = NREMint(1:4:end,:);
% delta = [1 4];
% theta = [6 10];
% spindle = [10 20];

ctxchannels = [spikegroups{probegroups{1}}];
%% STEP 1: Calculate Pairwise relationships (TBD) between cortical sites
% [gapowercorrmat,sites] = LFPPairMat('channels',[ctxchannels],...
%      'metric','powercorr','frange','gamma','numload',48,'int',NREMint);
% 
% [spcorrmat,sites] = LFPPairMat('channels',[ctxchannels],...
%     'metric','filtcorr','frange','spindles','numload',48,'int',NREMint);
% 
% [sppowercorrmat,sites] = LFPPairMat('channels',[ctxchannels],...
%     'metric','powercorr','frange','spindles','numload',48,'int',NREMint);
% 
% %% Figure
% figure
%     subplot(2,2,3)
%         imagesc(gapowercorrmat)
%         colorbar
%         title('Gamma Power Correlation')
%         xlabel('Channel (Sorted by Shank)');
%         ylabel('Channel (Sorted by Shank)');
%     subplot(2,2,2)
%         imagesc(spcorrmat)
%         title('Spindle-Filtered LFP Correlation')
%         colorbar
%         xlabel('Channel (Sorted by Shank)');
%         ylabel('Channel (Sorted by Shank)');
%     subplot(2,2,1)
%         imagesc(sppowercorrmat)
%         title('Spindle Power Correlation')
%         colorbar
%         xlabel('Channel (Sorted by Shank)');
%         ylabel('Channel (Sorted by Shank)');
%         
% saveas(gcf,[figloc,recname,'_spindlepowercorr'],'jpeg')

% %% Power-Power Correlation
% frange = [6 35];
% nfreqs = 200;
% ncyc = 8;
% 
% [freqs,t,spec] = WaveSpec(allLFP(:,2),frange,nfreqs,ncyc,1/sf_LFP,'log');
% spec = abs(spec);
% powercorr = corr((spec').^2);
% 
% %% FIgure
% figure
%     imagesc(log2(freqs),log2(freqs),powercorr)
%     axis xy
%     ColorbarWithAxis([0.25 1],'Corr.')
%     LogScale('xy',2)
%     xlabel('f (Hz)');ylabel('f (Hz)')
%     title('Power-Power Correlation')


%% Step 2a: Peak Method
%broadspband = [10 17];
broadspband = [8 20];
broadspband = [10 20];
broadspband = [7 20];
chanavg = false;
thresh = 2.5;
[ sppeaks,sppeakfreq,sppeakssign,sppeakmag ] = FindSpindlePeaks( ctxchannels,NREMint,broadspband,thresh,chanavg,figloc,recname );

[pSpindleInts_peak] = JoinSpindlePeaks(sppeaks,figloc,recname);


%% Step 2b: Identidy Spindle Intervals by amplitude envelope peak.
chanavg = true;
broadspband = [7 18]; %gardner
%broadspband = [8 18]; 
%broadspband = [10 15]; %peyrache
[pSpindleInts_env] = SpindleEnvelopeInts(ctxchannels,NREMint,broadspband,chanavg,figloc,recname);


%% Step 2c: Wavelet Method
%random site from each shank.
ctxprobechannel = cellfun(@(X) datasample(X,1),spikegroups(probegroups{1}));
[pSpindleInts_wav] = SpindleWaveletInts(ctxprobechannel,NREMint,figloc,recname);

%% Compare Spindle Ints
intnames = {'7-18 Env.','7-18 Wav'};
CompareSpindleIntSets(pSpindleInts_env,pSpindleInts_wav,20,intnames,figloc);

%% Adjust Start End Time to peaks, return spindle cycle normalized time
frange = [7 18];
[pSpindleInts,cycletimemap] = SetSpindleIntTime(pSpindleInts_env,frange,ctxchannels);


%% Delta Peak Times
[deltapeaks] = DeltaPeakTimes(ctxchannels,NREMint,figloc,recname);


%% pSpindle LFP and peaks (add peaks characterization from below)
SpindleIntLFP(pSpindleInts,cycletimemap,deltapeaks,ctxchannels,NREMint,figloc,recname)


%% Step 4: Characterize Properties of in-pSpindle Peaks
%pSpindleInts = pSpindleInts_wav;

chanavg = true;
thresh = 0.5;
broadspband = [7 18];
[sppeaks,sppeakfreq,sppeakssign,sppeakmag] = FindSpindlePeaks( ctxchannels,NREMint,broadspband,thresh,chanavg);

sppeaks_pos = sppeaks(sppeakssign>0);
sppeaks_neg = sppeaks(sppeakssign<0);
IPI_pos = sppeaks_pos(2:end)-sppeaks_pos(1:end-1);
IPI_neg = sppeaks_neg(2:end)-sppeaks_neg(1:end-1);
[~,posidx] = RestrictInts(sppeaks_pos,pSpindleInts);
[~,negidx] = RestrictInts(sppeaks_neg,pSpindleInts);
IPI_pos = IPI_pos(posidx(1:end-1));
IPI_neg = IPI_neg(negidx(1:end-1));
%%
support = 1;
pSpLen = pSpindleInts(:,2) - pSpindleInts(:,1);
[intts,sortindex,intts_sorted] = SortByIntTime(sppeaks,pSpindleInts+pSpLen*support*[-1 1],'norm');
intts = intts.*(2*support+1)-support;

binedges = linspace(-1,2,31);
[ bincenters,binmeans,binstd ] = BinDataTimes( sppeakfreq,intts,binedges );
%%
[freqhist,freqbins] = hist3([intts,sppeakfreq],[80 35]);

%%
negcolor = [0 0 0.8];
poscolor = [0.6 0 0];
figure
    subplot(2,2,1)
        plot(intts(sppeakssign==1),sppeakfreq(sppeakssign==1),'.','color',poscolor)
        hold on
        plot(intts(sppeakssign==-1),sppeakfreq(sppeakssign==-1),'.','color',negcolor)
       % ylim([6 18])
	subplot(2,2,2)
        imagesc(freqbins{1},freqbins{2},freqhist')
        hold on
        plot(bincenters,binmeans,'wo-')
        axis xy
        %ylim([6 18])
    subplot(2,2,3)
        plot(1./IPI_pos(1:end-1),1./IPI_pos(2:end),'.','color',poscolor)
        hold on
        plot(1./IPI_neg(1:end-1),1./IPI_neg(2:end),'.','color',negcolor)


%% SpindleCycle Normalized Time


downsamplefactor = 1;
allLFP = GetLFP_Down(ctxchannels,'intervals',NREMint,'downsample',downsamplefactor);
sf_LFP = 1250./downsamplefactor;

allLFP(:,2) = mean(allLFP(:,2:end),2);
allLFP(:,3:end) = [];


broadspband = [8 20];
LFPspindle = FilterLFP(allLFP,'passband',frange,'nyquist',sf_LFP./2);
[phase,amplitude,unwrapped] = Phase(samples,times)
PhaseToCycleTime(ttime,phase,ints)






%%

%%
% bins = [-3:0.0075:2];
% freqbins = [0:0.2:20];
bins = [-3:0.01:2];
freqbins = [0:0.33:20];
N_freq = hist3([1./(IPIs_same(1:end-1)),1./(IPIs_same(2:end))],{freqbins,freqbins});
N_same = hist3([log10(IPIs_same(1:end-1)),log10(IPIs_same(2:end))],{bins,bins});

 %%   Figure: IPIs

figure
subplot(2,1,1)
    [peakhist,histbins] = hist(log2(IPIs_same(IPIs_same>0.02 & IPIs_same<10& peaksign(1:end-1)==1)),250);
    [troughhist] = hist(log2(IPIs_same(IPIs_same>0.02 & IPIs_same<10 & peaksign(1:end-1)==-1)),histbins);
    bar(histbins,[troughhist' peakhist'],'stacked')
    hold on
   % bar(histbins,troughhist,'r')
    LogScale('x',2)
    xlabel('Inter-Peak Interval (s)')
    ylabel('Number of Peaks/Troughs')
    title('Inter-Peak Intervals: identified peaks/troughs on the same shank')
    
% create smaller axes in top right, and plot on it
axes('Position',[.55 .715 .3 .15])
box on
    hist((IPIs_same(IPIs_same>0.02 &IPIs_same<0.25)),250)
   % LogScale('x',10)
    xlabel('Inter-Peak Interval (s)')
    %ylabel('Number of Peaks/Troughs')
    title('Linear Scale for <250ms')
    
    
  subplot(2,1,2)  
    [peakhist,histbins] = hist((1./IPIs_same(IPIs_same>0.05 & IPIs_same<0.16& peaksign(1:end-1)==1)),30);
    [troughhist] = hist((1./IPIs_same(IPIs_same>0.05 & IPIs_same<0.16 & peaksign(1:end-1)==-1)),histbins);
    bar(histbins,[troughhist' peakhist'],'stacked')
    hold on
   % bar(histbins,troughhist,'r')
    %LogScale('x',2)
    xlabel('1/Inter-Peak Interval (Hz)')
    ylabel('Number of Peaks/Troughs')
    title({'Spindle-Peaks/Troughs from all sites,', 'Instantaneous Frequency'})
    legend('Inter-Trough Intervals','Inter-Peak Intervals')
    
saveas(gcf,[figloc,recname,'_IPIs'],'jpeg')    


%% InterPeak Temporal Statistics - summary figure
figure
    subplot(3,2,[1,3])
        plot(log10(IPIs_same(1:end-1)),log10(IPIs_same(2:end)),'k.','MarkerSize',0.1)
       LogScale('xy',10)
       xlabel('Inter-Peak Interval n (s)');
       ylabel('Inter-Peak Interval n+1 (s)');
        xlim([-1.5 1.5]);ylim([-1.5 2.5])
          
          hold on
    [peakhist,histbins] = hist(log10(IPIs_same(IPIs_same>0 & IPIs_same<50)),150);
    bar(histbins,([peakhist']./max(peakhist))+1.5,'BaseValue',1.5,'facecolor','k')
    
%   isimin = 0.055;
%   isimax = 0.16;    
%     subplot(2,3,2)
%         plot(log2(IPIs_same(1:end-1)),log2(IPIs_same(2:end)),'k.','MarkerSize',8)
%       LogScale('xy',2)
%        xlabel('Inter-Peak Interval n (s)');
%        ylabel('Inter-Peak Interval n+1 (s)');
%         xlim(log2([isimin isimax]));ylim(log2([isimin ((isimax-isimin)./log2(3))+isimax]))
%           
%           hold on
%     normmag = (log2(isimax)-log2(isimin))./3;
%     [peakhist,histbins] = hist(log2(IPIs_same(IPIs_same>isimin & IPIs_same<isimax)),25);
%     bar(histbins,([peakhist']./max(peakhist)).*normmag+log2(isimax),'BaseValue',log2(isimax),'facecolor','k')
%    
  fmin = 6;
  fmax = 18;
subplot(3,2,[2 4])
%     imagesc(freqbins(freqbins>fmin&freqbins<fmax),freqbins(freqbins>fmin&freqbins<fmax),...
%         log10(N_freq(freqbins>fmin&freqbins<fmax,freqbins>fmin&freqbins<fmax)))
plot(1./(IPIs_same(1:end-1)),1./(IPIs_same(2:end)),'k.','MarkerSize',8)
    axis xy
     xlim([fmin fmax])
   %LogScale('xy',10)
   xlabel('1/Inter-Peak Interval n (Hz)');
   %colorbar
  % caxis([0 100])
   ylabel('1/Inter-Peak Interval n+1 (Hz)');
   ylim([fmin ((fmax-fmin)./3)+fmax])
    
          hold on
    [peakhist,histbins] = hist(1./(IPIs_same(IPIs_same>0.05 & IPIs_same<0.16)),25);
    bar(histbins,(([peakhist']./max(peakhist))*(fmax-fmin)./3)+fmax,'BaseValue',fmax,'facecolor','k')
    
    
    saveas(gcf,[figloc,recname,'_peaktimingmasterfig'],'jpeg')    
    

%% IPI map
figure
% subplot(2,2,1)
%     plot(1./(IPIs_same(1:end-1)),1./(IPIs_same(2:end)),'.')
%    %LogScale('xy',10)
%    xlabel('1/Inter-Peak Interval n (Hz)');
%    ylabel('1/Inter-Peak Interval n+1 (Hz)');
% %   xlim([-3 1]);ylim([-3 1])
%    title('1/Inter-Peak Interval for peaks/troughs on the same shank')
   
subplot(2,2,4)
    imagesc(freqbins(freqbins>6&freqbins<19),freqbins(freqbins>6&freqbins<19),log10(N_freq(freqbins>6&freqbins<19,freqbins>6&freqbins<19)))
    axis xy
    %   xlim([-3 1]);ylim([-3 1])
   %LogScale('xy',10)
   xlabel('1/Inter-Peak Interval n (Hz)');
   %colorbar
  % caxis([0 100])
   ylabel('1/Inter-Peak Interval n+1 (Hz)');
subplot(2,2,1)
    plot(log10(IPIs_same(1:end-1)),log10(IPIs_same(2:end)),'.','MarkerSize',0.1)
   LogScale('xy',10)
   xlabel('Inter-Peak Interval n (s)');
   ylabel('Inter-Peak Interval n+1 (s)');
      xlim([-1.5 2]);ylim([-1.5 2])
      title('Spindle Peak/Trough Iterative Map')
   
subplot(2,2,3)
    imagesc(bins,bins,log10(N_same))
    axis xy
       xlim([-1.5 2]);ylim([-1.5 2])
   LogScale('xy',10)
   xlabel('Inter-Peak Interval n (s)');
   ylabel('Inter-Peak Interval n+1 (s)');
   
subplot(2,2,2)
    imagesc(bins./log10(2),bins./log10(2),log10(N_same))
    axis xy
       xlim(log2([0.05 0.16]));ylim(log2([0.05 0.16]))
   LogScale('xy',2)
   xlabel('Inter-Peak Interval n (s)');
   ylabel('Inter-Peak Interval n+1 (s)');


   saveas(gcf,[figloc,recname,'_ISImap'],'jpeg') 
   

%% Join Peaks/Troughs to Spindles
distthresh = 0.2;
[spints] = MergeSeparatedInts([joinedpeaks joinedpeaks],distthresh);
spdur = spints(:,2)-spints(:,1);

spmaxdur = 3;
spmindur = 0.3;

spints(spdur>spmaxdur | spdur<spmindur,:) = [];
spdur(spdur>spmaxdur | spdur<spmindur) = [];

ISpI = spints(2:end,1)-spints(1:end-1,2);

inspindlepeaks = joinedpeaks;
[inspindlepeaks,spindlepeakidx] = RestrictInts(joinedpeaks,spints);
inspindlepeaks = inspindlepeaks(:,1);
inspindlepeaksign = joinedpeakssign(spindlepeakidx);
speak_pos = inspindlepeaks(inspindlepeaksign==1);
speak_neg = inspindlepeaks(inspindlepeaksign==-1);
IspeakI_pos = speak_pos(2:end)-speak_pos(1:end-1);
IspeakI_neg = speak_neg(2:end)-speak_neg(1:end-1);
IspeakI_all = [IspeakI_pos;IspeakI_neg];



%% Figure: Spindle Stats
figure
    subplot(4,2,1)
        hist(spdur,20)
        hold on
        text(0.75,40,['Total # pSpindles: ',num2str(length(spdur))])
        xlabel('pSpindle Duration');
        ylabel('# pSpindles')
    subplot(4,2,3)
        hist(log10(ISpI),25)
        xlabel('Inter-pSpindle Interval');
        ylabel('# pSpindles')
        LogScale('x',10)
        xlim([-1 2])
    subplot(2,2,3)
        plot(log10(ISpI(1:end-1)),log10(ISpI(2:end)),'.')
        xlabel('Inter-pSpindle Interval n')
        ylabel('Inter-pSpindle Interval n+1')
        %title('Iterative Map: InterSpindle Intervals')
        LogScale('xy',10)
        xlim([-1 2]);ylim([-1 2])
    subplot(2,2,4)
        plot(log10(IspeakI_pos(1:end-1)),log10(IspeakI_pos(2:end)),'k.')
        hold on
        plot(log10(IspeakI_neg(1:end-1)),log10(IspeakI_neg(2:end)),'k.')
        xlabel('Inter-Peak Interval n')
        ylabel('Inter-Peak Interval n+1')
        LogScale('xy',10)
        xlim([-1.5 2]);ylim([-1.5 2])
        title('Iterative Map: Within-Spindle Peaks')
    subplot(4,2,2)
        hist(1./IspeakI_all,100)
        xlim([6 20])
        xlabel('Instantaneous Frequency (Hz)');
saveas(gcf,[figloc,recname,'_pspindlestats'],'jpeg')







   
%% Calculate Spectrogram for plot


exchannel = 2;

% Load downsampled LFP
specdownfactor = 10;
specLFP = GetLFP_Down(exchannel,'downsample',specdownfactor);
sf_specLFP = DATA.rates.lfp./specdownfactor;

freqlist = logspace(0,2,100);
%freqlist = linspace(0.5,55.5,111);
window = 10;
noverlap = 9;
window = window*sf_specLFP;
noverlap = noverlap*sf_specLFP;
[FFTspec,FFTfreqs,t_FFT] = spectrogram(specLFP(:,2),window,noverlap,freqlist,sf_specLFP);


FFTspec = abs(FFTspec);
[zFFTspec,mu,sig] = zscore(log10(FFTspec)');
%% Whole Recording
zoomwin = t_FFT([1 end]);
figure
    subplot(4,1,1)
        imagesc(t_FFT,log2(FFTfreqs),log10(FFTspec))
        StateScorePlot( {StateIntervals.WAKEstate,StateIntervals.SWSstate,StateIntervals.REMstate},...
            {'k','b','r'})
        axis xy
        caxis([min(mu)-2.5*max(sig) max(mu)+2.5*max(sig)])
       xlim([zoomwin])
       box(gca,'off')
       LogScale('y',2)


    subplot(4,1,2)
       plot(peakt(1:end-1),1./IPIs_same,'.','MarkerSize',0.1)
       ylim([0 20])
       xlabel('t (s)')
       ylabel('Instant. Freq. (Hz)')
       ylim([6 18])
       xlim([zoomwin])

saveas(gcf,[figloc,recname,'_wholerecsp'],'jpeg')


%% Sleep Episode
zoomwin = [2100 4150];
figure
    subplot(4,1,1)
        imagesc(t_FFT,log2(FFTfreqs),log10(FFTspec))
        StateScorePlot( {StateIntervals.WAKEstate,StateIntervals.SWSstate,StateIntervals.REMstate},...
            {'k','b','r'})
        axis xy
        caxis([min(mu)-2.5*max(sig) max(mu)+2.5*max(sig)])
       xlim([zoomwin])
       box(gca,'off')
       LogScale('y',2)


    subplot(4,1,2)
       plot(peakt(1:end-1),1./IPIs_same,'.','MarkerSize',0.1)
       ylim([0 20])
       xlabel('t (s)')
       ylabel('Instant. Freq. (Hz)')
       ylim([6 18])
       xlim([zoomwin])

saveas(gcf,[figloc,recname,'_sleepepisode'],'pdf')


%% Packets
zoomwin = [2150 2750];
figure
    subplot(3,1,1)
        imagesc(t_FFT,log2(FFTfreqs),log10(FFTspec))
        StateScorePlot( {StateIntervals.WAKEstate,StateIntervals.SWSstate,StateIntervals.REMstate},...
            {'k','b','r'})
        axis xy
        caxis([min(mu)-2.5*max(sig) max(mu)+2.5*max(sig)])
       xlim([zoomwin])
       box(gca,'off')
       LogScale('y',2)


    subplot(3,1,2)
       plot(peakt(1:end-1),1./IPIs_same,'.','MarkerSize',0.1)
       ylim([0 20])
       xlabel('t (s)')
       ylabel('Instant. Freq. (Hz)')
       ylim([6 18])
       xlim([zoomwin])

    subplot(6,1,5)
        plot(allLFP(:,1),allLFP(:,exchannel),'k')
       % plot(specLFP(:,1),specLFP(:,2),'k')
        
               xlim([zoomwin])
saveas(gcf,[figloc,recname,'_sleeppackets'],'pdf')

%% Few Spindles
%zoomwin = [3700 3727];
zoomwin = [3950 3973];
%zoomwin = [3630 3655];
%zoomwin = [2800 2815];
figure
    subplot(3,1,3)
        plot(peakt(peaksign==1),peakchan(peaksign==1),'b.')
        hold on
        plot(peakt(peaksign==-1),peakchan(peaksign==-1),'r.')
       xlabel('t (s)')
       ylabel('Site Number')
       title('Sites each peak was identified on')
       xlim([zoomwin])

    subplot(3,1,1)
       plot(peakt(1:end-1),1./IPIs_same,'.')%,'MarkerSize',0.1)
       ylim([0 20])
       xlabel('t (s)')
       ylabel('Instant. Freq. (Hz)')
       ylim([6 18])
       xlim([zoomwin])
       
LFP4plot = zscore(allLFP(:,exchannel));
    subplot(3,1,2)
        plot(allLFP(:,1),LFP4plot,'k')
       % plot(specLFP(:,1),specLFP(:,2),'k')
       % plot(peakt,specLFP(ismember(specLFP(:,1),peakt),2),'o')
        set(gca,'Yticklabel',[])
        ylabel('LFP')
               xlim([zoomwin])
saveas(gcf,[figloc,recname,'_sleepspindlesmeso'],'pdf')

%% Spindles
zoomwin = [3960 3964];
%zoomwin = [2800 2815];
figure
    subplot(3,1,3)
        plot(peakt(peaksign==1),peakchan(peaksign==1),'b.')
        hold on
        plot(peakt(peaksign==-1),peakchan(peaksign==-1),'r.')
       xlabel('t (s)')
       ylabel('Site Number')
       title('Sites each peak was identified on')
       xlim([zoomwin])

    subplot(3,1,1)
       plot(peakt(1:end-1),1./IPIs_same,'.')%,'MarkerSize',0.1)
       ylim([0 20])
       xlabel('t (s)')
       ylabel('Instant. Freq. (Hz)')
       ylim([6 18])
       xlim([zoomwin])
       
LFP4plot = zscore(allLFP(:,exchannel));
    subplot(3,1,2)
    plot(LFPspindle(:,1),LFPspindle(:,exchannel),'b','linewidth',0.2)
    hold on
        plot(allLFP(:,1),LFP4plot,'k','linewidth',1)
        hold on
        plot(unique(peakt),LFP4plot(ismember(allLFP(:,1),unique(peakt))),'.r')
       % plot(specLFP(:,1),specLFP(:,2),'k')
       % plot(peakt,specLFP(ismember(specLFP(:,1),peakt),2),'o')
        set(gca,'Yticklabel',[])
        ylabel('LFP')
               xlim([zoomwin])
saveas(gcf,[figloc,recname,'_sleepspindles'],'pdf')

%% Joining Peaks
minpeakdist = 0.02;

[joinedpeaks,joindex] = unique(peakt_sorted);
joinedpeakssign = peaksign_sorted(joindex);
joinedpeaksheight = peakheight_sorted(joindex);
pp = 1;
while pp < length(joinedpeaks)
    if joinedpeaks(pp+1)-joinedpeaks(pp)<=minpeakdist
        joinedpeaks(pp+1) = [];
        joinedpeakssign(pp+1) = [];
        joinedpeaksheight(pp) = max(joinedpeaksheight([pp pp+1]));
        joinedpeaksheight(pp+1) = [];
    else
        pp = pp+1;
    end
end

%%
peaks_pos = joinedpeaks(joinedpeakssign==1);
peaks_neg = joinedpeaks(joinedpeakssign==-1);
IPIs_pos = diff(peaks_pos);
IPIs_neg = diff(peaks_neg);
IPIs_joined = diff(joinedpeaks);


%%   Figure: IPI HIst

figure
subplot(2,3,1:2)
    [peakhist,histbins] = hist(log2(IPIs_pos(IPIs_pos<10)),125);
    [troughhist] = hist(log2(IPIs_neg( IPIs_neg<10)),histbins);
    bar(histbins,[troughhist' peakhist'],'stacked')
    hold on
   % bar(histbins,troughhist,'r')
    LogScale('x',2)
    xlim([-4.5 3])
    xlabel('Inter-Peak Interval (s)')
    ylabel('Number of Peaks/Troughs')
    title('Inter-Peak Intervals: Spindle Peaks/Troughs')

% create smaller axes in top right, and plot on it
axes('Position',[.3 .715 .3 .15])
box on
    hist([IPIs_pos(IPIs_pos<0.3); IPIs_neg(IPIs_neg<0.3)],50)
   % LogScale('x',10)
    xlabel('Inter-Peak Interval (s)')
    %ylabel('Number of Peaks/Troughs')
    title('Linear Scale for <250ms')
    set(gca,'Yticklabel',[])
    xlim([0.05 0.3])
    
    
  subplot(2,3,4:5)  
    [peakhist,histbins] = hist(1./(IPIs_pos(IPIs_pos<10)),125);
    [troughhist] = hist(1./(IPIs_neg(IPIs_neg<10)),histbins);
    bar(histbins,[troughhist' peakhist'],'stacked')
    hold on
   % bar(histbins,troughhist,'r')
    %LogScale('x',2)
    xlabel('1/Inter-Peak Interval (Hz)')
    ylabel('Number of Peaks/Troughs')
    title({'Spindle Peaks/Troughs: Instantaneous Frequency'})
    legend('Inter-Trough Intervals','Inter-Peak Intervals','location','northwest')
    xlim([6.5 18])
    
saveas(gcf,[figloc,recname,'_IPIsjoined'],'jpeg')  



%%
bins = [-3:0.01:2];
freqbins = [0:0.25:20];
N_freq = hist3([1./([IPIs_pos(1:end-1);IPIs_neg(1:end-1)]),1./([IPIs_pos(2:end);IPIs_neg(2:end)])],{freqbins,freqbins});
N_same = hist3([log10([IPIs_pos(1:end-1);IPIs_neg(1:end-1)]),log10([IPIs_pos(2:end);IPIs_neg(2:end)])],{bins,bins});
%% IPI map
figure
% subplot(2,2,1)
%     plot(1./(IPIs_same(1:end-1)),1./(IPIs_same(2:end)),'.')
%    %LogScale('xy',10)
%    xlabel('1/Inter-Peak Interval n (Hz)');
%    ylabel('1/Inter-Peak Interval n+1 (Hz)');
% %   xlim([-3 1]);ylim([-3 1])
%    title('1/Inter-Peak Interval for peaks/troughs on the same shank')
   
subplot(2,2,4)
    imagesc(freqbins(freqbins>6&freqbins<19),freqbins(freqbins>6&freqbins<19),log10(N_freq(freqbins>6&freqbins<19,freqbins>6&freqbins<19)))
    axis xy
    %   xlim([-3 1]);ylim([-3 1])
   %LogScale('xy',10)
   xlabel('1/Inter-Peak Interval n (Hz)');
   %colorbar
  % caxis([0 100])
   ylabel('1/Inter-Peak Interval n+1 (Hz)');
subplot(2,2,1)
    plot(log10(IPIs_pos(1:end-1)),log10(IPIs_pos(2:end)),'k.','MarkerSize',0.1)
    hold on
    plot(log10(IPIs_neg(1:end-1)),log10(IPIs_neg(2:end)),'k.','MarkerSize',0.1)
   LogScale('xy',10)
   xlabel('Inter-Peak Interval n (s)');
   ylabel('Inter-Peak Interval n+1 (s)');
      xlim([-1.5 2]);ylim([-1.5 2])
      title('Spindle Peak/Trough Iterative Map')
   
subplot(2,2,3)
    imagesc(bins,bins,log10(N_same))
    axis xy
       xlim([-1.5 2]);ylim([-1.5 2])
   LogScale('xy',10)
   xlabel('Inter-Peak Interval n (s)');
   ylabel('Inter-Peak Interval n+1 (s)');
   xlim([-1.25 1]);ylim([-1.25 1])
   
subplot(2,2,2)
    imagesc(bins./log10(2),bins./log10(2),log10(N_same))
    axis xy
       xlim(log2([0.05 0.16]));ylim(log2([0.05 0.16]))
   LogScale('xy',2)
   xlabel('Inter-Peak Interval n (s)');
   ylabel('Inter-Peak Interval n+1 (s)');


   saveas(gcf,[figloc,recname,'_ISImapjoined'],'jpeg') 


%% Peak Magnitude

figure
subplot(2,2,1)
    plot(1./([IPIs_pos; nan]),joinedpeaksheight(joinedpeakssign==1),'b.')
    hold on
    plot(1./([IPIs_neg; nan]),-joinedpeaksheight(joinedpeakssign==-1),'r.')
    xlim([6.5 17])
    xlabel('Instantaneous Frequency (Hz)');
    ylabel({'Peak Height', '(STDs above NREM Mean)'})
    
    
saveas(gcf,[figloc,recname,'_IPIspeakmag'],'jpeg')  

%%
figure
    subplot(2,1,1)
        hist(log10(IPIs_joined),150)
        LogScale('x',10)
        xlim([-1.7 1])
        ylabel('Number of Peaks/Toughs')
        xlabel('Inter-Peak/Trough Interval (s)')
        title({'Inter peak-or-trough interval histogram, for joining threshold'})
       % ylim([0 200])
        
  axes('Position',[.35 .7 .5 .15])
box on
        hist(log2(IPIs_joined),150)
        LogScale('x',2)
        xlim([-6 4])
        %ylabel('Number of Peaks/Toughs')
        set(gca,'ytick',[])
        xlabel('Inter-Peak/Trough Interval (s)')
        title('Zoom')
        ylim([0 150])   
        
saveas(gcf,[figloc,recname,'_IPIjoinhist'],'jpeg')  
        
       %% 
    subplot(2,3,4:5)
        hist((IPIs_joined(IPIs_joined<0.4)),60)
        xlabel('Inter-Peak/Trough Interval (s)')
        ylabel('Number of Peaks/Toughs')
        title('Linear Scale')
    %LogScale('x',10)
    %xlim([0 0.1])
    
    axes('Position',[.3 .25 .3 .15])
box on
    hist((IPIs_joined(IPIs_joined<0.4)),60)
   % LogScale('x',10)
    xlabel('Inter-Peak/Trough Interval (s)')
    
    %ylabel('Number of Peaks/Troughs')
    title('Zoom')
    set(gca,'Yticklabel',[])
    %xlim([0.05 0.3])
    ylim([0 100])

    saveas(gcf,[figloc,recname,'_IPIpeakandtrough'],'jpeg')  
%%

figure
    subplot(2,1,1)
        plot(LFPspindle(:,1),LFPspindle(:,2))
        
    subplot(2,1,2)
        plot(allLFP(:,1),allLFP(:,2))
        
%% Merge Spindles
mindists = [0.06:0.02:1];
numdists = length(mindists);
numbins = 100;
durbins = linspace(-2,0.5,numbins);
ISIbins = linspace(-1,2,numbins);

durhist = zeros(numdists,numbins);
ISIhist = zeros(numdists,numbins);

for dd = 1:numdists
%minspdist = 0.3;
    minspdist = mindists(dd);
    [spints] = MergeSeparatedInts([joinedpeaks joinedpeaks],minspdist);
    %%
    spdur = spints(:,2)-spints(:,1);
    ISpI = spints(2:end,1)-spints(1:end-1,2);
    
    durhist(dd,:) = hist(log10(spdur),durbins);
    durhist(dd,:) = durhist(dd,:)./sum(durhist(dd,:));
    ISIhist(dd,:) = hist(log10(ISpI),ISIbins);
    ISIhist(dd,:) = ISIhist(dd,:)./sum(ISIhist(dd,:));
end
%%

distthresh = 0.2;
%mindists==0.5





%% Histogram changes with joining threshold
figure
    subplot(2,2,3)
        bar(durbins,(durhist(mindists==distthresh,:))')
        hold on
        plot(log10([0.25 0.25]),get(gca,'ylim'),'r')
        plot(log10([2 2]),get(gca,'ylim'),'r')
        xlabel('Putative Spindle Duration (s)');
        ylabel('Proportion of Spindles');
        LogScale('x',10)
        xlim([-2.1 0.5])

    subplot(2,2,4)
        bar(ISIbins,(ISIhist(mindists==distthresh,:))')
        xlabel('Inter-pSpindle Interval (s)');
        ylabel('Proportion of Spindles');
        axis xy
        LogScale('x',10)
        xlim([-1 2])
    subplot(2,2,1)
        imagesc(mindists,durbins,(durhist)')
        hold on
        plot([distthresh distthresh],get(gca,'ylim'),'r')
        axis xy
        ylabel('Putative Spindle Duration (s)');
        xlabel('Minimum Spindle Separation');
        LogScale('y',10)
        %colorbar
        caxis([0 0.15])
    subplot(2,2,2)
        imagesc(mindists,ISIbins,(ISIhist)')
                hold on
        plot([distthresh distthresh],get(gca,'ylim'),'r')
        ylabel('Inter-pSpindle Interval (s)');
        xlabel('Minimum Spindle Separation');
        axis xy
        LogScale('y',10)

        saveas(gcf,[figloc,recname,'_joiningspindles'],'jpeg')   
%% Join Peaks/Troughs to Spindles
distthresh = 0.2;
[spints] = MergeSeparatedInts([joinedpeaks joinedpeaks],distthresh);
spdur = spints(:,2)-spints(:,1);

spmaxdur = 3;
spmindur = 0.3;

spints(spdur>spmaxdur | spdur<spmindur,:) = [];
spdur(spdur>spmaxdur | spdur<spmindur) = [];

ISpI = spints(2:end,1)-spints(1:end-1,2);

inspindlepeaks = joinedpeaks;
[inspindlepeaks,spindlepeakidx] = RestrictInts(joinedpeaks,spints);
inspindlepeaks = inspindlepeaks(:,1);
inspindlepeaksign = joinedpeakssign(spindlepeakidx);
speak_pos = inspindlepeaks(inspindlepeaksign==1);
speak_neg = inspindlepeaks(inspindlepeaksign==-1);
IspeakI_pos = speak_pos(2:end)-speak_pos(1:end-1);
IspeakI_neg = speak_neg(2:end)-speak_neg(1:end-1);
%%
figure
    subplot(2,2,1)
        hist(spdur,20)
        hold on
        text(0.75,40,['Total # pSpindles: ',num2str(length(spdur))])
        xlabel('pSpindle Duration');
        ylabel('# pSpindles')
    subplot(2,2,2)
        hist(log10(ISpI),20)
        xlabel('Inter-pSpindle Interval');
        ylabel('# pSpindles')
        LogScale('x',10)
    subplot(2,2,3)
        plot(log10(ISpI(1:end-1)),log10(ISpI(2:end)),'.')
        xlabel('Inter-pSpindle Interval n')
        ylabel('Inter-pSpindle Interval n+1')
        title('Iterative Map: InterSpindle Intervals')
        LogScale('xy',10)
    subplot(2,2,4)
        plot(log10(IspeakI_pos(1:end-1)),log10(IspeakI_pos(2:end)),'k.')
        hold on
        plot(log10(IspeakI_neg(1:end-1)),log10(IspeakI_neg(2:end)),'k.')
        xlabel('Inter-Peak Interval n')
        ylabel('Inter-Peak Interval n+1')
        LogScale('xy',10)
        title('Iterative Map: Within-Spindle Peaks')
saveas(gcf,[figloc,recname,'_pspindlestats'],'jpeg')



%% Sleep Episode
zoomwin = [2100 4150];
figure
    subplot(4,1,1)
        imagesc(t_FFT,log2(FFTfreqs),log10(FFTspec))
        StateScorePlot( {StateIntervals.WAKEstate,StateIntervals.SWSstate,StateIntervals.REMstate},...
            {'k','b','r'})
        axis xy
        caxis([min(mu)-2.5*max(sig) max(mu)+2.5*max(sig)])
       xlim([zoomwin])
       box(gca,'off')
       LogScale('y',2)


    subplot(4,1,2)
       plot(speak_pos(1:end-1),1./IspeakI_pos,'.','MarkerSize',0.5)
       hold on
       plot(speak_neg(1:end-1),1./IspeakI_neg,'.','MarkerSize',0.5)
       ylim([0 20])
       xlabel('t (s)')
       ylabel('Instant. Freq. (Hz)')
       ylim([6 18])
       xlim([zoomwin])

saveas(gcf,[figloc,recname,'_sleepepisode'],'jpeg')


%% Packets
zoomwin = [3550 4050];
figure
    subplot(3,1,1)
        imagesc(t_FFT,log2(FFTfreqs),log10(FFTspec))
        StateScorePlot( {StateIntervals.WAKEstate,StateIntervals.SWSstate,StateIntervals.REMstate},...
            {'k','b','r'})
        axis xy
        caxis([min(mu)-2.5*max(sig) max(mu)+2.5*max(sig)])
       xlim([zoomwin])
       box(gca,'off')
       LogScale('y',2)


    subplot(3,1,2)
       plot(speak_pos(1:end-1),1./IspeakI_pos,'.','MarkerSize',5)
       hold on
       plot(speak_neg(1:end-1),1./IspeakI_neg,'.','MarkerSize',5)
       ylim([0 20])
       xlabel('t (s)')
       ylabel('Instant. Freq. (Hz)')
       ylim([6 18])
       xlim([zoomwin])

saveas(gcf,[figloc,recname,'_sleeppackets'],'jpeg')



%% Spindles
figure
for exnum = 1:6
%exnum = 1;
speventnum = randi(length(spdur),1);
zoomwin = spints(speventnum,:)+[-1 1];
%zoomwin = [2800 2815];


       
LFP4plot = zscore(allLFP(:,exchannel));
    subplot(3,2,exnum)
    plot(LFPspindle(:,1),LFPspindle(:,exchannel),'b','linewidth',0.05)
    hold on
        plot(allLFP(:,1),LFP4plot,'k','linewidth',1)
        plot(spints',ones(size(spints))'.*4,'r')
        hold on
        plot(unique(joinedpeaks),LFP4plot(ismember(allLFP(:,1),unique(joinedpeaks))),'.r')
       % plot(specLFP(:,1),specLFP(:,2),'k')
       % plot(peakt,specLFP(ismember(specLFP(:,1),peakt),2),'o')
        set(gca,'Yticklabel',[])
        ylabel('LFP')
        ylim([-4.5 4.5])
               xlim([zoomwin])
end              
               
saveas(gcf,[figloc,recname,'_sleepspindle_example',num2str(exnum)],'pdf')
end

