function [SWchannum,THchannum] = PickSWTHChannel(datasetfolder,recname,figfolder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%TO DO
%   -Change from GetLFP to LoadBinary or readmulti
%% DEV
datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/DTData/';
recname = 'DT3_rLS_rCA1_20150927_298um_288um';
figfolder = '/Users/dlevenstein/Code Library/SleepScoreDevelopment/StateScoreFigures/';

% datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/GGData/';
% recname = 'Rat08-20130717';
xmlfilename = [datasetfolder,'/',recname,'/',recname,'.xml'];
eegfilename = [datasetfolder,'/',recname,'/',recname,'.eeg'];

%% FMA

SetCurrentSession(xmlfilename);
global DATA
numchans = DATA.nChannels;

%% Hist/Freqs Parms
numhistbins = 21;
numfreqs = 100;

%% Load LFP files from .lfp
%To do here: catch situations where LFP is bigger than RAM and/or
%downsample within GetLFP
allLFP = GetLFP('all');

% Load downsampled LFP
%allLFP = LoadBinary(eegfilename,'frequency',DATA.rates.lfp,'nchannels',nChannels,'channels',channels,'skip',skip);

%% Downsample the LFP to 250Hz
sf_LFP = 1/(allLFP(2,1)-allLFP(1,1));
if sf_LFP == 1250
    downsamplefactor = 5;
else
    display('sf not 1250... if only you made this able to set its own downsample...')
    downsamplefactor = 1;
end
allLFP = downsample(allLFP,downsamplefactor);
sf_LFP = sf_LFP./downsamplefactor;

%% For each channel, calculate the PC1 and check it
pc1hists = zeros(numhistbins,numchans);
THhist = zeros(numhistbins,numchans);
pc1coeff = zeros(numfreqs,numchans);
THmeanspec = zeros(numfreqs,numchans);
dipSW = zeros(numchans,1);
dipTH = zeros(numchans,1);
for cc = 1:numchans;
%cc = 40;
channum = cc+1;
    display(['Channel ',num2str(cc),' of ',num2str(numchans)])

    %Calcualte Spectrogram
    freqlist = logspace(0,2,numfreqs);
    window = 10;
    noverlap = 9;
    window = window*sf_LFP;
    noverlap = noverlap*sf_LFP;
    [FFTspec,FFTfreqs,t_FFT] = spectrogram(allLFP(:,channum),window,noverlap,freqlist,sf_LFP);
    FFTspec = abs(FFTspec);

    %% Find TRANSIENTS and set to 0 for PCA
    [zFFTspec,mu,sig] = zscore(log10(FFTspec)');
    totz = zscore(abs(sum(zFFTspec')));
    badtimes = find(totz>5);
    zFFTspec(badtimes,:) = 0;

    %% PCA
    smoothfact = 10; %si_FFT
    thsmoothfact = 15;
     [COEFF, SCORE, LATENT] = pca(zFFTspec);
     SCORE(:,1) = smooth(SCORE(:,1),smoothfact);
    SCORE(:,1) = (SCORE(:,1)-min(SCORE(:,1)))./max(SCORE(:,1)-min(SCORE(:,1)));

    %% Histogram and diptest of PC1
    histbins = linspace(0,1,numhistbins);
    [pcahist]= hist(SCORE(:,1),histbins);

    pc1hists(:,cc) = pcahist;
    pc1coeff(:,cc) = COEFF(:,1);
    
    dipSW(cc) = hartigansdiptest(sort(SCORE(:,1)));
    
    
    %% Calculate theta

    %NarrowbandTheta
    f_all = [2 16];
    %f_all = [2 20];
    f_theta = [5 10];
    thfreqlist = logspace(log10(f_all(1)),log10(f_all(2)),numfreqs);

    [thFFTspec,thFFTfreqs] = spectrogram(allLFP(:,channum),window,noverlap,thfreqlist,sf_LFP);
    thFFTspec = (abs(thFFTspec));
    allpower = sum(log10(thFFTspec),1);
    %Why log10? Does it matter?  Could not log transform make th stand out?
    thfreqs = find(thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
    thpower = sum(log10(thFFTspec(thfreqs,:)),1);

    thratio = thpower./allpower;    %Narrowband Theta
    thratio = smooth(thratio,thsmoothfact);
    thratio = (thratio-min(thratio))./max(thratio-min(thratio));
    
    %% Histogram and diptest of Theta
    THhist(:,cc) = hist(thratio,histbins);
    THmeanspec(:,cc) = mean(thFFTspec,2);
    
    dipTH(cc) = hartigansdiptest(sort(thratio));
end

%% Sort by dip and pick channels
[~,dipsortSW] = sort(dipSW);
[~,dipsortTH] = sort(dipTH);

goodSWidx = dipsortSW(end);
goodTHidx = dipsortTH(end);
SWchan = goodSWidx+1;
THchan = goodTHidx+1;

swLFP = allLFP(:,SWchan);
thLFP = allLFP(:,THchan);

SWchannum = goodSWidx-2;
THchannum = goodTHidx-2;



%% Find Inverted PC1s and flip them
invpc1 = mean(pc1coeff(freqlist<4,:))<0 & mean(pc1coeff(freqlist>50,:))>0;
pc1coeff(:,invpc1) = -pc1coeff(:,invpc1);
pc1hists(:,invpc1) = flipud(pc1hists(:,invpc1));
%% Test
%PC1 coefficients for NREM match
%Theta spectrum for isolated peak?


%% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PC1 Weights and Coefficients


swfig = figure;
    subplot(2,2,1)
        imagesc(log2(FFTfreqs),1:numchans,pc1coeff(:,dipsortSW)')
        ylabel('Channel #');xlabel('f (Hz)')
        LogScale('x',2)
        axis xy
        title('PC1 Frequency Coefficients: All Channels') 
    subplot(2,2,2)
        imagesc(histbins,1:numchans,pc1hists(:,dipsortSW)')
        ylabel('Channel #');xlabel('PC1 projection weight')
        title('PC1 Projection Histogram: All Channels')
        axis xy
    subplot(2,2,3)
        set(gca,'ColorOrder',RainbowColors(length(dipsortSW)))
        hold all
        plot(log2(FFTfreqs),pc1coeff')     
        plot(log2(FFTfreqs([1 end])),[0 0],'k')
        ylabel('PC1 Coefficient');xlabel('f (Hz)')
        LogScale('x',2)
        title('PC1 Frequency Coefficients: All Channels')
    subplot(2,2,4)
        set(gca,'ColorOrder',RainbowColors(length(dipsortSW)))
        hold all
        plot(histbins,pc1hists')
        ylabel('hist');xlabel('PC1 projection weight')
        title('PC1 Projection Histogram: All Channels')
        
saveas(swfig,[figfolder,recname,'_FindBestSW'],'jpeg')

%% Theta Hist and Coefficients

thfig = figure;
    subplot(2,2,1)
        imagesc(log2(thFFTfreqs),1:numchans,THmeanspec(:,dipsortTH)')
        ylabel('Channel #');xlabel('f (Hz)')
        LogScale('x',2)
        axis xy
        title('Spectrum: All Channels') 
    subplot(2,2,2)
        imagesc(histbins,1:numchans,THhist(:,dipsortTH)')
        ylabel('Channel #');xlabel('PC1 projection weight')
        title('Theta Ratio Histogram: All Channels')
        axis xy
    subplot(2,2,3)
        set(gca,'ColorOrder',RainbowColors(length(dipsortTH)))
        hold all
        plot(log2(thFFTfreqs),THmeanspec')     
        plot(log2(thFFTfreqs([1 end])),[0 0],'k')
        ylabel('Power');xlabel('f (Hz)')
        LogScale('x',2)
        title('Spectrum: All Channels')
    subplot(2,2,4)
        set(gca,'ColorOrder',RainbowColors(length(dipsortTH)))
        hold all
        plot(histbins,THhist')
        ylabel('hist');xlabel('PC1 projection weight')
        title('Theta Ratio Histogram: All Channels') 
        
saveas(thfig,[figfolder,recname,'_FindBestTH'],'jpeg')
%% Show Channels


    %Calculate PC1 for plot/return
    [FFTspec,FFTfreqs,t_FFT] = spectrogram(allLFP(:,SWchan),window,noverlap,freqlist,sf_LFP);
    FFTspec = abs(FFTspec);
    [zFFTspec,mu,sig] = zscore(log10(FFTspec)');

    totz = zscore(abs(sum(zFFTspec')));
    badtimes = find(totz>5);
    zFFTspec(badtimes,:) = 0;
    
    smoothfact = 10; %si_FFT
    thsmoothfact = 15;
     [COEFF, SCORE, LATENT] = pca(zFFTspec);
     SCORE(:,1) = smooth(SCORE(:,1),smoothfact);
    PC1power = (SCORE(:,1)-min(SCORE(:,1)))./max(SCORE(:,1)-min(SCORE(:,1)));

chanfig =figure;
	subplot(5,1,1:2)
        imagesc(t_FFT,log2(FFTfreqs),log10(FFTspec))
        axis xy
        LogScale('y',2)
        caxis([min(mu)-2.5*max(sig) max(mu)+2.5*max(sig)])
        ylim([log2(FFTfreqs(1)) log2(FFTfreqs(end))+0.2])
        xlim(t_FFT([1,end]))
        ylabel({'LFP - FFT','f (Hz)'})
        title('SW Channel');
        
    subplot(5,1,3)
        plot(t_FFT,PC1power,'k')
        xlim(t_FFT([1,end]))
     
    %Calculate Theta ratio for plot/return    
    [thFFTspec,thFFTfreqs,t_FFT] = spectrogram(allLFP(:,THchan),window,noverlap,thfreqlist,sf_LFP);
    thFFTspec = abs(thFFTspec);
    [zFFTspec,mu,sig] = zscore(log10(thFFTspec)');
        
    allpower = sum(log10(thFFTspec),1);
    %Why log10? Does it matter?  Could not log transform make th stand out?
    thfreqs = find(thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
    thpower = sum(log10(thFFTspec(thfreqs,:)),1);

    thratio = thpower./allpower;    %Narrowband Theta
    thratio = smooth(thratio,thsmoothfact);
    thratio = (thratio-min(thratio))./max(thratio-min(thratio));
    
subplot(5,1,4)
    plot(allLFP(:,1),allLFP(:,SWchan),'k')
    
        imagesc(t_FFT,log2(thFFTfreqs),log10(thFFTspec))
        hold on
        plot(t_FFT([1,end]),log2(f_theta([1,1])),'w')
        plot(t_FFT([1,end]),log2(f_theta([2,2])),'w')
        axis xy
        LogScale('y',2)
        caxis([min(mu)-2.5*max(sig) max(mu)+2.5*max(sig)])
        ylim([log2(thFFTfreqs(1)) log2(thFFTfreqs(end))+0.2])
        ylabel({'LFP - FFT','f (Hz)'})
        title('Theta Channel');
        xlim(t_FFT([1,end]))
        
subplot(5,1,5)
        plot(t_FFT,thratio,'k')
        xlim(t_FFT([1,end]))
        
saveas(chanfig,[figfolder,recname,'_SWTHChannels'],'jpeg')
end

