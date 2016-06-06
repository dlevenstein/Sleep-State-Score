function [SWchannum,THchannum,swLFP,thLFP,t_LFP] = PickSWTHChannel(datasetfolder,recordingname,figfolder,scoretime)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%TO DO
%   -Change from GetLFP to LoadBinary or readmulti
%% DEV
%datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/DTData/';
%recordingname = 'DT2_rPPC_rCCG_362um_218um_20160209_160209_183610';
% datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/BWData/~updated/Recordings (1)/';
% recordingname = 'c3po_160202';
% figfolder = '/Users/dlevenstein/Code Library/SleepScoreDevelopment/StateScoreFigures/';

%recname = 'c3po_160202';
%datasetfolder = '/Users/dlevenstein/Dropbox/Share Folders/Recordings/';

% datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/GGData/';
% recname = 'Rat08-20130717';
xmlfilename = [datasetfolder,'/',recordingname,'/',recordingname,'.xml'];
if exist (fullfile(datasetfolder,recordingname,[recordingname,'.lfp']),'file')
    rawlfppath = fullfile(datasetfolder,recordingname,[recordingname,'.lfp']);
elseif exist (fullfile(datasetfolder,recordingname,[recordingname,'.eeg']),'file')
    rawlfppath = fullfile(datasetfolder,recordingname,[recordingname,'.eeg']);
else 
    display('No .lfp file')
end

%% FMA
% 
% SetCurrentSession(xmlfilename);
% global DATA
%nChannels = DATA.nChannels;

Par = LoadPar_SleepScore(xmlfilename);
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

%% Hist/Freqs Parms
numhistbins = 21;
numfreqs = 100;


%% Pick channels to use
spkgroupchannels = [SpkGrps.Channels];

%Add reject channels here...
rejectchannels = [];


usechannels = setdiff(spkgroupchannels,rejectchannels);
numusedchannels = length(usechannels);

%% Load LFP files from .lfp
%To do here: catch situations where LFP is bigger than RAM and/or
%downsample within GetLFP
%allLFP = GetLFP('all');

% Load downsampled LFP
downsamplefactor = 10;
allLFP = LoadBinary_Down(rawlfppath,'frequency',Fs,...
    'nchannels',nChannels,'channels',usechannels+1,'downsample',downsamplefactor,...
    'start',scoretime(1),'duration',diff(scoretime));
Fs = Fs./downsamplefactor;

%% For each channel, calculate the PC1 and check it
pc1hists = zeros(numhistbins,numusedchannels);
THhist = zeros(numhistbins,numusedchannels);
pc1coeff = zeros(numfreqs,numusedchannels);
THmeanspec = zeros(numfreqs,numusedchannels);
dipSW = zeros(numusedchannels,1);
dipTH = zeros(numusedchannels,1);
%%
for chanidx = 1:numusedchannels;
%channum = 1;
    display(['Channel ',num2str(chanidx),' of ',num2str(numusedchannels)])

    %Calcualte Z-scored Spectrogram
    freqlist = logspace(0,2,numfreqs);
    window = 10;
    noverlap = 9;
    window = window*Fs;
    noverlap = noverlap*Fs;
    [FFTspec,FFTfreqs,t_FFT] = spectrogram(allLFP(:,chanidx),window,noverlap,freqlist,Fs);
    FFTspec = abs(FFTspec);
    [zFFTspec,mu,sig] = zscore(log10(FFTspec)');

    %% Remove transients before calculating SW histogram
    %this should be it's own whole section - removing/detecting transients
    totz = zscore(abs(sum(zFFTspec')));
    badtimes = find(totz>5);
    zFFTspec(badtimes,:) = 0;
    
    %% PCA for Broadband Slow Wave
    [COEFF, SCORE, LATENT] = pca(zFFTspec);
   % broadbandSlowWave = SCORE(:,1);
    
	%% Set Broadband filter weights for Slow Wave
    load('SWweights.mat')
    assert(isequal(freqlist,SWfreqlist), 'spectrogram freqs.  are not what they should be...')
    broadbandSlowWave = zFFTspec*SWweights';
    
    %% Smooth and 0-1 normalize
    smoothfact = 10; %units of si_FFT
    thsmoothfact = 15;
     
    broadbandSlowWave = smooth(broadbandSlowWave,smoothfact);
    broadbandSlowWave = (broadbandSlowWave-min(broadbandSlowWave))./max(broadbandSlowWave-min(broadbandSlowWave));

    %% Histogram and diptest of PC1
    histbins = linspace(0,1,numhistbins);
    [pcahist]= hist(broadbandSlowWave,histbins);

    pc1hists(:,chanidx) = pcahist;
    pc1coeff(:,chanidx) = COEFF(:,1);
    
    dipSW(chanidx) = hartigansdiptest(sort(broadbandSlowWave));
    
    
    %% Calculate theta

    %NarrowbandTheta
    f_all = [2 16];
    %f_all = [2 20];
    f_theta = [5 10];
    thfreqlist = logspace(log10(f_all(1)),log10(f_all(2)),numfreqs);

    [thFFTspec,thFFTfreqs] = spectrogram(allLFP(:,chanidx),window,noverlap,thfreqlist,Fs);
    thFFTspec = (abs(thFFTspec));
    allpower = sum(log10(thFFTspec),1);
    %Why log10? Does it matter?  Could not log transform make th stand out?
    thfreqs = find(thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
    thpower = sum(log10(thFFTspec(thfreqs,:)),1);

    thratio = thpower./allpower;    %Narrowband Theta
    thratio = smooth(thratio,thsmoothfact);
    thratio = (thratio-min(thratio))./max(thratio-min(thratio));
    
    %% Histogram and diptest of Theta
    THhist(:,chanidx) = hist(thratio,histbins);

    
    dipTH(chanidx) = hartigansdiptest(sort(thratio));
    
    %% Theta Peak in mean spectrum
    THmeanspec(:,chanidx) = mean(thFFTspec,2);
    meanthratio = sum(log10(THmeanspec(thfreqs,chanidx)))./sum(log10(THmeanspec(:,chanidx)));
    dipTH(chanidx) = meanthratio;
end

%% Sort by dip and pick channels
[~,dipsortSW] = sort(dipSW);
[~,dipsortTH] = sort(dipTH);

goodSWidx = dipsortSW(end);
goodTHidx = dipsortTH(end);

SWchannum = usechannels(goodSWidx);
THchannum = usechannels(goodTHidx);

swthLFP = LoadBinary_Down(rawlfppath,'frequency',Fs,...
    'nchannels',nChannels,'channels',[SWchannum+1,THchannum+1],...
    'start',scoretime(1),'duration',diff(scoretime));

swLFP = swthLFP(:,1);
thLFP = swthLFP(:,2);
t_LFP = [1:length(swLFP)]./Fs;

%% Find Inverted PC1s and flip them for plot
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
        imagesc(log2(FFTfreqs),1:numusedchannels,pc1coeff(:,dipsortSW)')
        ylabel('Channel #');xlabel('f (Hz)')
        LogScale('x',2)
        axis xy
        title('PC1 Frequency Coefficients: All Channels') 
    subplot(2,2,2)
        imagesc(histbins,1:numusedchannels,pc1hists(:,dipsortSW)')
        ylabel('Channel #');xlabel('PC1 projection weight')
        title('PC1 Projection Histogram: All Channels')
        axis xy
    subplot(2,2,3)
        set(gca,'ColorOrder',RainbowColors(length(dipsortSW)))
        hold all
        plot(log2(FFTfreqs),pc1coeff')  
        plot(log2(FFTfreqs),pc1coeff(:,goodSWidx)','k','LineWidth',1)
        plot(log2(FFTfreqs([1 end])),[0 0],'k')
        ylabel('PC1 Coefficient');xlabel('f (Hz)')
        LogScale('x',2)
        title('PC1 Frequency Coefficients: All Channels')
    subplot(2,2,4)
        set(gca,'ColorOrder',RainbowColors(length(dipsortSW)))
        hold all
        plot(histbins,pc1hists')
        plot(histbins,pc1hists(:,goodSWidx)','k','LineWidth',1)
        ylabel('hist');xlabel('PC1 projection weight')
        title('PC1 Projection Histogram: All Channels')
        
saveas(swfig,[figfolder,recordingname,'_FindBestSW'],'jpeg')

%% Theta Hist and Coefficients

thfig = figure;
    subplot(2,2,1)
        imagesc(log2(thFFTfreqs),1:numusedchannels,THmeanspec(:,dipsortTH)')
        ylabel('Channel #');xlabel('f (Hz)')
        LogScale('x',2)
        axis xy
        title('Spectrum: All Channels') 
    subplot(2,2,2)
        imagesc(histbins,1:numusedchannels,THhist(:,dipsortTH)')
        ylabel('Channel #');xlabel('PC1 projection weight')
        title('Theta Ratio Histogram: All Channels')
        axis xy
    subplot(2,2,3)
        set(gca,'ColorOrder',RainbowColors(length(dipsortTH)))
        hold all
        plot(log2(thFFTfreqs),THmeanspec')  
        plot(log2(thFFTfreqs),THmeanspec(:,goodTHidx)','k','LineWidth',1)
        ylabel('Power');xlabel('f (Hz)')
        LogScale('x',2)
        title('Spectrum: All Channels')
    subplot(2,2,4)
        set(gca,'ColorOrder',RainbowColors(length(dipsortTH)))
        hold all
        plot(histbins,THhist')
        plot(histbins,THhist(:,goodTHidx)','k','LineWidth',1)
        ylabel('hist');xlabel('PC1 projection weight')
        title('Theta Ratio Histogram: All Channels') 
        
saveas(thfig,[figfolder,recordingname,'_FindBestTH'],'jpeg')
%% Show Channels


    %Calculate PC1 for plot/return
    [FFTspec,FFTfreqs,t_FFT] = spectrogram(allLFP(:,goodSWidx),window,noverlap,freqlist,Fs);
    FFTspec = abs(FFTspec);
    [zFFTspec,mu,sig] = zscore(log10(FFTspec)');

    totz = zscore(abs(sum(zFFTspec')));
    badtimes = find(totz>5);
    zFFTspec(badtimes,:) = 0;
    
     %[COEFF, SCORE, LATENT] = pca(zFFTspec);
    %broadbandSlowWave = SCORE(:,1);
     broadbandSlowWave = zFFTspec*SWweights';
     broadbandSlowWave = smooth(broadbandSlowWave,smoothfact);
    broadbandSlowWave = (broadbandSlowWave-min(broadbandSlowWave))./max(broadbandSlowWave-min(broadbandSlowWave));

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
        plot(t_FFT,broadbandSlowWave,'k')
        xlim(t_FFT([1,end]))
     
    %Calculate Theta ratio for plot/return    
    [thFFTspec,thFFTfreqs,t_FFT] = spectrogram(allLFP(:,goodTHidx),window,noverlap,thfreqlist,Fs);
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
 %   plot(allLFP(:,1),allLFP(:,goodSWidx),'k')
    
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
        
saveas(chanfig,[figfolder,recordingname,'_SWTHChannels'],'jpeg')
end

