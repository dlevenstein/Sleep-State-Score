function [ output_args ] = PickSWChannel( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% DEV
datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/DTData/';
recname = 'DT3_rLS_rCA1_20150927_298um_288um';

% datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/GGData/';
% recname = 'Rat08-20130717';
filename = [datasetfolder,'/',recname,'/',recname,'.xml'];

%% FMA

SetCurrentSession(filename);
global DATA
numchans = DATA.nChannels;

%% Hist/Freqs Parms
numhistbins = 21;
numfreqs = 100;

%% Load LFP files from .lfp
allLFP = GetLFP('all');

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
pc1coeff = zeros(numfreqs,numchans);
dip = zeros(numchans,1);
for cc = 1:numchans;
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
    bins = linspace(0,1,numhistbins);
    [pcahist,histbins]= hist(SCORE(:,1),bins);

    pc1hists(:,cc) = pcahist;
    pc1coeff(:,cc) = COEFF(:,1);
    
    dip(cc) = hartigansdiptest(sort(SCORE(:,1)));
end

%% Sort by dip
[~,dipsort] = sort(dip);

%% Find Inverted PC1s
invpc1 = mean(pc1coeff(freqlist<4,:))<0 & mean(pc1coeff(freqlist>50,:))>0;
pc1coeff(:,invpc1) = -pc1coeff(:,invpc1);
pc1hists(:,invpc1) = flipud(pc1hists(:,invpc1));
%% Test
%PC1 histogram for bimodality
%PC1 coefficients for NREM and proper orientation
%Theta for bimodality
%Theta spectrum for isolated peak?

%% PC1 Weights and Coefficients


figure
    subplot(2,2,1)
        imagesc(log2(FFTfreqs),1:numchans,pc1coeff(:,dipsort)')
        ylabel('Channel #');xlabel('f (Hz)')
        LogScale('x',2)
        axis xy
        title('PC1 Frequency Coefficients: All Channels') 
    subplot(2,2,2)
        imagesc(histbins,1:numchans,pc1hists(:,dipsort)')
        ylabel('Channel #');xlabel('PC1 projection weight')
        title('PC1 Projection Histogram: All Channels')
        axis xy
    subplot(2,2,3)
        set(gca,'ColorOrder',RainbowColors(length(dipsort)))
        hold all
        plot(log2(FFTfreqs),pc1coeff')     
        plot(log2(FFTfreqs([1 end])),[0 0],'k')
        ylabel('PC1 Coefficient');xlabel('f (Hz)')
        LogScale('x',2)
        title('PC1 Frequency Coefficients: All Channels')
    subplot(2,2,4)
        set(gca,'ColorOrder',RainbowColors(length(dipsort)))
        hold all
        plot(histbins,pc1hists')
        ylabel('hist');xlabel('PC1 projection weight')
        title('PC1 Projection Histogram: All Channels')

        

%%
exchan = 60;
exchan = 20;
    [FFTspec,FFTfreqs,t_FFT] = spectrogram(allLFP(:,exchan),window,noverlap,freqlist,sf_LFP);
    FFTspec = abs(FFTspec);
    [zFFTspec,mu,sig] = zscore(log10(FFTspec)');

figure
	subplot(2,1,1)
        imagesc(t_FFT,log2(FFTfreqs),log10(FFTspec))
        axis xy
        set(gca,'YTick',(log2([1 2 4 8 16 32 64 128])))
        set(gca,'YTickLabel',{'1','2','4','8','16','32','64','128'})
        caxis([3.5 6.5])
        caxis([min(mu)-2.5*max(sig) max(mu)+2.5*max(sig)])
        %xlim(viewwin)
        %colorbar('east')
        ylim([log2(FFTfreqs(1)) log2(FFTfreqs(end))+0.2])
        ylabel({'LFP - FFT','f (Hz)'})
        %title([figloc,' Wake-Sleep ',num2str(ep)]);
subplot(2,1,2)
    plot(allLFP(:,1),allLFP(:,exchan),'k')
end

