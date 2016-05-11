function [ INT, IDX, t_IDX,PC1weights,PC1expvar ] = ClusterStates(LFP,thLFP,EMG,sf_LFP,sf_EMG,figloc,WSEpisodes)
%StateID(LFP,thLFP,EMG,sf_LFP,sf_EMG,figloc,WSEpisodes)
%   Detailed explanation goes here
%
%
%
%Dependencies: IDXtoINT, INTtoIDX
%
%Last Updated: 1/31/16
%DLevenstein




%% Min Win Parameters

minSWS = 2;
minWnexttoREM = 2;
minWinREM = 2;       
minREMinW = 2;
minREM = 2;
minWAKE = 2;


%% Downsample and filter
%Make Downsample to niquest frequency

if sf_LFP == 1250
    downsamplefactor = 5;
else
    display('sf not 1250... if only you made this able to set its own downsample...')
    downsamplefactor = 2;
end
LFP = downsample(LFP,downsamplefactor);
thLFP = downsample(thLFP,downsamplefactor);
sf_LFP = sf_LFP/downsamplefactor;


%filtbounds = [0.5 120];
%display(['Filtering ',num2str(filtbounds(1)),'-',num2str(filtbounds(2)),' Hz...']);
%LFP = FiltNPhase(LFP, filtbounds, sf_LFP );


%% Calculate Spectrogram
display('FFT Spectrum for Broadband LFP')

freqlist = logspace(0,2,100);
%freqlist = linspace(0.5,55.5,111);
window = 10;
noverlap = 9;
window = window*sf_LFP;
noverlap = noverlap*sf_LFP;
[FFTspec,FFTfreqs,t_FFT] = spectrogram(LFP,window,noverlap,freqlist,sf_LFP);


FFTspec = abs(FFTspec);

%% Find TRANSIENTS and set to 0 for PCA

[zFFTspec,mu,sig] = zscore(log10(FFTspec)');
totz = zscore(abs(sum(zFFTspec')));

badtimes = find(totz>5);
zFFTspec(badtimes,:) = 0;


%% PCA
smoothfact = 10; %si_FFT
thsmoothfact = 15;
 [COEFF, SCORE, ~, ~, EXPLAINED] = pca(zFFTspec);
 SCORE(:,1) = smooth(SCORE(:,1),smoothfact);
SCORE(:,1) = (SCORE(:,1)-min(SCORE(:,1)))./max(SCORE(:,1)-min(SCORE(:,1)));

PC1weights = COEFF(:,1);
PC1expvar = EXPLAINED(1);
 
%% Calculate theta
display('FFT Spectrum for Theta')

% %NarrowbandTheta
f_all = [2 16];
f_theta = [5 10];
freqlist = logspace(log10(f_all(1)),log10(f_all(2)),100);

%ThetaDelta
% f_delta = [0.5 4];
% f_theta = [5 10];
% freqlist = linspace((f_delta(1)),(f_theta(2)),50);

[thFFTspec,thFFTfreqs] = spectrogram(thLFP,window,noverlap,freqlist,sf_LFP);


thFFTspec = (abs(thFFTspec));


allpower = sum(log10(thFFTspec),1);

thfreqs = find(thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
thpower = sum(log10(thFFTspec(thfreqs,:)),1);

%Delta Power
% defreqs = find(thFFTfreqs>=f_delta(1) & thFFTfreqs<=f_delta(2));
% depower = sum(log10(thFFTspec(defreqs,:)),1);


thratio = thpower./allpower;    %Narrowband Theta
%thratio = thpower./depower;    %Theta/Delta
thratio = smooth(thratio,thsmoothfact);
thratio = (thratio-min(thratio))./max(thratio-min(thratio));
 
%% EMG
dtEMG = 1/sf_EMG;
t_EMG = (1:length(EMG))*dtEMG;
EMG = smooth(EMG,smoothfact/dtEMG);
EMG = (EMG-min(EMG))./max(EMG-min(EMG));

reclength = round(t_EMG(end));

%downsample to FFT time points;
[~,t_intersect] = intersect(t_EMG,t_FFT);
EMG = EMG(t_intersect);
t_EMG = t_EMG(t_intersect);



%% Divide PC1 for SWS
numbins = 10;
%numbins = 12; %for Poster...
[pcahist,histbins]= hist(SCORE(:,1),numbins);
 
[PKS,LOCS] = findpeaks(pcahist,'NPeaks',2,'SortStr','descend');
LOCS = sort(LOCS);

betweenpeaks = histbins(LOCS(1):LOCS(2));
[dip,diploc] = findpeaks(-pcahist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

thresh = betweenpeaks(diploc);

%Set transients to wake state
SCORE(badtimes,1)=histbins(LOCS(1));
 
 
%SWS time points
SWStimes = (SCORE(:,1) >thresh);


%% Then Divide EMG
numpeaks = 1;
numbins = 12; %for Poster...
while numpeaks ~=2
    [EMGhist,EMGhistbins]= hist(EMG(SWStimes==0),numbins);
    %[EMGhist,EMGhistbins]= hist(EMG,numbins);

    [PKS,LOCS] = findpeaks([0 EMGhist],'NPeaks',2);
    LOCS = sort(LOCS)-1;
    numbins = numbins+1;
    numpeaks = length(LOCS);
end

betweenpeaks = EMGhistbins(LOCS(1):LOCS(2));
[dip,diploc] = findpeaks(-EMGhist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

EMGthresh = betweenpeaks(diploc);

MOVtimes = (SCORE(:,1)<thresh & EMG>EMGthresh);


%% Then Divide Theta
numpeaks = 1;
numbins = 10;
numbins = 12; %for Poster...
while numpeaks ~=2 && numbins <=20
    %[THhist,THhistbins]= hist(thratio(SWStimes==0 & MOVtimes==0),numbins);
    [THhist,THhistbins]= hist(thratio(MOVtimes==0),numbins);

    [PKS,LOCS] = findpeaks(THhist,'NPeaks',2,'SortStr','descend');
    LOCS = sort(LOCS);
    numbins = numbins+1;
    numpeaks = length(LOCS);
end

numbins = 10;
numbins = 15; %for Poster...
while numpeaks ~=2 && numbins <=20
    [THhist,THhistbins]= hist(thratio(SWStimes==0 & MOVtimes==0),numbins);

    [PKS,LOCS] = findpeaks(THhist,'NPeaks',2,'SortStr','descend');
    LOCS = sort(LOCS);
    numbins = numbins+1;
    numpeaks = length(LOCS);
end

if length(PKS)==2
    betweenpeaks = THhistbins(LOCS(1):LOCS(2));
    [dip,diploc] = findpeaks(-THhist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

    THthresh = betweenpeaks(diploc);

    REMtimes = (SCORE(:,1)<thresh & EMG<EMGthresh & thratio>THthresh);
else
    THthresh = 0;
    REMtimes =(SCORE(:,1)<thresh & EMG<EMGthresh);
end

%%
%Index Vector: SWS=2, REM=3, MOV=6, NonMOV=1.   
%(Separate MOV for REM, then join later)
%IDX = SWStimes+2*REMtimes+5*MOVtimes+1;
%No separation of MOV and NonMOV WAKE
IDX = SWStimes+2*REMtimes+1;

%Start/end offset due to FFT


%% Minimum Interuptions
INT = IDXtoINT(IDX);


%Make the following repeated chunks of code into a single function.

%SWS  (to NonMOV)
Sints = INT{2};
Slengths = Sints(:,2)-Sints(:,1);
shortSints = {Sints(find(Slengths<=minSWS),:)};
shortSidx = INTtoIDX(shortSints,length(IDX));
%Change Short SWS to Wake
IDX(shortSidx==1) = 1;   
INT = IDXtoINT(IDX);

%NonMOV next to REM   (to REM)
Wints = INT{1};
trans = (diff(IDX)); %All State Transitions
WRtrans = find((trans)==-2)+1;  %Just transitions between WAKE, REM
%Convert to interval indices
[~,WRtransON] = intersect(Wints(:,1),WRtrans);
WRtrans = find((trans)==2);
[~,WRtransOFF] = intersect(Wints(:,2),WRtrans);
WRtrans = union(WRtransON,WRtransOFF); %On or offset are RW
%Find WAKE intervals that border REM and are less than min
Wlengths = Wints(:,2)-Wints(:,1);
shortWRints = find(Wlengths(WRtrans)<=minWnexttoREM);
shortWRints = WRtrans(shortWRints);
shortWRints = {Wints(shortWRints,:)};
shortWRidx = INTtoIDX(shortWRints,length(IDX));
%Convert wake to rem
IDX(shortWRidx==1) = 3;
INT = IDXtoINT(IDX);


%NonMOV in REM   (to REM)
Wints = INT{1};
trans = (diff(IDX)); %All State Transitions
WRtrans = find((trans)==-2)+1;  %Just transitions between WAKE, REM
%Convert to interval indices
[~,WRtransON] = intersect(Wints(:,1),WRtrans);
WRtrans = find((trans)==2);
[~,WRtransOFF] = intersect(Wints(:,2),WRtrans);
WRtrans = intersect(WRtransON,WRtransOFF); %Both onset and offset are RW
%Find WAKE intervals that border REM and are less than min
Wlengths = Wints(:,2)-Wints(:,1);
shortWRints = find(Wlengths(WRtrans)<=minWinREM);
shortWRints = WRtrans(shortWRints);
shortWRints = {Wints(shortWRints,:)};
shortWRidx = INTtoIDX(shortWRints,length(IDX));
%Convert wake to rem
IDX(shortWRidx==1) = 3;
IDX(IDX==6) = 1; %Convert NonMOV to WAKE
INT = IDXtoINT(IDX);


%REM in WAKE   (to WAKE)
Rints = INT{3};
trans = (diff(IDX)); %All State Transitions
WRtrans = find((trans)==2)+1;  %Just transitions between WAKE, REM
%Convert to interval indices
[~,WRtransON] = intersect(Rints(:,1),WRtrans);
WRtrans = find((trans)==-2);
[~,WRtransOFF] = intersect(Rints(:,2),WRtrans);
WRtrans = intersect(WRtransON,WRtransOFF); %Both onset and offset are RW
%Find WAKE intervals that border REM and are less than min
Rlengths = Rints(:,2)-Rints(:,1);
shortWRints = find(Rlengths(WRtrans)<=minREMinW);
shortWRints = WRtrans(shortWRints);
shortWRints = {Rints(shortWRints,:)};
shortWRidx = INTtoIDX(shortWRints,length(IDX));
%Convert REM to WAKE
IDX(shortWRidx==1) = 1;
INT = IDXtoINT(IDX);




%REM (only applies to REM in the middle of SWS)    (to WAKE)
Rints = INT{3};
Rlengths = Rints(:,2)-Rints(:,1);
shortRints = {Rints(find(Rlengths<=minREM),:)};
shortRidx = INTtoIDX(shortRints,length(IDX));

IDX(shortRidx==1) = 1;
INT = IDXtoINT(IDX);


%WAKE   (to SWS)     essentiall a minimum MA time
Wints = INT{1};
Wlengths = Wints(:,2)-Wints(:,1);
shortWints = {Wints(find(Wlengths<=minWAKE),:)};
shortWidx = INTtoIDX(shortWints,length(IDX));
IDX(shortWidx==1) = 2;

INT = IDXtoINT(IDX);

%SWS  (to NonMOV)
Sints = INT{2};
Slengths = Sints(:,2)-Sints(:,1);
shortSints = {Sints(find(Slengths<=minSWS),:)};
shortSidx = INTtoIDX(shortSints,length(IDX));
%Change Short SWS to Wake
IDX(shortSidx==1) = 1;   
INT = IDXtoINT(IDX);




%% Pad time to match recording time
offset = t_FFT(1)-1;

INT = cellfun(@(x) x+offset,INT,'UniformOutput',false);


 %% Figure
 
 if exist('figloc','var')
 viewwin  =[t_FFT(1) t_FFT(end)];
 %viewwin  =[32000 34000];
%viewwin=[9000 11000];
figure
	subplot(8,1,[1:2])
        imagesc(t_FFT,log2(FFTfreqs),log10(FFTspec))
        axis xy
        set(gca,'YTick',(log2([1 2 4 8 16 32 64 128])))
        set(gca,'YTickLabel',{'1','2','4','8','16','32','64','128'})
        caxis([3.5 6.5])
        caxis([min(mu)-2.5*max(sig) max(mu)+2.5*max(sig)])
        xlim(viewwin)
        colorbar('east')
        ylim([log2(FFTfreqs(1)) log2(FFTfreqs(end))+0.2])
        ylabel({'LFP - FFT','f (Hz)'})
        title(figloc);
	subplot(8,1,3)
        imagesc(t_FFT,log2(thFFTfreqs),log10(thFFTspec))
        axis xy
        set(gca,'YTick',(log2([1 2 4 8 16 32 64 128])))
        set(gca,'YTickLabel',{'1','2','4','8','16','32','64','128'})
        %caxis([3.5 6.5])
        caxis([min(mu)-2.5*max(sig) max(mu)+2.5*max(sig)])
        xlim(viewwin)
        %colorbar('east')
        ylim([log2(thFFTfreqs(1)) log2(thFFTfreqs(end))+0.2])
        ylabel({'thLFP - FFT','f (Hz)'})
        set(gca,'XTickLabel',{})
        
    subplot(8,1,4)
        %plot(t_FFT,-IDX,'LineWidth',2)
        hold on
        plot(INT{1}',-1*ones(size(INT{1}))','k','LineWidth',8)
        plot(INT{2}',-2*ones(size(INT{2}))','b','LineWidth',8)
        plot(INT{3}',-3*ones(size(INT{3}))','r','LineWidth',8)
        xlim(viewwin)
        ylim([-4 0])
        set(gca,'YTick',[-3:-1])
        set(gca,'YTickLabel',{'REM','SWS','Wake/MA'})
        set(gca,'XTickLabel',{})
        
   	subplot(6,1,4)
        hold on
        plot(t_FFT,SCORE(:,1),'k')
        %plot(synchtimes',thresh*ones(size(synchtimes))','r')
        ylabel('PC1')
        xlim([t_FFT(1) t_FFT(end)])
        xlim(viewwin)
        
   	subplot(6,1,5)
        hold on
        plot(t_FFT,thratio,'k')
        %plot(synchtimes',thresh*ones(size(synchtimes))','r')
        ylabel('Theta')
        xlim([t_FFT(1) t_FFT(end)])
        xlim(viewwin)
        
   	subplot(6,1,6)
        hold on
        plot(t_EMG,EMG,'k')
        %plot(synchtimes',thresh*ones(size(synchtimes))','r')
        ylabel('EMG')
        xlim([t_FFT(1) t_FFT(end)])
        xlim(viewwin)
        xlabel('t (s)')
        
	saveas(gcf,[figloc,'_ClusterResults'],'jpeg')
        
%%
if exist('WSEpisodes','var')
for ep = 1:length(WSEpisodes)
epstart = Start(WSEpisodes{ep},'s');
epend = End(WSEpisodes{ep},'s');

viewwin = [epstart(1) epend(2)];
%  viewwin  =[2150 2200];
%   viewwin  =[15000 20000];
 
figure
	subplot(8,1,[1:2])
        imagesc(t_FFT,log2(FFTfreqs),log10(FFTspec))
        axis xy
        set(gca,'YTick',(log2([1 2 4 8 16 32 64 128])))
        set(gca,'YTickLabel',{'1','2','4','8','16','32','64','128'})
        caxis([3.5 6.5])
        caxis([min(mu)-2.5*max(sig) max(mu)+2.5*max(sig)])
        xlim(viewwin)
        %colorbar('east')
        ylim([log2(FFTfreqs(1)) log2(FFTfreqs(end))+0.2])
        ylabel({'LFP - FFT','f (Hz)'})
        title([figloc,' Wake-Sleep ',num2str(ep)]);
	subplot(8,1,3)
        imagesc(t_FFT,log2(thFFTfreqs),log10(thFFTspec))
        axis xy
        set(gca,'YTick',(log2([1 2 4 8 16 32 64 128])))
        set(gca,'YTickLabel',{'1','2','4','8','16','32','64','128'})
        %caxis([3.5 6.5])
        caxis([min(mu)-2.5*max(sig) max(mu)+2.5*max(sig)])
        xlim(viewwin)
        %colorbar('east')
        ylim([log2(thFFTfreqs(1)) log2(thFFTfreqs(end))+0.2])
        ylabel({'thLFP - FFT','f (Hz)'})
        set(gca,'XTickLabel',{})
        
    subplot(8,1,4)
        %plot(t_FFT,-IDX,'LineWidth',2)
        hold on
        plot(INT{1}',-1*ones(size(INT{1}))','k','LineWidth',8)
        plot(INT{2}',-2*ones(size(INT{2}))','b','LineWidth',8)
        plot(INT{3}',-3*ones(size(INT{3}))','r','LineWidth',8)
        xlim(viewwin)
        ylim([-4 0])
        set(gca,'YTick',[-3:-1])
        set(gca,'YTickLabel',{'REM','SWS','Wake/MA'})
        set(gca,'XTickLabel',{})
        
   	subplot(6,1,4)
        hold on
        plot(t_FFT,SCORE(:,1),'k')
        %plot(synchtimes',thresh*ones(size(synchtimes))','r')
        ylabel('PC1')
        xlim([t_FFT(1) t_FFT(end)])
        xlim(viewwin)
        
   	subplot(6,1,5)
        hold on
        plot(t_FFT,thratio,'k')
        %plot(synchtimes',thresh*ones(size(synchtimes))','r')
        ylabel('Theta')
        xlim([t_FFT(1) t_FFT(end)])
        xlim(viewwin)
        
   	subplot(6,1,6)
        hold on
        plot(t_EMG,EMG,'k')
        %plot(synchtimes',thresh*ones(size(synchtimes))','r')
        ylabel('EMG')
        xlim([t_FFT(1) t_FFT(end)])
        xlim(viewwin)
        xlabel('t (s)')
        
    saveas(gcf,[figloc,'_WS',num2str(ep)],'jpeg')    
end        
end


        
%% 
figure

    subplot(2,3,1)
        scatter(SCORE(:,1),thratio,3,IDX,'filled')
        xlabel('Broadband PC1');ylabel('Narrowband Theta')
    subplot(2,3,2)
        scatter(SCORE(:,1),EMG,3,IDX,'filled')
        xlabel('Broadband PC1');ylabel('EMG')
    subplot(2,3,3)
        scatter(thratio,EMG,3,IDX,'filled')
        xlabel('Narrowband Theta');ylabel('EMG')

    subplot(2,3,4)
        scatter(SCORE(SWStimes==0,1),thratio(SWStimes==0),3,IDX(SWStimes==0),'filled')
        xlabel('Broadband PC1');ylabel('Narrowband Theta')
    subplot(2,3,5)
        scatter(SCORE(SWStimes==0,1),EMG(SWStimes==0),3,IDX(SWStimes==0),'filled')
        xlabel('Broadband PC1');ylabel('EMG')
        title('non-nonREM only')
    subplot(2,3,6)
        %scatter(thratio(SWStimes==0,1),EMG(SWStimes==0,1),3,IDX(SWStimes==0),'filled')
        plot(thratio(SWStimes==0 & IDX==1,1),EMG(SWStimes==0 & IDX==1,1),'k.')
        hold on
        plot(thratio(SWStimes==0 & IDX==3,1),EMG(SWStimes==0 & IDX==3,1),'r.')
        xlabel('Narrowband Theta');ylabel('EMG')

%% Figure: Split REM/Arousal  
figure
	subplot(3,2,1)
        hold on
        bar(histbins(histbins>thresh),pcahist(histbins>thresh),'FaceColor','b','barwidth',0.9,'linewidth',1)
        bar(histbins(histbins<=thresh),pcahist(histbins<=thresh),'FaceColor',0.9*[1 1 1],'barwidth',0.9,'linewidth',1)
        plot([thresh thresh],[0 max(pcahist)],'r','LineWidth',1)
        xlabel('PC 1')
        title('Step 1: PCA for SWS')
        

	subplot(3,2,3)
        hold on
        bar(EMGhistbins(EMGhistbins>EMGthresh),EMGhist(EMGhistbins>EMGthresh),'FaceColor','k','barwidth',0.9,'linewidth',1)
        bar(EMGhistbins(EMGhistbins<=EMGthresh),EMGhist(EMGhistbins<=EMGthresh),'FaceColor',0.9*[1 1 1],'barwidth',0.9,'linewidth',1)
        plot([EMGthresh EMGthresh],[0 max(EMGhist)],'r','LineWidth',1)
        xlabel('EMG')
        title('Step 2: EMG for Muscle Tone')
	subplot(3,2,5)
        hold on
        bar(THhistbins(THhistbins>=THthresh),THhist(THhistbins>=THthresh),'FaceColor','r','barwidth',0.9,'linewidth',1)
        bar(THhistbins(THhistbins<THthresh),THhist(THhistbins<THthresh),'FaceColor','k','barwidth',0.9,'linewidth',1)
        plot([THthresh THthresh],[0 max(THhist)],'r','LineWidth',1)
        xlabel('Theta')
        title('Step 3: Theta for REM')
        
        
    subplot(2,2,2)
        plot(SCORE(IDX==2,1),EMG(IDX==2),'b.')
        hold on
        plot(SCORE(EMG>EMGthresh & IDX==1,1),EMG(EMG>EMGthresh & IDX==1),'k.')
        plot(SCORE(EMG<EMGthresh & IDX==1|IDX==3,1),EMG(EMG<EMGthresh & IDX==1|IDX==3),'.','Color',0.8*[1 1 1])
        plot(thresh*[1 1],get(gca,'ylim'),'r','LineWidth',1)
        plot(thresh*[0 1],EMGthresh*[1 1],'r','LineWidth',1)
        xlabel('Broadband PC1');ylabel('EMG')
	subplot(2,2,4)
        %scatter(thratio(SWStimes==0,1),EMG(SWStimes==0,1),3,IDX(SWStimes==0),'filled')
        plot(thratio(SWStimes==0 & IDX==1,1),EMG(SWStimes==0 & IDX==1,1),'k.')
        hold on
        plot(thratio(SWStimes==0 & IDX==3,1),EMG(SWStimes==0 & IDX==3,1),'r.')
        xlabel('Narrowband Theta');ylabel('EMG')
        plot(THthresh*[1 1],EMGthresh*[0 1],'r','LineWidth',1)
        plot([0 1],EMGthresh*[1 1],'r','LineWidth',1)

saveas(gcf,[figloc,'_clust2'],'jpeg')
%saveas(gcf,['/Users/dlevenstein/Code Library/SleepScoreDevelopment/StateScoreFigures/','ThetaEMGExample'],'jpeg')
%% Figure: Clustering
colormat = [[0 0 0];[0 0 1];[1 0 0]];
coloridx = colormat(IDX,:);

figure
    subplot(1,3,[2,3])
        hold all
        scatter3(SCORE(:,1),thratio,EMG,2,coloridx,'filled')
        %rotate3d
        view(133.7,18.8);
        grid on
        xlabel('Broadband PC1');ylabel('Narrowband Theta');zlabel('EMG')
      
	subplot(3,3,1)
        hold on
        bar(histbins,pcahist,'FaceColor','none','barwidth',0.9,'linewidth',2)
        plot([thresh thresh],[0 max(pcahist)],'r')
        xlabel('PC 1')
        title('Step 1: PCA for SWS')
	subplot(3,3,4)
        hold on
        bar(EMGhistbins,EMGhist,'FaceColor','none','barwidth',0.9,'linewidth',2)
        plot([EMGthresh EMGthresh],[0 max(EMGhist)],'r')
        xlabel('EMG')
        title('Step 2: EMG for Muscle Tone')
	subplot(3,3,7)
        hold on
        bar(THhistbins,THhist,'FaceColor','none','barwidth',0.9,'linewidth',2)
        plot([THthresh THthresh],[0 max(THhist)],'r')
        xlabel('Theta')
        title('Step 3: Theta for REM')
        
	saveas(gcf,[figloc,'_clust'],'jpeg')
%saveas(gcf,['/Users/dlevenstein/Code Library/SleepScoreDevelopment/StateScoreFigures/','clust'],'jpeg')    
  %% Figure: Duration Distributions
  Wints = INT{1};
  Wlengths = Wints(:,2)-Wints(:,1);
  Sints = INT{2};
  Slengths = Sints(:,2)-Sints(:,1);
  Rints = INT{3};
  Rlengths = Rints(:,2)-Rints(:,1);
  
  figure
    subplot(2,3,1)
        hist(log10(Wlengths),10)
        set(gca,'XTick',0:3)
        set(gca,'XTickLabel',10.^[0:3])
        xlabel('Duration (s)')
        title('Wake Interval Durations')
    subplot(2,3,2)
        hist(log10(Slengths),10)
        set(gca,'XTick',0:3)
        set(gca,'XTickLabel',10.^[0:3])
        xlabel('Duration (s)')
        title('SWS Interval Durations')
    subplot(2,3,3)
        hist(log10(Rlengths),10)
        set(gca,'XTick',0:3)
        set(gca,'XTickLabel',10.^[0:3])
        xlabel('Duration (s)')
        title('REM Interval Durations')
    subplot(2,3,4)
        plot(log10(Wlengths(1:end-1)),log10(Wlengths(2:end)),'.')
        set(gca,'YTick',0:3)
        set(gca,'YTickLabel',10.^[0:3])
        set(gca,'XTick',0:3)
        set(gca,'XTickLabel',10.^[0:3])
        xlabel('Interval n Duration')
        ylabel('Interval n+1 Duration')
        title('Wake Interval Durations')
    subplot(2,3,5)
        plot(log10(Slengths(1:end-1)),log10(Slengths(2:end)),'.')
        set(gca,'YTick',0:3)
        set(gca,'YTickLabel',10.^[0:3])
        set(gca,'XTick',0:3)
        set(gca,'XTickLabel',10.^[0:3])
        xlabel('Interval n Duration')
        ylabel('Interval n+1 Duration')
        title('SWS Interval Durations')
    subplot(2,3,6)
        plot(log10(Rlengths(1:end-1)),log10(Rlengths(2:end)),'.')
        set(gca,'YTick',0:3)
        set(gca,'YTickLabel',10.^[0:3])
        set(gca,'XTick',0:3)
        set(gca,'XTickLabel',10.^[0:3])
        xlabel('Interval n Duration')
        ylabel('Interval n+1 Duration')
        title('REM Interval Durations')
        
        saveas(gcf,[figloc,'_intdur'],'jpeg')
        
    

 end
 
IDX = INTtoIDX(INT,reclength);
t_IDX = 1:length(IDX);
IDX = IDX';

end

