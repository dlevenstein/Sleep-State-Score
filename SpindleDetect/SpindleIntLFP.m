function [ SpindleStats ] = SpindleIntLFP( pSpindleInts, cycletimemap,deltapeaks,ctxchannels,NREMint,figfolder,recordingname )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%

savebool = 1;


%
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
     meanLFP = mean(allLFP(:,2:end),2);
    % allLFP(:,3:end) = [];
end


%% fparms
frange = [6 20];
nfreqs = 50;
ncyc = 8;

%% Wavelets and Frequency of Max Power
[freqs,t,spec] = WaveSpec(meanLFP,frange,nfreqs,ncyc,1/sf_LFP,'lin');
spec = zscore(abs(spec)'.^2);
%spec = zscore(abs(spec)');

%Max Frequency
    [maxval,maxfreqidx] = max(spec,[],2);
    maxfreq = freqs(maxfreqidx);

 %% Filter Gamma   
garange = [80 150];
LFPspindle = FilterLFP(allLFP,'passband',garange,'nyquist',sf_LFP./2);
[~,ga_amp] = Phase(LFPspindle);  
ga_amp = zscore(ga_amp);
ga_amp(:,2) = max(ga_amp(:,2:end),[],2);
ga_amp(:,3:end) = [];
% %% Gamma Frequency of Max Power
% [gafreqs,t,gaspec] = WaveSpec(allLFP(:,2),garange,nfreqs,ncyc,1/sf_LFP,'lin');
% gaspec = zscore(abs(gaspec)'.^2);
% %spec = zscore(abs(spec)');
% 
% %Max Frequency
%     [~,maxgafreqidx] = max(gaspec,[],2);
%     maxgafreq = gafreqs(maxgafreqidx);
% %%
% cycletimeall = cat(1,cycletimemap{:});
% figure
%     plot(allLFP(:,1),ga_amp(:,2),'r')
%     hold on
%     plot(allLFP(:,1),meanLFP,'k')
%     plot(cycletimeall(:,1),mod(cycletimeall(:,2),2*pi),'b')
% %% Find Gamma Peaks
%     [cpeak,cpeakt] = findpeaks(amplitude(:,cc+1),amplitude(:,1),...
%         'MinPeakHeight',peakthresh,'MinPeakDistance',minpeakdist);

%% Mean Wavelet Power and Gamma power phase for each pSpindle
nphases = 25;
spphasebins = linspace(0,2*pi,nphases+1);
spphasebins = spphasebins(1:end-1)+0.5*(spphasebins(2)-spphasebins(1));

numspindles = length(pSpindleInts(:,1));
spwavepower = zeros(numspindles,nfreqs);
spgammaphase = zeros(numspindles,nphases);

for ss = 1:numspindles
   % spint = pSpindleInts(ss,:);
    spindex =  interp1(allLFP(:,1),1:length(allLFP(:,1)),cycletimemap{ss}(:,1),'nearest');
    spwavepower(ss,:) = mean(spec(spindex,:),1);
    
    spphases = interp1(spphasebins,spphasebins,mod(cycletimemap{ss}(:,2),2*pi),'nearest');
    gapowers = ga_amp(spindex,2);
    spLFPs = meanLFP(spindex);
    for pp = 1:nphases
        spgammaphase(ss,pp) = mean(gapowers(spphases==spphasebins(pp)));
        spLFPphase(ss,pp) = mean(spLFPs(spphases==spphasebins(pp)));
    end
end
 
%% PCA for mean wavelet power
% [COEFF, SCORE,~,~,EXP] = pca(spwavepower);
% 
% [~,sortPC1] = sort(SCORE(:,1));
% [~,sortPC2] = sort(SCORE(:,2));

%% 
% figure
%     subplot(2,2,3)
%         imagesc(freqs,[1 numspindles],spwavepower(sortPC1,:))
%         colorbar
%         caxis([0 4])
%         xlabel('Frequency (Hz)');ylabel('Spindle (Sorted by PC1)')
%     subplot(2,2,4)
%         imagesc(freqs,[1 numspindles],spwavepower(sortPC2,:))
%         colorbar
%         caxis([0 4])
%         xlabel('Frequency (Hz)');ylabel('Spindle (Sorted by PC2)')
%     subplot(2,2,2)
%         plot(SCORE(:,1),SCORE(:,2),'.')
%         xlabel('PC1 Score');ylabel('PC2 Score')
%     subplot(4,2,3)
%         plot(freqs,COEFF(:,[1 2]))
%         xlabel('Frequency (Hz)');ylabel('Coefficient')
%         xlim(freqs([1 end]))
%         legend('PC1','PC2','location','southeast')
% 	subplot(4,2,1)
%         plot(EXP,'o-')
%         xlim([0 5])

        
        %% PCA for gamma phase
% [COEFF, SCORE,~,~,EXP] = pca(spgammaphase);
% 
% [~,sortPC1] = sort(SCORE(:,1));
% [~,sortPC2] = sort(SCORE(:,2));

%% 
% figure
%     subplot(2,2,3)
%         imagesc([0 4*pi],[1 numspindles],[spgammaphase(sortPC1,:) spgammaphase(sortPC1,:)])
%         colorbar
%         caxis([0 4])
%         xlabel('Spindle Phase');ylabel('Spindle (Sorted by PC1)')
%     subplot(2,2,4)
%         imagesc([0 4*pi],[1 numspindles],[spgammaphase(sortPC2,:) spgammaphase(sortPC2,:)])
%         colorbar
%         caxis([0 4])
%         xlabel('Spindle Phase');ylabel('Spindle (Sorted by PC2)')
%     subplot(2,2,2)
%         plot(SCORE(:,1),SCORE(:,2),'.')
%         xlabel('PC1 Score');ylabel('PC2 Score')
%     subplot(4,2,3)
%         plot(spphasebins,COEFF(:,[1 2]))
%         xlabel('Spindle Phase');ylabel('Coefficient')
%         xlim(spphasebins([1 end]))
%         legend('PC1','PC2','location','southeast')
% 	subplot(4,2,1)
%         plot(EXP,'o-')
%         xlim([1 5])
%         ylabel('Exp. Var.')
%         xlabel('PC #')
%         
        
%% pSpindle "Characteristic Frequency" and "Gamma Phase Peak"

[maxgammapower,maxgammaphase] = max(spgammaphase,[],2);
maxgammaphase = spphasebins(maxgammaphase);
[maxfreqpower,characteristicfreq] = max(spwavepower,[],2);
characteristicfreq = freqs(characteristicfreq);

[~,sortcharfreq] = sort(characteristicfreq);
[~,sortmaxpower] = sort(maxfreqpower);

[~,sortgammaphase] = sort(mod(maxgammaphase-pi,2*pi));
[~,sortgammapower] = sort(maxgammapower);

%% 
% figure
%     subplot(2,2,3)
%         imagesc(freqs,[1 numspindles],spwavepower(sortcharfreq,:))
%         colorbar
%         caxis([0 4])
%         xlabel('Frequency (Hz)');ylabel('Spindle (Sorted by PC1)')
%     subplot(2,2,4)
%         imagesc(freqs,[1 numspindles],spwavepower(sortmaxpower,:))
%         colorbar
%         caxis([0 4])
%         xlabel('Frequency (Hz)');ylabel('Spindle (Sorted by PC2)')
%     subplot(2,2,2)
%         plot(SCORE(:,1),SCORE(:,2),'.')
%         xlabel('PC1 Score');ylabel('PC2 Score')
%     subplot(4,2,3)
%         plot(freqs,COEFF(:,[1 2]))
%         xlabel('Frequency (Hz)');ylabel('Coefficient')
%         xlim(freqs([1 end]))
%         legend('PC1','PC2','location','southeast')
% 	subplot(4,2,1)
%         plot(EXP,'o-')
%         xlim([0 5])
% 
% 
% %% 
% figure
%     subplot(2,2,3)
%         imagesc([0 4*pi],[1 numspindles],[spgammaphase(sortgammaphase,:) spgammaphase(sortgammaphase,:)])
%         colorbar
%         caxis([0 4])
%         xlabel('Spindle Phase');ylabel('Spindle (Sorted by PC1)')
%     subplot(2,2,4)
%         imagesc([0 4*pi],[1 numspindles],[spgammaphase(sortgammapower,:) spgammaphase(sortgammapower,:)])
%         colorbar
%         caxis([0 4])
%         xlabel('Spindle Phase');ylabel('Spindle (Sorted by PC2)')
%     subplot(2,2,2)
%         plot(SCORE(:,1),SCORE(:,2),'.')
%         xlabel('PC1 Score');ylabel('PC2 Score')
%     subplot(4,2,3)
%         plot(spphasebins,COEFF(:,[1 2]))
%         xlabel('Spindle Phase');ylabel('Coefficient')
%         xlim(spphasebins([1 end]))
%         legend('PC1','PC2','location','southeast')
% 	subplot(4,2,1)
%         plot(EXP,'o-')
%         xlim([1 5])
%         ylabel('Exp. Var.')
%         xlabel('PC #') 
 
        
%% Delta Peaks
[delta_spnorm,~,~,overlaps] = SortByIntTime(deltapeaks,pSpindleInts,'norm');
delta_spnorm = [delta_spnorm' overlaps{:}];
[delta_spoffset,~,~,overlaps] = SortByIntTime(deltapeaks,[pSpindleInts(1:end-1,2) pSpindleInts(2:end,1)],'onset');
delta_spoffset = [delta_spoffset' overlaps{:}];
[delta_sponset,~,~,overlaps] = SortByIntTime(deltapeaks,[pSpindleInts(1:end-1,2) pSpindleInts(2:end,1)],'offset');
delta_sponset = [delta_sponset' overlaps{:}];

afterhist = hist(delta_spoffset(delta_spoffset>0 &delta_spoffset<1.5),15);
beforehist = hist(delta_sponset(delta_sponset>-1.5 &delta_sponset<0),15);
duringhist = hist(delta_spnorm); 


[sponset_beforedelta,~,~,overlaps] = SortByIntTime(pSpindleInts(:,1),[deltapeaks(1:end-1) deltapeaks(2:end)],'offset');
sponset_beforedelta = [sponset_beforedelta' overlaps{:}];
beforedeltahist_on = hist(sponset_beforedelta(sponset_beforedelta>-1.5 &sponset_beforedelta<0),15);     

[spoffset_beforedelta,~,~,overlaps] = SortByIntTime(pSpindleInts(:,2),[deltapeaks(1:end-1) deltapeaks(2:end)],'offset');
spoffset_beforedelta = [spoffset_beforedelta' overlaps{:}];
beforedeltahist_off = hist(spoffset_beforedelta(spoffset_beforedelta>-1.5 &spoffset_beforedelta<0),15);

[sponset_afterdelta,~,~,overlaps] = SortByIntTime(pSpindleInts(:,1),[deltapeaks(1:end-1) deltapeaks(2:end)],'onset');
sponset_afterdelta = [sponset_afterdelta' overlaps{:}];
afterdeltahist_on = hist(sponset_afterdelta(sponset_afterdelta>0 &sponset_afterdelta<1.5),15);     

[spoffset_afterdelta,~,~,overlaps] = SortByIntTime(pSpindleInts(:,2),[deltapeaks(1:end-1) deltapeaks(2:end)],'onset');
spoffset_afterdelta = [spoffset_afterdelta' overlaps{:}];
afterdeltahist_off = hist(spoffset_afterdelta(spoffset_afterdelta>0 &spoffset_afterdelta<1.5),15);



%% real time to indextime
t_idx = round(t*sf_LFP+1);
SpIntIdx = interp1(allLFP(:,1),t_idx,pSpindleInts,'nearest');

%%
% plotparms.xview = [0.5 1]*sf_LFP;
%  [epochs,intlengths] = IntervalPETH(meanLFP,SpIntIdx,100,1,plotparms)
%  
 
%%
% timemap = cat(1,cycletimemap{:});
% timemap = downsample(timemap,1);
% 
% timemap_LFP = interp1(allLFP(:,1),meanLFP,timemap(:,1));
% timemap_maxfreq = interp1(allLFP(:,1),maxfreq,timemap(:,1));
% timemap_gaamp = interp1(allLFP(:,1),ga_amp(:,2),timemap(:,1));
% %timemap_gafreq = interp1(allLFP(:,1),maxgafreq,timemap(:,1));


%% time map smooth

%[corrbin,bincenters] = PairMatHist(timemap_maxfreq,[timemap(:,2)./(2*pi) timemap(:,3)*10],150,[0 10]);
% [corrbin,bincenters] = PairMatHist(timemap_gaamp,[timemap(:,2)./(2*pi) timemap(:,3)*10],150,[0 10]);
%[corrbin,bincenters] = PairMatHist(timemap_LFP,[timemap(:,2)./(2*pi) timemap(:,3)*10],150,[0 10]);

%% Heatmap


% figure
%     imagesc(bincenters,bincenters./10,corrbin.mean')
%     peaklinex = [1:10; 1:10]';
%     peakliney = bsxfun(@(X,Y) X.*Y,get(gca,'ylim'),ones(size(peaklinex)));
%     hold on
%     plot(peaklinex',peakliney','w')
%     xlabel('Cycle-Normalized Time');ylabel('Spindle-Normalized Time')
%     title('Gamma Power')
% 
%     axis xy
%    colorbar
%   caxis([0 3])
 %% Trajectories
% figure
%     subplot(2,2,1)
%         scatter(timemap(:,2)./(2*pi),timemap(:,3),0.5*ones(size(timemap_LFP)),timemap_LFP)
%         xlabel('Cycle-Normalized Time');ylabel('Spindle-Normalized Time')
%         title('Raw LFP with Spindle Time')
%         xlim([0 10]);ylim([0 1])
%     subplot(2,2,2)
%         scatter(timemap(:,2)./(2*pi),timemap(:,3),0.2*ones(size(timemap_LFP)),timemap_maxfreq)
%         xlabel('Cycle-Normalized Time');ylabel('Spindle-Normalized Time')
%         title('Max Frequency with Spindle Time')
%         xlim([0 10]);ylim([0 1])
%       %  colorbar
%     subplot(2,2,3)
%         scatter(timemap(:,2)./(2*pi),timemap(:,3),1*ones(size(timemap_LFP)),timemap_gaamp)
%         xlabel('Cycle-Normalized Time');ylabel('Spindle-Normalized Time')
%         title('Gamma Amplitude with Spindle Time')
%         xlim([0 10]);ylim([0 1])
%        % colorbar
%         caxis([-0.5 1.5])
% 
%         
% %% Gamma Power with Spindle Phase
% [N,C] = hist3([mod(timemap(:,2),2*pi),timemap_gaamp],[100 100]);
% figure
%     subplot(2,2,1)
%         plot([mod(timemap(:,2),2*pi) mod(timemap(:,2),2*pi)+2*pi],...
%             [timemap_gaamp timemap_gaamp],'k.')
%     subplot(2,2,2)
%         imagesc([C{1} C{1}+2*pi],C{2},[N; N]')
%         xlabel('Spindle Phase');ylabel('Gamma Power (Z Score)')
%         axis xy
%         

%% Peaks
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

        
%% Durations: Absolute and Cycle time
spdur = pSpindleInts(:,2)-pSpindleInts(:,1);
spcycledur = cellfun(@(X) max(X(:,2)),cycletimemap)./(2*pi);
        
%% MASTER SUMMARY FIGURE
splfpfig = figure;
    subplot(3,3,1)
        imagesc(freqs,[1 numspindles],spwavepower(sortcharfreq,:))
        %colorbar
        caxis([0 4])
        %xlabel('Frequency (Hz)');
        ylabel('Sort: char. frequency')
        title({'All Spindles:', 'Spindle Power'})
    subplot(3,3,4)
        imagesc(freqs,[1 numspindles],spwavepower(sortmaxpower,:))
        %colorbar
        caxis([0 4])
        %xlabel('Frequency (Hz)');
        ylabel('Sort: tot. power')
	subplot(6,3,13)
        plot(freqs,mean(spwavepower,1))
        xlim(freqs([1 end]))
        xlabel('f (Hz)');ylabel('Mean Power')
    subplot(3,3,2)
        imagesc([0 4*pi],[1 numspindles],[spgammaphase(sortgammaphase,:) spgammaphase(sortgammaphase,:)])
        %colorbar
        caxis([0 4])
        %xlabel('Spindle Phase');
        ylabel('Sort: gamma phase')
        title({'All Spindles:', 'Gamma (80-150Hz) Power'})
    subplot(3,3,5)
        imagesc([0 4*pi],[1 numspindles],[spgammaphase(sortgammapower,:) spgammaphase(sortgammapower,:)])
        %colorbar
        caxis([0 4])
        %xlabel('Spindle Phase');
        ylabel('Sort: gamma power')
    subplot(6,3,14)      
        plot([spphasebins spphasebins+2*pi],[mean(spLFPphase,1), mean(spLFPphase,1)]./3+1,'k')
        hold on
        plot([spphasebins spphasebins+2*pi],[mean(spgammaphase,1), mean(spgammaphase,1)],'LineWidth',1)
        xlim([0 4*pi])
        xlabel('Spindle Phase');ylabel('Mean Power')     

	subplot(6,3,3)
        hist([1./IPI_pos; 1./IPI_neg()],50)
        title('Within-Spindle Peaks')
	
    subplot(6,3,[6 9])
        plot(1./IPI_pos(1:end-1),1./IPI_pos(2:end),'.','color',poscolor)
        hold on
        plot(1./IPI_neg(1:end-1),1./IPI_neg(2:end),'.','color',negcolor)
        
        xlabel('IPI n (Hz)');ylabel('IPI n+1 (Hz)')
        
	subplot(6,3,12)
        imagesc(freqbins{1},freqbins{2},freqhist')
        hold on
        ylabel('IPI (Hz)');xlabel('Norm. Spindle Time')
        %plot(bincenters,binmeans,'wo-')
        axis xy
        %ylim([6 18])  
        
    subplot(6,3,18)
        bar(linspace(-1.5,2.5,40),[beforehist,duringhist,afterhist])
        hold on
        plot([0 0],get(gca,'ylim'),'r','LineWidth',1)
        plot([1 1],get(gca,'ylim'),'r','LineWidth',1)
        set(gca,'XTick',[-1 0 1 2])
        set(gca,'XTickLabel',{'-1s','S','E','+1s'})
        xlabel('Spindle-Relative Time');ylabel('Delta Peaks')
        xlim([-1.5 2.5])

     maxcycles =ceil(max(spcycledur));
     mincycles =floor(min(spcycledur));
     
    subplot(6,3,17)
        bar(linspace(-1.5,1.5,30),[beforedeltahist_on,afterdeltahist_on;beforedeltahist_off,afterdeltahist_off]',...
            'stacked')
        hold on
        plot([0 0],get(gca,'ylim'),'r','LineWidth',1)
        xlabel('Offset from Delta Peak');ylabel('Spindle S/E')
        xlim([-1.5 1.5])

     maxcycles =ceil(max(spcycledur));
     mincycles =floor(min(spcycledur));
     
     
     
	subplot(6,6,29)
        hist(spdur,maxcycles-mincycles)
        xlabel('Spindle Duration (s)')

	subplot(6,6,30)
        hist(spcycledur,maxcycles-mincycles)
        xlim([mincycles maxcycles])
        xlabel('Spindle Duration (cyc)')
        
saveas(splfpfig,[figfolder,recordingname,'_SWSpindleProps'],'jpeg')

%%
SpindleStats.deltaCCGlag = linspace(-1.5,1.5,30);
SpindleStats.deltaCCG = [beforedeltahist_on,afterdeltahist_on;beforedeltahist_off,afterdeltahist_off]';
SpindleStats.spindledeltaCCGlag = linspace(-1.5,2.5,40),[beforehist,duringhist,afterhist];
SpindleStats.spindledeltaCCG = [beforehist,duringhist,afterhist];
SpindleStats.spphasebins = spphasebins;
SpindleStats.spLFPphase = spLFPphase;
SpindleStats.spgammaphase = spgammaphase;
SpindleStats.spfreqs = freqs;
SpindleStats.spwavepower = spwavepower;
SpindleStats.characteristicfreq = characteristicfreq;
SpindleStats.maxfreqpower = maxfreqpower;
SpindleStats.maxgammaphase = maxgammaphase;
SpindleStats.maxgammapower = maxgammapower;



end

