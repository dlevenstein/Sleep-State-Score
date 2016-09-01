function [ deltapeaks,peakheights ] = DetectDELTA( LFP,sf,t_LFP,NREMints )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%INPUT
%   LFP
%   sf
%   t_LFP (optional)
%
%OUTPUT
%   deltapeaks
%   peakheights
%   DOWNoffset (to be added)
%   DOWNonset   (to be added)
%%
%t_LFP = (1:length(LFP))./sf;

deltarange = [0.5 4];
%deltarange = [1 4];
deltaLFP = FiltNPhase(LFP,deltarange,sf);
deltaLFP = NormToInt(deltaLFP,NREMints,sf,'Z');
% gammarange = [75 120];
% [~,gammapower] = FiltNPhase(LFP,gammarange,sf);
% gammapower = zscore(gammapower);
% 
% deltagamma = deltaLFP.*gammapower;

%% Delta peaks
peakthresh = 1.5;
peakthresh = 1.3;
peakdist = 0.1;
[peakheights,deltapeaks] = findpeaks(deltaLFP,t_LFP,'MinPeakHeight',peakthresh,'MinPeakDistance',peakdist);

%%
% xwin = [4378 4385];
% xwin = [1250 1265];
% %xwin = [100 5000];
% figure
%     subplot(4,1,1)
%         plot(t_LFP,LFP)
%         hold on
%         plot(deltapeaks,peakheights,'o')
%         xlim(xwin)
%     subplot(4,1,2)
%         plot(t_LFP,deltaLFP)
%         xlim(xwin)
% %     subplot(4,1,3)
% %         plot(t_LFP,gammapower)
% %         xlim(xwin)
% %     subplot(4,1,4)
% %         plot(t_LFP,deltagamma)
% %         xlim(xwin)
% 
% %%
% figure
%     plot(log2(deltapeaks(2:end-1)-deltapeaks(1:end-2)),log2(deltapeaks(3:end)-deltapeaks(2:end-1)),'.')
%     LogScale_ss('xy',2)
end

