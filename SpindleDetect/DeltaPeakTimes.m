function [deltapeaks] = DeltaPeakTimes(ctxchannels,NREMint,figloc,recname);
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
%%
 [ deltapeaks,peakheights ] = DetectDELTA( allLFP(:,2),sf_LFP,allLFP(:,1) );

