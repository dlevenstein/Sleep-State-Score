%CompareRecordingPC1.m


datafolder = '/Users/dlevenstein/Dropbox/Research/Datasets/BWData/';
%figfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/TransitionsAndSpindles/AnalysisFigsMats/Spindle LFP/';
%Make List of all Recordings in data folder
recordings = dir(datafolder);
recordings(strncmpi('.',{recordings.name},1)) = [];   %Remove entries starting with .,~
recordings(strncmpi('~',{recordings.name},1)) = []; 

numrecs = length(recordings);


%%
for r = 1:numrecs;
 %%
 r = 1;
        close all
display(['Recording ',num2str(r),' of ',num2str(numrecs)])


LFPdata = [datafolder,recordings(r).name,'/',...
    recordings(r).name,'_ThetaLFP.mat'];
load(LFPdata)
thLFP = LFP;

LFPdata = [datafolder,recordings(r).name,'/',...
    recordings(r).name,'_LFP.mat'];
load(LFPdata)
sf_LFP = 1250;


EMGdata = [datafolder,recordings(r).name,'/',...
    recordings(r).name,'_EMGCorr.mat'];
load(EMGdata)
sf_EMG = 2;


%% Only take LFP, EMG from good sleep interval
% goodsleepmat = [datafolder,recordings(r).name,'/',...
%     recordings(r).name,'_WSWEpisodes.mat'];
% load(goodsleepmat)
% goodtime_s = [Start(GoodSleepInterval,'s') End(GoodSleepInterval,'s')];
% goodtime = round((goodtime_s*sf_LFP))+[1 0];
goodtime = [1 Inf]; %Entire Recording for all recordings

if goodtime(2)==Inf || goodtime(2)>length(LFP)
    goodtime(2) = length(LFP);
end

if length(goodtime) ~= 2
    display('Check goodtime')
    pause
end

LFP = LFP(goodtime(1):goodtime(2));
thLFP = thLFP(goodtime(1):goodtime(2));

%EMG = EMGCorr((EMGCorr(:,1)>goodtime_s(1) & EMGCorr(:,1)<goodtime_s(2)),2);
EMG = EMGCorr(:,2); %Entire Recording for all recordings

%%
%Save Figure Location
figloc = ['AllRecsStateID/',recordings(r).name];

[ INT, IDX, t_IDX ] = ClusterStates(LFP,thLFP,EMG,sf_LFP,sf_EMG)

states(states==3)=5;
states(states==2)=3;

% save([datafolder,recordings(r).name,'/',...
%     recordings(r).name,'_StateID.mat'],'stateintervals','states')
% save([datafolder,recordings(r).name,'/',...
%     recordings(r).name,'_StateIDNoMin.mat'],'stateintervals','states')
% save(['/Users/dlevenstein/Dropbox/dan&brendon/Sept15IntervalsDraftInf/MinImposed/',...
%     recordings(r).name,'_StateID.mat'],'stateintervals','states')
% save(['/Users/dlevenstein/Dropbox/dan&brendon/Sept8IntervalsDraft1/NoMins/',...
%     recordings(r).name,'_StateID.mat'],'stateintervals','states')


end

