function StateIntervals = StatesToFinalScoring(NREMints,WAKEints,REMints,MAints)
% Takes a series of state start/stop intervals (each interval is a single
% row with a start column and a stop column... all in a single nx2 matrix)
% and finds Microarousals, WAKE and Episodes using hard-coded duraiton
% max/mins.  
% * NOTE: if MAints is entered as an input, those will override MAInts as
% defined here... ie if users want to input after manual review in a way
% that mantains users having final say over MA.  Also NOTE, even if MA Ints
% is input, only 
%   
% Dan Levenstein and Brendon Watson 2016

minPacketDuration = 30;
maxMicroarousalDuration = 100;

maxEpisodeDuration = 40;
minSWSEpisodeDuration = 20;
minWAKEEpisodeDuration = 20;
minREMEpisodeDuration = 20;


SWSlengths = NREMints(:,2)-NREMints(:,1);
packetintervals = NREMints(SWSlengths>=minPacketDuration,:);

WAKElengths = WAKEints(:,2)-WAKEints(:,1);

if ~exist('MAints','var')
    MAIntervals = WAKEints(WAKElengths<=maxMicroarousalDuration,:);
    WAKEIntervals = WAKEints(WAKElengths>maxMicroarousalDuration,:);
else
    MAIntervals = MAints;
    WAKEIntervals = WAKEints;
end

[episodeintervals{2}] = IDStateEpisode(NREMints,maxEpisodeDuration,minSWSEpisodeDuration);
[episodeintervals{1}] = IDStateEpisode(WAKEints,maxEpisodeDuration,minWAKEEpisodeDuration);
[episodeintervals{3}] = IDStateEpisode(REMints,maxEpisodeDuration,minREMEpisodeDuration);


%% Save
StateIntervals.NREMstate = NREMints;
StateIntervals.REMstate = REMints;
StateIntervals.WAKEstate = WAKEIntervals;
StateIntervals.NREMepisode = episodeintervals{2};
StateIntervals.REMepisode = episodeintervals{3};
StateIntervals.WAKEeposode = episodeintervals{1};
StateIntervals.NREMpacket = packetintervals;
StateIntervals.MAstate = MAIntervals;

