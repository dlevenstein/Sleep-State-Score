function [ spints ] = JoinSpindlePeaks( joinedpeaks,figloc,recname )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Join Peaks/Troughs to Spindles
distthresh = 0.2;
[spints] = MergeSeparatedInts([joinedpeaks joinedpeaks],distthresh);
spdur = spints(:,2)-spints(:,1);

spmaxdur = 3;
spmindur = 0.3;

spints(spdur>spmaxdur | spdur<spmindur,:) = [];
spdur(spdur>spmaxdur | spdur<spmindur) = [];

ISpI = spints(2:end,1)-spints(1:end-1,2);

% inspindlepeaks = joinedpeaks;
% [inspindlepeaks,spindlepeakidx] = RestrictInts(joinedpeaks,spints);
% inspindlepeaks = inspindlepeaks(:,1);
% inspindlepeaksign = joinedpeakssign(spindlepeakidx);
% speak_pos = inspindlepeaks(inspindlepeaksign==1);
% speak_neg = inspindlepeaks(inspindlepeaksign==-1);
% IspeakI_pos = speak_pos(2:end)-speak_pos(1:end-1);
% IspeakI_neg = speak_neg(2:end)-speak_neg(1:end-1);
% IspeakI_all = [IspeakI_pos;IspeakI_neg];



%% Figure: Spindle Stats
figure
    subplot(4,2,1)
        hist(spdur,20)
        hold on
        text(0.75,40,['Total # pSpindles: ',num2str(length(spdur))])
        xlabel('pSpindle Duration (s)');
        ylabel('# pSpindles')
    subplot(4,2,3)
        hist(log10(ISpI),25)
        xlabel('Inter-pSpindle Interval (s)');
        ylabel('# pSpindles')
        LogScale('x',10)
       % xlim([-1 2])
    subplot(2,2,2)
        plot(log10(ISpI(1:end-1)),log10(ISpI(2:end)),'.')
        xlabel('Inter-pSpindle Interval n (s)')
        ylabel('Inter-pSpindle Interval n+1 (s)')
        %title('Iterative Map: InterSpindle Intervals')
        LogScale('xy',10)
        %xlim([-1 2]);ylim([-1 2])
%     subplot(2,2,4)
%         plot(log10(IspeakI_pos(1:end-1)),log10(IspeakI_pos(2:end)),'k.')
%         hold on
%         plot(log10(IspeakI_neg(1:end-1)),log10(IspeakI_neg(2:end)),'k.')
%         xlabel('Inter-Peak Interval n')
%         ylabel('Inter-Peak Interval n+1')
%         LogScale('xy',10)
%         xlim([-1.5 2]);ylim([-1.5 2])
%         title('Iterative Map: Within-Spindle Peaks')
%     subplot(4,2,2)
%         hist(1./IspeakI_all,100)
%         xlim([6 20])
%         xlabel('Instantaneous Frequency (Hz)');
saveas(gcf,[figloc,recname,'_pspindlestats'],'jpeg')





end

