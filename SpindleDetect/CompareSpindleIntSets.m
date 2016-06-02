function [  ] = CompareSpindleIntSets( intset1,intset2,exchans,intnames,figloc )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[ int1overlap, int2overlap ] = FindOverlappingInts( intset1,intset2 );



%%

only1 = sum(int1overlap==0);
only2 = sum(int2overlap==0);
overlap1 = sum(int1overlap);
overlap2 = sum(int2overlap);

only1ints = find(int1overlap==0);
only2ints = find(int2overlap==0);
overlapints1 = find(int1overlap);

%%
int1dur = intset1(:,2)-intset1(:,1);
int2dur = intset2(:,2)-intset2(:,1);

[durhist1,histbins] =hist(int1dur,20);
[durhist2] =hist(int2dur,histbins);


%% Merge ints
plotints = bsxfun(@(X,Y) X+Y,[intset1;intset2],[-1 1]);
[ plotints ] = MergeSeparatedInts( plotints,0 );

%% Load LFP to show examples and average
% Load downsampled LFP

%loadints = bsxfun(@(X,Y) X+Y,plotints,[-1 1]);

downsamplefactor = 2;
allLFP = GetLFP_Down(exchans,'intervals',plotints,'downsample',downsamplefactor);
sf_LFP = 1250./downsamplefactor;

allLFP(:,2:end) = zscore(allLFP(:,2:end));
allLFP(:,2) = mean(allLFP(:,2:end),2);
allLFP(:,3:end) = [];

%%
int1examples = [datasample(overlapints1,2);datasample(only1ints,2)];
int2examples = [datasample(only2ints,2)];

figure
    subplot(4,3,1:2)
        bar([overlap1,only1;overlap2,only2],'stacked')
        legend('Overlap','Nonoverlap','location','eastoutside')
        set(gca,'XTickLabel',intnames)
    subplot(4,3,3)
        plot(histbins,durhist1,'b');
        hold on
        plot(histbins,durhist2,'r');
        xlabel('pSpindle Duration');ylabel('# pSpindles');
        legend(intnames{1},intnames{2})
        
for ee = 1:4        
    subplot(4,2,ee+2)
        plot(intset1',3*ones(size(intset1))','b')
        hold on
        plot(intset2',3.2*ones(size(intset2))','r')
        plot(allLFP(:,1),allLFP(:,2),'k')


        xlim(intset1(int1examples(ee),:)+[-1 1])
end

for ee = 1:2        
    subplot(4,2,ee+6)
        plot(intset1',3*ones(size(intset1))','b')
        hold on
        plot(intset2',3.2*ones(size(intset2))','r')
        plot(allLFP(:,1),allLFP(:,2),'k')


        xlim(intset2(int2examples(ee),:)+[-1 1])
end
        

saveas(gcf,[figloc,'_pspindlecomparefig'],'jpeg')
end

