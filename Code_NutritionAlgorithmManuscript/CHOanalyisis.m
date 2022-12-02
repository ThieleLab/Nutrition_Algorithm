% CHO cell analysis accompanying "A Nutrition Algorithm to Optimize Feed
% and Medium Composition Using Genome-Scale Metabolic Models" 
% Bronson R Weston and Ines Thiele 2022

close all
results = readcell('ModelsAndData/Results from CHO experiments.xlsx');

%% Get statistics on randomize sets
[mAP2,stdAP2,popAP2,mBP2,stdBP2,popBP2]=getStatsOnRandomChoice(2,0.2,results);
[mAP3,stdAP3,popAP3,mBP3,stdBP3,popBP3]=getStatsOnRandomChoice(3,0.2,results);
[mAP4,stdAP4,popAP4,mBP4,stdBP4,popBP4]=getStatsOnRandomChoice(4,0.2,results);
[mAP5,stdAP5,popAP5,mBP5,stdBP5,popBP5]=getStatsOnRandomChoice(5,0.2,results);



%% Test case
distP=getCountDistribution({'thr_L', 'arachd'}, results); %get test case overall results
statsP=getDistStatistics(distP,results);
indexes1=find(distP==1)+1;
indexes2=find(distP==2)+1;
pPaAP=ranksum(cell2mat(results(indexes1,end)),cell2mat(results(indexes2,end))); %ranksum test, test case antibody
pPaBM=ranksum(cell2mat(results(indexes1,end-1)),cell2mat(results(indexes2,end-1))); %ranksum test, test case biomass
[popA,popB]=randomEquivalentSetAnalysis({'thr_L', 'arachd'},results);%Generate random equivalent set population
[HbTC,PbTC]=ttest(popB,mean(cell2mat(results(indexes2,end-1))),'Tail','left') %ttest, test case biomass
[HaTC,PaTC]=ttest(popA,mean(cell2mat(results(indexes2,end))),'Tail','left') %ttest, test case antibody
perbTC=getPercentile(popB,mean(cell2mat(results(indexes2,end-1))))  %percentile rank, test case biomass
peraTC=getPercentile(popA,mean(cell2mat(results(indexes2,end))))    %percentile rank, test case antibody

%% Algorithm Results -> Prediction1
dist1=getCountDistribution({'thr_L','lys'}, results);  %get predicton 1 overall results
stats1=getDistStatistics(dist1,results); 
indexes1=find(dist1==1)+1;
indexes2=find(dist1==2)+1;
p1aAP=ranksum(cell2mat(results(indexes1,end)),cell2mat(results(indexes2,end)));    %ranksum test, prediction 1 antibody
p1aBM=ranksum(cell2mat(results(indexes1,end-1)),cell2mat(results(indexes2,end-1))); %ranksum test, prediction 1 biomass
[popA,popB]=randomEquivalentSetAnalysis({'3mob or val_L', 'thr_L','lys'},results); %Generate random equivalent set population
[Hb1,Pb1]=ttest(popB,mean(cell2mat(results(indexes2,end-1))),'Tail','left'); %ttest, prediction 1 biomass
[Ha1,Pa1]=ttest(popA,mean(cell2mat(results(indexes2,end))),'Tail','left');   %ttest, prediction 1 antibody
perb1=getPercentile(popB,mean(cell2mat(results(indexes2,end-1)))); %percentile rank, prediction 1 biomass
pera1=getPercentile(popA,mean(cell2mat(results(indexes2,end))));  %percentile rank, prediction 1 antibody

%% Algorithm Results -> Prediction2
dist2=getCountDistribution({'3mob or val_L', 'thr_L','lys'}, results);  %get predicton 2 overall results
stats2=getDistStatistics(dist2,results); 
indexes2=find(dist2==2)+1;
indexes3=find(dist2==3)+1;
p1aAP=ranksum(cell2mat(results(indexes2,end)),cell2mat(results(indexes3,end)));    %ranksum test, prediction 2 antibody
p1aBM=ranksum(cell2mat(results(indexes2,end-1)),cell2mat(results(indexes3,end-1))); %ranksum test, prediction 2 biomass
[popA,popB]=randomEquivalentSetAnalysis({'3mob or val_L', 'thr_L','lys'},results); %Generate random equivalent set population
[Hb1,Pb1]=ttest(popB,mean(cell2mat(results(indexes3,end-1))),'Tail','left'); %ttest, prediction 2 biomass
[Ha1,Pa1]=ttest(popA,mean(cell2mat(results(indexes3,end))),'Tail','left');   %ttest, prediction 2 antibody
perb1=getPercentile(popB,mean(cell2mat(results(indexes3,end-1)))); %percentile rank, prediction 2 biomass
pera1=getPercentile(popA,mean(cell2mat(results(indexes3,end))));  %percentile rank, prediction 2 antibody

%% Algorithm Results -> Prediction 3
dist3=getCountDistribution({'arachd','3mob or val_L', 'thr_L', 'lys','trp_L'}, results); %get predicton 3 overall results
stats3=getDistStatistics(dist3,results);
indexes4=find(dist3==4)+1;
indexes5=find(dist3==5)+1;
p3aAP=ranksum(cell2mat(results(indexes4,end)),cell2mat(results(indexes5,end))); %ranksum test, prediction 3 antibody
p3aBM=ranksum(cell2mat(results(indexes4,end-1)),cell2mat(results(indexes5,end-1))); %ranksum test, prediction 3 biomass
[popA,popB]=randomEquivalentSetAnalysis({'arachd','3mob or val_L', 'thr_L', 'lys','trp_L'},results); %Generate random equivalent set population
[Hb3,Pb3]=ttest(popB,mean(cell2mat(results(indexes5,end-1))),'Tail','left'); %ttest, prediction 3 biomass
[Ha3,Pa3]=ttest(popA,mean(cell2mat(results(indexes5,end))),'Tail','left'); %ttest, prediction 3 antibody
perb3=getPercentile(popB,mean(cell2mat(results(indexes5,end-1)))); %percentile rank, prediction 3 biomass
pera3=getPercentile(popA,mean(cell2mat(results(indexes5,end)))); %percentile rank, prediction 3 antibody


%% Plot Results
figure()
set(gcf,'Position',[0 0 1250 660])

subplot(2,2,1)
model_series = [stats1{end-1,3},stats1{end,3}; stats2{end-1,3},stats2{end,3};stats3{end-1,3},stats3{end,3};  statsP{end-1,3},statsP{end,3}]; 
model_error = [stats1{end-1,4},stats1{end,4}; stats2{end-1,4},stats2{end,4}; stats3{end-1,4},stats3{end,4};  statsP{end-1,4},statsP{end,4}]; 
b = bar(model_series, 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',model_series,model_error,'k','linestyle','none');
leg=legend({['Missing One' newline 'Metabolite'],['All' newline 'Metabolites']},'Position',[0.5160    0.7925    0.1101    0.1342]);
xticklabels({'Prediction 1', 'Prediction 2', 'Prediction 3','Test Case'})
hold off
ylabel(['Integral Viable Cell Count' newline '(Million Cells)'])
set(gca, 'FontName', 'Times', 'FontSize', 14);

subplot(2,2,3)
model_series = [stats1{end-1,5},stats1{end,5}; stats2{end-1,5},stats2{end,5};stats3{end-1,5},stats3{end,5}; statsP{end-1,5},statsP{end,5}]; 
model_error = [stats1{end-1,6},stats1{end,6}; stats2{end-1,6},stats2{end,6};stats3{end-1,6},stats3{end,6};  statsP{end-1,6},statsP{end,6}]; 
b = bar(model_series, 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',model_series,model_error,'k','linestyle','none');
% legend({'Missing One Metabolite','All Metabolites'})
xticklabels({'Prediction 1', 'Prediction 2', 'Prediction 3','Test Case'})
hold off
ylabel(['Total mAb' newline 'Expression (\mug)'])
set(gca, 'FontName', 'Times', 'FontSize', 14);


%% Analysis Functions

function [popAP,popBP]=randomEquivalentSetAnalysis(predictionSet,results)
nMets=length(predictionSet);
Combinations=nchoosek(1:12,nMets);
meanAntibodyProduction=zeros(length(Combinations(:,1))-1,1);
meanBiomassProduction=zeros(length(Combinations(:,1))-1,1);
metaboliteList={'gln_L or asn_L','lys','trp_L','thr_L','3mob or val_L', 'his_L', 'thm', 'pydxn', 'thymd or dgsn', 'dcyt','retinol','arachd'};
[~,ia,~] = intersect(metaboliteList,predictionSet);
ia=ia.';
for i=1:length(Combinations(:,1))
    if isempty(setdiff(Combinations(i,:),ia))
        Combinations(i,:)=[];
        break
    end
end
for i=length(Combinations(:,1)):-1:1
    metabolites=metaboliteList(Combinations(i,:));
    dist= getCountDistribution(metabolites,results);
    distStats=getDistStatistics(dist,results);
    if distStats{end,1}~=nMets
        meanAntibodyProduction(i)=[];
        meanBiomassProduction(i)=[];
    else
        meanAntibodyProduction(i)=distStats{end,end-1};
        meanBiomassProduction(i)=distStats{end,end-3};
    end
end
popAP=meanAntibodyProduction;
popBP=meanBiomassProduction;
end


function [mAP,stdAP,popAP,mBP,stdBP,popBP]=getStatsOnRandomChoice(nMets,nSamples,results)
Combinations=nchoosek(1:12,nMets);
if nSamples<1
    nSamples=round(length(Combinations(:,1))*nSamples);
end
if nSamples>length(Combinations(:,1))
    error('nSamples greater than number of possible combinations')
end
ind=randsample(length(Combinations(:,1)),nSamples);
meanAntibodyProduction=zeros(nSamples,1);
meanBiomassProduction=zeros(nSamples,1);
metaboliteList={'gln_L or asn_L','lys','trp_L','thr_L','val_L or 3mob', 'his_L', 'thm', 'pydxn', 'thymd or dgsn', 'dcyt','retinol','arachd'};
for i=nSamples:-1:1
    metabolites=metaboliteList(Combinations(ind(i),:));
    dist= getCountDistribution(metabolites,results);
    distStats=getDistStatistics(dist,results);
    if distStats{end,1}~=nMets
        meanAntibodyProduction(i)=[];
        meanBiomassProduction(i)=[];
    else
        meanAntibodyProduction(i)=distStats{end,end-1};
        meanBiomassProduction(i)=distStats{end,end-3};
    end
end
popAP=meanAntibodyProduction;
mAP=mean(meanAntibodyProduction);
stdAP=std(meanAntibodyProduction);
mBP=mean(meanBiomassProduction);
stdBP=std(meanBiomassProduction);
popBP=meanBiomassProduction;
end

function [distribution]= getCountDistribution(metabolites,results)
distribution=zeros(length(results(:,1))-1,1);
for i=1:length(metabolites)
    strings=strsplit(metabolites{i},' or ');
    for t=1:length(results(:,1))-1
        for m=1:length(strings)
            f=find(strcmp(results(1,:),strings{m}));
            if results{t+1,f}==1
                distribution(t)=distribution(t)+1;
                break
            end
        end
    end
end
end

function [stats]=getDistStatistics(distribution,results)
counts=unique(distribution);
stats=cell(length(counts)+1,6);
stats(1,:)={'Count', 'N', 'Biomass Mean','Biomass Std','Antibody Mean','Antibody Std'};
for i=1:length(counts)
    stats{i+1,1}=counts(i);
    indexes=find(distribution==counts(i))+1;
    stats{i+1,2}=length(indexes);
    stats{i+1,3}=mean(cell2mat(results(indexes,end-1)));
    stats{i+1,4}=std(cell2mat(results(indexes,end-1)));
    stats{i+1,5}=mean(cell2mat(results(indexes,end)));
    stats{i+1,6}=std(cell2mat(results(indexes,end)));
end
end

function [percentile]=getPercentile(population,m)
n=0;
for i=1:length(population)
    if population(i)>=m
        n=n+1;
    end
end
percentile=n/length(population);
end
