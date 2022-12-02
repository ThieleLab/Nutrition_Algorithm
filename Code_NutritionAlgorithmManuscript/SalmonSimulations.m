% Salmon simulations accompanying "A Nutrition Algorithm to Optimize Feed
% and Medium Composition Using Genome-Scale Metabolic Models" 
% Bronson R Weston and Ines Thiele 2022

close all
clear options
load('ModelsAndData/salmonModel.mat');
LPSolver = 'ibm_cplex';
if ~exist('solverInstalled','var')
    [solverOK, solverInstalled] = changeCobraSolver(LPSolver, 'LP',0,1);
end


%% Testing if fish or soybean meal should be added to the diet (Table 1, Diet 1)
modelTest=salmonModel;
clear options
options.roiWeights=[10];

options.targetedDietRxns={'Diet_EX_FM[e]',0;'Diet_EX_SBM[e]',0};
modelTest=addElementTracker(modelTest,'N','Exit');
options.foodRemovedLimit=0;
modelTest.lb(718)=1;
[newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges]=nutritionAlgorithm(modelTest,{'N_Track'},{'min'},options);

%% Testing the flux on 1.3361 g of soybean meal (Table 1, Control)
modelTest=salmonModel;
clear options
modelTest=addElementTracker(modelTest,'N','Exit');
options.foodRemovedLimit=0;
options.foodAddedLimit=0;
modelTest.lb(find(contains(modelTest.rxns,'Diet_EX_SBM[e]')))=-1.3361;
modelTest.ub(find(contains(modelTest.rxns,'Diet_EX_SBM[e]')))=-1.3361;
[newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges]=nutritionAlgorithm(modelTest,{'N_Track'},{'min'},options);

%% Testing what amino acids should be added to 1.3361 g of fishmeal (Table 1, Diet 2)
modelTest=salmonModel;
AminoAcidsInDiet=find(contains(modelTest.rxns,'Diet_EX'));
AminoAcidsInDiet=AminoAcidsInDiet(contains(modelTest.rxns(AminoAcidsInDiet),'__L[e]'));
AminoAcidsInDiet=modelTest.rxns(AminoAcidsInDiet);

clear options
options.roiWeights=[10];

options.targetedDietRxns=[{'Diet_EX_FM[e]',1;'Diet_EX_SBM[e]',1};[AminoAcidsInDiet,num2cell(0.1*ones(length(AminoAcidsInDiet),1))]];
modelTest=addElementTracker(modelTest,'N','Exit');
options.foodRemovedLimit=0;
modelTest.lb(718)=1;
modelTest.lb(find(contains(modelTest.rxns,'Diet_EX_FM[e]')))=-1.3361;
modelTest.ub(find(contains(modelTest.rxns,'Diet_EX_FM[e]')))=-1.3361;
[newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges]=nutritionAlgorithm(modelTest,{'N_Track'},{'min'},options);


%% Maximizing net profit in salmon (Figure 2)

%setting up variables
modelTest=salmonModel;
clear options
options.OFS=0;
nweights=0:0.00005:0.005;
bsfSln=zeros(1,length(nweights));
fmSln=zeros(1,length(nweights));
sbmSln=zeros(1,length(nweights));
profit=zeros(1,length(nweights));
unchangedProfit=zeros(1,length(nweights));
minProfit=zeros(1,length(nweights));
maxProfit=zeros(1,length(nweights));
minN=zeros(1,length(nweights));
maxN=zeros(1,length(nweights));
biomass=zeros(1,length(nweights));
options.foodRemovedLimit=0;
options.display='off';
options.targetedDietRxns={'Diet_EX_FM[e]',1.3/1000;'Diet_EX_SBM[e]',0.4/1000;'Diet_EX_BSF[e]',3/1000};
modelTest=addElementTracker(modelTest,'N','Exit');
salmonValue=12.6/1000;

%find solution for each possible nitrogen cost in nweights
for i=1:length(nweights)
    options.roiWeights=[nweights(i),salmonValue];
    [newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges,detailedAnalysis]=nutritionAlgorithm(modelTest,{'N_Track','Biomass'},{'min','max'},options);
    minN(i)=detailedAnalysis.Rxn1.min.ND.f;
    maxN(i)=detailedAnalysis.Rxn1.max.ND.f;
    profit(i)=-1*pointsModelSln.f;
    bsfSln(i)=-1*(newDietModel.lb(find(contains(newDietModel.rxns,'EX_BSF')))+newDietModel.ub(find(contains(newDietModel.rxns,'EX_BSF'))))/2;
    fmSln(i)=-1*(newDietModel.lb(find(contains(newDietModel.rxns,'EX_FM')))+newDietModel.ub(find(contains(newDietModel.rxns,'EX_FM'))))/2;
    sbmSln(i)=-1*(newDietModel.lb(find(contains(newDietModel.rxns,'EX_SBM')))+newDietModel.ub(find(contains(newDietModel.rxns,'EX_SBM'))))/2;
    biomass(i)=detailedAnalysis.Rxn2.max.ND.f;
    minProfit(i)= salmonValue*biomass(i)-(1.9/1000)*fmSln(i)-(0.6/1000)*sbmSln(i)-(3/1000)*bsfSln(i)-nweights(i)*maxN(i);
    maxProfit(i)= salmonValue*biomass(i)-(1.9/1000)*fmSln(i)-(0.6/1000)*sbmSln(i)-(3/1000)*bsfSln(i)-nweights(i)*minN(i);
    unchangedProfit(i)=salmonValue*1.61348-(3*0.4/1000)-(nweights(i)*1.858);
end

%plot results
figure()
set(gcf,'Position', [ 573.0000  524.3333  670.0000  439.6667])
subplot(2,1,2)
plot(nweights,profit,'LineWidth',2)
hold on
plot(nweights,unchangedProfit,'LineWidth',2)
plot(nweights,zeros(1,length(nweights)),'k--')
xlabel(['Nitrogen Waste Cost (',char(8364),'/mmol)'],'FontSize',12)
ylabel(['Net Profit (', char(8364),'/day)'],'FontSize',12)
legend({'Optimal Feed','3g SBM'},'FontSize',10)
ylim([min(unchangedProfit)*0.9 max(profit)*1.1])
set(gca, 'FontName', 'Times');
subplot(2,1,1)
plot(nweights,sbmSln,'LineWidth',2,'Color',[0.6 0 0.9]);
hold on
plot(nweights,fmSln,'LineWidth',2,'Color',[0 0.6 0]);
plot(nweights,bsfSln,'LineWidth',2);
legend({'SBM','FM','BSF'},'FontSize',10)
xlabel(['Nitrogen Waste Cost (',char(8364),'/mmol)'],'FontSize',12)
ylabel({'Optimal Feed'; 'Composition (g/day)'},'FontSize',12)
set(gca, 'FontName', 'Times');
