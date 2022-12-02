%% Speed Test
solver= 'ibm_cplex';
changeCobraSolver(solver, 'LP', 0, -1);

%% Salmon
load('ModelsAndData/salmonModel.mat');
modelTest=salmonModel;
modelTest=addElementTracker(modelTest,'N','Exit');
clear options
options.slnType='Quick';
options.display='off';
modelTest.lb(find(contains(modelTest.rxns,'Diet_EX_FM[e]')))=-1.34;
modelTest.ub(find(contains(modelTest.rxns,'Diet_EX_FM[e]')))=-1.34;
for i=1:10
    options.roiWeights=10;
    tic
    [newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges]=nutritionAlgorithm(modelTest,{'Biomass'},{'max'},options);
    timeNA1(i)=toc;
    tic
    optimizeCbModel(modelTest);
    timeFBA(i)=toc;
    options.roiWeights=[5,10];
    tic
    [newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges]=nutritionAlgorithm(modelTest,{'N_Track','Biomass'},{'min','max'},options);
    timeNA2(i)=toc;
end

salmonTime1=mean(timeNA1);
salmonTimeStd1=std(timeNA1);
salmonTime2=mean(timeNA2);
salmonTimeStd2=std(timeNA2);
salmonTimeFBA=mean(timeFBA);
salmonTimeStdFBA=std(timeFBA);

%% Ecoli
load ecoli_core_model.mat;
clear options
options.slnType='Quick';
ecoliModel=model;
% model=addEnvironment2Model(model);
model=convert_EX_to_diet(model);
model.osenseStr='max';
options.roiWeights=30;
model.lb(36)=0;
model.ub(36)=0;
f=find(contains(model.rxns,'Diet_EX'));
[newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges] = nutritionAlgorithm(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2','max',options);
timeNA1=zeros(10,1);
timeNA2=timeNA1;
timeFBA=timeNA1;
options.display='off';

for i=1:10
    options.roiWeights=30;
    tic
    [~,~,~,~,~] = nutritionAlgorithm(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2','max',options);
    timeNA1(i)=toc;
    tic
    optimizeCbModel(model);
    timeFBA(i)=toc;
    options.roiWeights=[30,10];
    tic
    [~,~,~,~,~] = nutritionAlgorithm(model,{'Biomass_Ecoli_core_N(w/GAM)-Nmet2','Exit_EX_co2[e]'},{'max','max'},options);
    timeNA2(i)=toc;
end

ecTime1=mean(timeNA1);
ecTimeStd1=std(timeNA1);
ecTime2=mean(timeNA2);
ecTimeStd2=std(timeNA2);
ecTimeFBA=mean(timeFBA);
ecTimeStdFBA=std(timeFBA);

%% CHO cell
load('ModelsAndData/CHOmodelUnconstrained.mat')
choModel=convert_EX_to_diet(CHOmodel);
f=find(contains(choModel.rxns,'Diet'));
choModel.rxns(f)=regexprep(choModel.rxns(f),'\_e\_','\[e\]');
choModel=changeObjective(choModel,'biomass_cho');
choModel.ub(6619)=100;
choModel.lb(6619)=1;
scalar=2000;
default=-50;
T=table(choModel.rxns(f),choModel.lb(f),'VariableNames',{'rxns','flux'});
for i=1:length(T{:,1})
    ind=find(strcmp(choModel.rxns,T{i,1}));
    if T{i,2}==0
        n=default;
    elseif T{i,2}<-999
        n=T{i,2};
    else
        n=scalar*T{i,2};
    end
    choModel.ub(ind)=n;
    choModel.lb(ind)=n;
end
clear options
options.slnType='Quick';
options.display='off';
timeNA1=zeros(10,1);
timeNA2=timeNA1;
timeFBA1=timeNA1;
for i=1:10
    options.roiWeights=30;
    tic
    [~,~,~,~,~] = nutritionAlgorithm(choModel,{'biomass_cho'},{'max'},options);
    timeNA1(i)=toc;
    tic
    optimizeCbModel(choModel);
    timeFBA(i)=toc;
    options.roiWeights=[30,10];
    tic
    [~,~,~,~,~] = nutritionAlgorithm(choModel,{'DM_oms[c]','biomass_cho'},{'max','max'},options);
    timeNA2(i)=toc;
end
choTime1=mean(timeNA1);
choTimeStd1=std(timeNA1);
choTime2=mean(timeNA2);
choTimeStd2=std(timeNA2);
choTimeFBA=mean(timeFBA);
choTimeStdFBA=std(timeFBA);

%% Human
load('ModelsAndData/humanModel.mat');
% modelTest= setDietConstraints(modelHM);

obj='Whole_body_objective_rxn';
o=find(strcmp(humanModel.rxns,obj));
humanModel = changeObjective(humanModel,obj);
humanModel= setDietConstraints(humanModel);
humanModel.osenseStr='max';
humanModel.lb(o)=1;
humanModel.ub(o)=1;
options.slnType='Quick';

for i=1:10
    options.display='off';
    options.roiWeights=30;
    tic
    [newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges] = nutritionAlgorithm(humanModel,{'Brain_GABAVESSEC'},'max',options);
    timeNA1(i)=toc;
    tic
    optimizeCbModel(humanModel);
    timeFBA(i)=toc;
    options.roiWeights=[30,10];
    tic
    [newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges] = nutritionAlgorithm(humanModel,{'Brain_GABAVESSEC','Liver_EX_creat(e)_[bc]'},{'max','max'},options);
    timeNA2(i)=toc;
end
hmTime1=mean(timeNA1);
hmTimeStd1=std(timeNA1);
hmTime2=mean(timeNA2);
hmTimeStd2=std(timeNA2);
hmTimeFBA=mean(timeFBA);
hmTimeStdFBA=std(timeFBA);