% CHO cell simulations accompanying "A Nutrition Algorithm to Optimize Feed
% and Medium Composition Using Genome-Scale Metabolic Models" 
% Bronson R Weston and Ines Thiele 2022

solver= 'ibm_cplex';
changeCobraSolver(solver, 'LP', 0, -1);

%% Constrained Model Simulations

%Load and modify model as specified in materials and methods
load('ModelsAndData/CHOmodelConstrained.mat','CHOmodel')
choModel=convert_EX_to_diet(CHOmodel);
f=find(contains(choModel.rxns,'Diet'));
choModel.rxns(f)=regexprep(choModel.rxns(f),'\_e\_','\[e\]');
T=table(choModel.rxns(f),choModel.lb(f),'VariableNames',{'rxns','flux'});
default=-1e-2;
for i=1:length(T{:,1})
    ind=find(strcmp(choModel.rxns,T{i,1}));
    if T{i,2}==0
        n=default;
    else
        n=T{i,2};
    end
    choModel.ub(ind)=n;
    choModel.lb(ind)=n;
end

%Set up weights for the nutrition algorithm, scheme 1
value=0.01;
options.targetedDietRxns={'Diet_EX_gln_L[e]',value; 'Diet_EX_asn_L[e]',value; 'Diet_EX_lys_L[e]',value; 'Diet_EX_trp_L[e]', value; ...
    'Diet_EX_thr_L[e]',value; 'Diet_EX_val_L[e]',value; 'Diet_EX_his_L[e]',value; 'Diet_EX_thm[e]',value; 'Diet_EX_pydxn[e]',value; ...
    'Diet_EX_thymd[e]',value;  'Diet_EX_dcyt[e]', value; 'Diet_EX_3mob[e]',value; 'Diet_EX_dgsn[e]',value; 'Diet_EX_retinol[e]',value; ...
    'Diet_EX_arachd[e]', value};
options.roiWeights=[0.01];
[newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges] = nutritionAlgorithm(choModel,{'biomass_cho_producing'},{'max'},options);

%Set up weights for the nutrition algorithm, scheme 2
weights={0.019830863; 0.008553752; 0.00515338; 0.039269251; 0.010325727; 0.006543513; 0.018931955; 4.209422181; 2.68371286; 0.012808365; ...
    0.01441762; 0.006140063; 6.056063912; 0.908826559; 1.207456636};
options.targetedDietRxns(:,2)=weights;
[newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges] = nutritionAlgorithm(choModel,{'biomass_cho_producing'},{'max'},options);

%% Unonstrained Model Simulations

%Load and modify model as specified in materials and methods
load('ModelsAndData/CHOmodelUnconstrained.mat')
choModel=convert_EX_to_diet(CHOmodel);
f=find(contains(choModel.rxns,'Diet'));
choModel.rxns(f)=regexprep(choModel.rxns(f),'\_e\_','\[e\]');
choModel=changeObjective(choModel,'biomass_cho');
choModel.ub(6619)=100;
choModel.lb(6619)=1;
scalar=2000;
default=-50;
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

%Set up weights for the nutrition algorithm, scheme 1
options.roiWeights=[1e4,10];
options.foodAddedLimit=100000;
value=0.01;
options.targetedDietRxns={'All',value};
options.display='on';
options.targetedDietRxns={'Diet_EX_gln_L[e]',value; 'Diet_EX_asn_L[e]',value; 'Diet_EX_lys_L[e]',value; 'Diet_EX_trp_L[e]', value; ...
    'Diet_EX_thr_L[e]',value; 'Diet_EX_val_L[e]',value; 'Diet_EX_his_L[e]',value; 'Diet_EX_thm[e]',value; 'Diet_EX_pydxn[e]',value; ...
    'Diet_EX_thymd[e]',value;  'Diet_EX_dcyt[e]', value; 'Diet_EX_3mob[e]',value; 'Diet_EX_dgsn[e]',value; 'Diet_EX_retinol[e]',value; ...
    'Diet_EX_arachd[e]', value};
[newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges] = nutritionAlgorithm(choModel,{'DM_oms[c]','biomass_cho'},{'max','max'},options);

%Set up weights for the nutrition algorithm, scheme 2
weights={0.019830863; 0.008553752; 0.00515338; 0.039269251; 0.010325727; 0.006543513; 0.018931955; 4.209422181; 2.68371286; 0.012808365; ...
    0.01441762; 0.006140063; 6.056063912; 0.908826559; 1.207456636};
options.targetedDietRxns(:,2)=weights;
[newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges] = nutritionAlgorithm(choModel,{'DM_oms[c]','biomass_cho'},{'max','max'},options);



