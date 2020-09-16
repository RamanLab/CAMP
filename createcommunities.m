initCobraToolbox();
%get input models and generate all possible pairwise combinations
folder = 'LAB_models';
F = get_model_names('LAB_models');
v = {F};
C = nchoosek(v{1,1},2);

%generate community models and save them as .mat files
for iter = 1:length(C)
    %load the models
    modelpair = {C{iter,1},C{iter,2}};
    org1{iter,1} = readCbModel(modelpair{1,1});
    org2{iter,1} = readCbModel(modelpair{1,2});
    modelCell = {org1{iter,1} org2{iter,1}};
    
    %identify biomass and ATPM reactions
    spBm{1,iter} = org1{iter,1}.rxns(find(org1{iter,1}.c));
    spBm{2,iter} = org2{iter,1}.rxns(find(org2{iter,1}.c));
    spATPM{1,iter}= validatestring('DM_atp_c_',org1{iter,1}.rxns);
    spATPM{2,iter} = validatestring('DM_atp_c_',org2{iter,1}.rxns);
    minNorm = 1;
    %create options structure containing species biomass and ATPM fields
    options(iter).spBm= spBm(:,iter);
    options(iter).spATPM = spATPM(:,iter);
    options(iter).sepUtEx = false;
    options(iter).minNorm = minNorm;
    
    disp(iter);
    %use SteadyCom function createCommModel to create a community model
    [modelCom{iter,1},infoCom{iter,1},indCom{iter,1}] = createCommModel(modelCell, options(iter));
    modelCom{iter,1}.modelID = modelCom{iter,1}.description;
    modelCom{iter,1}.modelName = modelCom{iter,1}.description;
    modelComNew = modelCom{iter,1};
    modelComNew.subSystems = modelCom{iter,1}.subSystems;
    infoComNew = infoCom{iter,1};
    indComNew = indCom{iter,1};
    filename = (modelComNew.description);
    save(filename, 'modelComNew','infoComNew','indComNew','-mat');
end