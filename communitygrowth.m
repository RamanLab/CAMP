%In case of community-specific uptake nutrient state, set the glucose
%and xylose uptake rates of each individual model in the coculture
%GUP is an array of modelname and monoculture uptake of glucose at
%half-maximal growth
%XUP is array of modelname and monoculture uptake of xylose at half-maximal
%growth
%C is all possible combinations of communities obtained in
%createcommunities.m
glcup = cell2mat(cellfun(@(x) GUP{strmatch(x,GUP(:,1),'exact'),2},C,'UniformOutput', false));
xylup = cell2mat(cellfun(@(x) XUP{strmatch(x,XUP(:,1),'exact'),2},C,'UniformOutput', false));
%calculate total glucose and xylose uptake for each community model
for k = 1:length(glcup)
    glcsum(k,1) = glcup(k,1)+ glcup(k,2);
    xylsum(k,1) = xylup(k,1)+ xylup(k,2);
end
%load community models from folder
folder = 'LAB_communities';
Files = get_model_names('LAB_communities');
Files = Files';

for iter = 1:length(Files)
    disp(iter)
    modelCom{iter,1} = readCbModel(Files{iter,1});
    spBm{iter,1} = modelCom{iter,1}.rxns(find(modelCom{iter,1}.c));
    rxnNameList{iter,1} = {'EX_lac_D(u)'};
    minNorm = 1;
    options(iter).rxnNameList = rxnNameList{iter,1};  
    options(iter).minNorm = minNorm;
    %identify glucose and xylose exchange rxns in each organism
    glu_org1{iter,1} = findRxnIDs(modelCom{iter,1},'EX_glc_D(e)_org1');
    glu_org2{iter,1} = findRxnIDs(modelCom{iter,1},'EX_glc_D(e)_org2');
    xyl_org1{iter,1} = findRxnIDs(modelCom{iter,1},'EX_xyl_D(e)_org1');
    xyl_org2{iter,1} = findRxnIDs(modelCom{iter,1},'EX_xyl_D(e)_org2');
    %change bounds of glucose and xylose exchange based on community
    %specific nutrient uptakes (skip lines 37 to 50 in other nutrient
    %states)
    if glu_org1{iter,1} ~= 0
        modelCom{iter,1}=changeRxnBounds(modelCom{iter,1},{'EX_glc_D(e)_org1'},glcup(iter,1) ,'l');
    end
    if glu_org2{iter,1} ~= 0
        modelCom{iter,1}=changeRxnBounds(modelCom{iter,1},{'EX_glc_D(e)_org2'},glcup(iter,2) ,'l');
    end
    if xyl_org1{iter,1} ~=0
        modelCom{iter,1}=changeRxnBounds(modelCom{iter,1},{'EX_xyl_D(e)_org1'},xylup(iter,1) ,'l');
    end
    if xyl_org2{iter,1} ~=0
        modelCom{iter,1}=changeRxnBounds(modelCom{iter,1},{'EX_xyl_D(e)_org2'},xylup(iter,2) ,'l');
    end
    modelCom{iter,1}=changeRxnBounds(modelCom{iter,1},{'EX_xyl_D(u)'},xylsum(iter,1) ,'l');
    modelCom{iter,1}=changeRxnBounds(modelCom{iter,1},{'EX_glc_D(u)'},glcsum(iter,1) ,'l');
    %add MRS media components and ATPM constraints
    modelCom{iter,1}= communityuptake(modelCom{iter,1});
    ATPM_org1{iter,1} = findRxnIDs(modelCom{iter,1},'DM_atp_c__org1');
    ATPM_org2{iter,1} = findRxnIDs(modelCom{iter,1},'DM_atp_c__org2');
    modelCom{iter,1}.lb(ATPM_org1{iter,1}) = 0.36;
    modelCom{iter,1}.ub(ATPM_org1{iter,1}) = 0.36;
    modelCom{iter,1}.lb(ATPM_org2{iter,1}) = 0.36;
    modelCom{iter,1}.ub(ATPM_org2{iter,1}) = 0.36;
    
    %FBA with SteadyCom
    [sol{iter,1},result{iter,1}] = SteadyComCplex(modelCom{iter,1},options(iter));
    if (result{iter,1}.vBM(1)) >=0.01 && (result{iter,1}.vBM(2)) >= 0.01
        if isequal(result{iter,1}.stat,'optimal')
            if glu_org1{iter,1} ~= 0
                glcuptake_org1{iter,1} = result{iter,1}.flux(glu_org1{iter,1});
            end
            if glu_org2{iter,1} ~= 0
                glcuptake_org2{iter,1} = result{iter,1}.flux(glu_org2{iter,1});
            end
            if xyl_org1{iter,1} ~=0
                xyluptake_org1{iter,1} = result{iter,1}.flux(xyl_org1{iter,1});
            end
            if xyl_org2{iter,1} ~=0
                xyluptake_org2{iter,1} = result{iter,1}.flux(xyl_org2{iter,1});
            end
            %FVA with SteadyCom 
            [minFlux{iter,1},maxFlux{iter,1},minFD{iter,1},maxFD{iter,1}, GRvector{iter,1}, result2{iter,1},LP{iter,1}] = SteadyComFVACplex(modelCom{iter,1},options(iter));
            modelduo{iter,1} = modelCom{iter,1};
            commgrowth{iter,1} = result{iter,1}.GRmax;
            orggrowth{iter,1} = result{iter,1}.vBM;
            abundance{iter,1} = result{iter,1}.BM;
        end
    end
end
% identify only viable pairs list
A = horzcat(commgrowth, orggrowth, abundance);
B = num2cell(find(~cellfun(@isempty,A(:,2))));
for iter1=1:length(B)
    m = B{iter1,1};
    viable(iter1,:) = A(m,:);
    viablenames(iter1,:) = Files(m,1);
end

%extract maxflux lactate values for the viable pairs
MFviable = num2cell(find(~cellfun(@isempty,maxFlux(:,1))));
for iter2=1:length(MFviable)
    n = MFviable{iter2,1};
    viablelac{iter2,:} = maxFlux{n,:};
end

D = horzcat(viablenames, viable, viablelac);
Modelduo = modelduo(~cellfun('isempty',modelduo));
vresult = {};
for index = 1:length(Modelduo)
    for index2 = 1:length(modelCom)
        if strmatch(Modelduo{index,1}.modelID, modelCom{index2,1}.modelID, 'exact')
            vresult{index,1} = result2{index2,1};
            vindCom{index,1} = modelCom{index2,1}.indCom;
            vinfoCom{index,1} = modelCom{index2,1}.infoCom;
        end
    end
end

%identify cross-fed metabolites between each community and tabulate
EXch= {};
threshold = 0.01;
for iter3 = 1:length(Modelduo)
    [nEx(iter3,1),~] = size(vinfoCom{iter3,1}.EXcom);
    fed{iter3,1} = zeros(nEx(iter3,1),1);
    for iter4= 1:nEx(iter3,1)
        posflux = 0;
        negflux = 0;
        for k = 1:2
            if vindCom{iter3,1}.EXsp(iter4,k) ~= 0
                if vresult{iter3,1}.flux((vindCom{iter3,1}.EXsp(iter4,k))) >= threshold
                    posflux = 1;
                end
                if vresult{iter3,1}.flux((vindCom{iter3,1}.EXsp(iter4,k))) <= -(threshold)
                    negflux = 1;
                end
            end
        end
        if posflux == 1 & negflux == 1
            fed{iter3,1}(iter4)= 1;
        end
    end
    xchanged{iter3,1} = vinfoCom{iter3,1}.EXsp(fed{iter3,1}==1);
    xchangedflux{iter3,1} = vresult{iter3,1}.flux(vindCom{iter3,1}.EXsp((fed{iter3,1}==1),:));
    list{iter3,1}= horzcat(xchanged(iter3,1),xchangedflux(iter3,1));
    if ~isempty(xchanged{iter3,1})
        EXch = [EXch; xchanged{iter3,1}{:,1}];
    end
end
filename = 'communityresults.xlsx';
D1 = cell2table(D);
D1.Properties.VariableNames ={'Communitymodel';'CommunityGrowthRate';'OrganismGrowthRate';'OrganismAbundance';'ProductFlux_FVA'};
EX = cell2table(EXch);
EX.Properties.VariableNames ={'Cross_fed_metabolites'};
writetable(D1, filename,'Sheet',1);
writetable(EX, filename,'Sheet',2);
