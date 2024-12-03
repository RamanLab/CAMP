function communitygrowth(folder_name,RxnName,substrateRxns,reactionsAndBounds,growthrate_threshold,cross_fed_threshold)
% RxnName is the product/metaolite of interest reaction identifier in the
% model, such as 'EX_lac_D(u)'
% substrateRxn is a cell array with the carbon substrates which the community consumes such
% as {'EX_glc_D(e)';'EX_fru(e)'};
Files = get_model_names(folder_name);
Files = Files';
for iter = 1:length(Files)
    disp(iter)
    modelCom{iter,1} = readCbModel(Files{iter,1});
    spBm{iter,1} = modelCom{iter,1}.rxns(find(modelCom{iter,1}.c));
    rxnNameList{iter,1} = RxnName;
    minNorm = 1;
    options(iter).rxnNameList = rxnNameList{iter,1};
    options(iter).minNorm = minNorm;
    
    % Identify substrate rxns in each organism    
    substrateRxn_org1 = cell(size(substrateRxns));
    substrateRxn_org2 = cell(size(substrateRxns));
    substrate_org1 = cell(size(substrateRxns));
    substrate_org2 = cell(size(substrateRxns));

    for j = 1:length(substrateRxns)
        substrateRxn_org1{j} = [substrateRxns{j} '_org1'];
        substrateRxn_org2{j} = [substrateRxns{j} '_org2'];
        substrate_org1{iter,j} = findRxnIDs(modelCom{iter,1}, substrateRxn_org1{j});
        substrate_org2{iter,j} = findRxnIDs(modelCom{iter,1}, substrateRxn_org2{j});
    end
    
    %add media components and ATPM constraints
    modelCom{iter,1}= communityuptake(modelCom{iter,1},reactionsAndBounds);
    ATPM_org1{iter,1} = findRxnIDs(modelCom{iter,1},'DM_atp_c__org1');
    ATPM_org2{iter,1} = findRxnIDs(modelCom{iter,1},'DM_atp_c__org2');
    modelCom{iter,1}.lb(ATPM_org1{iter,1}) = 0.36;
    modelCom{iter,1}.ub(ATPM_org1{iter,1}) = 0.36;
    modelCom{iter,1}.lb(ATPM_org2{iter,1}) = 0.36;
    modelCom{iter,1}.ub(ATPM_org2{iter,1}) = 0.36;
    
    %FBA with SteadyCom
    try
        [sol{iter,1},result{iter,1}] = SteadyComCplex(modelCom{iter,1},options(iter));
    catch
        disp('Error in SteadyComCplex');
        result{iter,1}.GRmax = 0;
        result{iter,1}.vBM = [-1;-1];
        result{iter,1}.BM = [1;1];
        result{iter,1}.flux = zeros(size((modelCom{iter,1}.rxns),1),1);
        result{iter,1}.stat = 'error';
        
    end
    if (result{iter,1}.vBM(1)) >=growthrate_threshold && (result{iter,1}.vBM(2)) >= growthrate_threshold
        if isequal(result{iter,1}.stat,'optimal')
            for j = 1:length(substrateRxns)
                if substrate_org1{iter,j} ~= 0
                    substrate_org1{iter,j} = result{iter,1}.flux(substrate_org1{iter,j});
                end
                if substrate_org2{iter,j} ~= 0
                    substrate_org2{iter,j} = result{iter,1}.flux(substrate_org2{iter,j});
                end
            end
            %FVA with SteadyCom
            [minFlux{iter,1},maxFlux{iter,1},minFD{iter,1},maxFD{iter,1}, GRvector{iter,1}, result2{iter,1},LP{iter,1}] = SteadyComFVACplex(modelCom{iter,1},options(iter));
            modelduo{iter,1} = modelCom{iter,1};
            commgrowth{iter,1} = result{iter,1}.GRmax;
            orggrowth{iter,1} = result{iter,1}.vBM;
            abundance{iter,1} = result{iter,1}.BM;
            substrate_uptakeflux{iter,1} = [substrate_org1{iter,:}, substrate_org2{iter,:}]; 
        end
    end
end
% identify only viable pairs list
A = horzcat(commgrowth, orggrowth, abundance,substrate_uptakeflux);
B = num2cell(find(~cellfun(@isempty,A(:,2))));
for iter1=1:length(B)
    m = B{iter1,1};
    viable(iter1,:) = A(m,:);
    viablenames(iter1,:) = Files(m,1);
end

%extract maxflux values for the viable pairs
MFviable = num2cell(find(~cellfun(@isempty,maxFlux(:,1))));
for iter2=1:length(MFviable)
    n = MFviable{iter2,1};
    viableproduct{iter2,:} = maxFlux{n,:};
end

D = horzcat(viablenames, viable, viableproduct);
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
cross_fed_threshold = 0.0001;
for iter3 = 1:length(Modelduo)
    [nEx(iter3,1),~] = size(vinfoCom{iter3,1}.EXcom);
    fed{iter3,1} = zeros(nEx(iter3,1),1);
    for iter4= 1:nEx(iter3,1)
        posflux = 0;
        negflux = 0;
        for k = 1:2
            if vindCom{iter3,1}.EXsp(iter4,k) ~= 0
                if vresult{iter3,1}.flux((vindCom{iter3,1}.EXsp(iter4,k))) >= cross_fed_threshold
                    posflux = 1;
                end
                if vresult{iter3,1}.flux((vindCom{iter3,1}.EXsp(iter4,k))) <= -(cross_fed_threshold)
                    negflux = 1;
                end
            end
        end
        if posflux == 1 && negflux == 1
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
D1.Properties.VariableNames ={'Communitymodel';'CommunityGrowthRate';'OrganismGrowthRate';'OrganismAbundance';'SubstrateUptakeFluxesForOrg1And2';'ProductFlux_FVA'};
EX = cell2table(EXch);
EX.Properties.VariableNames ={'Cross_fed_metabolites'};
writetable(D1, filename,'Sheet',1);
writetable(EX, filename,'Sheet',2);

    function [file_names] = get_model_names(folder_name)
        %GET_MODEL_NAMES Returns an array of names of the .mat files inside
        %a given folder
        file_info = dir(strcat(strcat('./',folder_name),'/*.mat'));
        file_names = ({file_info.name});
    end
    function [modelCom] = communityuptake(modelCom, reactionsAndBounds)
        % INPUTS:
        %   modelCom               : community model
        %   reactionsAndBounds : Cell array with reaction names and corresponding bound values
        %                           e.g., {'EX_ala_L(u)', -1;
        %                           'EX_arg_L(u)', -1; ...} 'u' refers to
        %                           the community compartment
        % Apply bounds to reactions
        for i = 1:size(reactionsAndBounds, 1)
            reaction = reactionsAndBounds{i, 1};
            bound_value = reactionsAndBounds{i, 2};
            modelCom = changeRxnBounds(modelCom, {reaction}, bound_value, 'l');
        end
    end
end