%to find the glucose and xyl uptake fluxes needed for 50%biomass growth in each model 
folder = 'LAB_models';
Files = get_model_names('LAB_models');
Files = Files.';
for iter = 1:length(Files)
    disp(iter);
%read the models
    model{iter,1} = readCbModel(Files{iter,1});
    [selExc, ~] = findExcRxns(model{iter,1});
    xchange{iter,1} = find(selExc==1);
    model{iter,1}.lb(xchange{iter,1}) = 0;
    model{iter,1}.ub(xchange{iter,1}) = 1000;
%add nutrient constraints and ATPM constraints
    model{iter,1}= monocultureuptake(model{iter,1});
    ATPM{iter,1} = findRxnIDs(model{iter,1},'DM_atp_c_');
    model{iter,1}.lb(ATPM{iter,1}) = 0.36;
    model{iter,1}.ub(ATPM{iter,1}) = 0.36;
%FBA
    sol{iter,1} = optimizeCbModel(model{iter,1},'max');
    Growth{iter,1}= round(sol{iter,1}.f,4);
    growth50{iter,1} = (Growth{iter,1}/2);
end
%pinning biomass growth to 50% and identifying the glucose xylose uptake
%in monocultureuptake..have -1000 for EX_glc_D when computing glu uptake, and
%similar for xylose.
for iter = 1:length(Files)
    disp(iter);
    model{iter,1} = readCbModel(Files{iter,1});
    % glu_rxn{i,1} = findRxnIDs(model{i,1},'EX_glc_D(e)');
    xyl_rxn{iter,1} = findRxnIDs(model{iter,1},'EX_xyl_D(e)');
    if xyl_rxn{iter,1} ~= 0
        obj{iter,1} = find(model{iter,1}.c);
        model{iter,1}.c(obj{iter,1}) =0;
        % model{i,1}.c(glu_rxn{i,1}) = 1;
        model{iter,1}.c(xyl_rxn{iter,1}) = 1;
    end
    [selExc, ~] = findExcRxns(model{iter,1});
    xchange{iter,1} = find(selExc==1);
    model{iter,1}.lb(xchange{iter,1}) = 0;
    model{iter,1}.ub(xchange{iter,1}) = 1000;
    if growth50{iter,1} ~= 0
        model{iter,1}.lb(obj{iter,1}) = growth50{iter,1};
        model{iter,1}.ub(obj{iter,1}) = growth50{iter,1};
    end
    model{iter,1}= monocultureuptake(model{iter,1});
    ATPM{iter,1} = findRxnIDs(model{iter,1},'DM_atp_c_');
    model{iter,1}.lb(ATPM{iter,1}) = 0.36;
    model{iter,1}.ub(ATPM{iter,1}) = 0.36;
    sol{iter,1} = optimizeCbModel(model{iter,1},'max');
    %glucoseuptake{i,1} = round(sol{i,1}.f,4);
    xyloseuptake{iter,1} = round(sol{iter,1}.f,4);
end
%use these uptake fluxes of glucose and xylose obtained from above for community-specific
%nutrient condition