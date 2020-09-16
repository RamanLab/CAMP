function[Growth,fbaflux,maxflux]=monoculturegrowth(foldername,RxnName)
%This computes the monoculture growth rates of single organism models and
%the fluxes from FBA and FVA solutions for the particular reaction of interest
% INPUTS:
%   foldername :name of folder containing a list of COBRA models as
%               .mat files
%   RxnName  :exchange rxn name of the metabolite of interest in the model
%   e.g. 'EX_lac_D(e)'
% OUTPUTS:
%   Growth  : growth rates of each model with biomass rxn as objective
%   fbaflux : FBA flux value of input reaction
%   maxflux : FVA maximum flux value for the input reaction
%
%load input models
Files = get_model_names(foldername);
Files = Files.';
% computing growth rates for monoculture 
for iter = 1:length(Files)
    disp(iter);
    %read the models
    model{iter,1} = readCbModel(Files{iter,1});
    
    % add nutrient uptake bounds and set ATPM constraints, first set all exchange rxn bounds to zero
    [selExc, ~] = findExcRxns(model{iter,1});
    xchange{iter,1} = find(selExc==1);
    model{iter,1}.lb(xchange{iter,1}) = 0;
    model{iter,1}.ub(xchange{iter,1}) = 1000;
    
    model{iter,1}= monocultureuptake(model{iter,1});
    ATPM{iter,1} = findRxnIDs(model{iter,1},'DM_atp_c_');
    model{iter,1}.lb(ATPM{iter,1}) = 0.36;
    model{iter,1}.ub(ATPM{iter,1}) = 0.36;
    
    %identify glucose and xylose uptake rxns
    glu{iter,1} = findRxnIDs(model{iter,1},'EX_glc_D(e)');
    xyl{iter,1} = findRxnIDs(model{iter,1},'EX_xyl_D(e)');
    
    % for community-specific nutrient state, change bounds of glucose and
    % xylose accordingly (uncomment lines 40 to 45)
    %     if glu{iter,1} ~= 0
    %         model{iter,1} = changeRxnBounds(model{iter,1},{'EX_glc_D(e)'},glcup(iter,1),'l');
    %     end
    %     if xyl{iter,1} ~=0
    %         model{iter,1} = changeRxnBounds(model{iter,1},{'EX_xyl_D(e)'},xylup(iter,1),'l');
    %     end
    %FBA
    sol{iter,1} = optimizeCbModel(model{iter,1},'max');
    Growth{iter,1} = round(sol{iter,1}.f,4);
    
    %check uptake fluxes of glucose and xylose after optimization of
    %biomass
    if glu{iter,1} ~= 0
        glcuptake{iter,1} = sol{iter,1}.x(glu{iter,1});
    end
    if xyl{iter,1} ~=0
        xyluptake{iter,1} = sol{iter,1}.x(xyl{iter,1});
    end
    
    %FVA
    
    if sol{iter,1}.f > 1e-02
        [minFlux{iter,1}, maxFlux{iter,1}] = fluxVariability(model{iter,1},100,'max',RxnName);
        
        %product flux from FBA solution and maxFlux from FVA solution
        Productrxn{iter,1} = findRxnIDs(model{iter,1},RxnName);
        fbaflux{iter,1} = round(sol{iter,1}.x(Productrxn{iter,1}), 4);
        maxflux{iter,1} = round(maxFlux{iter,1}, 4);
    end
end