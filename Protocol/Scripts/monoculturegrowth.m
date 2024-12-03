function[Growth,fbaflux,maxflux,substrate_uptakeflux,maxproductyield_mono]=monoculturegrowth(foldername,RxnName,substrateRxns,reactionsAndBounds)
%This computes the monoculture growth rates of single organism models and
%the fluxes from FBA and FVA solutions for the particular reaction of interest
% INPUTS:
%   foldername :name of folder containing a list of COBRA models as
%               .mat files
%   RxnName  :exchange rxn name of the metabolite of interest in the model
%   e.g. 'EX_lac_D(e)'
%   reactionsAndBounds : Cell array with reaction names and corresponding bound values
%                           e.g., {'EX_ala_L(e)', -1; 'EX_arg_L(e)', -1; ...}
%   substrateRxns : substrates reaction names e.g. 'EX_glc_D(e)'
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
    model{iter,1}= monocultureuptake(model{iter,1},reactionsAndBounds);
    ATPM{iter,1} = findRxnIDs(model{iter,1},'DM_atp_c_');
    model{iter,1}.lb(ATPM{iter,1}) = 0.36;
    model{iter,1}.ub(ATPM{iter,1}) = 0.36;
    
    %identify substrate rxns
    substrateRxn = cell(size(substrateRxns));
    for j = 1:length(substrateRxn)
        substrateRxn{j} = substrateRxns{j};
        substrate{iter,j} = findRxnIDs(model{iter,1},substrateRxn{j});
    end
    
    %FBA
    sol{iter,1} = optimizeCbModel(model{iter,1},'max');
    Growth{iter,1} = round(sol{iter,1}.f,4);
    
    %check uptake fluxes of substrate rxn after optimization of
    %biomass
    
    for j = 1:length(substrateRxns)
        if substrate{iter,j} ~= 0
            substrate_uptakeflux{iter,j} = sol{iter,1}.x(substrate{iter,j});
        else
            substrate_uptakeflux{iter, j} = NaN;
        end
    end
    
    %FVA
    if sol{iter,1}.f > 1e-02
        [minFlux{iter,1}, maxFlux{iter,1}] = fluxVariability(model{iter,1},100,'max',RxnName);
        
        %product flux from FBA solution and maxFlux from FVA solution
        Productrxn{iter,1} = findRxnIDs(model{iter,1},RxnName);
        fbaflux{iter,1} = round(sol{iter,1}.x(Productrxn{iter,1}), 4);
        maxflux{iter,1} = round(maxFlux{iter,1}, 4);
        maxproductyield_mono{iter,1} = abs(maxflux{iter,1}/substrate_uptakeflux{iter,1});
    end
end

    function [file_names] = get_model_names(folder_name)
        %GET_MODEL_NAMES Returns an array of names of the .mat files inside
        %a given folder
        file_info = dir(strcat(strcat('./',folder_name),'/*.mat'));
        file_names = ({file_info.name});
    end
    function [model] = monocultureuptake(model, reactionsAndBounds)
        % INPUTS:
        %   model               : COBRA model to be adjusted
        %   reactions_with_bounds : Cell array with reaction names and corresponding bound values
        %                           e.g., {'EX_ala_L(e)', -1; 'EX_arg_L(e)', -1; ...}
        % Apply bounds to reactions
        for i = 1:size(reactionsAndBounds, 1)
            reaction = reactionsAndBounds{i, 1};
            bound_value = reactionsAndBounds{i, 2};
            model = changeRxnBounds(model, {reaction}, bound_value, 'l');
        end
    end
filename = 'monoculture_results.xlsx';
D = horzcat(Files,Growth,fbaflux,maxflux,substrate_uptakeflux,maxproductyield_mono);
D = cell2table(D);
D.Properties.VariableNames ={'Organism';'Growth';'fbaflux';'product_maxflux';'substrate_uptakeflux';'maxproductyield'};
writetable(D, filename,'Sheet',1);
end
