function [FESOFRxnList, KOtable] = FSEOFKnock(model,RxnName,biomassRxn)
%Predicts single or double knockouts in a community model for optimal production of target
%metabolite
% INPUTS:
%       model :community model 
%       RxnName :name of the target product Rxn
%       biomassRxn :cell array containing the RxnIDs of the biomass rxn of
%                   each organism
% OUTPUTS:
%       FSEOFRxnList :rxns whose fluxes decreases as product flux
%                     increases
%       KOtable :table of predicted knockouts 
% FBA and wild-type fluxes
fba = optimizeCbModel(model,'max','one');
productFluxWT = fba.x(strcmp(model.rxns, RxnName));
vbiomass = fba.f;
flux = fba.x;

%change objective to product rxn & compute mutant fluxes
modelnew = changeObjective(model, RxnName);
biomassorg1 = findRxnIDs(modelnew,biomassRxn{1,1});
biomassorg2 = findRxnIDs(modelnew,biomassRxn{2,1});
fbanew= optimizeCbModel(modelnew,'max','one');
vmaxpdt = (abs(fbanew.f));
vipdt = productFluxWT;
fluxnew = fbanew.x;
allrxns = length(model.rxns);
n = 10;
vepdt = [];
decflux =[];

%FSEOF
for k = 1:n
    vepdt(k) = vipdt + (k/n)*(vmaxpdt - vipdt);
    modelF = model;
    modelF = changeRxnBounds(modelF, RxnName, vepdt(k), 'b');
    fbaMT = optimizeCbModel(modelF, 'max','one');
    if fbaMT.stat ==1
        fluxMT(:,k) = fbaMT.x;
    end
end
% identify rxns with monotonically decreasing fluxes as product flux
% increases
for r = 1:allrxns
    if abs(flux(r)) > abs(fluxMT(r, 1:10))
        decflux = [decflux; r];
    end
end
decflux_fluxes = fluxMT(decflux,:);
decflux_rxns = model.rxns(decflux);
decflux_rxnnames = model.rxnNames(decflux);
decflux_rxnformula = printRxnFormula(model, model.rxns(decflux));

%identify exchange rxns and filter the list
[transRxns, nonTransRxns, transRxnsBool] = findTransRxns(model,'true');
FESOFRxnList = setdiff(decflux_rxns,transRxns);

%identify double or single rxn deletions from the selectedRxnList that
%increase product flux while retaining viability(>= 0.01 h-1 growth rate) of both orgs in the
%community
KO =[];
threshold = 30; %Can be changed depending on the number of rxns in the FSEOFRxnList

if length(FESOFRxnList) <= threshold
    C = nchoosek(FESOFRxnList,2);  %all possible combinations of 2 rxn deletions
    for q = 1:length(C)
        modelnew1 = model;
        rxnRemoveList1 = C{q,1}; 
        modelnew1 = removeRxns(modelnew1, rxnRemoveList1); %delete first rxn in the combination
        M1Biomassnew1=find(ismember(modelnew1.rxns, biomassRxn{1,1}));
        M2Biomassnew1=find(ismember(modelnew1.rxns, biomassRxn{2,1}));
        modelnew1.c(M1Biomassnew1)=1;
        modelnew1.c(M2Biomassnew1)=1;
        fbaMT1 = optimizeCbModel(modelnew1,'max','one');
        if fbaMT1.stat ==1
            productFluxMT1{q,1} = fbaMT1.x(strcmp(modelnew1.rxns, RxnName));
            growthRateMT1{q,1} = fbaMT1.f;
            growthOrg1_MT1{q,1} = fbaMT1.x(M1Biomassnew1);
            growthOrg2_MT1{q,1} = fbaMT1.x(M2Biomassnew1);
            if productFluxMT1{q,1} ~= 0 && round(productFluxMT1{q,1},4) > round(productFluxWT) && growthOrg1_MT1{q,1} >= 0.01 && growthOrg2_MT1{q,1} >= 0.01
                P1 = C(q,1);
                KOlist1{q,1} = horzcat(P1,productFluxMT1{q,1},growthRateMT1{q,1},growthOrg1_MT1{q,1},growthOrg2_MT1{q,1});
            else 
                KOlist1{q,1} = [];
            end
        end
        modelnew2 = model;
        rxnRemoveList2 = C{q,2};
        modelnew2 = removeRxns(modelnew2, rxnRemoveList2); %delete second rxn in the combination
        M1Biomassnew2=find(ismember(modelnew2.rxns, biomassRxn{1,1}));
        M2Biomassnew2=find(ismember(modelnew2.rxns, biomassRxn{2,1}));
        modelnew2.c(M1Biomassnew2)=1;
        modelnew2.c(M2Biomassnew2)=1;
        fbaMT2 = optimizeCbModel(modelnew2,'max','one');
        if fbaMT2.stat ==1
            productFluxMT2{q,1} = fbaMT2.x(strcmp(modelnew2.rxns, RxnName));
            growthRateMT2{q,1} = fbaMT2.f;
            growthOrg1_MT2{q,1} = fbaMT2.x(M1Biomassnew2);
            growthOrg2_MT2{q,1} = fbaMT2.x(M2Biomassnew2);
            if productFluxMT2{q,1} ~= 0 && round(productFluxMT2{q,1},4) > round(productFluxWT) && growthOrg1_MT2{q,1} >= 0.01 && growthOrg2_MT2{q,1} >= 0.01
                P2 = C(q,2);
                KOlist2{q,1} = horzcat(P2,productFluxMT2{q,1},growthRateMT2{q,1},growthOrg1_MT2{q,1},growthOrg2_MT2{q,1});
            else 
                KOlist2{q,1} = [];
            end
        end
        modelnew3 = model;
        rxnRemoveList3 = {C{q,1},C{q,2}};
        modelnew3 = removeRxns(modelnew3, rxnRemoveList3); %delete both rxns in the combination
        M1Biomassnew3=find(ismember(modelnew3.rxns, biomassRxn{1,1}));
        M2Biomassnew3=find(ismember(modelnew3.rxns, biomassRxn{2,1}));
        modelnew3.c(M1Biomassnew3)=1;
        modelnew3.c(M2Biomassnew3)=1;
        fbaMT3 = optimizeCbModel(modelnew3,'max','one');
        if fbaMT3.stat ==1
            productFluxMT3{q,1} = fbaMT3.x(strcmp(modelnew3.rxns, RxnName));
            growthRateMT3{q,1} = fbaMT3.f;
            growthOrg1_MT3{q,1} = fbaMT3.x(M1Biomassnew3);
            growthOrg2_MT3{q,1} = fbaMT3.x(M2Biomassnew3);
            if productFluxMT3{q,1} ~= 0 && round(productFluxMT3{q,1},4) > round(productFluxWT) && growthOrg1_MT3{q,1} >= 0.01 && growthOrg2_MT3{q,1} >= 0.01
                P3 = C(q,:);
                KOlist3{q,1} = horzcat(P3,productFluxMT3{q,1},growthRateMT3{q,1},growthOrg1_MT3{q,1},growthOrg2_MT3{q,1});
            else 
                KOlist3{q,1} = [];
            end
        end
% select only the deletion strategy that improves product flux, either
% deletion of first rxn or second rxn or both
        if round(productFluxMT1{q,1},4) > round(productFluxMT2{q,1},4) && round(productFluxMT1{q,1},4) > round(productFluxMT3{q,1},4)
            KOlist{q,1} = KOlist1{q,1};
        elseif round(productFluxMT2{q,1},4) > round(productFluxMT1{q,1},4) && round(productFluxMT2{q,1},4) > round(productFluxMT3{q,1},4)
            KOlist{q,1} = KOlist2{q,1};
        elseif round(productFluxMT3{q,1},4) > round(productFluxMT2{q,1},4) && round(productFluxMT3{q,1},4) > round(productFluxMT1{q,1},4)
            KOlist{q,1} = KOlist3{q,1};
        elseif round(productFluxMT1{q,1},4) == round(productFluxMT2{q,1},4) && round(productFluxMT3{q,1},4) >= round(productFluxMT2{q,1},4) 
            KOlist{q,1} = KOlist3{q,1};
        else
            KOlist{q,1} = [];
        end
    end
    KO= KOlist(~cellfun('isempty',KOlist));
    if isempty(KO)
        error('NO Knockouts found');  %Cannot identify Knockouts for the target metabolite of interest from this model
    else
        KOs={};
        for iter = 1:length(KO)
            KOs{iter} = cell(1,6);
            if length(KO{iter}) == 5
                KOs{iter}{1} = KO{iter}{1};
                KOs{iter}{2} = '';
                KOs{iter}(3:6) = KO{iter}(2:5);
            else
                KOs{iter} = KO{iter};
            end
        end
        KOs=KOs';
        KOtable=[];
        for iter2 = 1:length(KOs)
            KOtable = [KOtable; KOs{iter2,:}];
        end
        KOtable = unique(cell2table(KOtable),'rows');
        KOtable.Properties.VariableNames ={'Rxn1';'Rxn2';'ProductFlux';'CommunityGrowth';'GrowthOrg1';'GrowthOrg2'};
        KOtable = sortrows(KOtable,{'ProductFlux'},{'descend'});
    end
else
    for q = 1:length(FESOFRxnList) %perform single rxn deletions only
        modelnew = model;
        rxnRemoveList = FESOFRxnList(q);
        modelnew = removeRxns(modelnew, rxnRemoveList);
        M1Biomassnew=find(ismember(modelnew.rxns, biomassRxn{1,1}));
        M2Biomassnew=find(ismember(modelnew.rxns, biomassRxn{2,1}));
        modelnew.c(M1Biomassnew)=1;
        modelnew.c(M2Biomassnew)=1;
        fbaMT = optimizeCbModel(modelnew,'max','one');
        if fbaMT.stat ==1
            productFluxMT{q,1} = fbaMT.x(strcmp(modelnew.rxns, RxnName));
            growthRateMT{q,1} = fbaMT.f;
            growthOrg1{q,1} = fbaMT.x(M1Biomassnew);
            growthOrg2{q,1} = fbaMT.x(M2Biomassnew);
            if productFluxMT{q,1} ~= 0 && round(productFluxMT{q,1},4) > round(productFluxWT) && growthOrg1{q,1} >= 0.01 && growthOrg2{q,1} >= 0.01
                P = FESOFRxnList(q);
                KOlist{q,1} = horzcat(P,productFluxMT{q,1},growthRateMT{q,1},growthOrg1{q,1},growthOrg2{q,1});
            else
                KOlist{q,1} = [];
            end
        end
    end
    KO= KOlist(~cellfun('isempty',KOlist));
    if isempty(KO)
        error('NO Knockouts found');  %Cannot identify Knockouts for the target metabolite of interest from this model
    else
        KOtable={};
        for iter = 1:length(KO)
            KOtable = [KOtable; KO{iter,:}];
        end
        KOtable = unique(cell2table(KOtable),'rows');
        KOtable.Properties.VariableNames ={'RxnID';'ProductFlux';'CommunityGrowth';'GrowthOrg1';'GrowthOrg2'};
        KOtable = sortrows(KOtable,{'ProductFlux'},{'descend'});
    end
end