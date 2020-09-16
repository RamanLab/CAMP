%To generate cross-fes mets table for correlation analysis in crossfedmets_corr.R 
% array containing all metabolites cross-fed between all communities
M = {};
M= M(~cellfun('isempty',M));
%identify unique mets and their occurences
U= unique(M);
v = cellfun(@(x) sum(ismember(M,x)),U, 'UniformOutput', false);
%E is an array of org pairs and the cross-fed metabolites in each pair 
E={};
indexPre = zeros(length(U), length(E));
nonEmptyCells = ~cellfun('isempty',E);
for iter = 1:length(E)
    colmNo = nnz(nonEmptyCells(iter,:));
    for iter2 = 1:length(U)
        if ~isempty(strmatch(U{iter2},E(iter,3:colmNo)))
            indexPre(iter2,iter) = 1;
        end
    end
end
indexPre1 = indexPre';
metstable = array2table(indexPre1);
%save table as csv file 
filename = 'mets_corr.csv';
writetable(metstable, filename);
%change column headers as the metabolite names
