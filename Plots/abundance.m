%To generate abundance table used for plotting box plots of relative abundance of each organism in
%relativeabundance.R
%org is an array of all single model names
%orgabd is the names of each org in the community and their relative abundances (4 columns)
org = {};
orgabd = {};
for iter = 1:length(org)
        row = strmatch(org{iter,1},orgabd(:,1),'exact');
   if ~isempty(row)
        abd{iter,1} = orgabd(row,3);
   end
        row2 = strmatch(org{iter,1},orgabd(:,2),'exact');
   if ~isempty(row2)
        abd{iter,2} = orgabd(row2,4);
   end 
end
%concatenate the abundance entries for each org
for iter1 = 1:length(org)
   ABD{iter1,1} = cat(1,abd{iter1,1},abd{iter1,2});
end
%abundance table for each org
ABDtable=[];
for iter2 = 1:length(org)
    iter3 = 1;
    while iter3 <= length(ABD{iter2,1})
        ABDtable =[ABDtable;char2cell(org{iter2,1}), ABD{iter2,1}{iter3,1}];
        iter3 = iter3 + 1;
    end
end
ABDtable = cell2table(ABDtable);
ABDtable.Properties.VariableNames ={'Species','RelativeAbundance'};
%save table as xlsx file 
filename = 'GXminimal_abd.xlsx';
writetable(ABDtable, filename);