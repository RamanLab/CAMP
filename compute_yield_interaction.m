%To compute product yield in the community and interaction types
% list names of all single organisms 
orglist = {};
% monoculture growth rates of all orgs
g =[];
%monoculture max lactate yield yi (product flux/total substrate uptake) of all orgs
maxlacflux = {};

% list names of organism 1 and organism 2 in each community  
org1 = {};
org2 = {};
orgpair = horzcat(org1,org2);
%co-culture growth rates of each org as g1 and g2
g1 = {};
g2 = {};
%substrate uptake fluxes for each org in the pair (if 2 substrates, sum of glucose and
%xylose uptake fluxes in each org)
s1 = {};
s2 = {};

%calculate 10% higher and lower growth from monoculture growth rates for all org 
for n = 1 :length(orglist)
    per_high{n,1} = g(n) +(g(n)*0.1);    
    per_low{n,1} = g(n) -(g(n)*0.1);  
end
% concatenate arrays
org = horzcat(orglist,maxlacflux,per_low,per_high);

%lactate yield of each organism in monoculture as y1 and y2 
y1 = cell2mat(cellfun(@(x) org{strmatch(x,org(:,1),'exact'),2},org1,'UniformOutput', false));
y2 = cell2mat(cellfun(@(x) org{strmatch(x,org(:,1),'exact'),2},org2,'UniformOutput', false));

% lactate FVA maxflux value in the community divided by the total substrate
% uptake of community(s1+s2)
Moi= {};

%Adapting the ConYE model (Medlock et.al,2018)
%calculate expected yield of lactate from each community and compare with observed yield 
for m = 1:length(Moi)
    Mei{m,1} = (s1{m,1}*y1(m))+(s2{m,1}*y2(m));
    if Moi{m,1} > (10*Mei{m,1}) 
        fluxincrease(m,1) = cell2mat(Moi(m,1))-cell2mat(Mei(m,1));
        bestpairs{m,1} = {orgpair(m,:)};
    end
end

%extract candidate pairs which have 10*higher than expected lactate yield 
t = find(~cellfun(@isempty,bestpairs));
for iter=1:length(t)
    j = t(iter,1);
    candidate(iter,:) = orgpair(j,:);
end

%compare co-culture growth rates of each org to 10%higher or 10%lower growth
%of monoculture
l1 = cell2mat(cellfun(@(x) org{strmatch(x,org(:,1),'exact'),3},org1,'UniformOutput', false));
l2 = cell2mat(cellfun(@(x) org{strmatch(x,org(:,1),'exact'),3},org2,'UniformOutput', false));
h1 = cell2mat(cellfun(@(x) org{strmatch(x,org(:,1),'exact'),4},org1,'UniformOutput', false));
h2 = cell2mat(cellfun(@(x) org{strmatch(x,org(:,1),'exact'),4},org2,'UniformOutput', false));
for p = 1:length(g1)
    if g1{p,1} < l1(p)
        c1{p,1} = -1;
    elseif g1{p,1} > h1(p)
        c1{p,1} = 1;
    elseif g1{p,1} == l1(p) || g1{p,1} == h1(p)
        c1{p,1} = 0;
    else
        c1{p,1} = 0;
    end  
end
for b = 1:length(g2)    
    if g2{b,1} < l2(b)
        c2{b,1} = -1;
    elseif g2{b,1} > h2(b)
        c2{b,1} = 1;
    elseif g2{b,1} == l2(b) || g2{b,1} == h2(b)
        c2{b,1} = 0;
    else
        c2{b,1} = 0;
    end
end
% identify interaction type in each pair 
z = horzcat(c1,c2);
for q = 1: length(z)
    if z{q,1} == -1 && z{q,2} == -1
        Interaction{q,1} = 'Competition';
    end
    if z{q,1} == 1 &&  z{q,2} == -1
        Interaction{q,1} = 'Parasitism';
    end
    if z{q,1} == -1 && z{q,2} == 1
        Interaction{q,1} = 'Parasitism';
    end
    if z{q,1} == -1 && z{q,2} == 0
        Interaction{q,1} = 'Amensalism';
    end
    if z{q,1} == 0 && z{q,2} == -1
        Interaction{q,1} = 'Amensalism';
    end
    if z{q,1} == 1 && z{q,2} == 1
        Interaction{q,1} = 'Mutualism';
    end
    if z{q,1} == 1 && z{q,2} == 0
        Interaction{q,1} = 'Commensalism';
    end
    if z{q,1} == 0 && z{q,2} == 1
        Interaction{q,1} = 'Commensalism';
    end
    if z{q,1} == 0 && z{q,2} == 0
        Interaction{q,1} = 'Neutral';
    end
end
% count no. of pairs with particular interaction type
A = count(Interaction, 'Amensalism');
Amensal= find(A==1);
M = count(Interaction, 'Mutualism');
Mutualistic = find(M==1);
CM = count(Interaction, 'Competition');
Competition = find(CM==1);
P = count(Interaction, 'Parasitism');
Parasitic = find(P ==1);
CO = count(Interaction, 'Commensalism');
Commensal = find(CO==1);
N = count(Interaction, 'Neutral');
Neutral = find(N==1);

%viable pairs and their interaction type
for iter2=1:length(Amensal)
    r = Amensal(iter2,1);
    Amensalpair(iter2,:) = orgpair(r,:);
end
for iter3=1:length(Parasitic)
    s = Parasitic(iter3,1);
    Parasiticpair(iter3,:) = orgpair(s,:);
end
for iter4=1:length(Competition)
    u = Competition(iter4,1);
    Competitionpair(iter4,:) = orgpair(u,:);
end
for iter5=1:length(Neutral)
    v = Neutral(iter5,1);
    Neutralpair(iter5,:) = orgpair(v,:);
end
for iter6=1:length(Commensal)
    x = Commensal(iter6,1);
    Commensalpair(iter6,:) = orgpair(x,:);
end
for iter7=1:length(Mutualistic)
    e = Mutualistic(iter7,1);
    Mutualisticpair(iter7,:) = orgpair(e,:);
end
% extract candidate pairs and their interaction types
candidate_Amensal = {};
for index = 1:length(candidate)
    rowNos = strmatch(candidate{index,1},Amensalpair(:,1),'exact');
    if ~isempty(rowNos)
        if strmatch(candidate{index,2},Amensalpair(rowNos,2),'exact')
            candidate_Amensal = [candidate_Amensal;candidate(index,:)];
        end
    end
end

candidate_Parasitic = {};
for index1 = 1:length(candidate)
    rowNos1 = strmatch(candidate{index1,1},Parasiticpair(:,1),'exact');
    if ~isempty(rowNos1)
        if strmatch(candidate{index1,2},Parasiticpair(rowNos1,2),'exact')
            candidate_Parasitic = [candidate_Parasitic;candidate(index1,:)];
        end
    end
end
candidate_Competition = {};
for index2 = 1:length(candidate)
    rowNos2 = strmatch(candidate{index2,1},Competitionpair(:,1),'exact');
    if ~isempty(rowNos2)
        if strmatch(candidate{index2,2},Competitionpair(rowNos2,2),'exact')
            candidate_Competition = [candidate_Competition;candidate(index2,:)];
        end
    end
end

candidate_Neutral = {};
for index3 = 1:length(candidate)
    rowNos3 = strmatch(candidate{index3,1},Neutralpair(:,1),'exact');
    if ~isempty(rowNos3)
        if strmatch(candidate{index3,2},Neutralpair(rowNos3,2),'exact')
            candidate_Neutral = [candidate_Neutral;candidate(index3,:)];
        end
    end
end
candidate_Commensal = {};
for index4 = 1:length(candidate)
    rowNos4 = strmatch(candidate{index4,1},Commensalpair(:,1),'exact');
    if ~isempty(rowNos4)
        if strmatch(candidate{index4,2},Commensalpair(rowNos4,2),'exact')
            candidate_Commensal = [candidate_Commensal;candidate(index4,:)];
        end
    end
end
candidate_Mutualistic = {};
for index5 = 1:length(candidate)
    rowNos5 = strmatch(candidate{index5,1},Mutualisticpair(:,1),'exact');
    if ~isempty(rowNos5)
        if strmatch(candidate{index5,2},Mutualisticpair(rowNos5,2),'exact')
            candidate_Mutualistic = [candidate_Mutualistic;candidate(index5,:)];
        end
    end
end
