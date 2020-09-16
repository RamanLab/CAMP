%To compute interaction type of each community based on growth change
% list names of all single organisms 
orglist = {};
% list names of organism 1 and organism 2 in each community  
org1 = {};
org2 = {};
orgpair = horzcat(org1,org2);
% monoculture growth rates of all orgs
g =[];
%co-culture growth rates of each org as g1 and g2
g1 = {};
g2 = {};
%calculate 10% higher and 10% lower growth rate from monoculture growth rates for all orgs 
for n = 1 :length(orglist)
    per_high{n,1} = g(n) +(g(n)*0.1);    
    per_low{n,1} = g(n) -(g(n)*0.1);  
end
% concatenate arrays
orginfo = horzcat(orglist,per_low,per_high);
%compare co-culture growth rates of each org to 10%higher or 10%lower growth
%of monoculture
l1 = cell2mat(cellfun(@(x) orginfo{strmatch(x,orginfo(:,1),'exact'),2},org1,'UniformOutput', false));
l2 = cell2mat(cellfun(@(x) orginfo{strmatch(x,orginfo(:,1),'exact'),2},org2,'UniformOutput', false));
h1 = cell2mat(cellfun(@(x) orginfo{strmatch(x,orginfo(:,1),'exact'),3},org1,'UniformOutput', false));
h2 = cell2mat(cellfun(@(x) orginfo{strmatch(x,orginfo(:,1),'exact'),3},org2,'UniformOutput', false));
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
% define interaction type in each pair based on the following criteria
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
    Competitivepair(iter4,:) = orgpair(u,:);
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
