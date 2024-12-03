function [candidate_Amensal, candidate_Parasitic, candidate_Competition, candidate_Neutral, candidate_Commensal, candidate_Mutualistic, candidate] = ComputeCommunityInteraction(orglist_mono, growth_mono, maxproductyield_mono, org1, org2, g1, g2, s1, s2, comm_FVAmaxflux)

% orglist_mono is a nx1 cell array that lists names of n organisms
% growth_mono is a nx1 double with monoculture growth rates of n organisms
% maxproductyield_mono is the monoculture max product yield yi (product flux/total substrate uptake) of n orgs
% org1 and org2 are cell arrays that have names of organism 1 and organism 2 in each community
% g1 and g2 are cell arrays with growth rates of each org in the community simulations obtained in communitygrowth.m
% s1 and s2 are cell arrays with substrate uptake fluxes for each org in the pair (if 2 substrates, e.g. sum of glucose and
% xylose uptake fluxes in each org)
% comm_FVAmaxflux is the FVA product maxflux value in each community


orgpair = horzcat(org1,org2);
% Calculate 10% higher and lower growth from monoculture growth rates for all org
for n = 1:length(orglist_mono)
    per_high{n,1} = growth_mono(n) + (growth_mono(n) * 0.1);
    per_low{n,1} = growth_mono(n) - (growth_mono(n) * 0.1);
end

% Concatenate arrays
org = horzcat(orglist_mono, maxproductyield_mono, per_low, per_high);

% Product yield of each organism in monoculture
y1 = cell2mat(cellfun(@(x) org{strcmp(x, org(:,1)), 2}, org1, 'UniformOutput', false));
y2 = cell2mat(cellfun(@(x) org{strcmp(x, org(:,1)), 2}, org2, 'UniformOutput', false));

% Calculate Moi for each community 
% Moi is FVA product maxflux value in the community divided by the total substrate
% uptake of community(s1+s2)
for k= 1 :size(s1,1)
    totalsubstrate_uptake{k} = {s1{k}+s2{k}};
    Moi{k,1} = comm_FVAmaxflux(k)/cell2mat(totalsubstrate_uptake{k});
end


% Adapting the ConYE model (Medlock et.al,2018)
% Calculate expected yield of lactate from each community and compare with observed yield
foldchange = 2;
for m = 1:length(Moi)
    Mei{m,1} = (s1{m,1} * y1(m)) + (s2{m,1} * y2(m));
    if abs(Moi{m,1}) >= (foldchange * abs(Mei{m,1}))
        fluxincrease(m,1) = Moi{m,1} - Mei{m,1};
        bestpairs{m,1} = orgpair(m,:);
    end
end

% Extract candidate pairs which have 10*higher than expected lactate yield
t = find(~cellfun(@isempty, bestpairs));
candidate = orgpair(t,:);

% Compare co-culture growth rates of each org to 10% higher or 10% lower growth of monoculture
l1 = cell2mat(cellfun(@(x) org{strcmp(x, org(:,1)), 3}, org1, 'UniformOutput', false));
l2 = cell2mat(cellfun(@(x) org{strcmp(x, org(:,1)), 3}, org2, 'UniformOutput', false));
h1 = cell2mat(cellfun(@(x) org{strcmp(x, org(:,1)), 4}, org1, 'UniformOutput', false));
h2 = cell2mat(cellfun(@(x) org{strcmp(x, org(:,1)), 4}, org2, 'UniformOutput', false));

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

% Identify interaction type in each pair
z = horzcat(c1,c2);
Interaction = cell(size(z));
for q = 1:size(z,1)
    if z{q,1} == -1 && z{q,2} == -1
        Interaction{q,1} = 'Competition';
    elseif z{q,1} == 1 &&  z{q,2} == -1
        Interaction{q,1} = 'Parasitism';
    elseif z{q,1} == -1 && z{q,2} == 1
        Interaction{q,1} = 'Parasitism';
    elseif z{q,1} == -1 && z{q,2} == 0
        Interaction{q,1} = 'Amensalism';
    elseif z{q,1} == 0 && z{q,2} == -1
        Interaction{q,1} = 'Amensalism';
    elseif z{q,1} == 1 && z{q,2} == 1
        Interaction{q,1} = 'Mutualism';
    elseif z{q,1} == 1 && z{q,2} == 0
        Interaction{q,1} = 'Commensalism';
    elseif z{q,1} == 0 && z{q,2} == 1
        Interaction{q,1} = 'Commensalism';
    elseif z{q,1} == 0 && z{q,2} == 0
        Interaction{q,1} = 'Neutral';
    end
end

% Count no. of pairs with particular interaction type
% Interaction_unique = unique(Interaction);
% for i = 1:length(Interaction_unique)
%     count_interaction(i) = sum(strcmp(Interaction, Interaction_unique{i}));
% end

% Extract candidate pairs for each interaction type
candidate_Amensal = extract_candidate(orgpair, candidate, Interaction, 'Amensalism');
candidate_Parasitic = extract_candidate(orgpair, candidate, Interaction, 'Parasitism');
candidate_Competition = extract_candidate(orgpair, candidate, Interaction, 'Competition');
candidate_Neutral = extract_candidate(orgpair, candidate, Interaction, 'Neutral');
candidate_Commensal = extract_candidate(orgpair, candidate, Interaction, 'Commensalism');
candidate_Mutualistic = extract_candidate(orgpair, candidate, Interaction, 'Mutualism');
end

function candidates = extract_candidate(orgpair,candidate,Interaction,type)
indices = strcmp(Interaction, type);
candidates = orgpair(indices,:);
[~,ia,~] = unique(candidates, 'rows', 'stable');
if ia ~= 0
    candidates = candidates(ia(1),:);
end
end
