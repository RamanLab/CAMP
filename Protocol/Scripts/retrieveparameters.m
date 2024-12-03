function[orglist_mono, growth_mono, maxproductyield_mono, org1, org2, g1, g2, s1, s2, comm_FVAmaxflux] = retrieveparameters(monoculture_results,communityresults)
% function that retrieves data required for ComputeCommunityInteraction
% using the monoculture and community results .xlsx files that are outputs
% of monoculturegrowth and communitygrowth functions

% monoculture_results = 'monoculture_results.xlsx';
mono_data = readtable(monoculture_results, 'Sheet',1);
orglist_mono = mono_data.Organism;
growth_mono = mono_data.Growth;
maxproductyield_mono  = num2cell(mono_data.maxproductyield);

% communityresults = 'communityresults.xlsx';
comm_data = readtable(communityresults, 'Sheet', 1);
Communitymodel = comm_data.Communitymodel;
orgnames = cell(length(Communitymodel), 2);
org1 = cell(length(Communitymodel), 1);
org2 = cell(length(Communitymodel), 1);
for k = 1:length(Communitymodel)
    orgnames(k,1:2) = regexp(Communitymodel{k,1}, '&', 'split');
    org1{k,1} = orgnames{k,1};
    org2{k,1} = orgnames{k,2};
end
g1 = num2cell(comm_data.OrganismGrowthRate_1);
g2 = num2cell(comm_data.OrganismGrowthRate_2);
s1 = num2cell(comm_data.SubstrateUptakeFluxesForOrg1And2_1);
s2 = num2cell(comm_data.SubstrateUptakeFluxesForOrg1And2_2);
comm_FVAmaxflux = comm_data.ProductFlux_FVA;
end

