function [modelCom] = communityuptake(modelCom,reactionsAndBounds)
% INPUTS:
%   modelCom               : community model
%   reactionsAndBounds : Cell array with reaction names and corresponding bound values
%                           e.g., {'EX_ala_L(u)', -1;
%                           'EX_arg_L(u)', -1; ...} 'u' refers to
%                           the community compartment
% Apply bounds to reactions
for i = 1:size(reactionsAndBounds, 1)
    reaction = reactionsAndBounds{i, 1};
    bound_value = reactionsAndBounds{i, 2};
    modelCom = changeRxnBounds(modelCom, {reaction}, bound_value, 'l');
end
end
