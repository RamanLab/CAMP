function [model] = monocultureuptake(model, reactions_with_bounds)
% INPUTS:
%   model               : COBRA model to be adjusted
%   reactions_with_bounds : Cell array with reaction names and corresponding bound values
%                           e.g., {'EX_ala_L(e)', -1; 'EX_arg_L(e)', -1; ...}
% Apply bounds to reactions
for i = 1:size(reactions_with_bounds, 1)
    reaction = reactions_with_bounds{i, 1};
    bound_value = reactions_with_bounds{i, 2};
    model = changeRxnBounds(model, {reaction}, bound_value, 'l');
end
end
