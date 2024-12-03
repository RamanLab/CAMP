 function create_community_models(input_folder, output_folder)
% This function is used to create all possible 2-species community models
% using a list of model files
%  input_folder :name of folder containing a list of COBRA models as
%               .mat files  
%  output_folder :name of folder which will be used to store the community
%  model .mat files

    % Get model names from the input folder
    v = get_model_names(input_folder);
    C = nchoosek(v, 2);

    % Initialize variables to store models, options, and results
    org1 = cell(length(C), 1);
    org2 = cell(length(C), 1);
    spBm = cell(2, length(C));
    spATPM = cell(2, length(C));
    options = struct('spBm', {}, 'spATPM', {}, 'sepUtEx', {}, 'minNorm', {});

    % Generate community models and save them as .mat files
    for iter = 1:length(C)
        % Load the models
        modelpair = {C{iter, 1}, C{iter, 2}};
        org1{iter, 1} = readCbModel(fullfile(input_folder, modelpair{1}));
        org2{iter, 1} = readCbModel(fullfile(input_folder, modelpair{2}));
        modelCell = {org1{iter, 1}, org2{iter, 1}};

        % Identify biomass and ATPM reactions
        spBm{1, iter} = org1{iter, 1}.rxns(find(org1{iter, 1}.c));
        spBm{2, iter} = org2{iter, 1}.rxns(find(org2{iter, 1}.c));
        spATPM{1, iter} = validatestring('DM_atp_c_', org1{iter, 1}.rxns);
        spATPM{2, iter} = validatestring('DM_atp_c_', org2{iter, 1}.rxns);
        minNorm = 1;

        % Create options structure containing species biomass and ATPM fields
        options(iter).spBm = spBm(:, iter);
        options(iter).spATPM = spATPM(:, iter);
        options(iter).sepUtEx = false;
        options(iter).minNorm = minNorm;

        disp(iter);

        % Use SteadyCom function createCommModel to create a community model
        [modelCom{iter, 1}, infoCom{iter, 1}, indCom{iter, 1}] = createCommModel(modelCell, options(iter));
        %modelCom{iter, 1}.description = strrep(modelCom{iter, 1}.description, '.mat', '');
        modelCom{iter, 1}.modelID = modelCom{iter, 1}.description;
        modelCom{iter, 1}.modelName = modelCom{iter, 1}.description;
        modelComNew = modelCom{iter, 1};
        modelComNew.subSystems = modelCom{iter, 1}.subSystems;
        infoComNew = infoCom{iter, 1};
        indComNew = indCom{iter, 1};
        positions = strfind(modelComNew.description, '.mat');
        if length(positions) > 1
            modelComNew.description = eraseBetween(modelComNew.description, positions(1), positions(end - 1) + 3);
        end
        %modelComNew.description = strrep(modelComNew.description,'.mat','');
        filename = fullfile(output_folder, modelComNew.description);
        save(filename, 'modelComNew', 'infoComNew', 'indComNew', '-mat');
    end
end

function [file_names] = get_model_names(folder_name)
    % GET_MODEL_NAMES Returns an array of names of the .mat files inside
    % a given folder
    file_info = dir(fullfile(folder_name, '*.mat'));
    file_names = {file_info.name};
end
