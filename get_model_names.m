function [file_names] = get_model_names(folder_name)
%GET_MODEL_NAMES Returns an array of names of the .mat files inside
%a given folder
file_info = dir(strcat(strcat('./',folder_name),'/*.mat'));
file_names = ({file_info.name});
end
