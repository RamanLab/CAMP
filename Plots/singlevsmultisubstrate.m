
%A,B,C etc are the six interaction types ..no.of such interactions in each of the nine
%conditions(glucose-excess, glucose-minimal, glucose-community-specific etc)
A = [214,55,5;215,64,19;23,1,2];
B = [100,38,19;93,42,26;1,2,3];
C = [4,261,226;0,224,235;78,68,95];
D = [30,5,31;30,2,35;0,0,0];
E = [0,45,44;0,62,19;70,31,34];
F = [0,110,65;0,98,20;1,89,71];
Z = cat(3,A,B,C,D,E,F);  %create a 3-D matrix

stackData = Z; 
groupLabels = {'Glucose(E,M,CS)' 'Glucose+Xylose(E,M,CS)' 'Xylose(E,M,CS)'}; 
h = plotBarStackGroups(stackData, groupLabels);
% Change the colors of each bar segment
colors = jet(size(h,2)); %or define your own color order; 1 for each m segments
colors = repelem(colors,size(h,1),1); 
colors = mat2cell(colors,ones(size(colors,1),1),3);
set(h,{'FaceColor'},colors)

