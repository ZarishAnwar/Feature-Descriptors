% Folder=strcat('','D:\MATLAB\New Folder\CVPR09-ScSPM\ScSPM\data\KTH_TIPS','')
% for i=1:10
% C=[i];
% FileList=dir(fullfile(Folder,'*.mat'));
% FileList.name
% for iFile = 1:numel(FileList)  
%  Data= load(fullfile(Folder,FileList(iFile).name))
% histo=(Data.feaSet.feaArr);
% B=imhist(histo);
% B=B.';
% B=horzcat(B, C);
% Thisto=[Thisto ; B];
% end
% 
% end

Thisto=[];
for i=1:10
% Folder='/Users/uettaxila/Desktop/TestRP/Brod111//',i;
Folder=strcat('','D:\MATLAB\New Folder\CVPR09-ScSPM\ScSPM\data\KTH_TIPS\',num2str(i))
% histo=[];

C=[i]
FileList = dir(fullfile(Folder, '*.mat'));  % List of all MAT files
FileList.name
for iFile = 1:numel(FileList)               % Loop over found files
  Data   = load(fullfile(Folder, FileList(iFile).name));
  histo=(Data.feaSet.feaArr);
%      B = histo.';
 B=imhist(histo);
 B=B.';
     B =horzcat(B , C);
  Thisto = [Thisto ; B];
  
end
end

