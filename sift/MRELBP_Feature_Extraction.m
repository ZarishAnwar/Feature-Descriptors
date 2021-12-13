clear all;

global imgPath;
global sampleNames;
global sampleIdxClass;
global allTrainSamples;
global texDatabaseName;
global outexProb;
global lbpRadius;
global lbpPoints;
global numLBPbins;
global lbpMethod;
global blockSize;

blockSize = 5; % Parametr Initialization
numLBPbins = 0; % Parametr Initialization
lbpPoints = 0; % Parametr Initialization

% This is the default parameter setting used in our TIP paper
% It gives fairly good results.
samplingScheme = 'EightSch1';
lbpRadiusSet = [2 4 6 8];
lbpMethodAll{1} = 'MELBP';
lbpMethodAll{1} = strcat(lbpMethodAll{1},samplingScheme);


if strcmp(samplingScheme,'EightSch1')
    % This is the default parameter setting used in our TIP paper
    % It gives fairly good results.
    lbpPointsSet = [8 8 8 8]; 
elseif strcmp(samplingScheme,'Vary')
    lbpPointsSet = [8 16 24 24];
else
    error('Not Such Case!');
end
% *************************************************************************

for idxLBPmethod = 1
    lbpMethod = lbpMethodAll{idxLBPmethod};
    filePathSaved = 'outex_dataset.txt';
    texDatabaseName = 'outex';
    outexProbSet = get_outexProbSet();
    for idxLbpRadius = 1 : length(lbpRadiusSet)
        for idxNumProbs = 3 : 4 
            outexProb = outexProbSet{idxNumProbs};
            imgPath = get_img_path(filePathSaved,idxNumProbs);
            [subPath,allTrainSamples] = get_subPath(outexProb);
            sampleNamesWithClass = textread(strcat(imgPath,subPath),'%s',allTrainSamples*2+1);
            sampleNames = sampleNamesWithClass(2:2:end);
            sampleIdxClassStr = sampleNamesWithClass(3:2:end);
            sampleIdxClass = zeros(1,length(sampleIdxClassStr));
            for i = 1 : length(sampleIdxClassStr)
                sampleIdxClass(i) = str2double(sampleIdxClassStr{i});
            end
            % ==================================================
            lbpRadius = lbpRadiusSet(idxLbpRadius);
            lbpPoints = lbpPointsSet(idxLbpRadius);
            mapping = get_mapping_info(lbpRadius,lbpPoints);
            if idxLbpRadius > 1
                lbpRadiusPre = lbpRadiusSet(idxLbpRadius-1);
            else
                lbpRadiusPre = 0;
            end

            if strcmp(lbpMethod,'MELBPEightSch1')
                MRELBP_Outex(mapping,lbpRadiusPre);
            end
            
        end % loop idxLbpRadius
    end % loop idxOutexProb
end % loop idxLBPmethod




