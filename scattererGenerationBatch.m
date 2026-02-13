%% Generate multiple scatterer locations

totalRealizations = 100;
scatterersNum = [3, 6, 9, 12];

% define the street parameters
streetParams.streetWidth = 5;
streetParams.streetLength = 70;
streetParams.sideWalkWidth = 3;

% define the rcs parameters
rcsParams.rcsMean = 5;
rcsParams.rcsVar = 2;

% generate all the realizations for each number of scatterers
for i = scatterersNum
    scattererParams = cell(totalRealizations,1);
    for j = 1:totalRealizations
        scattererParams{j} = getScattererPositions(i, streetParams, rcsParams);
    end
    % write to file
    filename = strcat('scattererParams', int2str(i), '.mat');
    save(filename, 'scattererParams');
end