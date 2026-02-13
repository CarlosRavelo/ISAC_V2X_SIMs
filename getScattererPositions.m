function scattererParams = getScattererPositions(numScatterers, streetParams, rcsParams)
    % obtain the positions of the scatterers on a street scenario and the
    % rcs of each scatterer
    % Input Parameters:
    % numScatterers             :           number of scatterers in the
    %                                       scenario
    % streetParams              :           structure with three fields, street width
    %                                       street length and sidewalk width 
    % rcsParams                 :           structure with two fields, rcs
    %                                       mean and rcs variance
    % Output Parameters:
    % scattererParams           :           structure with fields of
    %                                       position and rcs

    % process the position
    scattererPos = zeros(2, numScatterers);
    
    streetWidth = streetParams.streetWidth;
    sideWalkWidth = streetParams.sideWalkWidth;
    streetLength = streetParams.streetLength;

    scattererPos(1, :) = streetWidth/2 + sideWalkWidth*rand(numScatterers, 1);
    scattererPos(1, floor(numScatterers/2)+1:end) = -scattererPos(1, floor(numScatterers/2)+1:end);
    scattererPos(2, :) = randi(streetLength, 1, numScatterers);

    scattererParams.Positions = scattererPos;

    % process the rcs
    rcsScatterers = rcsParams.rcsVar*randn(1, numScatterers) + rcsParams.rcsMean;

    scattererParams.RCSs = rcsScatterers;
end