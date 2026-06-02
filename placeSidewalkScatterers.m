function scatterers = placeSidewalkScatterers(numScatterers, streetParams, rcsParams)
    % PLACESIDEWALKSCATTERERS places scatterers along the sidewalk of a
    % street
    % scatterers = PLACESSIDEWALKSCATTERERS(numScatterers, streetParams,
    %                                       rcsParams)
    % Generates scattererers for simulations. Example of random scatterer
    % generation. Feel free to replace with your own model.
    % Inputs:
    %       numScatterers       : number of scatterers
    %       streetParams        : structure witht the street parameters.
    %                             Expected fields are streetWidth, sidewalkWidth 
    %                             and streetLength
    %       rcsParams           : structure with the rcs parameters.
    %                             Assumes the rcs of scatterers will be modeled 
    %                             as a normal variable. Expected fields are 
    %                             rcsMean and rcsVar
    % Output:
    %       scatterers          : structure with fields positions, velocities and RCSs.
    %                             positions is a matrix of 2 x
    %                             numScatterers with [x, y] of each
    %                             scatterer. RCSs is an array of 1 x
    %                             numScatterers
    scattererPos = zeros(numScatterers, 2);
    
    streetWidth = streetParams.streetWidth;
    sidewalkWidth = streetParams.sidewalkWidth;
    streetLength = streetParams.streetLength;

    scattererPos(:, 1) = streetWidth/2 + sidewalkWidth*rand(numScatterers, 1);
    scattererPos(floor(numScatterers/2)+1:end, 1) = -scattererPos(floor(numScatterers/2)+1:end, 1);
    scattererPos(:, 2) = randi(streetLength, numScatterers, 1);

    scatterers.positions = scattererPos;
    % static scatterers
    scatterers.velocities = zeros(numScatterers, 2);

    rcsScatterers = rcsParams.rcsVar*randn(1, numScatterers) + rcsParams.rcsMean;

    scatterers.RCSs = rcsScatterers;    
end