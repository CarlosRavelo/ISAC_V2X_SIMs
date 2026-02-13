function velRad = getRadialVelocity(posTx, posScatterer, velTx, velScatterer)
    % 1. obtain the joining vector
    joiningVector = posScatterer - posTx;
    joiningVector = joiningVector/norm(joiningVector);
    % 2. project source velocity into joining vector
    vRelTx = dot(joiningVector, velTx);   
    % 3. project scatterer velocity into joining vector
    velRelScat = dot(joiningVector, velScatterer);
    % 4. obtain radial velocity
    velRad = vRelTx - velRelScat;
end