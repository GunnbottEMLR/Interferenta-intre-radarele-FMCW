function vpaths = helperGenerateTgtPaths(tgtPoses,actProf,lambda)
tgtPos = reshape([zeros(1,3) tgtPoses.Position],3,[]);
tgtVel = reshape([zeros(1,3) tgtPoses.Velocity],3,[]);

% Position point targets at half of each target's height
tgtPos(3,:) = tgtPos(3,:)+0.5*[actProf.Height];
% Number of Targets
Ntargets = length(tgtPos(:,1))-1;

vpaths = repmat(struct(...
    'PathLength', zeros(1, 1), ...
    'PathLoss', zeros(1, 1), ...
    'ReflectionCoefficient', zeros(1,1), ...
    'AngleOfDeparture', zeros(2, 1), ...
    'AngleOfArrival', zeros(2, 1), ...
    'DopplerShift', zeros(1, 1)),...
    1,Ntargets);

for m = 1:Ntargets
    % Each target is already in the path length, angle
    [plength,tang] = rangeangle(tgtPos(:,m+1),tgtPos(:,1),eye(3));

    % path loss
    ploss = fspl(plength,lambda);

    % reflection gain
    rgain = aperture2gain(actProf(m+1).RCSPattern(1),lambda);

    % Doppler
    dop = speed2dop(2*radialspeed(tgtPos(:,m+1),tgtVel(:,m+1),tgtPos(:,1),tgtVel(:,1)),...
        lambda);

    % Target paths to victim radar
    vpaths(m).PathLength = 2*plength;
    vpaths(m).PathLoss = 2*ploss;
    vpaths(m).ReflectionCoefficient = db2mag(rgain);
    vpaths(m).AngleOfDeparture = tang;
    vpaths(m).AngleOfArrival = tang;
    vpaths(m).DopplerShift = dop;
end
end
