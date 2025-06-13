function ipaths = helperGenerateIntPaths(tgtPoses,actProf,lambda)
tgtPos = reshape([zeros(1,3) tgtPoses.Position],3,[]);
tgtVel = reshape([zeros(1,3) tgtPoses.Velocity],3,[]);

% Position point targets at half of each target's height
tgtPos(3,:) = tgtPos(3,:)+0.5*[actProf.Height];
% Number of Targets
Ntargets = length(tgtPos(:,1))-1;

% Interference paths to victim radar
ipaths = repmat(struct(...
    'PathLength', zeros(1, 1), ...
    'PathLoss', zeros(1, 1), ...
    'ReflectionCoefficient', zeros(1,1), ...
    'AngleOfDeparture', zeros(2, 1), ...
    'AngleOfArrival', zeros(2, 1), ...
    'DopplerShift', zeros(1, 1)),...
    1,Ntargets);

% Each target is already in the path length, angle
[plengthicar2vcar,tangicar2vcar] = rangeangle(tgtPos(:,3),tgtPos(:,1),eye(3));
[plengthicar2tgt,tangicar2tgt] = rangeangle(tgtPos(:,3),tgtPos(:,2),eye(3));
[plengthtgt2vcar,tangtgt2vcar] = rangeangle(tgtPos(:,2),tgtPos(:,1),eye(3));

% path loss
plossicar2vcar = fspl(plengthicar2vcar,lambda);  % from interfering radar to victim radar
plossicar2tgt = fspl(plengthicar2tgt,lambda);    % from interfering radar to target
plosstgt2vcar = fspl(plengthtgt2vcar,lambda);    % from target to victim radar

% reflection gain
rgainicar2vcar = aperture2gain(0.01,lambda);
rgainicar2tgt = aperture2gain(actProf(2).RCSPattern(1),lambda);

% Doppler
dopicar2vcar = speed2dop(2*radialspeed(tgtPos(:,3),tgtVel(:,3),tgtPos(:,1),tgtVel(:,1)),...
    lambda);

% Generate interference path from interfering radar to victim radar
ipaths(1).PathLength = plengthicar2vcar;
ipaths(1).PathLoss = plossicar2vcar;
ipaths(1).ReflectionCoefficient = db2mag(rgainicar2vcar);
ipaths(1).AngleOfDeparture = tangicar2vcar;
ipaths(1).AngleOfArrival = tangicar2vcar;
ipaths(1).DopplerShift = dopicar2vcar;

% Generate interference path from interfering radar to target and from target to victim radar
ipaths(2).PathLength = plengthicar2tgt+plengthtgt2vcar;
ipaths(2).PathLoss = plossicar2tgt+plosstgt2vcar;
ipaths(2).ReflectionCoefficient = db2mag(rgainicar2tgt);
ipaths(2).AngleOfDeparture = tangicar2tgt;
ipaths(2).AngleOfArrival = tangtgt2vcar;
ipaths(2).DopplerShift = dopicar2vcar;
end
