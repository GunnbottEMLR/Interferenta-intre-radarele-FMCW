function fmcwwav = helperFMCWWaveform(fc,maxrange,rangeres,maxvel)
c = physconst('LightSpeed');                  % viteza luminii (m/s)
lambda = freq2wavelen(fc,c);                  % lungimea de unda (m)

%  chirp duration =5*Round-Trip Travel Time,adica cat ii ia ca radarul FMCW
%  sa faca un sweep adica timpul de transmisie si receprie=5*timpul de
%  transmisie si inapoi
tm = 5*range2time(maxrange,c);                % Chirp duration (s)

% determinarea largimii benzii
bw = rangeres2bw(rangeres,c);                 % largimea benzii (Hz)

% Set the sampling rate to satisfy both range and velocity requirements
sweepSlope = bw/tm;                           % FMCW sweep slope=cet de repede creste frecv pe secunda (Hz/s)
fbeatMax = range2beat(maxrange,sweepSlope,c); % Maximum beat frequency=distanta maxima intre transmitator si reciever (Hz)
fdopMax = speed2dop(2*maxvel,lambda);         % Maximum Doppler shift=schimbul de frecventa datorita miscarii obiectului (Hz)
fifMax = fbeatMax+fdopMax;                    % Maximum received IF=cea mai mare frecventa posibila receptionata de radar (Hz)
fs = max(2*fifMax,bw);                        % Sampling rate=cat de repede colecteaza radarul date (Hz)

%facerea FCMW waveform
fmcwwav = phased.FMCWWaveform('SweepTime',tm,'SweepBandwidth',bw,'SampleRate',fs,'SweepDirection','Up');
end

