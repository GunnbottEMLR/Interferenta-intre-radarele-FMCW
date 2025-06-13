%Define Victim Baseband FMCW Radar Waveform.
%For illustration purposes, in this example, configure the radar to a maximum range of 150 m. The victim radar operates at a 77 GHz frequency. The radar is required to resolve objects that are at least 1 meter apart. Because this is a forward-facing radar application, the maximum relative speed of vehicles on highway is specified as 230 km/h.

%Radarul Victima Colaterala

%FMCW-Frequency-Modulated Continuous Wave(modulare in frecventa unda-
%-continua.

rng(2023);             % cand folosesc functii random asta face ca functia mereu sa dea acelasi rezultat
fc = 77e9;             % frecventa centrala
c = physconst('LightSpeed');             % viteza luminii (m/s)
lambda = freq2wavelen(fc,c);             % lambda/lungimea undei (m)
rangeMax = 150;        % distanta maxima la care poate radaru sa vada (m)
rangeRes = 1;          % distanta minima la care obiectele trebuie sa fie separate ca radarul sa le rezolve (m)
vMax = 230*1000/3600;  % viteza maxima relativa a masinilor (m/s)
fmcwwav1 = helperFMCWWaveform(fc,rangeMax,rangeRes,vMax);% call la functia helperFMCWWaveform
global sig;
sig = fmcwwav1(); %stocheaza in sig
Ns = numel(sig);%cate sampleluri avem





%%
%Define Interfering Baseband FMCW Radar Waveform

%Radarul care interfereaza

fcRdr2 = 77e9;         % frecv centrala (Hz)
lambdaRdr2 = freq2wavelen(fcRdr2,c);              % lungimea de unda (m)
rangeMaxRdr2 = 100;    %  distanta maxima la care poate radaru sa vada(m)
rangeResRdr2 = 0.8;     % distanta minima la care obiectele trebuie sa fie separate ca radarul sa le rezolve (m)
vMaxRdr2 = vMax; % viteza maxima relativa a masinilor  (m/s)
global fmcwav2;
 fmcwwav2 = helperFMCWWaveform(fcRdr2,rangeMaxRdr2,rangeResRdr2,vMaxRdr2);% call la functia helperFMCWWaveform
sigRdr2 = fmcwwav2();%stocheaza in sig

%%
%Plotarea formei de unda de baza a radarului victima in primele 3 pulsatii

 figure %deschide o fereastra noua pentru a pune chestii pe ea
 pspectrum(repmat(sig,3,1),fmcwwav1.SampleRate,'spectrogram', ...
     'Reassign',true,'FrequencyResolution',10e6)
 axis([0 15 -80 80]); title('Victim Radar'); colorbar off

%calculeaza si afiseaza spectograma semnalului sig,ce e in repmat face ca
%semnalul sig sa se repete de 3 ori

%%
%Plotarea formei de unda de baza a radarului interferenta in primele 3 pulsatii

figure %deschide o fereastra noua pentru a pune chestii pe ea
pspectrum(repmat(sigRdr2,3,1),fmcwwav1.SampleRate,'spectrogram',...
    'Reassign',true,'FrequencyResolution',10e6,'MinThreshold',-13)
axis([0 15 -80 80]); title('Interfering Radar'); colorbar off;
%calculeaza si afiseaza spectograma semnalului sigRdr,ce e in repmat face ca
%semnalul sigRdr sa se repete de 3 ori


%%
% Modelul transceiverelor radar MIMO: victimă și interferent
% Se consideră că ambele radare utilizează un array liniar uniform (ULA) atât pentru transmisie, 
% cât și pentru recepția undelor radar.

% Utilizarea unui array de recepție permite radarului să estimeze direcția azimutală a energiei 
% reflectate de la ținte potențiale.

% Utilizarea unui array de transmisie permite radarului să formeze un array virtual mare la recepție, 
% îmbunătățind astfel rezoluția unghiului azimutal.

% Pentru a recepționa semnalele de la țintă și zgomotul la radarul victimă, modelați transceiverul
% monostatic al radarului victimă folosind radarTransceiver și specificând array-urile ULA de transmisie 
% și recepție ale radarului victimă în proprietățile sale.



% Setez parametrii de baza ai transceiverului radar
antAperture = 6.06e-4;                          % Deschiderea antenei (m^2)
antGain = aperture2gain(antAperture,lambda);    % Castigul antenei (dB)
txPkPower = db2pow(13)*1e-3;                    % Puterea de varf a emitătorului (W)
rxNF = 4.5;                                     % Figura de zgomot a receptorului (dB)

% Construiește array-ul de recepție pentru radarul victima
Nvr =  16;     % Numarul de elemente de recepție pentru radarul victima
vrxEleSpacing = lambda/2;        % Distanța intre elementele de receptie pentru radarul victima
antElmnt = phased.IsotropicAntennaElement('BackBaffled',false); %Creează un element de antenă izotrop.
% O antenă izotropă este o antenă ideală care radiază în mod uniform în toate direcțiile. În realitate, 
% astfel de antene nu există, dar ele sunt utile în modelarea teoretică a radiației antenei.
vrxArray = phased.ULA('Element',antElmnt,'NumElements',Nvr,...
    'ElementSpacing',vrxEleSpacing); % Creează un Uniform Linear Array (ULA), adică o aranjare liniară 
% de antene. Un ULA este o configurație de antene în linie dreaptă, folosită frecvent în radare și sisteme 
% de comunicații pentru a capta semnale din diferite direcții.

% Construiește array-ul de transmisie pentru radarul victimă
Nvt = 2;     % Numărul de elemente de transmisie pentru radarul victimă
vtxEleSpacing = Nvr*vrxEleSpacing;     % Distanța între elementele de transmisie pentru radarul victimă
vtxArray = phased.ULA('Element',antElmnt,'NumElements',Nvt,...
    'ElementSpacing',vtxEleSpacing);% Creează un Uniform Linear Array (ULA), adică o aranjare liniară de antene. 
% Un ULA este o configurație de antene în linie dreaptă, folosită frecvent în radare și sisteme de comunicații 
% pentru a capta semnale din diferite direcții.

% Modelează transceiverul monostatic al radarului victimă pentru recepționarea semnalelor de țintă și zgomot
vradar = radarTransceiver("Waveform",fmcwwav1,'ElectronicScanMode','Custom'); %Creează un obiect radarTransceiver
% care reprezintă un transceiver radar.  'ElectronicScanMode', Configurează modul de scanare electronică al radarului
% ca fiind personalizat, indicând că radarul poate să-și modifice directivitatea antenei pentru a urmări ținta.
vradar.Transmitter = phased.Transmitter('PeakPower',txPkPower,'Gain',antGain);%Creează obiectul Transmitter al 
% radarului, care reprezintă emițătorul radarului.
vradar.TransmitAntenna = phased.Radiator('Sensor',vtxArray,'OperatingFrequency',fc,'WeightsInputPort',true);%Creează 
% obiectul TransmitAntenna care reprezintă antena de transmisie a radarului. 'Sensor', Asociază antena de transmisie
% cu array-ul de antene de transmisie vtxArray, care a fost definit anterior. 'WeightsInputPort', true Permite aplicarea
% de cântare pentru elementele array-ului antenei, care sunt necesare pentru a direcționa semnalul într-o anumită direcție.
vradar.ReceiveAntenna = phased.Collector('Sensor',vrxArray,'OperatingFrequency',fc); % Creează obiectul ReceiveAntenna,
% care reprezintă antena de recepție a radarului. 'Sensor',Asociază antena de transmisie cu array-ul de antene de transmisie
% vtxArray, care a fost definit anterior.
vradar.Receiver = phased.ReceiverPreamp('Gain',antGain,'NoiseFigure',rxNF,'SampleRate',fmcwwav1.SampleRate);%Creează obiectul
% Receiver, care reprezintă receptorul radarului. 'NoiseFigure' Setează figura de zgomot a receptorului, care indică nivelul
% de zgomot pe care receptorul îl introduce în semnalul de recepție. rxNF este o valoare dată (în dB).



%%

%În diferit față de semnalele de țintă, interferența este transmisă de la ULA-ul de transmisie al radarului de 
% interferență și recepționată de ULA-ul de recepție al radarului victimă. Modelează transceiverul de interferență
% folosind funcția helperInterferenceTransceiver, care include ULA-ul de transmisie al radarului de interferență și ULA-ul
% de recepție al radarului victimă ca intrări.


% Construiește array-ul de transmisie pentru radarul de interferență
Nit = 3;     % Numărul de elemente de transmisie pentru radarul de interferență
Nir = 4;     % Numărul de elemente de recepție pentru radarul de interferență
irxEleSpacing = lambda/2;        % Distanța între elementele de recepție pentru radarul de interferență
itxEleSpacing = Nir*irxEleSpacing;     % Distanța între elementele de transmisie pentru radarul de interferență
itxArray = phased.ULA('Element',antElmnt,'NumElements',Nit,...
    'ElementSpacing',itxEleSpacing);%Creează un Uniform Linear Array (ULA), adică o configurație de antene plasate pe o
% linie dreaptă. Un ULA este folosit pentru a direcționa semnalele de transmisie sau recepție într-o anumită direcție și 
% pentru a colecta semnalele de la diverse unghiuri. În acest caz, ULA-ul va fi utilizat pentru transmisia semnalului de 
% interferență. 'Element' Acesta specifică tipul de elemente de antenă folosite în array. antElmnt este un obiect definit
% anterior în cod (probabil un element de antenă izotropă sau o altă configurație de antenă).

% Modelează transceiverul de interferență pentru recepționarea semnalului de interferență la radarul victimă
iradar = helperInterferenceTransceiver(fmcwwav2,txPkPower,antGain,fc,itxArray,vrxArray);

%%
%Simulate Driving Scenario




% Create driving scenario
[scenario, egoCar] = helperAutoDrivingScenario;
%activare/afisare
drivingScenarioDesigner(scenario)
%%

% Inițializează pozițiile țintelor în cadrul de referință al vehiculului ego
tgtPoses = targetPoses(egoCar);

% Distanța și unghiul vehiculului țintă relativ la radarul victimă
[tarRange, tarAngle] = rangeangle(tgtPoses(1).Position');

% Viteza radială a vehiculului țintă relativ la radarul victimă (semnul negativ
% este din cauza axelor de referință diferite între radialspeed și drivingScenarioDesigner)
tarVelocity = -radialspeed(tgtPoses(1).Position', tgtPoses(1).Velocity');

% Distanța și unghiul vehiculului de interferență relativ la radarul victimă
[intRange, intAngle] = rangeangle(tgtPoses(2).Position');

% Viteza radială a vehiculului de interferență relativ la radarul victimă
intVelocity = -radialspeed(tgtPoses(2).Position', tgtPoses(2).Velocity');

%%

%Transceiver Signal Processing




%% Obtain Interference-free and Interfered Data Cubes
% Definește eșantioanele de timp rapid și timp lent
Nft = round(fmcwwav1.SweepTime*fmcwwav1.SampleRate);  % Numărul de eșantioane
% de timp rapid pentru radarul victimă
Nsweep = 192;                                         % Numărul de eșantioane
% de timp lent pentru radarul victimă
Nift = round(fmcwwav2.SweepTime*fmcwwav2.SampleRate); % Numărul de eșantioane 
% de timp rapid pentru radarul de interferență
Nisweep = ceil(Nft*Nsweep/Nift);                      % Numărul de eșantioane 
% de timp lent pentru radarul de interferență

% Generează codul TDM-MIMO pentru elementele antenei de interferență
wi = helperTDMMIMOEncoder(Nit, Nisweep);

% Asamblează trenul de pulsiuni de interferență la timpul de scenariul radarului
% de interferență
rxIntTrain = zeros(Nift*Nisweep, Nvr);

% Inițializează timpul de scenariul
time = scenario.SimulationTime;

% Inițializează profilele actorilor
actProf = actorProfiles(scenario);

% Obține trenul de pulsiuni de interferență la timpul de scenariul radarului 
% de interferență
for l = 1:Nisweep
    % Generează toate drumurile către radarul victimă
    ipaths = helperGenerateIntPaths(tgtPoses, actProf, lambdaRdr2);

    % Obține semnalul primit doar cu interferență
    rxInt = iradar(ipaths, time, wi(:, l));
    rxIntTrain((l-1)*Nift+1:l*Nift, :) = rxInt;

    % Obține timpul curent al scenariului
    time = time + fmcwwav2.SweepTime;

    % Mută țintele înainte în timp pentru următorul sweep
    tgtPoses(1).Position = [tgtPoses(1).Position] + [tgtPoses(1).Velocity] * fmcwwav2.SweepTime;
    tgtPoses(2).Position = [tgtPoses(2).Position] + [tgtPoses(2).Velocity] * fmcwwav2.SweepTime;
end

% Trunchiază trenul de pulsiuni de interferență la lungimea trenului de pulsiuni 
% de țintă
rxIntTrain = rxIntTrain(1:Nft*Nsweep, :);

% Asamblează trenul de pulsiuni de interferență la timpul de scenariul radarului 
% victimă
rxIntvTrain = permute(reshape(rxIntTrain, Nft, Nsweep, Nvr), [1, 3, 2]);


%%

% Dechirpează semnalul țintei la radarul victimă pentru a obține cubul de date 
% fără interferență, 
% care este folosit pentru a trasa răspunsul țintei ca referință ideală de performanță. 
% De asemenea, dechirpează semnalul combinat de țintă și interferență pentru a 
% obține cubul de date cu interferență, 
% care este folosit pentru a trasa efectele interferenței asupra detectării țintei.

% Generează codul TDM-MIMO pentru elementele de antenă ale radarului victimă
wv = helperTDMMIMOEncoder(Nvt, Nsweep);

% Asamblează cubul de date fără interferență
XcubeTgt = zeros(Nft, Nvr, Nsweep);

% Asamblează cubul de date cu interferență
Xcube = zeros(Nft, Nvr, Nsweep);

% Inițializează timpul de scenariul
time = scenario.SimulationTime;

% Inițializează pozițiile țintelor în cadrul de referință al vehiculului ego
tgtPoses = targetPoses(egoCar);

% Inițializează profilele actorilor
actProf = actorProfiles(scenario);

% Obține cubul de date la timpul de scenariul radarului victimă
for l = 1:Nsweep
    % Generează drumurile țintelor către radarul victimă
    vpaths = helperGenerateTgtPaths(tgtPoses, actProf, lambda);

    % Obține semnalul primit fără interferență
    rxVictim = vradar(vpaths, time, wv(:, l));

    % Dechirpează semnalul țintei fără interferență
    rxVsig = dechirp(rxVictim, sig);

    % Salvează sweep-ul în cubul de date
    XcubeTgt(:,:,l) = rxVsig;

    % Obține semnalul de interferență primit
    rxInt = rxIntvTrain(:,:,l);

    % Dechirpează semnalul cu interferență
    rx = dechirp(rxInt + rxVictim, sig);
    Xcube(:,:,l) = rx;

    % Obține timpul curent al scenariului
    time = time + fmcwwav1.SweepTime;

    % Mută țintele înainte în timp pentru următorul sweep
    tgtPoses(1).Position = [tgtPoses(1).Position] + [tgtPoses(1).Velocity] * fmcwwav1.SweepTime;
    tgtPoses(2).Position = [tgtPoses(2).Position] + [tgtPoses(2).Velocity] * fmcwwav1.SweepTime;
end



%%
% Calculate number of range samples
Nrange = 2^nextpow2(Nft);

% Define range response
rngresp = phased.RangeResponse('RangeMethod','FFT', ...
    'SweepSlope',fmcwwav1.SweepBandwidth/fmcwwav1.SweepTime, ...
    'RangeFFTLengthSource','Property','RangeFFTLength',Nrange, ...
    'RangeWindow','Hann','SampleRate',fmcwwav1.SampleRate);

% Calculate the range response of interference-free data cube
XrngTgt = rngresp(XcubeTgt);
%%
% Decode TDM-MIMO waveform
[XdecTgt,Nsweep] = helperTDMMIMODecoder(XrngTgt,wv);

% Number of Doppler samples
Ndoppler = 2^nextpow2(Nsweep);

% Size of virtual array
Nv = Nvr*Nvt;

% Doppler FFT with Hann window
XrngdopTgt = helperVelocityResp(XdecTgt,Nrange,Nv,Ndoppler);




%% dopler map

% Pulse repetition interval for TDM-MIMO radar
tpri = fmcwwav1.SweepTime*Nvt;

% Maximum unambiguous velocity
vmaxunambg = lambda/4/tpri;

% Plot range-Doppler map
helperRangeVelocityMapPlot(XrngdopTgt,rngresp.SweepSlope,...
    fmcwwav1.SampleRate,tpri,rangeMax,vmaxunambg,Nrange,Ndoppler,lambda);

%% angle map

% Number of angle samples
Nangle = 2^nextpow2(Nv);

% Angle FFT with Hann window
XrngangdopTgt = helperAngleResp(XrngdopTgt,Nrange,Nangle,Ndoppler);

% Target Doppler bin
tarDopplerBin = ceil(tarVelocity/lambda*2*tpri*Ndoppler+Ndoppler/2);

% Plot range-angle map
helperRangeAngleMapPlot(XrngangdopTgt,rngresp.SweepSlope,...
    fmcwwav1.SampleRate,rangeMax,vrxEleSpacing,Nrange,Nangle,tarDopplerBin,lambda);
