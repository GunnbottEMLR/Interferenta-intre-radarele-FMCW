function main_gui
    % Creează fereastra principală
    fig = uifigure('Name', 'PROIECT RC', 'Position', [100 200 1350 550]);
 global rangeRes
    rangeRes = 1; % Valoare inițială
ax = uiaxes(fig, ...
        'Position', [370 -100 500 550]);  % poziționează unde vrei în GUI
img = imread('1.png'); 
    imshow(img, 'Parent', ax);
ax = uiaxes(fig, ...
        'Position', [870 -100 500 1050]);  % poziționează unde vrei în GUI
img = imread('3.png'); 
    imshow(img, 'Parent', ax);
ax = uiaxes(fig, ...
        'Position', [870 -100 500 550]);  % poziționează unde vrei în GUI
img = imread('4.png'); 
    imshow(img, 'Parent', ax);

    lbl = uilabel(fig, ...
        'Text', 'Rezoluția radarului (m):', ...
        'Position', [100 160 200 22]);
lbl2 = uilabel(fig, ...
        'Text', ' Interferenta FMCW(unda continuu modulate prin frecventa)', ...
        'Position', [400 460 400 122]);
lbl3 = uilabel(fig, ...
        'Text', '  intre radarele automotive', ...
        'Position', [500 460 400 102]);
lbl4 = uilabel(fig, ...
        'Text', ' Radarele astea sunt folosite în mașinile moderne pentru chestii ca:', ...
        'Position', [400 460 400 52]);
lbl5 = uilabel(fig, ...
        'Text', ' --asistență la frânare', ...
        'Position', [400 460 400 32]);
lbl6 = uilabel(fig, ...
        'Text', ' --detectare obstacole, pietoni, alte mașini etc.', ...
        'Position', [400 360 400 202]);
lbl7 = uilabel(fig, ...
        'Text', ' În loc să trimită „pulse” scurte de semnal (ca unele radare militare)', ...
        'Position', [400 360 400 162]);
lbl8 = uilabel(fig, ...
        'Text', ', radarul FMCW trimite un semnal continuu care își schimbă frecvența ', ...
        'Position', [400 360 400 142]);
lbl9 = uilabel(fig, ...
        'Text', ' în timp (crește sau scade – asta se numește „chirp”).', ...
        'Position', [400 360 400 122]);
lbl9 = uilabel(fig, ...
        'Text', 'Imaginează-ți două mașini una lângă alta, fiecare cu radarul ei.', ...
        'Position', [400 360 400 82]);
lbl9 = uilabel(fig, ...
        'Text', '  Ambele emit semnale FMCW.Problema e că semnalele se pot amesteca.', ...
        'Position', [400 360 400 62]);
lbl9 = uilabel(fig, ...
        'Text', ' Radarul tău poate primi semnalul emis de altă mașină, crezând ', ...
        'Position', [400 360 400 42]);
lbl9 = uilabel(fig, ...
        'Text', 'că e al lui întors de la un obstacol.', ...
        'Position', [400 360 400 22]);

lbl9 = uilabel(fig, ...
        'Text', 'Producătorii de radare și mașini încearcă tot felul de soluții, gen:', ...
        'Position', [400 260 400 182]);
lbl9 = uilabel(fig, ...
        'Text', 'Frecvențe ușor diferite pentru fiecare radar – ca să nu se calce', ...
        'Position', [400 260 400 162]);
lbl9 = uilabel(fig, ...
        'Text', ' pe bătătură.Rezolutia radarului mai mica.Numarul de sweepuri Doppler', ...
        'Position', [400 260 400 142]);
lbl9 = uilabel(fig, ...
        'Text', ' mai multe', ...
        'Position', [400 260 400 122]);

    sld = uislider(fig, ...
        'Limits', [1 2], ...
        'Value', rangeRes, ...
        'Position', [100 140 200 3], ...
        'ValueChangedFcn', @(sld,event) updateRangeRes(sld));

    global dopplerSweeps
    dopplerSweeps = 192; % Valoare default

    uilabel(fig, ...
        'Text', 'Număr sweep-uri Doppler:', ...
        'Position', [100 420 200 22]);

    sldDoppler = uislider(fig, ...
        'Limits', [64 512], ...
        'MajorTicks', [64, 128, 192, 256, 384, 512], ...
        'Value', dopplerSweeps, ...
        'Position', [100 400 200 3], ...
        'ValueChangedFcn', @(sldDoppler, event) updateDopplerSweeps(sldDoppler));

 % Slider pentru offset de frecvență între radare
    global freqOffsetFactor
    freqOffsetFactor = 1; % valoare inițială = 1x (adică fără offset)

    uilabel(fig, ...
        'Text', 'Offset frecvență radar 2:', ...
        'Position', [100 280 200 22]);

    sldOffset = uislider(fig, ...
        'Limits', [1 100], ...
        'Value', freqOffsetFactor, ...
        'Position', [100 250 200 3], ...
        'ValueChangedFcn', @(sldOffset, event) updateFreqOffset(sldOffset));


    % Creează butonul
    btn = uibutton(fig, ...
        'Text', 'Start Simulare FMCW radar', ...
        'Position', [100 60 250 40], ...
        'ButtonPushedFcn', @(btn,event) start_project(fig));
end
function updateRangeRes(sld)
    global rangeRes
   rangeRes = round(sld.Value, 1);
    
end
function updateFreqOffset(sldOffset)
    global freqOffsetFactor
    freqOffsetFactor = sldOffset.Value;
end
function updateDopplerSweeps(sldDoppler)
    global dopplerSweeps
    dopplerSweeps = round(sldDoppler.Value);
end
function start_project(mainFig)
%Define Victim Baseband FMCW Radar Waveform.
%For illustration purposes, in this example, configure the radar to a maximum range of 150 m. The victim radar operates at a 77 GHz frequency. The radar is required to resolve objects that are at least 1 meter apart. Because this is a forward-facing radar application, the maximum relative speed of vehicles on highway is specified as 230 km/h.

%Radarul Victima Colaterala

%FMCW-Frequency-Modulated Continuous Wave(modulare in frecventa unda-
%-continua.

rng(2023);             % cand folosesc functii random asta face ca functia mereu sa dea acelasi rezultat
fc = 77e9;    % frecventa centrala
c = physconst('LightSpeed');             % viteza luminii (m/s)
global lambda
lambda = freq2wavelen(fc,c);             % lambda/lungimea undei (m)
 global rangeMax
rangeMax = 150;        % distanta maxima la care poate radaru sa vada (m)
global rangeRes  % distanta minima la care obiectele trebuie sa fie separate ca radarul sa le rezolve (m)
 
vMax = 230*1000/3600;  % viteza maxima relativa a masinilor (m/s)
global fmcwwav1;
fmcwwav1 = helperFMCWWaveform(fc,rangeMax,rangeRes,vMax);% call la functia helperFMCWWaveform
global sig;
sig = fmcwwav1(); %stocheaza in sig
Ns = numel(sig);%cate sampleluri avem





%%
%Define Interfering Baseband FMCW Radar Waveform

%Radarul care interfereaza

global freqOffsetFactor
fcRdr2 = 77e9 + (freqOffsetFactor - 1)*200e6;  % între 77.0 GHz și 77.2 GHz          % frecv centrala (Hz)
lambdaRdr2 = freq2wavelen(fcRdr2,c);              % lungimea de unda (m)
rangeMaxRdr2 = 100;    %  distanta maxima la care poate radaru sa vada(m)
rangeResRdr2 = 0.8;     % distanta minima la care obiectele trebuie sa fie separate ca radarul sa le rezolve (m)
vMaxRdr2 = vMax; % viteza maxima relativa a masinilor  (m/s)
global sigRdr2;
global fmcwwav2;

 fmcwwav2 = helperFMCWWaveform(fcRdr2,rangeMaxRdr2,rangeResRdr2,vMaxRdr2);% call la functia helperFMCWWaveform
sigRdr2 = fmcwwav2();%stocheaza in sig

%%
%%
% Modelul transceiverelor radar MIMO: victimă și interferent
% Se consideră că ambele radare utilizează un array liniar uniform (ULA) atât pentru transmisie, cât și pentru recepția undelor radar.

% Utilizarea unui array de recepție permite radarului să estimeze direcția azimutală a energiei reflectate de la ținte potențiale.

% Utilizarea unui array de transmisie permite radarului să formeze un array virtual mare la recepție, îmbunătățind astfel rezoluția unghiului azimutal.

% Pentru a recepționa semnalele de la țintă și zgomotul la radarul victimă, modelați transceiverul monostatic al radarului victimă folosind radarTransceiver și specificând array-urile ULA de transmisie și recepție ale radarului victimă în proprietățile sale.



% Setez parametrii de baza ai transceiverului radar
antAperture = 6.06e-4;                          % Deschiderea antenei (m^2)
antGain = aperture2gain(antAperture,lambda);    % Castigul antenei (dB)
txPkPower = db2pow(13)*1e-3;                    % Puterea de varf a emitătorului (W)
rxNF = 4.5;                                     % Figura de zgomot a receptorului (dB)

% Construiește array-ul de recepție pentru radarul victima
Nvr =  16;     % Numarul de elemente de recepție pentru radarul victima
global vrxEleSpacing
vrxEleSpacing = lambda/2;        % Distanța intre elementele de receptie pentru radarul victima
antElmnt = phased.IsotropicAntennaElement('BackBaffled',false); %Creează un element de antenă izotrop. O antenă izotropă este o antenă ideală care radiază în mod uniform în toate direcțiile. În realitate, astfel de antene nu există, dar ele sunt utile în modelarea teoretică a radiației antenei.
vrxArray = phased.ULA('Element',antElmnt,'NumElements',Nvr,...
    'ElementSpacing',vrxEleSpacing); % Creează un Uniform Linear Array (ULA), adică o aranjare liniară de antene. Un ULA este o configurație de antene în linie dreaptă, folosită frecvent în radare și sisteme de comunicații pentru a capta semnale din diferite direcții.

% Construiește array-ul de transmisie pentru radarul victimă
Nvt = 2;     % Numărul de elemente de transmisie pentru radarul victimă
vtxEleSpacing = Nvr*vrxEleSpacing;     % Distanța între elementele de transmisie pentru radarul victimă
vtxArray = phased.ULA('Element',antElmnt,'NumElements',Nvt,...
    'ElementSpacing',vtxEleSpacing);% Creează un Uniform Linear Array (ULA), adică o aranjare liniară de antene. Un ULA este o configurație de antene în linie dreaptă, folosită frecvent în radare și sisteme de comunicații pentru a capta semnale din diferite direcții.

% Modelează transceiverul monostatic al radarului victimă pentru recepționarea semnalelor de țintă și zgomot
vradar = radarTransceiver("Waveform",fmcwwav1,'ElectronicScanMode','Custom'); %Creează un obiect radarTransceiver care reprezintă un transceiver radar.  'ElectronicScanMode', Configurează modul de scanare electronică al radarului ca fiind personalizat, indicând că radarul poate să-și modifice directivitatea antenei pentru a urmări ținta.
vradar.Transmitter = phased.Transmitter('PeakPower',txPkPower,'Gain',antGain);%Creează obiectul Transmitter al radarului, care reprezintă emițătorul radarului.
vradar.TransmitAntenna = phased.Radiator('Sensor',vtxArray,'OperatingFrequency',fc,'WeightsInputPort',true);%Creează obiectul TransmitAntenna care reprezintă antena de transmisie a radarului. 'Sensor', Asociază antena de transmisie cu array-ul de antene de transmisie vtxArray, care a fost definit anterior. 'WeightsInputPort', true Permite aplicarea de cântare pentru elementele array-ului antenei, care sunt necesare pentru a direcționa semnalul într-o anumită direcție.
vradar.ReceiveAntenna = phased.Collector('Sensor',vrxArray,'OperatingFrequency',fc); % Creează obiectul ReceiveAntenna, care reprezintă antena de recepție a radarului. 'Sensor',Asociază antena de transmisie cu array-ul de antene de transmisie vtxArray, care a fost definit anterior.
vradar.Receiver = phased.ReceiverPreamp('Gain',antGain,'NoiseFigure',rxNF,'SampleRate',fmcwwav1.SampleRate);%Creează obiectul Receiver, care reprezintă receptorul radarului. 'NoiseFigure' Setează figura de zgomot a receptorului, care indică nivelul de zgomot pe care receptorul îl introduce în semnalul de recepție. rxNF este o valoare dată (în dB).



%%

%În diferit față de semnalele de țintă, interferența este transmisă de la ULA-ul de transmisie al radarului de interferență și recepționată de ULA-ul de recepție al radarului victimă. Modelează transceiverul de interferență folosind funcția helperInterferenceTransceiver, care include ULA-ul de transmisie al radarului de interferență și ULA-ul de recepție al radarului victimă ca intrări.


% Construiește array-ul de transmisie pentru radarul de interferență
Nit = 3;     % Numărul de elemente de transmisie pentru radarul de interferență
Nir = 4;     % Numărul de elemente de recepție pentru radarul de interferență
irxEleSpacing = lambda/2;        % Distanța între elementele de recepție pentru radarul de interferență
itxEleSpacing = Nir*irxEleSpacing;     % Distanța între elementele de transmisie pentru radarul de interferență
itxArray = phased.ULA('Element',antElmnt,'NumElements',Nit,...
    'ElementSpacing',itxEleSpacing);%Creează un Uniform Linear Array (ULA), adică o configurație de antene plasate pe o linie dreaptă. Un ULA este folosit pentru a direcționa semnalele de transmisie sau recepție într-o anumită direcție și pentru a colecta semnalele de la diverse unghiuri. În acest caz, ULA-ul va fi utilizat pentru transmisia semnalului de interferență. 'Element' Acesta specifică tipul de elemente de antenă folosite în array. antElmnt este un obiect definit anterior în cod (probabil un element de antenă izotropă sau o altă configurație de antenă).

% Modelează transceiverul de interferență pentru recepționarea semnalului de interferență la radarul victimă
iradar = helperInterferenceTransceiver(fmcwwav2,txPkPower,antGain,fc,itxArray,vrxArray);

%%
% Create driving scenario
global scenario
global egoCar
[scenario, egoCar] = helperAutoDrivingScenario;

%%

% Inițializează pozițiile țintelor în cadrul de referință al vehiculului ego
tgtPoses = targetPoses(egoCar);

% Distanța și unghiul vehiculului țintă relativ la radarul victimă
[tarRange, tarAngle] = rangeangle(tgtPoses(1).Position');

% Viteza radială a vehiculului țintă relativ la radarul victimă (semnul negativ este din cauza axelor de referință diferite între radialspeed și drivingScenarioDesigner)
tarVelocity = -radialspeed(tgtPoses(1).Position', tgtPoses(1).Velocity');

% Distanța și unghiul vehiculului de interferență relativ la radarul victimă
[intRange, intAngle] = rangeangle(tgtPoses(2).Position');

% Viteza radială a vehiculului de interferență relativ la radarul victimă
intVelocity = -radialspeed(tgtPoses(2).Position', tgtPoses(2).Velocity');

%%

%Transceiver Signal Processing




%% Obtain Interference-free and Interfered Data Cubes
% Definește eșantioanele de timp rapid și timp lent
Nft = round(fmcwwav1.SweepTime*fmcwwav1.SampleRate);  % Numărul de eșantioane de timp rapid pentru radarul victimă
global dopplerSweeps
    Nsweep = dopplerSweeps;                                      % Numărul de eșantioane de timp lent pentru radarul victimă
Nift = round(fmcwwav2.SweepTime*fmcwwav2.SampleRate); % Numărul de eșantioane de timp rapid pentru radarul de interferență
Nisweep = ceil(Nft*Nsweep/Nift);                      % Numărul de eșantioane de timp lent pentru radarul de interferență

% Generează codul TDM-MIMO pentru elementele antenei de interferență
wi = helperTDMMIMOEncoder(Nit, Nisweep);

% Asamblează trenul de pulsiuni de interferență la timpul de scenariul radarului de interferență
rxIntTrain = zeros(Nift*Nisweep, Nvr);

% Inițializează timpul de scenariul
time = scenario.SimulationTime;

% Inițializează profilele actorilor
actProf = actorProfiles(scenario);

% Obține trenul de pulsiuni de interferență la timpul de scenariul radarului de interferență
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

% Trunchiază trenul de pulsiuni de interferență la lungimea trenului de pulsiuni de țintă
rxIntTrain = rxIntTrain(1:Nft*Nsweep, :);

% Asamblează trenul de pulsiuni de interferență la timpul de scenariul radarului victimă
rxIntvTrain = permute(reshape(rxIntTrain, Nft, Nsweep, Nvr), [1, 3, 2]);


%%

% Dechirpează semnalul țintei la radarul victimă pentru a obține cubul de date fără interferență, 
% care este folosit pentru a trasa răspunsul țintei ca referință ideală de performanță. 
% De asemenea, dechirpează semnalul combinat de țintă și interferență pentru a obține cubul de date cu interferență, 
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
global Nrange
Nrange = 2^nextpow2(Nft);

% Define range response

global rngresp

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
global Ndoppler
Ndoppler = 2^nextpow2(Nsweep);

% Size of virtual array
Nv = Nvr*Nvt;

% Doppler FFT with Hann window
global XrngdopTgt


XrngdopTgt = helperVelocityResp(XdecTgt,Nrange,Nv,Ndoppler);

%% dopler map

% Pulse repetition interval for TDM-MIMO radar
global tpri
 tpri = fmcwwav1.SweepTime*Nvt;
global vmaxunambg
% Maximum unambiguous velocity
vmaxunambg = lambda/4/tpri;
%% angle map

global XrngangdopTgt

global Nangle
global tarDopplerBin
% Number of angle samples
Nangle = 2^nextpow2(Nv);

% Angle FFT with Hann window
XrngangdopTgt = helperAngleResp(XrngdopTgt,Nrange,Nangle,Ndoppler);

% Target Doppler bin
tarDopplerBin = ceil(tarVelocity/lambda*2*tpri*Ndoppler+Ndoppler/2);
%%
    % Închide fereastra principală
    close(mainFig);

    % Deschide fereastra secundară
    secondFig = uifigure('Name', 'Meniu Principal', 'Position', [200 200 450 350]);

    % Buton 1
    btn1 = uibutton(secondFig, ...
        'Text', 'Radarul Victima', ...
        'Position', [125 230 100 40], ...
        'ButtonPushedFcn', @(btn,event) functie1());

    % Buton 2
    btn2 = uibutton(secondFig, ...
        'Text', 'Radarul care interfereaza', ...
        'Position', [125 180 150 40], ...
        'ButtonPushedFcn', @(btn,event) functie2());

    % Buton 3
    btn3 = uibutton(secondFig, ...
        'Text', 'Simulare Scenariu Trafic', ...
        'Position', [125 130 150 40], ...
        'ButtonPushedFcn', @(btn,event) functie3());

    % Buton 4
    btn4 = uibutton(secondFig, ...
        'Text', 'Doppler Map', ...
        'Position', [125 70 100 40], ...
        'ButtonPushedFcn', @(btn,event) functie4());

% Buton 5
    btn5 = uibutton(secondFig, ...
        'Text', 'Angle Map', ...
        'Position', [125 20 100 40], ...
        'ButtonPushedFcn', @(btn,event) functie5());
btn6 = uibutton(secondFig, ...
        'Text', 'INF:AmbeleRadare', ...
        'Position', [325 130 100 40], ...
        'ButtonPushedFcn', @(btn,event) functie6());
btn7 = uibutton(secondFig, ...
        'Text', 'INF:DrivingScenario', ...
        'Position', [325 70 100 40], ...
        'ButtonPushedFcn', @(btn,event) functie7());
btn8 = uibutton(secondFig, ...
        'Text', 'INF:AngleDoplerMap', ...
        'Position', [325 20 100 40], ...
        'ButtonPushedFcn', @(btn,event) functie8());
end

function functie1()
global sig;
global fmcwwav1;
   figure %deschide o fereastra noua pentru a pune chestii pe ea
pspectrum(repmat(sig,3,1),fmcwwav1.SampleRate,'spectrogram', ...
    'Reassign',true,'FrequencyResolution',10e6)
axis([0 15 -80 80]); title('Victim Radar'); colorbar off
end

function functie2()

global sigRdr2;
global fmcwwav1;
figure %deschide o fereastra noua pentru a pune chestii pe ea
pspectrum(repmat(sigRdr2,3,1),fmcwwav1.SampleRate,'spectrogram',...
    'Reassign',true,'FrequencyResolution',10e6,'MinThreshold',-13)
axis([0 15 -80 80]); title('Interfering Radar'); colorbar off;
end

function functie3()
   %%
%Simulate Driving Scenario

global scenario
global egoCar



%activare/afisare
drivingScenarioDesigner(scenario)
end
function functie4()
global tpri
global fmcwwav1
global rangeMax
global vmaxunambg
global Nrange
global Ndoppler
global lambda
global XrngdopTgt
global rngresp
% Plot range-Doppler map
helperRangeVelocityMapPlot(XrngdopTgt,rngresp.SweepSlope,...
    fmcwwav1.SampleRate,tpri,rangeMax,vmaxunambg,Nrange,Ndoppler,lambda);
end
function functie5()
global tpri
global fmcwwav1
global rangeMax
global vmaxunambg
global Nrange
global Ndoppler
global lambda
global XrngdopTgt
global rngresp

global XrngangdopTgt
global vrxEleSpacing
global Nangle
global tarDopplerBin
% Plot range-angle map
helperRangeAngleMapPlot(XrngangdopTgt,rngresp.SweepSlope,...
    fmcwwav1.SampleRate,rangeMax,vrxEleSpacing,Nrange,Nangle,tarDopplerBin,lambda);
end
function functie6()
    % Calea către fișierul PDF
    filePath = 'ambeleradare.pdf';

    % Deschide fișierul PDF cu aplicația asociată (Adobe, browser etc.)
    open(filePath);
end
function functie7()
    % Calea către fișierul PDF
    filePath = 'drivingscenario.pdf';

    % Deschide fișierul PDF cu aplicația asociată (Adobe, browser etc.)
    open(filePath);
end
function functie8()
    % Calea către fișierul PDF
    filePath = 'angledoplermap.pdf';

    % Deschide fișierul PDF cu aplicația asociată (Adobe, browser etc.)
    open(filePath);
end