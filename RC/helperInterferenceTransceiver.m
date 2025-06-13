function iradar = helperInterferenceTransceiver(waveform, txPkPower, antGain, fc, itxArray, vrxArray)
% Refolosește radarTransceiver pentru a transmite semnalul de interferență la
% radarul de interferență și pentru a recepționa semnalul de interferență la radarul victimă
iradar = radarTransceiver("Waveform", waveform, 'ElectronicScanMode', 'Custom');

% Specifică emițătorul radarului de interferență
iradar.Transmitter = phased.Transmitter('PeakPower', txPkPower, 'Gain', antGain);

% Specifică array-ul de transmisie ca fiind array-ul de transmisie al radarului de interferență
iradar.TransmitAntenna = phased.Radiator('Sensor', itxArray, 'OperatingFrequency', fc, 'WeightsInputPort', true);

% Specifică array-ul de recepție ca fiind array-ul de recepție al radarului victimă
iradar.ReceiveAntenna = phased.Collector('Sensor', vrxArray, 'OperatingFrequency', fc);

% Specifică receptorul radarului victimă fără zgomot
iradar.Receiver = phased.ReceiverPreamp('Gain', antGain, 'NoiseMethod', 'Noise power', 'NoisePower', 0);
end
%Verificare pls
