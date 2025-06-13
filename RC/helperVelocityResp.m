function Xrngdop = helperVelocityResp(Xdec,Nrange,Nv,Ndoppler)
Xrngdop = complex(zeros(Nrange,Nv,Ndoppler));
for n = 1: Nrange
    for m = 1:Nv
        % Data in different pulses
        XdecPulse = squeeze(Xdec(n,m,:));

        % Add Hann window on data
        XdecHannPulse = hanning(length(XdecPulse)).*XdecPulse;

        % Doppler FFT
        Xrngdop(n,m,:) = fftshift(fft(XdecHannPulse,Ndoppler));
    end
end
end