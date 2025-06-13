function Xrngangdop = helperAngleResp(Xrngdop,Nrange,Nangle,Ndoppler)
Xrngangdop = complex(zeros(Nrange,Nangle,Ndoppler));
for n = 1: Nrange
    for l = 1: Ndoppler
        % Data in different virtual array elements
        XrngdopArray = squeeze(Xrngdop(n,:,l)');

        % Add Hann window on data
        XrngdopHannArray = hanning(length(XrngdopArray)).*XrngdopArray;

        % Angle FFT
        Xrngangdop(n,:,l) = fftshift(fft(XrngdopHannArray,Nangle));
    end
end
end