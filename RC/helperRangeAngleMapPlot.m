function [angGrid,rangeGrid] = helperRangeAngleMapPlot(Xrngangdop,sweepSlope,sampleRate,rangeMax,vrxEleSpacing,Nrange,Nangle,tarDopplerBin,lambda)
% Power of range-angle-Doppler response
XrngangdopPow = abs(Xrngangdop).^2;

% Beat frequency grid
fbeat = sampleRate*(-Nrange/2:Nrange/2-1)'/Nrange; 

% Convert beat frequency to range
rangeGrid = beat2range(fbeat,sweepSlope);

% Angle grid
angGrid = asind(lambda/vrxEleSpacing*(-(Nangle/2):(Nangle/2)-1)'/Nangle);

% Obtain range-angle meshgrid
[Range, Angle] = meshgrid(rangeGrid,angGrid);

% Plot range-angle map
figure
surf(Range, Angle, pow2db(squeeze(XrngangdopPow(:,:,tarDopplerBin)).')) 
view(0,90);
xlim([0, rangeMax]); ylim([-60, 60]);
xlabel('Range (m)'); ylabel('Angle (degrees)'); zlabel('Power (dB)')
title('Range-Angle Map','FontSize',12); colormap(jet); 
cb = colorbar; cb.Label.String = 'Power (dB)'; cb.FontSize = 10; 
axis square; shading interp 
end