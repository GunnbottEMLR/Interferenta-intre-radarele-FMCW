function helperRangeVelocityMapPlot(Xrngdop,sweepSlope,sampleRate,tpri,rangeMax,vmaxunambg,Nrange,Ndoppler,lambda)
% Noncoherent integration over all virtual array elements
XrngdopPowIncoInt = 2*pow2db(squeeze(pulsint(Xrngdop)));

% Beat frequency grid
fbeat = sampleRate*(-Nrange/2:Nrange/2-1)'/Nrange; 

% Convert beat frequency to range
rangeGrid = beat2range(fbeat,sweepSlope);

% Doppler frequency grid
fDoppler = 1/tpri*(-Ndoppler/2:Ndoppler/2-1)'/Ndoppler; 

% Convert Doppler frequency to velocity for two-way propagation
velocityGrid = dop2speed(fDoppler,lambda)/2;

% Obtain range-Doppler meshgrid
[Range, Velocity] = meshgrid(rangeGrid,velocityGrid);

% Plot range-Doppler map
figure
surf(Range,Velocity,XrngdopPowIncoInt.') 
view(0,90);
xlim([0, rangeMax]); ylim([-vmaxunambg,vmaxunambg])
xlabel('Range (m)'); ylabel('Velocity (m/s)'); zlabel('Power (dB)')
title('Range-Doppler Map','FontSize',12); colormap(jet); 
cb = colorbar; cb.Label.String = 'Power (dB)'; cb.FontSize = 10;
axis square; shading interp 
end