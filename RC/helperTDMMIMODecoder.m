function [Xdec,Nsweep] = helperTDMMIMODecoder(X,w)
[Nrange,Nr,Nsweep] = size(X);
Nt = size(w,1);

% Reduce number of sweeps after TDM decoding
Nsweep = floor(Nsweep/Nt);

% TDM decoding
X = X(:,:,1:Nsweep*Nt);
Xdec = reshape(X,Nrange,Nt*Nr,Nsweep);
end
