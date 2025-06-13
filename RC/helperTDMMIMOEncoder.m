function w = helperTDMMIMOEncoder(Nt,Nsweep)
% Time division multiplexing (TDM) code applied on each antenna element
% at each sweep
w = zeros(Nt,Nsweep);
for l = 1:Nsweep
    wl = int8((1:Nt)' == mod(l-1,Nt)+1);
    w(:,l) = wl;
end
end