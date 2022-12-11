function [G] = MskBocSpectrum(l, freqs, ft)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
G = zeros(1, length(freqs));
if mod(l,2) == 0
    G = 8/(l^2*pi^2*ft)*(sin(pi*freqs/ft)./ (1-(2*freqs/(l*ft)).^2)).^2;
else
    G = 8/(l^2*pi^2*ft)*(cos(pi*freqs/ft)./ (1-(2*freqs/(l*ft)).^2)).^2;
end
end

