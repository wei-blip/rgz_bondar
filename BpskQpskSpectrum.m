function [G] = BpskQpskSpectrum(freqs, ft)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
G = (1/ft)*(sin(pi*freqs/ft)./ (pi*freqs/ft)).^2;
end

