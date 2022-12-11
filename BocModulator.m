function [boc] = BocModulator(sequence, samples_in_symbol, order)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
hperiod = samples_in_symbol/order;
pulse = [1 ones(1,hperiod) -ones(1,hperiod)];
L = length(sequence);
pulse_signal = zeros(1, L);
for i=1:L/(2*hperiod)
    begin = 1 + 2*hperiod*(i-1);
    endl = 2*hperiod*i + 1;
    pulse_signal(1,begin:endl) = pulse;
end
boc = pulse_signal .* sequence;
end

