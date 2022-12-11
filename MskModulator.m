function [I_msk, Q_msk] = MskModulator(sequence, samples_in_symbol, l)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
L = length(sequence);
Nbe = floor(L/2);
if mod(L,2) == 1
    Nbo = Nbe+1;
else
    Nbo = Nbe;
end

o = zeros(1, Nbo);
e = zeros(1, Nbe);
o_ind = 1;
e_ind = 1;
for i=1:L
    if mod(i,2) == 1
        o(o_ind) = sequence(i);
        o_ind = o_ind + 1;
    else
        e(e_ind) = sequence(i);
        e_ind = e_ind + 1;
    end
end

one_chip = ones(1, l*samples_in_symbol);
O_sig = zeros(1, Nbo*l*samples_in_symbol);
E_sig = zeros(1, Nbe*l*samples_in_symbol);
for i=1:Nbo
    begin = 1+(i-1)*l*samples_in_symbol;
    endl = begin+l*samples_in_symbol - 1;
    O_sig(1,begin:endl) = sign(o(i)) * one_chip;
end

for i=1:Nbe
    begin = 1+(i-1)*l*samples_in_symbol;
    endl = begin+l*samples_in_symbol - 1;
    E_sig(1,begin:endl) = sign(e(i)) * one_chip;
end

t = 0:l*samples_in_symbol*L;
s = sin(2*pi*1/(2*samples_in_symbol)*t);
figure
plot(s)

figure
subplot(2,1,1)
plot(O_sig)
subplot(2,1,2)
plot(E_sig)

OI_sig = O_sig .* s(1,1:length(O_sig));
EQ_sig = E_sig .* s(1,1:length(E_sig));
figure
subplot(2,1,1)
plot(OI_sig)
subplot(2,1,2)
plot(EQ_sig)

I_msk = OI_sig;
Q_msk = [zeros(1,samples_in_symbol/2) EQ_sig];
figure
subplot(2,1,1)
plot(I_msk)
title('Составляющая I')
grid minor
xlabel('Время, нс')
ylabel('Значение сигнала, В')
xlim([1 length(I_msk)])
subplot(2,1,2)
plot(Q_msk)
grid minor
title('Составляющая Q')
xlabel('Время, нс')
ylabel('Значение сигнала, В')
xlim([1 length(Q_msk)])
end

