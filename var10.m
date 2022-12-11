clc
clear
close all
l = 2;
N = 31;
ft = 1;
gamma = 0.9;
q = 10;
sigma_tau_ns = 20;

%% Пункт 1
close all
a = ones(1, 5);
d = ones(1, N);
for i=1:N
    if a(4) == 0
       d(i) = -1;
    end
    b = xor(a(4), a(1));
    a = circshift(a, [1, 1]);
    a(1) = b;
end

samples_in_symbol = 1000;
[I_msk, Q_msk] = MskModulator(d, samples_in_symbol, 2);
len = min(length(I_msk), length(Q_msk));

f = 0.003;
t = 0:1:len-1;
msk = I_msk(1,1:len) .* cos(2*f*pi*t) - Q_msk(1,1:len) .* sin(2*f*pi*t);
% figure
% plot(t, d_signal)
% xlim([min(t) max(t)])
% grid minor
% xlabel('Время, мкс')
% ylabel('Значение сигнала, B')

figure
plot(t, msk, 'LineWidth', 1)
xlabel('Время, нс')
ylabel('Значение сигнала, B')
grid minor
xlim([min(t) max(t)])

% figure
% subplot(2,1,1)
% plot(t, I_msk(1,1:end-samples_in_symbol))
% title('Cоставляющая I')
% xlabel('Время, мкс')
% ylabel('Значение сигнала, B')
% grid minor
% xlim([min(t) max(t)])
% subplot(2,1,2)
% plot(t, Q_msk(1,1:end-samples_in_symbol))
% title('Cоставляющая Q')
% xlabel('Время, мкс')
% ylabel('Значение сигнала, B')
% grid minor
% xlim([min(t) max(t)])

%% Пункт 2
close all
df = 0.09;
f = 0:df:10*ft;
P = zeros(1, length(f));
G = MskBocSpectrum(l, f, 1);
P(1) = G(1)*df;
for i = 2:length(f)
    P(i) = P(i-1)+G(i)*df;
end
P = 2*P;
g = gamma*ones(1, length(f));
figure
plot(f, P, f, g)
xlabel('Нормированная частота')
ylabel('Мощность сигнала, Вт')
grid minor
ylim([min(P) 1.1*max(P)])

%% Пункт 3
close all
% Из предыдущего пункта получили что Fв = 1.35
F_high = 1.35;
df = 0.0011;
f = df:df:F_high; 
f = [f F_high];
G = MskBocSpectrum(l, f, ft);
Fe = sqrt((trapz(f,f.*f.*G)/trapz(f,G))); 

F_big_high = 100; % Значение частоты при котором гамма квадрат равно 0.9999
f0 = df:df:F_big_high;
f0 = [f0 F_big_high];
G = MskBocSpectrum(1, f0, ft);
Fee = sqrt((trapz(f0,f0.*f0.*G)/trapz(f0,G))); 

loseEnergy = Fee/Fe;
loseEnergydB = 20*log10(loseEnergy);
%% Пункт 4
close all
% Из предыдущего раздела получили что Fe = 0.8654
% sigma_tau = 1/(2*pi*F_effective*db2mag(q));
Fe = 0.8654;
F_e_needed_MHz = 1/(2*pi*sigma_tau_ns/1000*gamma*db2mag(q));
f_t_new = F_e_needed_MHz/Fe;


%% Пункт 5
% В предыдущих разделах получили что f_t_new = 3.2310, F_high
close all
F_high = 1.35;
f_t_new = 3.2310;
T_t = 1/f_t_new;
df = 0.0011;
f = df:df:F_high*f_t_new;
Fd = f_t_new*100;
tau = -1.5*T_t:0.01*T_t:1.5*T_t;
G = MskBocSpectrum(l, f, f_t_new);
R = zeros(1, length(tau));
for i = 1:length(tau)
    R(i) = 2*trapz(f, G.*cos(2*pi*f*tau(i))); 
end

figure
plot(tau,R)
grid minor 
xlabel('Временной сдвиг, мкс')
ylabel('Значение АКФ')
title('АКФ')
xlim([min(tau) max(tau)])

figure
plot(f, 10*log10(G))
grid minor
xlabel('Частота, МГц')
ylabel('Значение спектра, Вт')
title('Энергетический спектр')
ylim([-60 0])
%% Пункт 6
% Из предыдущих разделов получили что Fe_needed = 2.7961
close all
Fe_needed = 2.7961;
q_dB = 0:20;
q = db2mag(q_dB);
sigma_tau = 1./(2*pi*Fe_needed*gamma*q);
sigma_phi = 1./(gamma*q);
figure
plot(q_dB, sigma_tau)
grid minor
xlabel('ОСШ, дБ')
ylabel('Значение ошибки, мкс')
title('Зависимость ошибки кодовой синхронизации от ОСШ')
figure
plot(q_dB, sigma_phi)
grid minor
title('Зависимость ошибки фазовой синхронизации от ОСШ')
xlabel('ОСШ, дБ')
ylabel('Значение ошибки, рад.')
