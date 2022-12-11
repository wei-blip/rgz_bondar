clc
clear
close all
l = 2;
N = 31;
ft = 1;
gamma = 0.99;
q = 8;
sigma_tau_ns = 22;

%% ����� 1
close all
a = ones(1, 5);
d = ones(1, N);
for i=1:N
    if a(2) == 0
       d(i) = -1;
    end
    b = xor(a(3), a(1));
    a = circshift(a, [1, 1]);
    a(1) = b;
end

samples_in_symbol = 1000;
[I_msk, Q_msk] = MskModulator(d, samples_in_symbol, 3);
len = min(length(I_msk), length(Q_msk));

d_signal = zeros(1,length(d)*samples_in_symbol);
for i=1:length(d)
    d_signal(1,1+samples_in_symbol*(i-1):samples_in_symbol*i) = d(i);
end

t = 0:1/samples_in_symbol:N-1/samples_in_symbol;
figure
plot(t, d_signal)
xlim([min(t) max(t)])
grid minor
xlabel('�����, ���')
ylabel('�������� �������, B')

f = 0.003;
t = 0:1:len-1;
msk = I_msk(1,1:len) .* cos(2*f*pi*t) - Q_msk(1,1:len) .* sin(2*f*pi*t);
% figure
% plot(t, d_signal)
% xlim([min(t) max(t)])
% grid minor
% xlabel('�����, ���')
% ylabel('�������� �������, B')

figure
plot(t, msk, 'LineWidth', 1)
xlabel('�����, ��')
ylabel('�������� �������, B')
grid minor
xlim([min(t) max(t)])

% figure
% subplot(2,1,1)
% plot(t, I_msk(1,1:end-samples_in_symbol))
% title('C����������� I')
% xlabel('�����, ���')
% ylabel('�������� �������, B')
% grid minor
% xlim([min(t) max(t)])
% subplot(2,1,2)
% plot(t, Q_msk(1,1:end-samples_in_symbol))
% title('C����������� Q')
% xlabel('�����, ���')
% ylabel('�������� �������, B')
% grid minor
% xlim([min(t) max(t)])

%% ����� 2
close all
df = 0.0011;
f = 0.0011:df:15*ft;
P = zeros(1, length(f));
G = MskBocSpectrum(3, f, 1);
P(1) = G(1)*df;
for i = 2:length(f)
    P(i) = P(i-1)+G(i)*df;
end
P = 2*P;
g = gamma*ones(1, length(f));
figure
plot(f, P, f, g)
xlabel('������������� �������')
ylabel('�������� �������, ��')
grid minor
ylim([min(P) 1.05*max(P)])
xlim([0 15])

%% ����� 3
close all
% �� ����������� ������ �������� ��� F� = 2.936
F_high = 2.936;
df = 0.0011;
f = df:df:F_high; 
f = [f F_high];
G = MskBocSpectrum(3, f, ft);
Fe = sqrt((trapz(f,f.*f.*G)/trapz(f,G)));

f_ = df:df:14.53;
G_ = MskBocSpectrum(3, f_, ft);
Fe_ = sqrt((trapz(f_,f_.*f_.*G_)/trapz(f_,G_)));
loseEnergy = Fe_/Fe;
loseEnergydB = 20*log10(loseEnergy);
%% ����� 4
close all
% �� ����������� ������� �������� ��� Fe = 1.4393
% sigma_tau = 1/(2*pi*F_effective*db2mag(q));
Fe = 1.4393;
F_e_needed_MHz = 1/(2*pi*sigma_tau_ns/1000*gamma*db2mag(q));
f_t_new = F_e_needed_MHz/Fe;


%% ����� 5
% � ���������� �������� �������� ��� f_t_new = 3.0922, F_high = 1.8
close all
F_high = 2.936;
f_t_new = 3.0212;
T_t = 1/f_t_new;
df = 0.0011;
f = df:df:F_high*f_t_new;
Fd = f_t_new*100;
tau = -1.5*T_t:0.01*T_t:1.5*T_t;
G = MskBocSpectrum(3, f, f_t_new);
R = zeros(1, length(tau));
for i = 1:length(tau)
    R(i) = 2*trapz(f, G.*cos(2*pi*f*tau(i))); 
end

figure
plot(tau,R)
grid minor 
xlabel('��������� �����, ���')
ylabel('�������� ���')
title('���')
xlim([min(tau) max(tau)])

figure
plot(f, 10*log10(G))
grid minor
xlabel('�������, ���')
ylabel('�������� �������, ��')
title('�������������� ������')
ylim([-60 0])
%% ����� 6
% �� ���������� �������� �������� ��� Fe_needed = 2.9091
close all
Fe_needed = 2.9091;
q_dB = 0:20;
q = db2mag(q_dB);
sigma_tau = 1./(2*pi*Fe_needed*gamma*q);
sigma_phi = 1./(gamma*q);
figure
plot(q_dB, sigma_tau)
grid minor
xlabel('���, ��')
ylabel('�������� ������, ���')
title('����������� ������ ������� ������������� �� ���')
figure
plot(q_dB, sigma_phi)
grid minor
title('����������� ������ ������� ������������� �� ���')
xlabel('���, ��')
ylabel('�������� ������, ���.')
