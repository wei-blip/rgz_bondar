close all
clc 
clear
gamma = 0.99;   %Внутриполосная мощность 
q = 6;  %Отношение сигнал/шум, дБ
sigma_tau = 0.025;  %СКО ошибки 
n=1;
N = 31; %Длина кодовой ПСП 
ft = 1*10^6;   %Тактовая частота 
Fe = 1/(2*pi*sigma_tau*sqrt(gamma)*10^(q/20));  %Эффек-тивная ширина спектра, 
%где exp(q/20) q в разах
l = Fe/ft;    %Значение нормированной эффективной шири-ны спектра
f0 = 20*l*ft; %Несущая частота
fd = 10*f0;   %Частота  дискретизации 
T = 1/ft; %Длительность элемента ШПС
Tc = N*T; %Определяющая длительность ШПС
A = 1;
t1 = 0:1/fd:T;
t = 0:1/fd:N*T;
f11 = 1/T;
phi = 0;
%% Формирование ПСП длины N=31
a = ones(1, 5);
N = 31;
f_prs = ones(1, N);
for i=1:N
if a(5) == 0
   f_prs(i) = -1; 
end
b0 = xor(a(5), a(4));
b1 = xor(a(3), b0);
b = xor(a(1),b1);
a = circshift(a, [1, 1]);
a(1) = b;
end
 
F = 100;
w = 2*pi;
t = 0:1/F:N;                         
f_bpsk = ones(1, length(t));         
                                      
for i = 1:N
    for j = F*(i-1)+1:F*i
        if f_prs(i) == -1
            f_bpsk(j) = -1;
        end 
    end
end   
 
S_cos_t = sin(w*t);
% S_cos_t = 1;
bpsk = S_cos_t.*f_bpsk;
figure
plot(t, bpsk);
xlabel('Время, мкс')
ylabel('Значение сигнала, В')
grid on
figure
plot((1:length(t))/100, f_bpsk, 'Linewidth', 2);
xlabel('Время, мкс')
ylabel('Значение сигнала, В')
grid on
%% 2.Расчёт внутриполосной мощности сигнала 
Fv = zeros(1,101);
Fv(1) = 1;
P_F = zeros(1,101);
f = 0.001:0.001:10;
ft1 = 1;
G_f = 1/ft1*(sin(pi*f/ft1)./(pi*f/ft1)).^2;    %Энерге-тический спектр сигнала
P_F = cumtrapz(f,G_f);
P_F_full = trapz(f,G_f);
P_F_otn=P_F./P_F_full;
figure
plot(f,G_f)
grid on;
figure
plot(f,P_F_otn)
Fv_vector=find(P_F_otn>0.99);
Fv2=f(Fv_vector(1));
grid on;
xlabel 'f/f_т'
ylabel 'P(F_в)/P_c' 
%% 3.Определение эффективной ширины спектра 
f = 0.001:0.001:5.345;
G_f = 1/ft1*(sin(pi*f/ft1)./(pi*f/ft1)).^2;
l = (trapz(f,f.*f.*G_f)/trapz(f,G_f))^(1/2); 
%% 4.Определение эквивалентных энергетических потерь из-за ограничения спектра 
% по предыдущим расчетам определена Fв(при гамма в квадрате, равной 0.99)=5.31
f = 0.001:0.001:9.729;
G_f=1/ft1*(sin(pi*f/ft1)./(pi*f/ft1)).^2;
Fe_0_99=(trapz(f,f.*f.*G_f)/trapz(f,G_f))^(1/2);
Lose_energy = 20*log10(Fe_0_99/l); 
%% 5.Выбор длительности элемента ШПС 
Fe_trebuemaya = 1/(2*pi*sigma_tau*sqrt(gamma)*10^(q/20));
ft=Fe_trebuemaya/l;
T = 1/(ft); 
%% 6.Построение графика энергетического спектра
f = 0.1:0.1:10*ft;
G_f = 1/ft*(sin(pi*f/ft)./(pi*f/ft)).^2;
Cdb=10*log(G_f);
figure;
plot(f, Cdb);
xlabel 'f/f_т'
ylabel 'G(f)'
grid on; 
%% 7. Расчёт АКФ сигнала  
v=9.41; %значение нормированной частоты п.2
Fv = v*ft;
f = 0.01:0.01:Fv;
G_f = 1/ft*(sin(pi*f/ft)./(pi*f/ft)).^2;
tau = -1.5*T:0.01*T:1.5*T;
R = zeros(1,length(tau));
for i=1:length(tau)
    R(i)=2*trapz(f,G_f.*cos(2*pi*f*tau(i)));
end
figure;
plot(tau, R)
ylabel 'R(\tau)'
xlabel '\tau, мкс'
grid on; 
xlim([min(tau) max(tau)])
%% 8. Построение графика зависимости СКО ошибки от от-ношения сигнал/шум
q_7 = 0:20;
qr=10.^(q_7./20);
sigma_tau = 1./(2*pi*Fe_trebuemaya*sqrt(gamma)*qr);
figure;
plot(q_7, sigma_tau);
ylabel '\sigma_{\tau}, мкс'
xlabel 'q, дБ'
grid on; 
%% 9.Построение графика зависимости СКО ошибки от отно-шения сигнал/шум
sigma_phi = 1./(sqrt(gamma)*qr);
figure;
plot(q_7, sigma_phi);
ylabel '\sigma_{\phi}'
xlabel 'q, дБ'
grid on; 
%% 10. Расчёт вероятности ошибки на бит
q15=0:15;
q_dB = sqrt(gamma)*10.^(q15./20);
phi1 = zeros(1, length(q_dB));
tau1 = zeros(1, length(q_dB));
R_10 = zeros(1, length(q_dB));
P_b = zeros(length(q_dB), 4); 
for j = 1:4
    phi = 1./(sqrt(gamma)*q_dB);
    tau = 1./(2*pi*Fe_trebuemaya*sqrt(gamma)*q_dB);
    if j == 1
        phi = phi1;
        tau = tau1;
    end
    if j == 2
        phi = phi1;
    end
    if j == 3
        tau = tau1;
    end
    for i = 1:length(q_dB)
    R_10(i) = 2*trapz(f,G_f.*cos(2*pi*f*tau(i)));
    x = -10:0.1:q_dB(i)*R_10(i)*cos(phi(i));
    P_b(i,j) = 1-1/sqrt(2*pi)*trapz(x,exp(-x.*x/2));
    end 
end
figure;
semilogy((0:15), P_b(:,1),'--');
hold on 
semilogy((0:15), P_b(:,2),'-.');
semilogy((0:15), P_b(:,3),':','Linewidth',2);
semilogy((0:15), P_b(:,4));
legend('\sigma_{\tau} = \sigma_{\phi} = 0','\sigma_{\phi} = 0','\sigma_{\tau} = 0','\sigma_{\tau} и \sigma_{\phi} не 0')
ylabel 'P_b'
xlabel 'q_b, дБ'
grid on;
