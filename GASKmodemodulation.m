clear all;
close all;
clc;

%% 0. Modulation Parmeter Initialization

% Gaussian parameter
samples_per_bit = 64;
n = [-samples_per_bit/2:samples_per_bit/2];
B = 0.5;
k1 = sqrt(2*pi/reallog(2));
h = k1*B*exp(-2*k1*k1*pi.*B.*B*(1/samples_per_bit.*n).^2);

% Modualtion Parameter
SampleRate = 44100;
t = 0:1/SampleRate:5;
fc = 18000;

%% 1. Modulation Process 

bit = [0,0,0,0,0,1,0,0,1,1,0,1,0,1,1,0];
% bit = [0 bit 0]; 
bit_raw =[];
percentil = 0;
for i = 1:length(bit)
    bit_raw = [ bit_raw zeros(1,samples_per_bit*percentil) ones(1,samples_per_bit*(1-percentil))*bit(i) zeros(1,samples_per_bit*percentil)];
end

%pre padding
bit_raw = [ ones(1,samples_per_bit/2)*bit(1)  bit_raw];

%post padding
bit_raw = [ bit_raw  ones(1,samples_per_bit/2)*bit(length(bit))];

% Gaussian Filtering
s = filter(h,1,bit_raw);
s = s(samples_per_bit/2+1 : length(s) - samples_per_bit/2);

figure;
    subplot(3,1,1);
    hold on;
    plot(s);
    xlim([0,1200]);
    hold off;
    subplot(3,1,2);
    hold on;
    plot(s,'r');
    xlim([0,1200]);
    hold off;

x = cos(2*pi*fc*t);
x = x(1:length(s)).*s;
max_v = max(x);

subplot(3,1,3);
z = x/max_v;
plot(z,'Color',[0.313725501298904 0.313725501298904 0.313725501298904]);
ylim([-1,1]);
xlim([0,1200]);

% Generate modulated signal 
waves=z;
z = z.*0.99;
filename = 'GASK_gen_signal.wav';
audiowrite(filename,z,SampleRate);




%% 2. Demodulation Process


% 2.1 Read the modulated signal
[z,Fs] = audioread('GASK_test_signal.wav');
figure;
plot(z);
% Focusing the region of interesting
z=z(59500:60524);
z= z';
% t = (0:length(z)-1)/SampleRate;
% de_x0 = z(1:length(z)).*cos(2*pi*(fc).*t(1:length(t)));

%% 2.2 Band Pass Filter
Wn = [15000 19000];
[b,a] = butter(5,Wn/(SampleRate/2),'bandpass');
de_x0 = filter(b,a,z);

%% 2.2 Low Pass Filter
[b,a] = butter(4,500/(SampleRate/2),'low');
de_x0 = filter(b,a,abs(de_x0));

%% 2.3 normalization
de_x0 = de_x0./max(de_x0);

%% 2.4 Decision for extracting the data (sigmoid func)
de_x0 = 1./(1+exp(-20*(de_x0-0.5)))-0.5;
de_x = de_x0;

figure;
subplot(4,1,1);
    plot(z,'Color',[0.313725501298904 0.313725501298904 0.313725501298904]);
subplot(4,1,2);
    hold on;
    plot(de_x(54:length(de_x))+0.5);
    xlim([0,1200]);
    ylim([0,1.01]);
    subsample = de_x(54:length(de_x))+0.5;
    stem(54:64:length(subsample),subsample(54:64:length(subsample)));
    hold off;

dem = (de_x);

subplot(4,1,3);
    hold on;
    plot(dem(54:length(dem)),'r');
    plot(1:length(de_x)-1:length(de_x),zeros(1,2),'k');
    xlabel('Sample');
    
    for i = 1 : length(dem)
        if(dem(i) < 0)  
            dem(i) = 0 ;
        else
            dem(i) = 1 ;
        end;
    end
    hold off;


%% bit seperator
bit_offset = 33;
cnt = 0;
j = 1;
for i = bit_offset:length(dem)
    cnt = cnt + dem(i);
    if mod(i-bit_offset+1,samples_per_bit) == 0
        cnt_v(j) = (cnt > 0);
        cnt = 0;
        j = j + 1;
    end
end

cnt_v(j) = (cnt>samples_per_bit/2);
cnt_v;

for i = 1:length(de_x)/samples_per_bit
    hold on;
    plot(i*samples_per_bit*ones(2)+bit_offset,-0.5:0.5,'k');    %each bits vertical line
    hold off;
    ylim([-1 1]);
end

subplot(4,1,4)
stem(subsample(54:64:length(subsample)));
round(subsample(54:64:length(subsample)),0)
xlim([0, 18.75])
xlabel('bit stream');
