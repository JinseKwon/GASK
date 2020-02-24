### GASK (Gaussian Amplitude Shift Keying)
---
#### Introduction

GASK provides an optimized modulation scheme for Ultrasonic indoor localization on mobile devices.

We have written the example code in matlab R2018b and provide an audio file and an .m file.

------

#### Modualtion Process
![](https://github.com/JinseKwon/GASK/tree/master/image/GASK_modulation.png)

Input return to zero(RZ) bits are convolved with a Gaussian filter to create a smooth bit stream.
```
// Gaussian parameter
samples_per_bit = 64;
n = [-samples_per_bit/2:samples_per_bit/2];
B = 0.5;
k1 = sqrt(2*pi/reallog(2));
h = k1*B*exp(-2*k1*k1*pi.*B.*B*(1/samples_per_bit.*n).^2);

...

// Gaussian Filtering
s = filter(h,1,bit_raw);
```

------

#### Demodualtion Process
![](https://github.com/JinseKwon/GASK/tree/master/image/GASK_demodulation.png)

In the demodulation process, the carrier signal is processed by a band pass filter, and then a low pass filter is used to remove noise.
```
// 2.1 Read the modulated signal
[z,Fs] = audioread('GASK_test_signal.wav');

// 2.2 Band Pass Filter
Wn = [15000 19000];
[b,a] = butter(5,Wn/(SampleRate/2),'bandpass');
de_x0 = filter(b,a,z);

// 2.2 Low Pass Filter
[b,a] = butter(4,500/(SampleRate/2),'low');
de_x0 = filter(b,a,abs(de_x0));

// 2.3 normalization
de_x0 = de_x0./max(de_x0);
```
At the end of the GASK demodulation process, the received bits are extracted through the slicer for decion.

------

#### Demonstration
Additional supplementary videos are available for further understanding. <https://youtu.be/tn9yxw23ZNY>
