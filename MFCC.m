[wav, Fs] = audioread('a.wav');
t = 0 : 1/Fs : length(wav)/Fs;
figure(1)
plot(t(1:length(t)-1)*1000, wav); xlabel('Time [ms]'); ylabel('Amplitude');

center = fix(length(wav) / 2);
cuttime = 0.04;
wavdata = wav(center-fix(cuttime/2*Fs) : center+fix(cuttime/2*Fs));
time = t(center-fix(cuttime/2*Fs) : center+fix(cuttime/2*Fs));
han_window = 0.5 - 0.5 * cos(2 * pi * [0 : 1/length(wavdata) : 1]);
wavdata = han_window(1:length(wavdata))' .* wavdata;

fftsize = 2048;
dft = fft(wavdata, fftsize);
Pdft = (real(dft) .^ 2) + (imag(dft) .^ 2);
Adft = sqrt(Pdft);
fscale = linspace(0, Fs, fftsize);

Angdft = atan2(imag(dft), real(dft));

Adft_log = log10(abs(dft));
Pdft_log = log10(abs(dft) .^ 2);

figure(2)
cps = real(ifft(Adft_log))
quefrency = linspace(0, cuttime, cuttime * Fs);
plot(quefrency(1: fftsize / 2) * 1000, cps(1: fftsize / 2));
xlabel('Quefrency [ms]'); ylabel('Logarithmic Amplitude [dB]');
saveas(gcf, 'cepstrum.png');

cps_lif = cps;
cps_lif(120: length(cps) - 118) = 0;

figure(3)
dftSpc = fft(cps, fftsize);
AdftSpc = abs(10 .^ dftSpc);
PdftSpc = AdftSpc .^ 2;
subplot(2, 1, 1); plot(fscale(1: fftsize / 2), PdftSpc(1: fftsize / 2));
xlabel('Frequency [Hz]'); xlim([0, 5000])

dftSpc = fft(cps_lif, fftsize);
AdftSpc = abs(10 .^ dftSpc);
PdftSpc = AdftSpc .^ 2;
subplot(2, 1, 2); plot(fscale(1: fftsize / 2), PdftSpc(1: fftsize / 2));
xlabel('Frequency [Hz]'); xlim([0, 5000]);
saveas(gcf, 'liftering.png');
