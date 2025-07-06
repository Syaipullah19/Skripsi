% Analisis Pengaruh Noise terhadap Akurasi Pengukuran Partial Discharge
% Bbtv busbar atas
% Berbasis Analisis Frekuensi dan High Pass Filter
% Membaca data waveform dari file .csv di folder repo

clc;

% --- BACA DATA WAVEFORM DARI CSV ---
data = readmatrix('Bbtv_busbar.csv');
if size(data,2) > 1
    waveforms = data(:,2); % Ambil kolom kedua sebagai amplitudo
else
    waveforms = data(:);   % Jika hanya 1 kolom, langsung pakai semua data
end
N = length(waveforms);

dt = 0.0023e-6;             % Interval waktu per data (detik)
fs = 1/dt;                  % Sampling frequency (Hz)
t = (0:N-1)*dt*1e6;         % Waktu dalam mikrodetik

SNR_values = [30, 20, 10, 5, 0]; % dB

fc = 3e6;                 % High pass filter cutoff 3 MHz
order = 4;
[b, a] = butter(order, fc/(fs/2), 'high');

results = struct;           % Untuk menyimpan hasil analisis

figure('Position',[100 100 1200 900]);
for i = 1:length(SNR_values)
    % Tambahkan noise Gaussian
    noisy_signal = awgn(waveforms, SNR_values(i), 'measured');
    % Filtering
    filtered_signal = filtfilt(b, a, noisy_signal);
    filtered_signal = filtered_signal(:); % Pastikan vektor 1D

    % FFT
    Y_before = abs(fft(noisy_signal));
    Y_after = abs(fft(filtered_signal));
    f = (0:N-1)*(fs/N)/1e6;

    % Plot FFT sebelum filtering
    subplot(length(SNR_values),2,2*i-1);
    plot(f, Y_before, 'b');
    title(['FFT Noisy Sinyal, SNR = ', num2str(SNR_values(i)), ' dB']);
    xlabel('Frekuensi (MHz)');
    ylabel('Magnitudo');
    xlim([0 fs/(2*1e6)]);
    grid on;

    % Plot FFT setelah filtering
    subplot(length(SNR_values),2,2*i);
    plot(f, Y_after, 'r');
    title(['FFT Setelah High Pass, SNR = ', num2str(SNR_values(i)), ' dB']);
    xlabel('Frekuensi (MHz)');
    ylabel('Magnitudo');
    xlim([0 fs/(2*1e6)]);
    grid on;

    % --- PD Point Detection ---
    [pks, locs] = findpeaks(filtered_signal, 'MinPeakHeight', mean(filtered_signal) + std(filtered_signal));
    pd_times = t(locs);   % Waktu terjadinya PD

    % --- Perhitungan dB dan ppc ---
    ppc = max(filtered_signal) - min(filtered_signal);
    noise_level = std(noisy_signal(1:round(N/5))); % Ambil baseline dari 1/5 awal data
    if noise_level == 0, noise_level = 1e-12; end % Hindari log(0)
    dB = 20*log10(max(abs(filtered_signal))/noise_level);

    % --- Penentuan jenis PD sederhana ---
    if ppc < 0.5 && dB == 0 && dB <= 9
        Tingkat_pd = "Tidak memerlukan perhatian";
    elseif ppc >= 0.5 && ppc <= 6 && dB >= 10 && dB <= 19
        Tingkat_pd = "Kemungkinan PD tingkat Rendah";
    elseif ppc >= 0.5 && ppc <= 6 && dB >= 20 && dB <= 29
        Tingkat_pd = "Kemungkinan PD tingkat menengah";
    elseif ppc >= 0.5 && ppc <= 6 && dB > 30
        Tingkat_pd = "Kemungkinan PD tingkat tinggi";    
    elseif ppc > 6 && ppc <= 30 && dB == 0 && dB <= 9
        Tingkat_pd = "Tidak memerlukan perhatian";
    elseif ppc > 6 && ppc <= 30 && dB >= 10 && dB <= 19
        Tingkat_pd = "Kemungkinan Discharge Permukaan_Cek Ultrasound";
    elseif ppc > 6 && ppc <= 30 && dB >= 20 && dB <= 29
        Tingkat_pd = "Kemungkinan Discharge Permukaan_Cek Ultrasound";
    elseif ppc > 6 && ppc <= 30 && dB > 30
        Tingkat_pd = "Kemungkinan logam mengambang/koneksi jelek/kendor";
    elseif ppc >= 30 && dB == 0 && dB <= 9
        Tingkat_pd = "Tidak memerlukan perhatian";
    elseif ppc >= 30 && dB >= 10 && dB <= 29
        Tingkat_pd = "Noise";
    elseif ppc > 30 && dB > 30
        Tingkat_pd = "Kemungkinan logam mengambang/koneksi jelek/kendor";
    else
        Tingkat_pd = "Noise";
    end

    % Simpan hasil
    results(i).SNR = SNR_values(i);
    results(i).pd_times = pd_times;
    results(i).pk_values = pks;
    results(i).dB = dB;
    results(i).ppc = ppc;
    results(i).jenis_pd = Tingkat_pd;

    % --- Plot waktu & PD point ---
    figure(100+i);
    plot(t, filtered_signal, 'b', 'LineWidth', 1.2); hold on;
    plot(pd_times, pks, 'ro', 'MarkerFaceColor','r');
    xlabel('Waktu (\mus)');
    ylabel('Amplitudo (setelah filter)');
    title(sprintf('Sinyal PD setelah Filtering (SNR = %d dB) \nTitik PD merah, Jenis: %s, dB: %.2f, ppc: %.2f', ...
        SNR_values(i), Tingkat_pd, dB, ppc));
    legend('Filtered Signal', 'PD Points');
    grid on;
end

sgtitle('Analisis Pengaruh Noise pada Pengukuran Partial Discharge (FFT & High Pass Filter)');

% Tabel Hasil
fprintf('\nRekapitulasi Analisis Partial Discharge:\n');
fprintf('SNR(dB)\tdB nilai\tppc\tJenis PD\t# Titik PD\n');
for i = 1:length(SNR_values)
    fprintf('%d\t%.2f\t%.2f\t%s\t\t%d\n', ...
        results(i).SNR, results(i).dB, results(i).ppc, string(results(i).jenis_pd), length(results(i).pd_times));
end

% Visualisasi sinyal pada domain waktu (opsional)
figure(200);
plot(t, waveforms, 'k', 'DisplayName', 'Original');
hold on
plot(t, awgn(waveforms, SNR_values(end), 'measured'), 'r', 'DisplayName', ['Noisy (SNR=', num2str(SNR_values(end)), ' dB)']);
xlabel('Waktu (\mus)');
ylabel('Amplitudo');
title('Sinyal Partial Discharge (Domain Waktu)');
legend show
grid on