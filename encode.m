
block_size = 8
[data_in, sample_rate] = load_input('snare-808');
[bass, noise_env, noise_color] = prepare(data_in, sample_rate, block_size);


% Load wav-file
function [data, sample_rate] = load_input(stub)
	filename = ['in/'  stub  '.wav'];
	[data, sample_rate] = audioread(filename);
end

% Create butterworth filter at given frequency
function [b, a] = butter_norm(sample_rate, order, cutoff_freq, type)
	nyquist = sample_rate / 2.0;
	[b, a] = butter(order, cutoff_freq / nyquist, type);
end

% Prepare data for compression
function [bass, noise_env, noise_col_b, noise_col_a] = prepare(data_in, sample_rate, block_size)
	nyquist = sample_rate/2.0;
	separation_freq = nyquist/block_size/2.0;

	% prepare separation filters
	[bass_b, bass_a] = butter_norm(sample_rate, 6, separation_freq, 'low');
	[noise_b, noise_a] = butter_norm(sample_rate, 6, separation_freq, 'high');

	% Run separation filters
	bass = filter(bass_b, bass_a, data_in);
	noise = filter(noise_b, noise_a, data_in);
	
	% Create noise volume envelope and fit it with piecewise linear approximation
	noise_env = abs(hilbert(noise));
	noise_env = filter(bass_b, bass_a, noise_env);
    
    % extract noise color
    num_freqs = 1024;
    noise_col_h = fft(noise, num_freqs);
    noise_col_h = noise_col_h(1:num_freqs/2+1);
    noise_col_w = (0:(num_freqs/2)) * (2*pi)/num_freqs;
    [noise_col_b, noise_col_a] = invfreqz(noise_col_h, noise_col_w, 3, 3);
    
    noise_col_scaler = max(abs([noise_col_a(2:end) noise_col_b]));
    noise_col_a = noise_col_a / noise_col_scaler
    noise_col_b = noise_col_b / noise_col_scaler

    % plot
    [noise_col_h2, noise_col_w2] = freqz(noise_col_b, noise_col_a, num_freqs);
    figure;
    subplot(2,1, 1);
    plot(noise_col_w*sample_rate/(2*pi), abs(noise_col_h)/num_freqs);
    subplot(2,1, 2);
    plot(noise_col_w2*22100.0/(2*pi), abs(noise_col_h2)/num_freqs)

    % Quantize this shit
    [noise_env_d, noise_env_d_scaler] = digitize_noise_env(noise_env);
    [bass_d, bass_d_scaler] = digitize_bass(bass);

	% plots
	figure;
	subplot(2,2, 1);
	plot(1:size(data_in), bass);
	subplot(2,2, 2);
	plot(1:size(data_in), noise_env);
	subplot(2,2, 3);
	plot(1:size(data_in), bass_d);
	subplot(2,2, 4);
	plot(1:size(data_in), noise_env_d);
        
end

% bass is in [-1:1]
function [res, scaler] = digitize_bass(val)
    bits = 4;

    % normalize
    scaler = max(abs(val));
    res = val / scaler;
    
    % digitize
    limit = 2 .^ bits;
    res = (res + 1.) * ((limit-1) / 2.);
    res = res + rand(size(res));
    res = round(res);
    res = min(limit-1, max(0.0, res));
end

% noise is in [0:1]
function [res, scaler] = digitize_noise_env(val)
    bits = 4;

    % normalize
    scaler = max(abs(val));
    res = val / scaler;
    
    % digitize
    limit = 2 .^ bits;
    res = res * (limit-1);
    res = res + rand(size(res)) - 0.5;
    res = round(res);
    res = min(limit-1, max(0.0, res));
end



