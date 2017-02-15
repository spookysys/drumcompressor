block_size = 8;
[data_in, sample_rate] = load_input('kick-808');
[bass, noise_env, noise_color] = prepare(data_in, sample_rate, block_size);


% Load wav-file
function [data, sample_rate] = load_input(stub)
	filename = ['in/'  stub  '.wav'];
	[data, sample_rate] = audioread(filename);
    scaler = max(abs(data)); % normalize
    data = data / scaler; % normalize
end

% Create butterworth filter at given frequency
function [b, a] = butter_norm(sample_rate, order, cutoff_freq, type)
	nyquist = sample_rate / 2.0;
	[b, a] = butter(order, cutoff_freq / nyquist, type);
end

% Prepare data for compression
function [bass, bass_env_fit, noise_env_fit, noise_col_b, noise_col_a] = prepare(data_in, sample_rate, block_size)
	nyquist = sample_rate/2.0;
	separation_freq = nyquist/block_size/2.0;
    data_range = (0:size(data_in)-1)';

	% prepare separation filters
	[bass_b, bass_a] = butter_norm(sample_rate, 6, separation_freq, 'low');
	[noise_b, noise_a] = butter_norm(sample_rate, 6, separation_freq, 'high');

	% Run separation filters
	bass = filter(bass_b, bass_a, data_in);
	noise = filter(noise_b, noise_a, data_in);
	
	% Create bass volume envelope
	bass_env = abs(hilbert(bass));
	% bass_env = filter(bass_b, bass_a, bass_env);
    
	% Create noise volume envelope
	noise_env = abs(hilbert(noise));
	noise_env = filter(bass_b, bass_a, noise_env);
    
    % extract noise color
    num_freqs = 1024;
    noise_col_h = fft(noise, num_freqs);
    noise_col_h = noise_col_h(1:num_freqs/2+1);
    noise_col_w = (0:(num_freqs/2)) * (2*pi)/num_freqs;
    [noise_col_b, noise_col_a] = invfreqz(noise_col_h, noise_col_w, 2, 2);
    
    % normalize noise color 
    noise_col_scaler = max(abs([noise_col_a(2:end) noise_col_b]));
    noise_col_a = noise_col_a / noise_col_scaler;
    noise_col_b = noise_col_b / noise_col_scaler;
    [noise_col_h_res, noise_col_w_res] = freqz(noise_col_b, noise_col_a, num_freqs);

    % create polynomials approximating the falloff of bass and noise
    bass_env_fit = fit(data_range, bass_env, 'exp2');
    noise_env_fit = fit(data_range, noise_env, 'exp2');
    bass_env_fit_vals = arrayfun(@(x) feval(bass_env_fit, x), data_range);
    noise_env_fit_vals = arrayfun(@(x) feval(noise_env_fit, x), data_range);
        
    % remove falloff from bass and noise
    bass_norm = bass ./ bass_env_fit_vals;
    noise_norm = noise ./ noise_env_fit_vals;
    
    bass_reduced = bass(1 : block_size : end);
    bass_norm_reduced = bass_norm(1 : block_size : end);
    bass_env_fit_vals_reduced = bass_env_fit_vals(1 : block_size : end);
    disp(['Finding palette for bass'])
    [bass_palette, bass_err] = find_palette(bass_reduced, ones(length(bass_reduced)));
    disp(['Finding palette for bass_norm'])
    [bass_norm_palette, bass_norm_err] = find_palette(bass_norm_reduced, bass_env_fit_vals_reduced);
    
    disp(['Error for bass:', num2str(bass_err), ' palette: ', num2str(bass_palette)]);
    disp(['Error for bass_norm:', num2str(bass_norm_err), ' palette:', num2str(bass_norm_palette)]);
    
    % Quantize this shit
    % [noise_env_d, noise_env_d_scaler] = digitize_noise_env(noise_env);
    % [bass_d, bass_d_scaler] = digitize_bass(bass);

    figure;
	subplot(4,1, 1);
	plot(1:size(data_in), bass);
    subplot(4,1, 2);
	plot(1:size(data_in), bass_norm);
    ylim([-4, 4])
    subplot(4,1, 3);
	plot(1:size(data_in), bass_env);
    subplot(4,1, 4);
	plot(1:size(data_in), bass_env_fit_vals);

    figure;
    subplot(4,1, 1);
	plot(1:size(data_in), noise);
    subplot(4,1, 2);
	plot(1:size(data_in), noise_norm);
    ylim([-4, 4])
    subplot(4,1, 3);
	plot(1:size(data_in), noise_env);
    subplot(4,1, 4);
	plot(1:size(data_in), noise_env_fit_vals);
    
    figure;
    plot(noise_col_w*sample_rate/(2*pi), abs(noise_col_h)/num_freqs/50); hold on;
    plot(noise_col_w_res*sample_rate/(2*pi), abs(noise_col_h_res)/num_freqs); hold off;
	    
end



% bass is in [-1:1]
function [res, scaler] = digitize_bass(val)
    bits = 4;

    % normalize
    scaler = max(abs(val));
    res = val / scaler;
    
    % digitize
    limit = 2 ^ bits;
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
    limit = 2 ^ bits;
    res = res * (limit-1);
    res = res + rand(size(res)) - 0.5;
    res = round(res);
    res = min(limit-1, max(0.0, res));
end


% Run a simplified compression with given palette and return total error
function total_error = get_compression_error(palette, data_in, weights)
    total_error = 0;
    recon_1 = 0;
    recon_2 = 0;
    for i = 1:length(data_in)
        val = data_in(i);
        recon_slope = recon_1 - recon_2;
        prediction = recon_1 + recon_slope;
        if (prediction >= 0)
            palette_adj = -palette;
        else
            palette_adj = palette;
        end
        recon_opts = prediction + palette_adj;
        error_opts = ((recon_opts - val) * weights(i)) .^ 2;
        [error, index] = min(error_opts);
        total_error = total_error + error;
        recon_2 = recon_1;
        recon_1 = recon_opts(index);
    end
end


function [palette, error] = find_palette(data_in, weights)
    palette = [0, 0, 0, 0];
    error = get_compression_error(palette, data_in, weights);
    disp(['Initial: Error: ', num2str(error), ' Palette: ', num2str(palette)]);
 
    %boost = 2;
    %hi = ones(1, boost*256)/2.;
    %mid = ones(1, boost*256)/64.;
    %lo = ones(1, boost*256)/256.;
    %for temperature = [hi mid lo]
        
    t = linspace(.5, -.25, 256*6);
    for temperature = t
        temperature_clamped = max(1./256, temperature);
        test_palette = palette + randn(size(palette)) * temperature_clamped;
        test_palette = max(-1.0, min(1.0, test_palette));
        test_palette = sort(test_palette);
        test_error = get_compression_error(test_palette, data_in, weights);
        if (test_error < error)
            error = test_error;
            palette = test_palette;
            disp(['Temperature: ', num2str(temperature), ' Error: ', num2str(error), ' Palette: ', num2str(palette)]);
        end
    end
end



