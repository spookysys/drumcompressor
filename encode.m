
block_size = 8
[data_in, sample_rate] = load_input('openhat-acoustic01');
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
function [bass_out, bass_env_poly, noise_env_poly, noise_col_b, noise_col_a] = prepare(data_in, sample_rate, block_size)
	nyquist = sample_rate/2.0;
	separation_freq = nyquist/block_size/2.0;
    data_range = (0:size(data_in)-1)';
    env_poly_deg = 3

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
    bass_env_poly = polyfit(data_range, bass_env, env_poly_deg);
    noise_env_poly = polyfit(data_range, noise_env, env_poly_deg);
    bass_env_poly_vals = polyval(bass_env_poly, data_range);
    noise_env_poly_vals = polyval(noise_env_poly, data_range);
    
    % remove falloff from bass and noise
    bass_out = bass ./ bass_env_poly_vals;
    noise_out = noise ./ noise_env_poly_vals;
    
    
    % Quantize this shit
    % [noise_env_d, noise_env_d_scaler] = digitize_noise_env(noise_env);
    % [bass_d, bass_d_scaler] = digitize_bass(bass);

    figure;
	subplot(4,2, 1);
	plot(1:size(data_in), bass);
	subplot(4,2, 2);
	plot(1:size(data_in), noise);
    subplot(4,2, 3);
	plot(1:size(data_in), bass_out);
    ylim([-4, 4])
    subplot(4,2, 4);
	plot(1:size(data_in), noise_out);
    ylim([-4, 4])
    subplot(4,2, 5);
	plot(1:size(data_in), bass_env);
    subplot(4,2, 6);
	plot(1:size(data_in), noise_env);
    subplot(4,2, 7);
	plot(1:size(data_in), bass_env_poly_vals);
    subplot(4,2, 8);
	plot(1:size(data_in), noise_env_poly_vals);
    
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



function ret = constrained_polyfit(x, y, x0, y0, n)
    x = x(:); % reshape the data into a column vector
    y = y(:);

    % 'C' is the Vandermonde matrix for 'x'
    V(:,n+1) = ones(length(x), 1, class(x));
    for j = n:-1:1
         V(:,j) = x.*V(:,j+1);
    end
    C = V;

    % 'd' is the vector of target values, 'y'.
    d = y;
    
    % There are no inequality constraints in this case, i.e., 
    A = [];
    b = [];
    
    % We use linear equality constraints to force the curve to hit the required point. In
    % this case, 'Aeq' is the Vandermoonde matrix for 'x0'
    Aeq = x0.^(n:-1:0);
    % and 'beq' is the value the curve should take at that point
    beq = y0;
    
    p = lsqlin( C, d, A, b, Aeq, beq )
    
    % We can then use POLYVAL to evaluate the fitted curve
    yhat = polyval( p, x );
    
    % Plot original data
    plot(x,y,'.b-') 
    hold on

    % Plot point to go through
    plot(x0,y0,'gx','linewidth',4) 
    
    % Plot fitted data
    plot(x,yhat,'r','linewidth',2) 
    hold off
    
    ret = yhat
end
