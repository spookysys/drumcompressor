block_size = 8;
adpcm_bits = 2;

files = dir('in/*.wav');
for file = files'
    filename = file.name;
    disp(['Processing ' filename]);
    main(['in/' filename], ['out/' filename], ['out_check/' filename '/'], block_size, adpcm_bits);
end


function main(infile, outfile, checkdir, block_size, adpcm_bits)
    samples_per_byte = 8 / adpcm_bits;
    
    mkdir(checkdir);
    [data_in, sample_rate] = load_input(infile, block_size);

    save_output([checkdir 'orig.wav'], data_in, sample_rate);

    [bass_env_fit, bass_sig_norm] = analyze_bass(data_in, sample_rate, block_size, samples_per_byte, checkdir);
    [treble_env_fit, treble_color_b, treble_color_a, len] = analyze_treble(data_in, length(bass_sig_norm)*block_size, sample_rate, block_size, checkdir);
    
    % write stats
    fid = fopen([checkdir 'stats.txt'],'w');
    size1 = 4 + 4 + 3 + 2 + 2 + 4; % bass_fit, treble_fit, treble_b, treble_a, bass_length, palette
    size2 = length(bass_sig_norm) / samples_per_byte; % encoded bass size in bytes (one sample in this array is one block, and 2 bits per block)
    size = size1 + size2;
    fprintf(fid, [ 'Original: ' num2str(length(data_in)) '\r\n']);
    fprintf(fid, [ 'Cropped: ' num2str(len) '\r\n']);
    fprintf(fid, [ 'Bass-Crop: ' num2str(length(bass_sig_norm)*block_size) '\r\n']);
    fprintf(fid, [ 'Final Size: ' num2str(size1) ' + ' num2str(size2) ' = ' num2str(size) '\r\n']);
    fprintf(fid, [ 'Kbps: ' num2str(size*8 / (len/sample_rate) / 1024) '\r\n']);
    fprintf(fid, [ 'Compression Ratio: ' num2str(len/size) '\r\n']);
    fprintf(fid, [ 'Bits Per Sample: ' num2str(size*8 / len) '\r\n']);
    fclose(fid);

    treble_recon = reconstruct_treble(len, treble_env_fit, treble_color_b, treble_color_a, checkdir, sample_rate);
    save_output([checkdir 'treble_recon.wav'], treble_recon, sample_rate);

    adpcm_weights = max(1./256, feval(bass_env_fit, (0:length(bass_sig_norm)-1)'));
    [bass_adpcm_data, bass_adpcm_palette, volume_adj] = compress_adpcm(bass_sig_norm, adpcm_weights, adpcm_bits, 8);
    bass_env_fit.a = bass_env_fit.a / volume_adj;
    bass_env_fit.c = bass_env_fit.c / volume_adj;
    clear volume_adj;
    clear adpcm_weights;
    
    bass_env_fit_values = feval(bass_env_fit, (0:length(bass_sig_norm)-1)');
    bass_norm_recon = decompress_adpcm(bass_adpcm_data, bass_adpcm_palette, 8);
    bass_recon = bass_norm_recon .* bass_env_fit_values;
    if (~isempty(bass_recon))
        save_output([checkdir 'bass_norm_recon.wav'], resample(bass_norm_recon, block_size, 1), sample_rate);
        save_output([checkdir 'bass_recon.wav'], resample(bass_recon, block_size, 1), sample_rate);
        bass_recon_upsampled = resample(bass_recon, block_size, 1);
    else
        bass_recon_upsampled = [];
    end

    % mix and output
    bass_recon_padded = pad_to(bass_recon_upsampled, len, 0);
    treble_recon_padded = pad_to(treble_recon, len, 0);
    mix = bass_recon_padded + treble_recon_padded;
    save_output(outfile, mix, sample_rate);
    save_output([checkdir 'result.wav'], mix, sample_rate);
end


function data_out = clip_sound(data_in, bits)
    max_val = 2^(bits-1)-1;
    min_val = -2^(bits-1);
    data_out = max(min(data_in, max_val), min_val);
end

% Load wav-file
function [data, sample_rate] = load_input(filename, block_size)
	[data, sample_rate] = audioread(filename);
    data = sum(data, 2); % combine channels
    data = pad(data, block_size); % pad
    scaler = max(abs(data)); % normalize
    data = data / scaler; % normalize
end


% Save wav-file
function save_output(filename, data, sample_rate)
    audiowrite(filename, data, sample_rate);
end

% Prepare data for compression
function [env_fit, color_b, color_a, cut_point] = analyze_treble(data_in, min_len, sample_rate, block_size, checkdir)
    nyquist = sample_rate/2.0;
	separation_freq = nyquist/block_size/2.0;

    % Separate treble
	[sep_b, sep_a] = butter(6, separation_freq / nyquist, 'high');
	treble = filter(sep_b, sep_a, data_in);
    save_output([checkdir 'treble.wav'], treble, sample_rate);

	% Find volume envelopes
	env = abs(hilbert(treble));
    env_smooth = smooth(env, 50);
    clear env;
    
    % Fit curve to volume envelope
    env_fit_x = (0:length(data_in)-1)';
    env_fit = fit(env_fit_x, env_smooth, 'exp2');
    
    % find where volume drops below 1/256
    cut_point = 0;
    limit = 1. / 2^8;
    for i = 1:length(treble)
        vol_real = env_smooth(i);
        vol_fit = feval(env_fit, i-1);
        if (vol_real > limit && vol_fit > limit)
            cut_point = i;
        end
    end
    cut_point = max(cut_point, min_len);
    treble = pad_to(treble, cut_point, 0);
    env_smooth = pad_to(env_smooth, cut_point, 'clamp');
    
    % Normalize treble signal
    treble_norm = treble ./ env_smooth;
    if (length(treble_norm)>0)
        save_output([checkdir 'treble_cut.wav'], treble, sample_rate);
        save_output([checkdir 'treble_norm.wav'], treble_norm, sample_rate);
    end
    
    % Extract treble color
    num_freqs = 2^18;
    color_h = fft(treble, num_freqs);
    color_h = mps(color_h);
    color_h = color_h(1:num_freqs/2+1);
    color_w = (0:(num_freqs/2)) * (2*pi)/num_freqs;
    [color_b, color_a] = invfreqz(color_h, color_w, 2, 2, [], 100, 0.01);
    %figure; zplane(color_b, color_a);
    
    % Adjust volume empirically
    white = rand(length(data_in), 1) * 2 - 1;
    colored = filter(color_b, color_a, white);
    volume_adjustment = mean(abs(hilbert(treble_norm))) / mean(abs(hilbert(colored)));
    color_b = color_b * volume_adjustment;
 end


% Prepare data for compression
function [env_fit, bass_norm] = analyze_bass(data_in, sample_rate, block_size, samples_per_byte, checkdir)
	% extract and resample bass
    bass = resample(data_in, 1, block_size);
    save_output([checkdir 'bass.wav'], resample(bass, block_size, 1), sample_rate);
    
    % calculate bass volume envelope
    env = abs(hilbert(bass));
    env_smooth = smooth(env, 50);
    
    % fit volume envelope
    env_fit_x = (0:length(bass)-1)';
    env_fit = fit(env_fit_x, env, 'exp2');
    clear env;
    
    % cut at the point where volume drops below 1/256
    cut_point = 0;
    limit = 1. / 2^8;
    for i = 1:length(bass)
        vol_real = env_smooth(i);
        vol_fit = feval(env_fit, i-1);
        if (abs(vol_real) > limit && abs(vol_fit) > limit)
            cut_point = i;
        end
    end
    cut_point = ceil(cut_point/samples_per_byte) * samples_per_byte;
    bass = pad_to(bass, cut_point, 0);
    env_smooth = pad_to(env_smooth, cut_point, 'clamp');
    
    
    % normalize bass to fit curve
    bass_norm = bass ./ env_smooth;

    if (length(bass_norm)>0)
        save_output([checkdir 'bass_cut.wav'], resample(bass, block_size, 1), sample_rate);
        save_output([checkdir 'bass_norm.wav'], resample(bass_norm, block_size, 1), sample_rate);
    end
end





function data_out = decompress_adpcm(data_in, palette, bitdepth)
    data_out = zeros(length(data_in), 1);
    recon_1 = 0;
    recon_2 = 0;
    for i = 1:length(data_in)
        index = data_in(i);
        palette_val = palette(index);
        recon_slope = recon_1 - recon_2;
        prediction = recon_1 + recon_slope;
        if (recon_1 < 0)
            palette_val = -palette_val;
        end
        recon = prediction + palette_val;
        data_out(i) = recon;
        recon_2 = recon_1;
        recon_1 = recon;
    end
    data_out = data_out ./ 2^(bitdepth-1);
end


function data_out = reconstruct_treble(num_samples, env_fit, color_b, color_a, checkdir, sample_rate)
    white = rand(num_samples, 1) * 2 - 1;
    colored = filter(color_b, color_a, white);
    save_output([checkdir 'treble_norm_recon.wav'], colored, sample_rate);
    env_fit_x = (0:num_samples-1)';
    env_values = feval(env_fit, env_fit_x);
    data_out = colored .* env_values;
end


function data = pad(data_in, block_size)
    num = ceil(length(data_in) / block_size) * block_size - length(data_in);
    padding = zeros(num, 1);
    data = [data_in; padding];
end

function data = pad_to(data_in, len, val)
    if (len > length(data_in))
        if (val == 'clamp') 
            val = data_in(end);
        end
        padding = ones(len - length(data_in), 1) .* val;
        data = [data_in; padding];
    else
        data = data_in(1:len);
    end
end


