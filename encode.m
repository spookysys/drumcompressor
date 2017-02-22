block_size = 8;
adpcm_bits = 2;

files = dir('in/*.wav');
for file = files'
    filename = file.name;
    disp(['Processing ' filename]);
    main(['in/' filename], ['out/' filename], ['out_check/' filename '/'], block_size, adpcm_bits);
end


function main(infile, outfile, checkdir, block_size, adpcm_bits)
    mkdir(checkdir);
    [data_in, sample_rate] = load_input(infile, block_size);

    audiowrite([checkdir 'orig.wav'], data_in, sample_rate);

    [treble_env_fit, treble_color_b, treble_color_a, treble_cut_point] = analyze_treble(data_in, sample_rate, block_size, checkdir);
    [bass_env_fit, bass_sig_norm] = analyze_bass(data_in, sample_rate, block_size, checkdir);
    len = max(length(bass_sig_norm)*block_size, ceil(treble_cut_point/block_size)*block_size);
    
    % write stats
    fid = fopen([checkdir 'stats.txt'],'w');
    size1 = 4 + 4 + 3 + 2 + 2 + 4; % bass_fit, treble_fit, treble_b, treble_a, bass_length, palette
    size2 = ceil(length(bass_sig_norm) * 2./8); % encoded bass size in bytes (one sample in this array is one block, and 2 bits per block)
    size = size1 + size2;
    fprintf(fid, [ 'Original: ' num2str(length(data_in)) '\r\n']);
    fprintf(fid, [ 'Cropped: ' num2str(len) '\r\n']);
    fprintf(fid, [ 'Blocks: ' num2str(len/block_size) '\r\n']);
    fprintf(fid, [ 'Bass-Blocks: ' num2str(length(bass_sig_norm)) '\r\n']);
    fprintf(fid, [ 'Final Size: ' num2str(size1) ' + ' num2str(size2) ' = ' num2str(size) '\r\n']);
    fprintf(fid, [ 'Kbps: ' num2str(size*8 / (len/sample_rate) / 1024) '\r\n']);
    fprintf(fid, [ 'Compression Ratio: ' num2str(len/size) '\r\n']);
    fprintf(fid, [ 'Bits Per Sample: ' num2str(size*8 / len) '\r\n']);
    fclose(fid);

    treble_recon = reconstruct_treble(len, treble_env_fit, treble_color_b, treble_color_a);
    audiowrite([checkdir 'treble_recon.wav'], treble_recon, sample_rate);

    bass_env_fit_x = (0:length(bass_sig_norm)-1)';
    bass_env_fit_values = feval(bass_env_fit, bass_env_fit_x);
    [bass_norm_adpcm, bass_norm_palette] = compress_adpcm(bass_sig_norm, bass_env_fit_values, adpcm_bits);
    bass_norm_recon = decompress_adpcm(bass_norm_adpcm, bass_norm_palette);
    bass_recon = bass_norm_recon .* bass_env_fit_values;
    if (length(bass_recon)>0)
        audiowrite([checkdir 'bass_recon.wav'], resample(bass_recon, block_size, 1), sample_rate);
    end

    % mix and output
    bass_recon_upsampled = resample(bass_recon, block_size, 1);
    bass_recon_padded = pad_to(bass_recon_upsampled, len);
    treble_recon_padded = pad_to(treble_recon, len);
    mix = bass_recon_padded + treble_recon_padded;
    audiowrite(outfile, mix, sample_rate);
    audiowrite([checkdir 'result.wav'], mix, sample_rate);
end


% Load wav-file
function [data, sample_rate] = load_input(filename, block_size)
	[data, sample_rate] = audioread(filename);
    data = sum(data, 2); % combine channels
    data = pad(data, block_size); % pad
    data = normalize(data); % normalize
end

% Prepare data for compression
function [env_fit, color_b, color_a, cut_point] = analyze_treble(data_in, sample_rate, block_size, checkdir)
    nyquist = sample_rate/2.0;
	separation_freq = nyquist/block_size/2.0;

    % Separate treble
	[sep_b, sep_a] = butter(6, separation_freq / nyquist, 'high');
	treble = filter(sep_b, sep_a, data_in);
    audiowrite([checkdir 'treble.wav'], treble, sample_rate);

	% Find volume envelopes
	env = abs(hilbert(treble));
    env_smooth = smooth(env, 50);
    clear env;
    
    % Fit curve to volume envelope
    env_fit_x = (0:length(data_in)-1)';
    env_fit = fit(env_fit_x, env_smooth, 'exp2');
    
    % find where volume drops below 1/512
    cut_point = 0;
    limit = 1. / 2^8;
    for i = 1:length(treble)
        vol_real = env_smooth(i);
        vol_fit = feval(env_fit, i-1);
        if (vol_real > limit && vol_fit > limit)
            cut_point = i;
        end
    end
    treble = treble(1:cut_point);
    env_smooth = env_smooth(1:cut_point);
    
    % Normalize treble signal
    treble_norm = treble ./ env_smooth;
    if (length(treble_norm)>0)
        audiowrite([checkdir 'treble_cut.wav'], treble, sample_rate);
        audiowrite([checkdir 'treble_norm.wav'], treble_norm, sample_rate);
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
function [env_fit, bass_norm] = analyze_bass(data_in, sample_rate, block_size, checkdir)
	% extract and resample bass
    bass = resample(data_in, 1, block_size);
    audiowrite([checkdir 'bass.wav'], resample(bass, block_size, 1), sample_rate);
    
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
    bass = bass(1:cut_point);
    env_smooth = env_smooth(1:cut_point);
    
    % normalize bass to fit curve
    bass_norm = bass ./ env_smooth;

    if (length(bass_norm)>0)
        audiowrite([checkdir 'bass_cut.wav'], resample(bass, block_size, 1), sample_rate);
        audiowrite([checkdir 'bass_norm.wav'], resample(bass_norm, block_size, 1), sample_rate);
    end
end





function data_out = decompress_adpcm(data_in, palette)
    data_out = zeros(length(data_in), 1);
    recon_1 = 0;
    recon_2 = 0;
    for i = 1:length(data_in)
        index = data_in(i);
        palette_val = palette(index+1);
        recon_slope = recon_1 - recon_2;
        prediction = recon_1 + recon_slope;
        if (prediction >= 0)
            palette_val = -palette_val;
        end
        recon = prediction + palette_val;
        data_out(i) = recon;
        recon_2 = recon_1;
        recon_1 = recon;
    end
end


function data_out = reconstruct_treble(num_samples, env_fit, color_b, color_a)
    white = rand(num_samples, 1) * 2 - 1;
    colored = filter(color_b, color_a, white);
    env_fit_x = (0:num_samples-1)';
    env_values = feval(env_fit, env_fit_x);
    data_out = colored .* env_values;
end


function [data_out, scaler] = normalize(data_in) 
    scaler = max(abs(data_in));
    data_out = data_in / scaler;
end

function data = pad(data_in, block_size)
    num = ceil(length(data_in) / block_size) * block_size - length(data_in);
    padding = zeros(num, 1);
    data = [data_in; padding];
end

function data = pad_to(data_in, len)
    if (len > length(data_in))
        padding = zeros(len - length(data_in), 1);
        data = [data_in; padding];
    else
        data = data_in(1:len);
    end
end


