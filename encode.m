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

    [bass, bass_norm, bass_env_fit_values, treble, treble_env_fit, treble_b, treble_a] = separate(data_in, sample_rate, block_size);

    % write stats
    fid = fopen([checkdir 'stats.txt'],'w');
    size1 = 4 + 4 + 3 + 2 + 1; % bass_fit, treble_fit, treble_b, treble_a, bass_length
    size2 = length(bass_norm) * 2./8; % encoded bass size in bytes (one sample in this array is one block, and 2 bits per block)
    size = size1 + size2;
    fprintf(fid, [ 'Samples: ' num2str(length(data_in)) '\n']);
    fprintf(fid, [ 'Blocks: ' num2str(length(data_in)/block_size) '\n']);
    fprintf(fid, [ 'Bass-blocks: ' num2str(length(bass_norm)) '\n']);
    fprintf(fid, [ 'Final size: ' num2str(size1) ' + ' num2str(size2) ' = ' num2str(size) '\n']);
    fprintf(fid, [ 'kbps: ' num2str(size*8 / (length(data_in)/sample_rate) / 1024) '\n']);
    fprintf(fid, [ 'Compression ratio: ' num2str(length(data_in)/size) '\n']);
    fprintf(fid, [ 'bits per sample: ' num2str(size*8 / length(data_in)) '\n']);
    fclose(fid);
    
    audiowrite([checkdir 'bass.wav'], bass, round(sample_rate/block_size));
    audiowrite([checkdir 'treble.wav'], treble, sample_rate);
        
    treble_recon = reconstruct_treble(length(data_in), treble_env_fit, treble_b, treble_a, sample_rate);
    audiowrite([checkdir 'treble_recon.wav'], treble_recon, sample_rate);

    bass_env_identity = ones(length(bass), 1);
    [bass_norm_adpcm, bass_norm_palette] = compress_adpcm(bass_norm, bass_env_identity, adpcm_bits);
    bass_norm_recon = decompress_adpcm(bass_norm_adpcm, bass_norm_palette);
    bass_recon = bass_norm_recon .* bass_env_fit_values;
    audiowrite([checkdir 'bass_recon.wav'], bass_recon, round(sample_rate/block_size));

    % mix and output
    padding = zeros(length(data_in)/block_size - length(bass), 1);
    bass_recon_padded = [bass_recon; padding];
    bass_recon_upsampled = resample(bass_recon_padded, block_size, 1);
    mix = bass_recon_upsampled + treble_recon;
    audiowrite(outfile, mix, sample_rate);
    audiowrite([checkdir 'result.wav'], mix, sample_rate);
end


% Load wav-file
function [data, sample_rate] = load_input(filename, block_size)
	[data, sample_rate] = audioread(filename);
    padding = ceil(length(data) / block_size) * block_size - length(data); % pad
    pad = zeros(padding, 1);
    data = [data; pad]; % pad
    scaler = max(abs(data)); % normalize
    data = data / scaler; % normalize
end

% Prepare data for compression
function [bass, bass_norm, bass_env_fit_values, treble, treble_env_fit, treble_b, treble_a] = separate(data_in, sample_rate, block_size)
	nyquist = sample_rate/2.0;
	separation_freq = nyquist/block_size/2.0;

    % extract and resample bass
    bass = resample(data_in, 1, block_size);
    %bass = normalize(bass);

    % calculate bass volume envelope
    bass_env_real_values = abs(hilbert(bass));
    
    % fit bass volume envelope
    bass_env_x = (1:length(bass))';    
    bass_env_fit = fit(bass_env_x, bass_env_real_values, 'exp2');
    bass_env_fit_values = feval(bass_env_fit, bass_env_x);
    
    % cut bass at the point where fitted volume drops below 1/256
    bass_cut_point = 0;
    limit = 1. / 2^8;
    for i = 1:length(bass)
        if (bass_env_fit_values(i) > limit && bass_env_real_values(i) > limit)
            bass_cut_point = i;
        end
    end
    bass = bass(1:bass_cut_point);
    bass_env_real_values = bass_env_real_values(1:bass_cut_point);
    bass_env_fit_values = bass_env_fit_values(1:bass_cut_point);

    % normalize bass to fit curve
    bass_norm = bass ./ bass_env_fit_values;
    
	% separate treble
	[sep_treble_b, sep_treble_a] = butter(6, separation_freq / nyquist, 'high');
	treble = filter(sep_treble_b, sep_treble_a, data_in);
    treble = normalize(treble);
	    
	% find volume envelopes
	treble_env = abs(hilbert(treble));
    treble_env_fit_x = (1:length(data_in))';
    treble_env_fit = fit(treble_env_fit_x, treble_env, 'exp2');

    % extract treble color
    num_freqs = 1024;
    treble_h = fft(treble, num_freqs);
    treble_h = treble_h(1:num_freqs/2+1);
    treble_w = (0:(num_freqs/2)) * (2*pi)/num_freqs;
    [treble_b, treble_a] = invfreqz(treble_h, treble_w, 2, 2);
    treble_ab_scaler = max(abs([treble_a(2:end) treble_b]));
    treble_a = treble_a / treble_ab_scaler;
    treble_b = treble_b / treble_ab_scaler;    
end


% Run a simplified compression with given palette and return total error
function total_error = get_compression_error(data_in, palette, weights)
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
        error_opts = (recon_opts - val) .^ 2;
        [error, index] = min(error_opts);
        error = error * weights(i);
        total_error = total_error + error;
        recon_2 = recon_1;
        recon_1 = recon_opts(index);
    end
end

% Run a simplified (for now) compression and return data
function data_out = compress_with_palette(data_in, palette)
    data_out = zeros(length(data_in), 1);
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
        error_opts = (recon_opts - val) .^ 2;
        [~, index] = min(error_opts);
        data_out(i) = index - 1;
        recon_2 = recon_1;
        recon_1 = recon_opts(index);
    end
end


function [palette, error] = find_palette(data_in, weights, bits)
    palette = zeros(1, 2^bits);
    error = get_compression_error(data_in, palette, weights);
    disp(['Initial: Error: ', num2str(error), ' Palette: ', num2str(palette)]);

    t = linspace(.5, -.25, 1000000 / length(data_in));
    for temperature = t
        temperature_clamped = max(1./256, temperature);
        test_palette = palette + randn(1, length(palette)) * temperature_clamped;
        test_palette = max(-1.0, min(1.0, test_palette));
        test_palette = sort(test_palette);
        test_error = get_compression_error(data_in, test_palette, weights);
        if (test_error < error)
            error = test_error;
            palette = test_palette;
            disp(['Temperature: ', num2str(temperature), ' Error: ', num2str(error), ' Palette: ', num2str(palette)]);
        end
    end
end


function [data_out, palette_out] = compress_adpcm(data_in, weights, bits)
    [palette_out, ~] = find_palette(data_in, weights, bits);
    data_out = compress_with_palette(data_in, palette_out);
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


function data_out = reconstruct_treble(num_samples, env_fit, treble_b, treble_a, sample_rate)
    white = rand(num_samples, 1) * 2 - 1;

    colored = filter(treble_b, treble_a, white);
    colored = normalize(colored);

    env_fit_x = (1:num_samples);
    env_values = feval(env_fit, env_fit_x); 
    data_out = colored .* env_values;
end


function [data_out, scaler] = normalize(data_in) 
    scaler = max(abs(data_in));
    data_out = data_in / scaler;
end