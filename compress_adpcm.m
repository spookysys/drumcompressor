



function [indices, palette, volume_adj] = compress_adpcm(data_in, weights, bits_per_sample, bitdepth)
    if (length(data_in)==0)
        indices = [];
        palette = [];
        volume_adj = 1;
    else
        % produce a reasonable volume for histogram algorithm to work
        scale = 2^(bitdepth-1);
        minval = -scale;
        maxval =  scale-1;
        data = data_in; clear data_in;
        
        % signal should be in signed byte range
        signal_adj = [minval-.49 maxval+.49] ./ [min(data) max(data)]; % scale smallest data value to -128 and biggest data value to +127
        signal_adj = min(signal_adj(signal_adj > 0));
        data = data .* signal_adj;

        % starting point
        palette = [0 0 0 0];
       
        % Optimize palette
        effort = 300;
        [error, palette] = optimize_palette(data, weights, palette, 2, 100, effort * min(1, 2000/length(data)), bits_per_sample, bitdepth, 1.5);
        

        % both palette and signal should be in signed byte range
        palette_adj = [minval-.49 maxval+.49] ./ [min(palette) max(palette)]; % scale most negative delta to -128 and most positive delta to +127
        palette_adj = min(palette_adj(palette_adj > 0));
        palette_adj = min(1, palette_adj); % never increase volume as this would bring signal out of byte range
        
        palette = round(palette .* palette_adj);
        data = data .* palette_adj;
        
        disp(['max data: ' num2str(max(data)) ' min data: ' num2str(min(data))]);
        disp(['Adjusted palette: ' num2str(palette)]);
                
        % Optimize palette again
        effort = 300;
        [error, palette] = optimize_palette(data, weights, palette, 2, 3, effort * min(1, 2000/length(data)), bits_per_sample, bitdepth, 1);
        
        % output the thing
        disp('Compressing');
        [error, indices] = compress_with_palette(data, weights, palette, 6);
        disp(['Final error: ', num2str(error), ' Palette: ', num2str(palette)]);
        disp('Done');
        
        % output result of volume adjustments made
        volume_adj = signal_adj * palette_adj;
    end
end




% run test compressions to try to optimize palette
function [error, palette] = optimize_palette(data_in, weights, palette, lookahead, temperature, iterations, bits_per_sample, bitdepth, clamping)
    disp(['Lookahead: ' num2str(lookahead) ' Temperature: ' num2str(temperature) ' Iterations: ' num2str(iterations)]);
    scale = 2^(bitdepth-1); % produce a reasonable volume for histogram algorithm to work
    minpal = (-scale) * clamping;
    maxpal = (scale-1) * clamping;
    [error, ~] = compress_with_palette(data_in, weights, palette, lookahead);
    disp(['Iteration 0 Error: ' num2str(error) ' Palette: ' num2str(palette)]);
    for i = 1:iterations
        diff = zeros(1, length(palette));
        while (sum(abs(diff))==0)
            diff = round(randn(1, length(palette)) * temperature);
        end
        test_palette = palette + diff;
        test_palette = min(maxpal, max(minpal, test_palette));
        [test_error, ~] = compress_with_palette(data_in, weights, test_palette, lookahead);
        if (test_error <= error)
            test_palette = sort(test_palette);
            error = test_error;
            palette = test_palette;
            disp(['Iteration ' num2str(i) ' Error: ', num2str(error), ' Palette: ', num2str(palette)]);
        end
    end
end


% Simulated annealing with histogram-based error estimation
function [error, palette] = first_estimate(data, weights, bits_per_sample)
    [hist_deltas, hist_weights] = gen_hist_map(data, weights); % generate histogram for quick palette check        
    maxpal = max(hist_deltas);
    minpal = min(hist_deltas);
    disp('first_estimate');
    disp(['max: ' num2str(maxpal) ' min: ' num2str(minpal)]);
    palette = ones(1, 2^bits_per_sample) * ((maxpal+minpal)/2);
    error = quick_palette_test(hist_deltas, hist_weights, palette);
    disp(['Initial: Error: ', num2str(error), ' Palette: ', num2str(palette)]);
    for temperature = linspace((maxpal-minpal)/2, 0, 100000000 / length(data))
        test_palette = palette + randn(1, length(palette)) * temperature;
        test_palette = min(maxpal, max(minpal, test_palette));
        test_error = quick_palette_test(hist_deltas, hist_weights, test_palette);
        if (test_error <= error)
            test_palette = sort(test_palette);
            error = test_error;
            palette = test_palette;
            disp(['Temperature: ', num2str(temperature), ' Error: ', num2str(error), ' Palette: ', num2str(palette)]);
        end
    end
    disp('end');
end



% prepare for quick_palette_test
function [hist_deltas, hist_weights] = gen_hist_map(data_in, weights)
    hist_map = containers.Map('KeyType','int32', 'ValueType','double');
    recon_1 = 0;
    recon_2 = 0;
    for i = 1:length(data_in)
        val = data_in(i);
        recon_slope = recon_1 - recon_2;
        prediction = recon_1 + recon_slope;
        delta = val - prediction;
        if (recon_1 < 0)
            delta = -delta;
        end
        key = round(delta);
        inc = weights(i);
        if isKey(hist_map, key)
            hist_map(key) = hist_map(key) + inc;
        else
            hist_map(key) = inc;
        end
        if (recon_1 < 0)
            estimate = prediction - delta;
        else
            estimate = prediction + delta;
        end
        assert(abs(estimate-val)<0.1);
        recon_2 = recon_1;
        recon_1 = val;
    end
    hist_deltas = double(cell2mat(keys(hist_map)));
    hist_weights = cell2mat(values(hist_map));
end



% todo: optimize
function [total_error] = quick_palette_test(hist_deltas, hist_weights, palette)
    error_matrix = abs(palette' - hist_deltas);
    column_errors = min(error_matrix);
    column_errors_weighted = column_errors .* hist_weights;
    total_error = sum(column_errors_weighted);
end




% Run a simplified (for now) compression and return data
function [total_error, data_out] = compress_with_palette(data_in, weights, palette, lookahead)
    data_out = zeros(length(data_in), 1);
    recon_1 = 0;
    recon_2 = 0;
    total_error = 0;
    for i = 1:length(data_in)
        [local_error, ~, recon, index] = compress_with_palette_step(data_in, weights, palette, i, recon_1, recon_2, lookahead);
        total_error = total_error + local_error;
        data_out(i) = index;
        recon_2 = recon_1;
        recon_1 = recon;
    end
end


function [local_error, total_error, recon, index] = compress_with_palette_step(data_in, weights, palette, pos, recon_1, recon_2, lookahead)
    val = round(data_in(pos));
    recon_slope = recon_1 - recon_2;
    prediction = recon_1 + recon_slope;
    if (recon_1 < 0)
        palette_adj = -palette;
    else
        palette_adj = palette;
    end
    recon_opts = prediction + palette_adj;
    local_error_opts = abs(recon_opts - val) * weights(pos);
    if lookahead>0 && pos<length(data_in)
        tail_error_opts = zeros(1, length(palette));
        for i = 1:length(palette)
            recon = recon_opts(i);
            [~, tail_error, ~, ~] = compress_with_palette_step(data_in, weights, palette, pos+1, recon, recon_1, lookahead-1);
            tail_error_opts(i) = tail_error;
        end
        total_error_opts = local_error_opts + tail_error_opts;
    else
        total_error_opts = local_error_opts;
    end
    [total_error, index] = min(total_error_opts);
    recon = recon_opts(index);
    local_error = local_error_opts(index);
end



