



function [indices, palette, volume_adj] = compress_adpcm(data_in, weights, bits_per_sample, bitdepth)
    if (length(data_in)==0)
        indices = [];
        palette = [];
        volume_adj = 1;
    else
        scale = 2^(bitdepth-1);
        minval = -scale;
        maxval =  scale-1;
        data = data_in .* scale;
        clear data_in;
        
        % first estimate using approximate algorithm
        [error, palette] = first_estimate(data, weights, bits_per_sample);
        
        % palette must fit in one byte per value, so we must rescale
        volume_adj = minval / min(palette);
        palette = palette .* volume_adj;
        assert(max(palette) <= maxval);
        data = data .* volume_adj;
        
        % output the thing
        indices = compress_with_palette(data, palette);
    end
end


% Simulated annealing with histogram-based error estimation
function [error, palette] = first_estimate(data, weights, bits_per_sample)
    [hist_deltas, hist_weights] = gen_hist_map(data, weights); % generate histogram for quick palette check        
    maxpal = max(hist_deltas);
    minpal = min(hist_deltas);
    disp('first_estimate');
    disp(['max: ' num2str(maxpal) ' min: ' num2str(minpal)]);
    error = -1;
    palette = ones(1, 2^bits_per_sample) * ((maxpal+minpal)/2);
    for temperature = linspace((maxpal-minpal)/2, 0, 10000000 / length(data))
        test_palette = palette + randn(1, length(palette)) * temperature;
        test_error = quick_palette_test(hist_deltas, hist_weights, test_palette);
        if (error==-1 || test_error <= error)
            test_palette = min(maxpal, max(minpal, test_palette));
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
function data_out = compress_with_palette(data_in, palette)
    data_out = zeros(length(data_in), 1);
    recon_1 = 0;
    recon_2 = 0;
    for i = 1:length(data_in)
        val = round(data_in(i));
        recon_slope = recon_1 - recon_2;
        prediction = recon_1 + recon_slope;
        if (recon_1 < 0)
            palette_adj = -palette;
        else
            palette_adj = palette;
        end
        recon_opts = prediction + palette_adj;
        error_opts = abs(recon_opts - val);
        [~, index] = min(error_opts);
        data_out(i) = index - 1;
        recon_2 = recon_1;
        recon_1 = recon_opts(index);
    end
end


