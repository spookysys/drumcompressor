
function [data_out, palette_out] = compress_adpcm(data_in, weights, bits)
    if (length(data_in)==0)
        palette_out = [];
        data_out = [];
    else
        [palette_out, ~] = find_palette(data_in, weights, bits);
        data_out = compress_with_palette(data_in, palette_out);       
    end
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



