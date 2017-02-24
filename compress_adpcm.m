



function [indices, palette] = compress_adpcm(data_in, weights, bits_per_sample, bitdepth)
    if (length(data_in)==0)
        indices = [];
        palette = [];
    else
        % Determine optimal palette
        palette = [0 0 0 0];
        error = Inf;
        iteration = 0;
        for i = 1:600
            [error, iteration, palette] = optimize_palette(error, iteration, data_in, weights, palette, 2);
        end
        
        % output the thing
        disp('Compressing');
        [error, indices] = compress_with_palette(data_in, weights, palette, 6);
        disp(['Final error: ', num2str(error), ' Palette: ', num2str(palette)]);
        disp('Done');
    end
end




% run test compressions to try to optimize palette
function [error, iteration, palette] = optimize_palette(error, iteration, data_in, weights, palette, lookahead)
    test_palette = palette;
    while test_palette == palette
        mask = randi([0 2], 1, length(palette));
        diffs1 = round(randn(1, length(palette)) * 2);
        diffs2 = round(randn(1, length(palette)) * 4);
        test_palette(mask==0) = palette(mask==0);
        test_palette(mask==1) = palette(mask==1) + diffs1(mask==1);
        test_palette(mask==2) = palette(mask==2) + diffs2(mask==2);
        test_palette(mask==3) = floor(rand() * 256) - 128;
    end
    [test_error, ~] = compress_with_palette(data_in, weights, test_palette, lookahead);
    if (test_error <= error)
        palette = sort(test_palette);
        error = test_error;
        disp(['Iteration ' num2str(iteration) ' Error: ', num2str(error), ' Palette: ', num2str(palette)]);
    end
    iteration = iteration + 1;
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
    prediction = recon_1 + recon_1 - recon_2;
    if (recon_1 < 0)
        recon_opts = wrap8(prediction - palette);
    else
        recon_opts = wrap8(prediction + palette);
    end
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


function [val] = clamp8(val)
    val = double(int8(val));
end

function [val] = wrap8(val)
    val = mod(val+128, 256) - 128;
end

