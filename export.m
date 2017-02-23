function [ output_args ] = export(wav_filename, bass_data, bass_palette, bass_env_fit, treble_b, treble_a, treble_env_fit)
    [data_fid, struct_fid, name] = start(wav_filename);
    
    % output bass_data
    bass_data_reshaped = reshape(bass_data, [4, length(bass_data)/4]);
    bass_data_packed = zeros(1, length(bass_data)/4);
    for i = 1:length(bass_data_reshaped)
        tmp = bass_data_reshaped(:,i)';
        tmp = tmp .* [2^0 2^2 2^4 2^6];
        tmp = sum(tmp);
        bass_data_packed(i) = tmp;
    end
    fprintf(data_fid, '%.0f,' , bass_data_packed);
    clear bass_data_reshaped;
    clear bass_data_packed;
    
    % output struct
    bass_env_coeff = coeffvalues(bass_env_fit);
    treble_env_coeff = coeffvalues(treble_env_fit);
    
    % a and c values should be scaled to -128:+127 to give max volume
    ac = [bass_env_coeff(1), bass_env_coeff(3), treble_env_coeff(1), treble_env_coeff(3)];
    boost_candidates = [127.49/max(ac) -128.49/min(ac)];
    boost = min(boost_candidates(boost_candidates > 0));
    ac = ac * boost;
    ac = round(ac);
    bass_env_coeff(1) = ac(1);
    bass_env_coeff(3) = ac(2);
    treble_env_coeff(1) = ac(3);
    treble_env_coeff(3) = ac(4);
    
    fclose(data_fid);
    fclose(struct_fid);
end

function [data_fid, struct_fid, name] = start(wav_filename)
     [~,~,~] = mkdir(fullfile('export', 'data'));
     [~,~,~] = mkdir(fullfile('export', 'struct'));
     [~, name, ~] = fileparts(wav_filename);
     name = lower(name);
     
     fid = fopen(fullfile('export', 'datas.inc'), 'a');
     fprintf(fid, [ '#ifdef INC_' upper(name) '\n']);
     fprintf(fid, [ 'const unsigned char ' name '_data[] = {\n']);
     fprintf(fid, [ '#include \"export/data/' name '.inc\"\n']);
     fprintf(fid, [ '};\n']);
     fprintf(fid, [ '#endif\n']);
     fprintf(fid, [ '\n']);
     fclose(fid);
     clear fid;
     
     fid = fopen(fullfile('export', 'structs.inc'), 'a');
     fprintf(fid, [ '#ifdef INC_' upper(name) '\n']);
     fprintf(fid, [ '{\n']);
     fprintf(fid, [ '#include \"export/struct/' name '.inc\"\n']);
     fprintf(fid, [ '};\n']);
     fprintf(fid, [ '#endif\n']);
     fprintf(fid, [ '\n']);
     fclose(fid);
     clear fid;     
          
     fid = fopen(fullfile('export', 'enums.inc'), 'a');
     fprintf(fid, [ '#ifdef INC_' upper(name) '\n']);
     fprintf(fid, [ upper(name) ',\n']);
     fprintf(fid, [ '#endif\n']);
     fclose(fid);
     clear fid;
     
     data_fid = fopen(fullfile('export', 'data', [name '.inc']), 'w');
     struct_fid = fopen(fullfile('export', 'struct', [name '.inc']), 'w');
end
