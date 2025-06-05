%% PROCESSAMENT DE TOTS ELS FITXERS .CDF D'UNA CARPETA

clc;
clear all;
close all;

% === DEFINICIÓ DE CARPETES ===
input_folder = "/Users/marina/Desktop/TFG BIOMED FINAL/ORINA";  % Carpeta d'entrada amb arxius .cdf
output_folder = "/Users/marina/Desktop/TFG BIOMED FINAL/ORINA_CENTROIDED";  % Carpeta de sortida

% Crea la carpeta de sortida si no existeix
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Obté tots els fitxers .cdf dins la carpeta
files = dir(fullfile(input_folder, "*.cdf"));

% === BUCLE PRINCIPAL PER PROCESSAR CADA FITXER ===
for f = 1:length(files)
    % Rutes d'entrada i sortida
    original_file = fullfile(input_folder, files(f).name);
    new_file = fullfile(output_folder, strrep(files(f).name, ".cdf", "_modified.cdf"));
    
    disp(['Processant: ', files(f).name]);

    % Llegeix el fitxer original
    sample = mzcdfread(original_file);

    % Extracció de dades rellevants
    mz_values = round(sample.mass_values);  % Arrodoneix valors m/z
    intensity_values = sample.intensity_values;
    scan_index = sample.scan_index;

    % Inicialització de vectors comprimits
    mz_compressed = [];
    intensity_compressed = [];
    scan_index_compressed = [0];

    % Comprimir cada escaneig: centroidatge bàsic (suma intensitats per m/z únics)
    for i = 1:size(scan_index, 1)
        if i < size(scan_index, 1)
            mz_current = mz_values(scan_index(i)+1 : scan_index(i+1));
            intensity_current = intensity_values(scan_index(i)+1 : scan_index(i+1));
        else
            mz_current = mz_values(scan_index(i)+1 : end);
            intensity_current = intensity_values(scan_index(i)+1 : end);
        end

        % Agrupa per m/z únic mantenint l'ordre i suma intensitats
        [mz_unique, ~, idx] = unique(mz_current, 'stable');
        intensity_unique = accumarray(idx(:), intensity_current(:));

        % Desa valors al vector comprimit
        mz_compressed = [mz_compressed mz_unique'];
        intensity_compressed = [intensity_compressed intensity_unique'];
        scan_index_compressed = [scan_index_compressed scan_index_compressed(end) + length(mz_unique)];
    end

    % Filtra intensitats > 0
    valid_idx = intensity_compressed > 0;
    mz_centroided = mz_compressed(valid_idx);
    intensity_centroided = intensity_compressed(valid_idx);

    % Recalcula scan_index per dades filtrades
    scan_index_centroided = [0];
    current_index = 1;
    for i = 2:length(scan_index_compressed)
        num_valid = sum(valid_idx(current_index : scan_index_compressed(i)));
        scan_index_centroided = [scan_index_centroided, scan_index_centroided(end) + num_valid];
        current_index = scan_index_compressed(i) + 1;
    end

    % Transposa vectors (per compatibilitat amb NetCDF)
    mz_centroided = transpose(mz_centroided);
    intensity_centroided = transpose(intensity_centroided);
    scan_index_centroided = transpose(scan_index_centroided(1:end-1));

    % Recalcul de point_count per escaneig (nombre de punts per espectre)
    point_count_centroided = diff([scan_index_centroided; length(mz_centroided)]);
    point_count_centroided = transpose(point_count_centroided);

    % Obté informació de metadades del fitxer original
    file_info = ncinfo(original_file);

    % Elimina fitxer existent amb el mateix nom
    if exist(new_file, 'file')
        delete(new_file);
    end

    % Torna a llegir el fitxer (per accedir a valors originals)
    sample = mzcdfread(original_file);

    % === ESCRIPTURA DE VARIABLES AL NOU FITXER ===
    for i = 1:length(file_info.Variables)
        var_name = file_info.Variables(i).Name;
        
        switch var_name
            case 'mass_values'
                var_data = mz_centroided;
                var_dim = struct('Name', 'point_number','Length', length(mz_centroided),'Unlimited', 1);
                var_type = file_info.Variables(i).Datatype;
            case 'intensity_values'
                var_data = intensity_centroided;
                var_dim = struct('Name', 'point_number','Length', length(mz_centroided),'Unlimited', 1);
                var_type = 'double';
            case 'scan_index'
                var_data = scan_index_centroided;
                var_dim = file_info.Variables(i).Dimensions;
                var_type = file_info.Variables(i).Datatype;
            case 'point_count'
                var_data = point_count_centroided;
                var_dim = file_info.Variables(i).Dimensions;
                var_type = file_info.Variables(i).Datatype;
            otherwise
                var_data = sample.(var_name);
                var_dim = file_info.Variables(i).Dimensions;
                var_type = file_info.Variables(i).Datatype;
        end

        % Conversió de tipus de dades si cal
        if ischar(var_data)
            var_data = char(var_data);
            var_type = 'char';
        elseif isnumeric(var_data)
            if strcmp(var_type, 'int32')
                var_data = int32(var_data);
            elseif strcmp(var_type, 'double')
                var_data = double(var_data);
            end
        end

        % Crea la variable al nou fitxer NetCDF
        dim_names = {var_dim.Name};
        dim_sizes = num2cell([var_dim.Length]);
        nccreate(new_file, var_name, "Datatype", var_type, "Dimensions", [dim_names; dim_sizes], "Format", "classic");

        % Escriu les dades si existeixen
        if ~isempty(var_data)
            ncwrite(new_file, var_name, var_data);
        else
            disp(['Variable buida, no escrita: ', var_name]);
        end

        % Copia els atributs de la variable
        for a = 1:length(file_info.Variables(i).Attributes)
            ncwriteatt(new_file, var_name, ...
                file_info.Variables(i).Attributes(a).Name, ...
                file_info.Variables(i).Attributes(a).Value);
        end
    end

    % Copia els atributs globals del fitxer original
    for a = 1:length(file_info.Attributes)
        ncwriteatt(new_file, '/', ...
            file_info.Attributes(a).Name, ...
            file_info.Attributes(a).Value);
    end

    disp(['Fitxer processat: ', files(f).name]);
end

disp("Tots els fitxers s'han processat correctament.");

% Mostra informació del darrer fitxer processat
ncdisp(new_file);
