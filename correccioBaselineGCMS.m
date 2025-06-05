% =========================================================================
% SCRIPT: Correcció de baseline en fitxers .CDF de GC-MS (TFG BIOMED)
%
% Aquest script llegeix tots els fitxers .CDF d’una carpeta, reconstrueix 
% els perfils d’ions extrets (EICs) per a cada valor m/z, aplica una 
% correcció de línia base per eliminar senyals de fons mitjançant `movmin`, 
% i reconstrueix el fitxer amb les dades corregides. Es calcula també la 
% variable `point_count` per mantenir la compatibilitat amb l’estructura 
% NetCDF esperada.
%
% Output: fitxers .CDF corregits guardats en una nova carpeta.
% =========================================================================

clc; 
clear all;
close all;

% === PARÀMETRES ===
input_folder = '/Users/marina/Desktop/TFG BIOMED FINAL/aaa';  % Carpeta d'origen
output_folder = [input_folder, ' Baseline Corrected'];         % Carpeta de sortida

% === CREAR CARPETA DE SORTIDA SI NO EXISTEIX ===
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% === LLISTAR FITXERS .CDF ===
files = dir(fullfile(input_folder, '*.cdf'));

% === PROCESSAR CADA FITXER ===
for f = 1:length(files)
    original_file = fullfile(input_folder, files(f).name);
    output_file = fullfile(output_folder, files(f).name);

    disp(['Processant: ', files(f).name]);

    % === LLEGIR FITXER .CDF ===
    sample = mzcdfread(original_file);
    mz_vals = round(sample.mass_values);
    intensities = sample.intensity_values;
    scan_index = sample.scan_index;
    n_scans = length(scan_index);
    
    % === GENERAR EICs PER CADA m/z ÚNIC ===
    mz_unique = unique(mz_vals);
    mz_eics = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
    for i = 1:length(mz_unique)
        mz = mz_unique(i);
        eic = zeros(1, n_scans);
        for s = 1:n_scans
            if s < n_scans
                idx = scan_index(s)+1 : scan_index(s+1);
            else
                idx = scan_index(s)+1 : length(mz_vals);
            end
            mz_scan = mz_vals(idx);
            inten_scan = intensities(idx);
            eic(s) = sum(inten_scan(mz_scan == mz));
        end
        mz_eics(mz) = eic;
    end

    % === CORREGIR BASELINE AMB MOVMIN PER CADA m/z ===
    for i = 1:length(mz_unique)
        mz = mz_unique(i);
        y = mz_eics(mz);
        if all(y == 0), continue; end
        baseline = movmin(y, 50, 'Endpoints', 'shrink');
        y_corr = max(y - baseline, 0);
        mz_eics(mz) = y_corr;
    end

    % === RECONSTRUIR DADES CORREGIDES ===
    mz_corr = []; int_corr = []; new_scan_index = [0];
    for s = 1:n_scans
        mzs = []; ints = [];
        for i = 1:length(mz_unique)
            mz = mz_unique(i);
            inten = mz_eics(mz);
            if inten(s) > 0
                mzs(end+1) = mz;
                ints(end+1) = inten(s);
            end
        end
        mz_corr = [mz_corr mzs];
        int_corr = [int_corr ints];
        new_scan_index(end+1) = new_scan_index(end) + length(mzs);
    end

    mz_corr = mz_corr(:);
    int_corr = int_corr(:);
    new_scan_index = new_scan_index(1:end-1)';

    % === CALCULAR point_count ===
    point_count = diff([new_scan_index; length(mz_corr)]);
    point_count = point_count(:);

    % === COPIAR METADADES I GUARDAR FITXER CORREGIT ===
    file_info = ncinfo(original_file);
    if exist(output_file, 'file')
        delete(output_file);
    end

    sample = mzcdfread(original_file);
    for i = 1:length(file_info.Variables)
        name = file_info.Variables(i).Name;
        switch name
            case 'mass_values'
                data = mz_corr;
                dim = struct('Name', 'point_number','Length', length(mz_corr),'Unlimited', 1);
                dtype = file_info.Variables(i).Datatype;
            case 'intensity_values'
                data = int_corr;
                dim = struct('Name', 'point_number','Length', length(mz_corr),'Unlimited', 1);
                dtype = 'double';
            case 'scan_index'
                data = new_scan_index;
                dim = file_info.Variables(i).Dimensions;
                dtype = file_info.Variables(i).Datatype;
            case 'point_count'
                data = point_count;
                dim = file_info.Variables(i).Dimensions;
                dtype = file_info.Variables(i).Datatype;
            otherwise
                if isfield(sample, name)
                    data = sample.(name);
                else
                    data = [];
                end
                dim = file_info.Variables(i).Dimensions;
                dtype = file_info.Variables(i).Datatype;
        end

        % Crear i escriure variable només si no és buida
        if ~isempty(data)
            nccreate(output_file, name, 'Datatype', dtype, 'Dimensions', {dim.Name, dim.Length}, 'Format', 'classic');
            ncwrite(output_file, name, data);
        else
            disp(['Variable buida, no escrita: ', name]);
        end

        % Copiar atributs
        for a = 1:length(file_info.Variables(i).Attributes)
            att = file_info.Variables(i).Attributes(a);
            ncwriteatt(output_file, name, att.Name, att.Value);
        end
    end

    % === COPIAR ATRIBUTS GLOBALS ===
    for a = 1:length(file_info.Attributes)
        attr = file_info.Attributes(a);
        val = attr.Value;
        ncwriteatt(output_file, '/', attr.Name, val);
    end

    disp(['Fitxer corregit i guardat: ', files(f).name]);
end

disp("Tots els fitxers han estat processats amb correcció de baseline.");