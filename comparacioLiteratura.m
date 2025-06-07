close all
clear all
%% ==== CÀRREGA DE DADES ====
load('biomarcadors_complets.mat'); % Conté 'biomarcadors_matriu'
bio_names = biomarcadors_matriu.ids;
bio_RIs = biomarcadors_matriu.ri;
bio_compounds = biomarcadors_matriu.compounds;
bio_spectra = biomarcadors_matriu.spectra;

lib_table = readtable('libs_ri.xlsx');
lib_names = lib_table{:,1};
lib_RIs = lib_table{:,2};

load('nist_spectra.mat'); % Conté 'spectra' [n_compounds x m/z bins]

%% ==== LLISTA DE COMPOSTOS D'INTERÈS ====
target_names = {
    'Testosterone', 'Androsterone', '5.alpha.-androstan-3.alpha.,17.beta.-diol', ...
    '4-Hydroxytestosterone', '7.beta.-Hydroxydehydroepiandrosterone', ...
    'Dehydroepiandrosterone, (E)-, TMS derivative', ...
    'Dehydroepiandrosterone, (E)-, TBDMS derivative', ...
    'Etiocholanolone, TBDMS derivative', ...
    '5-Hydroxyindoleacetic acid, 3TMS derivative', ...
    'Dehydroepiandrosterone, (3.beta.)-, TMS derivative', ...
    'Vanillyl alcohol, 2TMS derivative', ...
    'Etiocholanolone, TBDMS derivative', ...
    '7.beta.-Hydroxydehydroepiandrosterone', ...
    'Pentasiloxane, 1,1,3,3,5,5,7,7,9,9-decamethyl-', ...
    'Trisiloxane, 1,1,1,5,5,5-hexamethyl-3,3-bis[(trimethylsilyl)oxy]-', ...
    'Phthalic acid', '4-Nitro-4''-chlorodiphenylsulfoxide', ...
    '1-Propylpentachlorotriphosphazene', ...
    '2,3,5,6-Detetrahydrocyclohexanone, 2,6-di-t-butyl-4-hydroxymethylene-', ...
    'Pentanal', '2,6-Dimethylbenzaldehyde', '2,5-Dimethylbenzaldehyde-dnph', ...
    'Ethylbenzene', '3-Heptanone', 'Methyl benzoate, imine', ...
    'Benzaldehyde, 3-methyl-', 'Butyrolactone', 'Methyl vinyl ketone', ...
    'Methylamine', 'N-Ethylformamide', 'Acetonitrile, (dimethylamino)-', ...
    'Pyridine', 'N‐methylformamide', 'Acetaldehyde', 'Acetamide', ...
    '2-Methylpiperidine', '1-Piperidineacetonitrile', 'Dimethylamine', ...
    'Pyrrole', 'Methacrolein', 'N,N-dimethylformamide-D7', ...
    '2-Octanone', '2-n-Butylacrolein', 'Disulfide, methyl propyl'
};

% Filtrar referències d'interès
found_idx = find(ismember(lib_names, target_names));
ref_names = lib_names(found_idx);
ref_RIs = lib_RIs(found_idx);
ref_spectra = spectra(found_idx, 35:400); % Rang m/z d'interès

% Convertir RI a numèric
ref_RIs_num = NaN(size(ref_RIs));
for k = 1:length(ref_RIs)
    if ischar(ref_RIs{k}) && ~isempty(ref_RIs{k})
        ref_RIs_num(k) = str2double(ref_RIs{k});
    end
end

%% ==== FUNCIÓ DE COSINE SIMILARITY AMB MÀSCARA DE m/z ====
function cs = masked_cosine_similarity(ref, query)
    mask = ref > 0;
    ref_masked = ref(mask);
    query_masked = query(mask);
    if any(ref_masked) && any(query_masked)
        cs = dot(ref_masked, query_masked) / ...
             (norm(ref_masked) * norm(query_masked));
    else
        cs = 0;
    end
end

%% ==== CÀLCUL DE COINCIDÈNCIES ====
n_features = size(bio_spectra, 1);
detailed_results = {};

for i = 1:n_features
    query_id = bio_names(i);
    query_name = bio_compounds{i};
    query_spectrum = bio_spectra(i, :);
    query_RI = bio_RIs(i);

    % Cosine similarity clàssic
    cosine_sim = (ref_spectra * query_spectrum') ./ (vecnorm(ref_spectra, 2, 2) .* norm(query_spectrum));

    % Reverse similarity personalitzada (només m/z amb intensitat al ref)
    reverse_sim = zeros(size(ref_spectra, 1), 1);
    for j = 1:size(ref_spectra, 1)
        reverse_sim(j) = masked_cosine_similarity(ref_spectra(j,:), query_spectrum);
    end

    % Score combinat
    combined_score = (cosine_sim + reverse_sim) / 2;

    % TOP 1: millor score amb RI ±50
    mask_ri50 = abs(ref_RIs_num - query_RI) <= 50;
    combined_score_ri50 = combined_score;
    combined_score_ri50(~mask_ri50) = -Inf;
    [~, idx_top1] = max(combined_score_ri50);

    % TOP 2: millor score espectral sense restricció
    [~, idx_top2] = max(combined_score);

    % TOP 3: millor score amb RI ±30
    mask_ri30 = abs(ref_RIs_num - query_RI) <= 30;
    combined_score_ri30 = combined_score;
    combined_score_ri30(~mask_ri30) = -Inf;
    [~, idx_top3] = max(combined_score_ri30);

    % Guardar resultats
    detailed_results{end+1,1} = query_id;
    detailed_results{end,2} = query_RI;
    detailed_results{end,3} = query_name;

    % Top 1
    detailed_results{end,4} = ref_names(idx_top1);
    detailed_results{end,5} = ref_RIs_num(idx_top1);
    detailed_results{end,6} = cosine_sim(idx_top1);
    detailed_results{end,7} = reverse_sim(idx_top1);
    detailed_results{end,8} = combined_score(idx_top1);

    % Top 2
    detailed_results{end,9} = ref_names(idx_top2);
    detailed_results{end,10} = ref_RIs_num(idx_top2);
    detailed_results{end,11} = cosine_sim(idx_top2);
    detailed_results{end,12} = reverse_sim(idx_top2);
    detailed_results{end,13} = combined_score(idx_top2);

    % Top 3
    detailed_results{end,14} = ref_names(idx_top3);
    detailed_results{end,15} = ref_RIs_num(idx_top3);
    detailed_results{end,16} = cosine_sim(idx_top3);
    detailed_results{end,17} = reverse_sim(idx_top3);
    detailed_results{end,18} = combined_score(idx_top3);
end

%% ==== EXPORTACIÓ ====
detailed_table = cell2table(detailed_results, ...
    'VariableNames', {
        'FeatureID', 'FeatureRI', 'FeatureCompound', ...
        'Top1_Name', 'Top1_RI', 'Top1_CosSim', 'Top1_RevSim', 'Top1_CombSim', ...
        'Top2_Name', 'Top2_RI', 'Top2_CosSim', 'Top2_RevSim', 'Top2_CombSim', ...
        'Top3_Name', 'Top3_RI', 'Top3_CosSim', 'Top3_RevSim', 'Top3_CombSim'
    });

writetable(detailed_table, 'resultats_top3.xlsx');







