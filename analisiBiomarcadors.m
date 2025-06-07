close all 
clear all 

%% PART 1: Preprocessament de dades GC-MS

% === Càrrega de dades ===
% S'importa la taula amb les dades GC-MS
features_merged = readtable('/Users/marina/features_merged_urine.xlsx');

% Es crea la matriu numèrica amb les intensitats (columnes 6 fins a la tercera per la cua)
data_raw = table2array(features_merged(:, 6:(width(features_merged)-3)));
[numFeatures, numSamples] = size(data_raw);

% Es recuperen metadades i noms d’identificació
headers = features_merged.Properties.VariableNames;
sampleNames = headers(6:end-3);                 % noms de les mostres
featureNames = features_merged.id;              % identificadors de les característiques (features)
compoundNames = features_merged.Compound;       % noms dels compostos associats

% === Etiquetatge de mostres segons prefix del nom ===
labels = cell(size(sampleNames));
for i = 1:length(sampleNames)
    if startsWith(sampleNames{i}, 'CO')
        labels{i} = 'CO';        % Control
    elseif startsWith(sampleNames{i}, 'CP')
        labels{i} = 'CP';        % Pacient amb càncer de pròstata
    elseif startsWith(sampleNames{i}, 'BL') || startsWith(sampleNames{i}, 'BLK')
        labels{i} = 'BLK';       % Mostra blanc (blank)
    else
        labels{i} = 'ALTRE';     % Altres casos no identificats
    end
end

% === Índexs lògics per a cada grup de mostres ===
idxCO = strcmp(labels, 'CO');
idxCP = strcmp(labels, 'CP');
idxBLK = strcmp(labels, 'BLK');
nCO = sum(idxCO); nCP = sum(idxCP); nBLK = sum(idxBLK);

% === Separació de dades entre blancs i mostres biològiques ===
dataBLK = data_raw(:, idxBLK);        % dades de blancs
dataBiol = data_raw(:, ~idxBLK);      % dades de mostres biològiques
labelsBiol = labels(~idxBLK);
featureNamesBiol = featureNames;

% === 1. Filtrat per presència mínima (mínim 60% de presència en CO o CP) ===
dataCO = dataBiol(:, strcmp(labelsBiol, 'CO'));
dataCP = dataBiol(:, strcmp(labelsBiol, 'CP'));
presentCO = sum(dataCO > 0, 2) / size(dataCO, 2) >= 0.8;
presentCP = sum(dataCP > 0, 2) / size(dataCP, 2) >= 0.8;
keepPresence = presentCO | presentCP;

fprintf('Filtrat per presència 60%%: %d de %d característiques retingudes (%d eliminades)\n', ...
    sum(keepPresence), numFeatures, numFeatures - sum(keepPresence));

% === 2. Eliminació de features afectades per senyal en blanc (si la mitjana del blanc > 10% del senyal biològic) ===
meanBlank = mean(dataBLK, 2);
meanBiol = mean([dataCO, dataCP], 2);
ratioBlank = meanBlank ./ (meanBiol + eps);    % eps evita divisió per zero
keepBlank = ratioBlank <= 0.1;

fprintf('Filtrat per blancs (>10%% en mitjana): %d de %d característiques eliminades\n', ...
    sum(~keepBlank), numFeatures);

% === Aplicació dels filtres combinats ===
keepFeatures = keepPresence & keepBlank;
filteredData = dataBiol(keepFeatures, :);
filteredFeatureNames = featureNames(keepFeatures);
filteredcompoundNames = compoundNames(keepFeatures);

% === 3. Imputació de valors zero (0) amb 1/5 de la mitjana dels valors no nuls de la mateixa feature ===
imputedData = filteredData;
for i = 1:size(filteredData, 1)
    featureVals = filteredData(i, :);
    nonZeroVals = featureVals(featureVals > 0);
    if ~isempty(nonZeroVals)
        imputVal = mean(nonZeroVals) / 5;
        featureVals(featureVals == 0) = imputVal;
        imputedData(i, :) = featureVals;
    end
end
fprintf('Zeros imputats com a mitjana/5 per feature (ignorant zeros en el càlcul)\n');

% === 4. Normalització per TIC (Total Ion Current): es divideix cada columna pel seu total ===
normdata = imputedData ./ sum(imputedData, 1);

% === 5. Autoscaling per z-score (normalització de cada fila a mitjana 0 i desviació típica 1) ===
scaleddata = zscore(normdata, 0, 2);

% === Etiquetes finals (només mostres biològiques) ===
labelsFinal = labels(~idxBLK);

fprintf('--- Preprocessament completat ---\n');
fprintf('Features finals: %d\n', size(scaleddata, 1));
fprintf('Mostres finals: %d (CO: %d, CP: %d)\n', ...
    length(labelsFinal), sum(strcmp(labelsFinal, 'CO')), sum(strcmp(labelsFinal, 'CP')));
%% PART 2: Anàlisi estadística univariant i selecció preliminar de biomarcadors
fprintf('\n--- PART 2: Anàlisi univariant de biomarcadors ---\n');

% === Eliminació de mostres BLANK ===
data = data_raw(:, ~idxBLK);           % Es descarten els blancs
sampleNames = sampleNames(~idxBLK);    % Es filtren els noms de mostra
labels = labelsBiol;                   % Es mantenen les etiquetes de mostres biològiques

% Índexs per a controls i pacients
idxCO = strcmp(labels, 'CO');
idxCP = strcmp(labels, 'CP');

% === Divisió en conjunt d'entrenament i test (70% / 30%) ===
rng(1);  % Assegura resultats reproduïbles
cv = cvpartition(strcmp(labels, 'CP'), 'HoldOut', 0.3);
trainIdx = training(cv);
testIdx = test(cv);

% === Separació de dades segons el conjunt ===
Xtrain_scaled = scaleddata(:, trainIdx);
Xtest_scaled = scaleddata(:, testIdx);
Xtrain_norm   = normdata(:, trainIdx);
labels_train  = labels(trainIdx);
sampleNames_train = sampleNames(trainIdx);
idxCO_train   = strcmp(labels_train, 'CO');
idxCP_train   = strcmp(labels_train, 'CP');
labels_bin_train = [zeros(1,sum(idxCO_train)), ones(1,sum(idxCP_train))];

% Inicialització de vectors per resultats estadístics
numFeaturesFiltered = size(Xtrain_scaled, 1);
pvals = zeros(numFeaturesFiltered,1);
fc = zeros(numFeaturesFiltered,1);
testType = strings(numFeaturesFiltered,1);
auc_values = zeros(numFeaturesFiltered,1);
normalitat = false(numFeaturesFiltered,1);

% === Cicle d'anàlisi estadística per cada feature ===
for i = 1:numFeaturesFiltered
    x1 = Xtrain_scaled(i, idxCO_train);
    x2 = Xtrain_scaled(i, idxCP_train);

    % Test de normalitat (Lilliefors)
    [h1,~] = lillietest(x1');
    [h2,~] = lillietest(x2');
    normal = (h1 == 0) && (h2 == 0);
    normalitat(i) = normal;

    fprintf('Feature %d (%d): ', i, filteredFeatureNames(i));
    if normal
        fprintf('NORMAL\n');
        [~, p] = ttest2(x1, x2, 'Vartype','unequal');
        testType(i) = "t-test";
    else
        fprintf('NO normal\n');
        p = ranksum(x1, x2);
        testType(i) = "Mann–Whitney";
    end
    pvals(i) = p;

    % Càlcul de fold change (FC)
    fc(i) = mean(Xtrain_norm(i, idxCP_train)) / mean(Xtrain_norm(i, idxCO_train));
    
    % Càlcul de l'AUC (àrea sota la corba ROC)
    scores = [x1, x2];
    [~, ~, ~, AUC] = perfcurve(labels_bin_train, scores, 1);
    auc_values(i) = AUC;
end
%%
%% PCA abans del filtratge

% Dades sense filtrar (amb BLKs eliminats)
rawBiolData = data_raw(:, ~idxBLK);  % sense blanks
labelsBiolAll = labels;

% Normalització per TIC
rawNorm = rawBiolData ./ sum(rawBiolData, 1);

% Imputació simple per evitar zeros (opcional, per evitar errors amb zscore)
rawNorm(rawNorm == 0) = min(rawNorm(rawNorm > 0)) / 5;

% Autoscaling
rawScaled = zscore(rawNorm, 0, 2);

% PCA
[coeff_raw, score_raw, ~, ~, expl_raw] = pca(rawScaled');

% Representació
figure;
gscatter(score_raw(:,1), score_raw(:,2), labelsBiolAll', 'br', 'xo', 8);
xlabel(['PC1 (' num2str(expl_raw(1), '%.1f') '%)']);
ylabel(['PC2 (' num2str(expl_raw(2), '%.1f') '%)']);
title('PCA abans del filtratge');
legend('CO','CP');
grid on;
box on;

%% PCA després del filtratge

[coeff_filt, score_filt, ~, ~, expl_filt] = pca(scaleddata');

figure;
gscatter(score_filt(:,1), score_filt(:,2), labelsFinal', 'br', 'xo', 8);
xlabel(['PC1 (' num2str(expl_filt(1), '%.1f') '%)']);
ylabel(['PC2 (' num2str(expl_filt(2), '%.1f') '%)']);
title('PCA després del filtratge');
legend('CO','CP');
grid on;
box on;

%% Correcció FDR i selecció
[~, ~, ~, adjP] = fdr_bh(pvals);
alpha = 0.05;
fcThreshold = 1;  
biomarcador_idx = (adjP < alpha) & (abs(log2(fc)) > fcThreshold);
[~, sorted_idx] = sort((auc_values(biomarcador_idx)), 'descend');

%% === Volcano Plot ===
fprintf('\n--- Volcano Plot ---\n');

log2FC = log2(fc);
negLog10P = -log10(adjP);
fc_thresh = log2(2);
p_thresh = 0.05;

significant = (adjP < p_thresh) & (abs(log2FC) > fc_thresh);
colorVec = repmat([0.6 0.6 0.6], numFeaturesFiltered, 1);
colorVec(significant, :) = repmat([1 0 0], sum(significant), 1);

figure;
scatter(log2FC, negLog10P, 40, colorVec, 'filled');
hold on;
yline(-log10(p_thresh), '--k', 'FDR = 0.05');
xline(fc_thresh, '--k', 'log2FC = 1');
xline(-fc_thresh, '--k', 'log2FC = -1');
xlabel('log_2 Fold Change');
ylabel('-log_{10} FDR');
title('Volcano Plot');
grid on; box on;

%% === Taula de resultats ===
biomarcadors = table(filteredFeatureNames(biomarcador_idx), ...
    filteredcompoundNames(biomarcador_idx), ...
    pvals(biomarcador_idx), adjP(biomarcador_idx), ...
    log2(fc(biomarcador_idx)), auc_values(biomarcador_idx), ...
    testType(biomarcador_idx), normalitat(biomarcador_idx), ...
    'VariableNames', {'Feature','Compound','pval','adjP','log2FC','AUC','Test','Normalitat'});
biomarcadors = biomarcadors(sorted_idx, :);

% Exportar resultats
writetable(biomarcadors, 'biomarcadors_significatius.xlsx');
save('biomarcadors_complets.mat', 'biomarcadors');

fprintf('Biomarcadors seleccionats (adjP < %.2f: %d\n', ...
    alpha, height(biomarcadors));

%% === Exportació de biomarcadors amb espectres ===
% Índexs absoluts respecte a la matriu original
biomarcador_idx_abs = find(keepFeatures);
biomarcador_idx_abs = biomarcador_idx_abs(biomarcador_idx);

% Informació de cada feature
ids = filteredFeatureNames(biomarcador_idx);

ri = features_merged.ri(biomarcador_idx_abs);
compounds = filteredcompoundNames(biomarcador_idx);

% Càrrega de pseudoespectres
spectramatrix_merged = load("/Users/marina/spectramatrix_merged.mat");
spectramatrix_merged = spectramatrix_merged.spectramatrix_merged;
spectramatrix_merged_filtered = spectramatrix_merged(keepFeatures,:);
spectra = spectramatrix_merged(biomarcador_idx_abs, :);

% Atributs estadístics
pval = biomarcadors.pval;
adjP = biomarcadors.adjP;
log2FC = biomarcadors.log2FC;
auc = biomarcadors.AUC;

% Taula final
biomarcadors_matriu = table(ids, ri, compounds, spectra, pval, adjP, log2FC, auc);
biomarcadors_matriu=biomarcadors_matriu(sorted_idx, :);
save('biomarcadors_complets.mat', 'biomarcadors_matriu');

%%

% === Histogrames dels top 10 biomarcadors ===
nPlot = min(10, height(biomarcadors));
selected_idx = find(biomarcador_idx);
top_features = selected_idx(sorted_idx(1:nPlot));

for i = 1:nPlot
    idx = top_features(i);
    x1 = Xtrain_scaled(idx, idxCO_train);
    x2 = Xtrain_scaled(idx, idxCP_train);
    figure;
    histogram(x1, 'Normalization','probability', 'FaceColor','b');
    hold on;
    histogram(x2, 'Normalization','probability', 'FaceColor','r');
    title(sprintf('Distribució Feature: %d', filteredFeatureNames(idx)));
    xlabel('Intensitat (autoscalada)');
    ylabel('Freqüència relativa');
    legend('CO','CP');
    grid on;
end


%%
% === BOXPLOTS dels Top 10 biomarcadors ===
nPlot = min(10, height(biomarcadors));
selected_idx = find(biomarcador_idx);
top_features = selected_idx(sorted_idx(:));
figure('Name','Boxplots dels Top 10 Biomarcadors');
for i = 1:nPlot
    idx = top_features(i);
    dataCO_i = Xtrain_scaled(idx, idxCO_train)';
    dataCP_i = Xtrain_scaled(idx, idxCP_train)';
    
    % Vector de valors
    allData = [dataCO_i; dataCP_i];
    
    % Vector de grups
    groups = [repmat({'CO'}, length(dataCO_i), 1); repmat({'CP'}, length(dataCP_i), 1)];

    % Subplot i boxplot
    subplot(2,5,i);
    boxplot(allData, groups);
    title(sprintf('Feature %d', filteredFeatureNames(idx)));
    ylabel('Intensitat (z-score)');
    grid on;
end

%% === ROC i AUC dels Top 10 ===
figure('Name','ROC dels Top 10 Biomarcadors');
for i = 1:nPlot
    idx = top_features(i);
    scores = Xtrain_scaled(idx, :)';
    labels_bin = strcmp(labels_train, 'CP');  % 0 = CO, 1 = CP
    [X,Y,~,AUC] = perfcurve(labels_bin, scores, 1);
    
    subplot(2,5,i);
    plot(X, Y, 'b', 'LineWidth', 1.5);
    hold on; plot([0 1], [0 1], 'k--');
    title(sprintf('ROC Feature %d\nAUC = %.2f', filteredFeatureNames(idx), AUC));
    xlabel('1 - Specificity');
    ylabel('Sensitivity');
    grid on;
end

%% === PCA dels biomarcadors seleccionats ===
fprintf('\n--- PCA amb biomarcadors seleccionats ---\n');

% Extracció de les dades corresponents als biomarcadors
selected_scaled = Xtrain_scaled(biomarcador_idx, :);
selected_scaled = selected_scaled(sorted_idx, :);  % ordenats per AUC

%% PCA
[coeff, score, ~, ~, expl] = pca(selected_scaled');

% Representació
figure;
gscatter(score(:,1), score(:,2), labels_train', 'br', 'xo');
xlabel(sprintf('PC1 (%.1f%%)', expl(1)));
ylabel(sprintf('PC2 (%.1f%%)', expl(2)));
title('PCA - Biomarcadors seleccionats (univariant)');
legend('CO','CP');
grid on;
box on;
sampleNamesTrain=sampleNames(trainIdx);

%% PART 2.5: Filtrat per correlació entre biomarcadors significatius
fprintf('\n--- PART 2.5: Filtrat per correlació entre biomarcadors ---\n');

% Dades dels biomarcadors significatius ordenats per AUC
selected_data = Xtrain_scaled(selected_idx(sorted_idx), :);  % [features x samples]
featureNames_selected = filteredFeatureNames(selected_idx(sorted_idx));

% Matriu de correlació (Pearson)
R = corr(selected_data');

% Representació gràfica
figure;
imagesc(R);
colorbar;
title('Correlació entre biomarcadors seleccionats');
xlabel('Biomarcadors');
ylabel('Biomarcadors');
xticks(1:length(featureNames_selected));
yticks(1:length(featureNames_selected));
xticklabels(featureNames_selected);
yticklabels(featureNames_selected);
xtickangle(90);
set(gca, 'FontSize', 6);

% Filtrat: elimina un de cada parell amb r > llindar
correlationThreshold = 0.9;
toRemove = false(size(R,1),1);
for i = 1:size(R,1)
    for j = i+1:size(R,1)
        if abs(R(i,j)) > correlationThreshold
            toRemove(j) = true;
        end
    end
end

% Missatge de resum
numBefore = size(R,1);
numAfter = sum(~toRemove);
fprintf('Shan eliminat %d de %d biomarcadors per correlació > %.2f\n', ...
    numBefore - numAfter, numBefore, correlationThreshold);

% Actualitza les variables per a la Part 3
filtered_idx_corr = selected_idx(sorted_idx(~toRemove));
filtered_names_corr = filteredFeatureNames(filtered_idx_corr);


%% PART 3: Selecció òptima de biomarcadors amb GA + SVM
fprintf('\n--- PART 3: Genetic Algorithm amb SVM ---\n');

% === Dades: només els biomarcadors seleccionats ===
X = Xtrain_scaled(filtered_idx_corr, :)';  % després de filtre per correlació
X_test=Xtest_scaled (filtered_idx_corr, :)';
y = strcmp(labels_train, 'CP');  % 0 = CO, 1 = CP

% === Paràmetres del GA ===
numFeaturesGA = size(X,2);
opts = optimoptions('ga',...
    'PopulationSize', 50,...
    'MaxGenerations', 40,...
    'EliteCount', 2,...
    'CrossoverFraction', 0.8,...
    'MutationFcn', @mutationuniform,...
    'UseParallel', false,...
    'Display', 'iter');

% === Funció objectiu: error de classificació mitjà amb SVM (k-fold CV) ===
fitnessFcn = @(chrom) meanError_CV_SVM(X, y, chrom, 5);  % 5-fold CV

% === Llançament del GA ===
fprintf('Executant GA amb SVM...\n');
[bestChromosome, bestErr] = ga(fitnessFcn, numFeaturesGA, [], [], [], [], ...
                               zeros(1,numFeaturesGA), ones(1,numFeaturesGA), [], opts);

fprintf('\nMillor subconjunt trobat (error classificació = %.3f):\n', bestErr);
disp(find(bestChromosome));

%% === Subconjunt seleccionat final ===
X_opt = X(:, logical(bestChromosome));
X_test_opt = X_test(:, logical(bestChromosome));
%% === PCA del subconjunt final ===
[~, scoreGA, ~, ~, explGA] = pca(X_opt);
figure;
gscatter(scoreGA(:,1), scoreGA(:,2), y, 'br', 'xo');
xlabel(['PC1 (' num2str(explGA(1), '%.1f') '%)']);
ylabel(['PC2 (' num2str(explGA(2), '%.1f') '%)']);
title('PCA - Subconjunt seleccionat per GA + SVM');
grid on;

%% --- Funció objectiu per al GA amb SVM ---
function error = meanError_CV_SVM(X, y, chrom, k)
    selected = logical(chrom);
    if sum(selected) < 1
        error = 1.0;  % penalització
        return;
    end
    Xsel = X(:, selected);
    cv = cvpartition(y, 'KFold', k);
    errors = zeros(k,1);
    for i = 1:k
        trainIdx = training(cv,i);
        testIdx = test(cv,i);

        Xtrain = Xsel(trainIdx,:);
        ytrain = y(trainIdx);
        Xtest  = Xsel(testIdx,:);
        ytest  = y(testIdx);

        model = fitcsvm(Xtrain, ytrain, 'KernelFunction', 'linear', 'Standardize', true);
        ypred = predict(model, Xtest);
        errors(i) = mean(ypred ~= ytest');
    end
    error = mean(errors);  % error mitjà
end



%% Model final amb features del GA

fprintf('\n--- Validació del model amb features seleccionats pel GA ---\n');

% Labels binaris
ytest_bin = strcmp(labels(testIdx), 'CP');  % test binari

% Entrenament amb les dades del train amb les features seleccionades pel GA
model_GA = fitcsvm(X_opt, y, 'KernelFunction', 'linear', 'Standardize', true);

% Predicció sobre test amb les mateixes features
ypred_GA = predict(model_GA, X_test_opt);

% Avaluació
acc_GA = mean(ypred_GA == ytest_bin');
fprintf('Accuracy sobre el test set (model GA): %.2f%%\n', acc_GA * 100);

% Matriu de confusió
cm_GA = confusionmat(ytest_bin, ypred_GA);
disp('Matriu de confusió (model GA):');
disp(array2table(cm_GA, 'VariableNames', {'Pred_CO','Pred_CP'}, ...
                          'RowNames', {'True_CO','True_CP'}));

% Càlcul AUC
[~, scores_GA] = predict(model_GA, X_test_opt);
[~,~,~,AUC_GA] = perfcurve(ytest_bin, scores_GA(:,2), 1);
fprintf('AUC (Area Under Curve) - Model GA: %.2f\n', AUC_GA);

% Corba ROC
figure;
[Xroc_GA, Yroc_GA] = perfcurve(ytest_bin, scores_GA(:,2), 1);
plot(Xroc_GA, Yroc_GA, 'b-', 'LineWidth', 1.5);
hold on; plot([0 1], [0 1], 'k--');
xlabel('1 - Specificity'); ylabel('Sensitivity');
title(sprintf('ROC - Model GA (AUC = %.2f)', AUC_GA));
grid on;
%% PART 4: Forward Feature Selection (FFS) amb SVM
fprintf('\n--- PART 4: Forward Feature Selection amb SVM ---\n');

% === Dades: només biomarcadors seleccionats (com a X i y) ===
Xffs = Xtrain_scaled(selected_idx(sorted_idx), :)';  % [samples x features]
yffs = strcmp(labels_train, 'CP');  % 0 = CO, 1 = CP
yffs=yffs';
nFeatures = size(Xffs,2);
selected = false(1, nFeatures);
remaining = true(1, nFeatures);
featureOrder = [];
accuracyHistory = [];

% === FFS Loop ===
for k = 1:nFeatures
    bestAcc = 0;
    bestIdx = -1;
    
    for i = find(remaining)
        tempSelected = selected;
        tempSelected(i) = true;
        
        % Subconjunt de features seleccionades temporalment
        X_temp = Xffs(:, tempSelected);
        
        % Avaluació mitjançant cross-validation
        acc = 1 - crossval('mcr', X_temp, yffs, ...
            'Predfun', @(xtrain,ytrain,xtest) predict(fitcsvm(xtrain,ytrain,'KernelFunction','linear','Standardize',true), xtest), ...
            'KFold', 5);
        
        if acc > bestAcc
            bestAcc = acc;
            bestIdx = i;
        end
    end
    
    if bestIdx == -1
        break;  % no hi ha millora
    end
    if ~isempty(accuracyHistory) && accuracyHistory(end) == 1 && length(accuracyHistory)>1
        fprintf('Aturat: Accuracy màxima assolit al pas %d\n', k-1);
        break;
    end
    
    selected(bestIdx) = true;
    remaining(bestIdx) = false;
    featureOrder(end+1) = bestIdx;
    accuracyHistory(end+1) = bestAcc;
    
    fprintf('Pas %d: Afegit Feature #%d (Accuracy = %.2f%%)\n', ...
        k, bestIdx, bestAcc*100);
end

% === Resultats finals ===
X_ff = Xffs(:, selected);
selectedNames = filteredFeatureNames(selected_idx(sorted_idx(selected)));
fprintf('\nTotal features seleccionades per FFS: %d\n', sum(selected));
disp('Llista de biomarcadors seleccionats (ordenats per incorporació):');
disp(selectedNames);

% === PCA per visualitzar els seleccionats amb FFS ===
[~, scoreFFS, ~, ~, explFFS] = pca(X_ff);
figure;
gscatter(scoreFFS(:,1), scoreFFS(:,2), yffs, 'br', 'xo');
xlabel(['PC1 (' num2str(explFFS(1), '%.1f') '%)']);
ylabel(['PC2 (' num2str(explFFS(2), '%.1f') '%)']);
title('PCA - Biomarcadors seleccionats per FFS');
grid on;
%%
% === PCA per visualitzar els seleccionats amb FFS ===
[~, scoreFFS, ~, ~, explFFS] = pca(X_ff);

figure;
gscatter(scoreFFS(:,1), scoreFFS(:,2), yffs, 'br', 'xo');
xlabel(['PC1 (' num2str(explFFS(1), '%.1f') '%)']);
ylabel(['PC2 (' num2str(explFFS(2), '%.1f') '%)']);
title('PCA - Biomarcadors seleccionats per FFS');
grid on;

% === Afegir els noms de mostra al costat de cada punt ===
sampleNamesTrain=sampleNames(trainIdx);
for i = 1:size(sampleNamesTrain,2)
    sampleName=sampleNamesTrain{i}(1:6);
    text(scoreFFS(i,1)+0.05, scoreFFS(i,2), sampleNamesTrain{i}(1:6), 'FontSize',8, 'Interpreter', 'none');
end

%%
%% PART 5: Entrenament final del model i test
fprintf('\n--- PART 5: Entrenament i test final del model ---\n');

% === Partició Train/Test (70% entrenament / 30% test) ===
labels_bin_all = strcmp(labels, 'CP')';  % 0 = CO, 1 = CP

% === Dades originals escalades (totes les features filtrades) ===
X_all = scaleddata(selected_idx(sorted_idx), :)';  % [samples x features]

% === Partició de les dades ===
Xtrain_all = X_all(trainIdx, :);  % dades de train amb totes les features
Xtest_all  = X_all(testIdx, :);   % dades de test amb totes les features
ytrain = labels_bin_all(trainIdx);
ytest  = labels_bin_all(testIdx);
labels_train = labels(trainIdx);  % només per a etiquetes textuals



% === Aplicar selecció de features ===
Xtrain = Xtrain_all(:, selected);
Xtest  = Xtest_all(:, selected);

fprintf('Total features seleccionades al TRAIN: %d\n', sum(selected));

% === Entrenament del model SVM ===
model = fitcsvm(Xtrain, ytrain, 'KernelFunction', 'linear', 'Standardize', true);

% === Predicció i mètriques ===
ypred = predict(model, Xtest);
acc = mean(ypred == ytest);
fprintf('Accuracy sobre el test set: %.2f%%\n', acc*100);

% === Matriu de confusió ===
cm = confusionmat(ytest, ypred);
disp('Matriu de confusió:');
disp(array2table(cm, 'VariableNames', {'Pred_CO','Pred_CP'}, ...
                      'RowNames', {'True_CO','True_CP'}));

% === Càlcul AUC ===
[~, scores] = predict(model, Xtest);
[~,~,~,AUC] = perfcurve(ytest, scores(:,2), 1);
fprintf('AUC (Area Under Curve): %.2f\n', AUC);

% === ROC Curve ===
figure;
[Xroc, Yroc] = perfcurve(ytest, scores(:,2), 1);
plot(Xroc, Yroc, 'b-', 'LineWidth', 1.5);
hold on; plot([0 1], [0 1], 'k--');
xlabel('1 - Specificity'); ylabel('Sensitivity');
title(sprintf('ROC - Model final (AUC = %.2f)', AUC));
grid on;



