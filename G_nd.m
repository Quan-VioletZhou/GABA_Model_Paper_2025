

T = readtable('CombinedMasterGLX_new.csv');
% ND_RSA = readtable('distinctiveness_measures_combined_all_mindo_0520_2.csv');
% ND_RSA.Subject = erase(ND_RSA.Subject, 'p');

T_filtered = T(T.AgeCategory_Young1 == 2 & T.Subgroup == 1,:);
% T_filtered = T(T.AgeCategory_Young1 == 2 & T.Subgroup == 1 & T.Wave_Num == 1, :);

% idx_wave_1 = find(T_filtered.Wave_Num == 1);

% cd('D:\GABA-AD\results_G_ATC_a.2_S_R_only03_30_15.30\')
% cd('D:\GABA-AD\results_b.203_03_11.41_ALL');
% cd('D:\GABA-AD\results_b.103_03_17.14_ALL');
% cd('D:\GABA-AD\results_G_UNC_b.103_07_16.36')
% cd('D:\GABA_AD\CCN_conference\results_case_c_stimulation_final04_09_16.19');
% cd('D:\GABA-AD\RVV\results_case_c_stimulation_final05_20_09.09');

% Paper Results:
% cd('D:\GABA-AD\Paper_results\Final\results_case_a_stimulation_final05_29_13.22');
% cd('D:\GABA-AD\Paper_results\Final\results_case_b_stimulation_final05_29_18.09');
cd('D:\GABA-AD\Paper_results\Final\results_case_c_stimulation_final04_09_16.19')


% Using RVV results
% cd('D:\GABA-AD\RVV\results_case_c_stimulation_final05_19_13.56')


ND_all = readtable('cosine_similarity_matrix.csv');
Com_ND = mean(ND_all{:, 1700:2000}, 2);

T_filtered.Com_ND = Com_ND;

T_wave_1 = T_filtered(T_filtered.Wave_Num == 1, :);
% combined_table = innerjoin(ND_RSA, T_wave_1, 'Keys','Subject');
%% Now use the ND_RSA measure to run regression:
ALL_avg_within = mean([combined_table.All_Within_Face, combined_table.All_Within_House], 2);

ALL_RSA_Distinctiveness = ALL_avg_within ./ combined_table.All_Between;
combined_table.ALL_RSA_Distinctiveness = ALL_RSA_Distinctiveness;

combined_table_clean = combined_table(combined_table.D_vis_sm_sm2000_P >= 0, :);

lm_a_w1_rsa = fitlm(combined_table_clean, 'Com_ND ~ G_ATC_RVV');
lm_b_w1_rsa = fitlm(combined_table_clean, 'All_Overall_Distinctiveness ~ Com_ND + G_ATC_RVV');
lm_c_w1_rsa = fitlm(combined_table_clean, 'All_Overall_Distinctiveness ~ G_ATC_RVV');
lm_d_w1_rsa = fitlm(combined_table_clean, 'All_Overall_Distinctiveness ~ Com_ND');


%% run the model with account for individual slope and intercept

T_filtered_clean = T_filtered(T_filtered.D_vis_sm_sm2000_P >= 0, :);

% T_filtered_clean = T_filtered;

lme_a_1 = fitlme(T_filtered_clean, 'Com_ND ~ G_ATC_VV + (1 + G_ATC_VV|SubNum)');
lme_b_1 = fitlme(T_filtered_clean, 'D_vis_sm_sm2000_P ~ Com_ND + G_ATC_VV + (Com_ND + G_ATC_VV|SubNum)');
lme_c_1 = fitlme(T_filtered_clean, 'D_vis_sm_sm2000_P ~ G_ATC_VV + (1 + G_ATC_VV|SubNum)');

lme_d_1 = fitlme(T_filtered_clean, ...
    'D_vis_sm_sm2000_P ~ Com_ND + (1 + Com_ND | SubNum)');

lme_va_g = fitlme(T_filtered_clean, 'NIH_VA_USS ~ G_ATC_VV + (1 + G_ATC_VV|SubNum)');
%% 

Com_ND(Com_ND == 0) = NaN;
% VA(VA<30) = NaN;

T_filtered.Com_ND = Com_ND;
% T_filtered.NIH_VA_USS_clean = VA;

% lme_va_ND = fitlme(T_filtered, 'NIH_VA_USS_clean ~ D_vis_sm_sm2000_P + (1 + D_vis_sm_sm2000_P|SubNum)');
% lme_va_g = fitlme(T_filtered_clean, 'NIH_VA_USS_clean ~ G_ATC_VV + (1 + G_ATC_VV|SubNum)');


% a_coef = lme_a.Coefficients{'G_ATC_VV','Estimate'};
% b_coef = lme_b.Coefficients{'Com_ND','Estimate'};
% indirect_effect = a_coef * b_coef;

%% filter only the wave 1 subject data:
T_filtered_wave_1 = T_filtered(T_filtered.Wave_Num == 1, :);

T_wave_1_cleaned = T_wave_1(T_wave_1.D_vis_sm_sm2000_P > 0, :);
T_wave_1_cleaned.Com_ND(T_wave_1_cleaned.Com_ND == 0) = NaN;

% T_wave_1_cleaned(134,:) = []; 


lm_a_w1 = fitlm(T_wave_1_cleaned, 'Com_ND ~ G_ATC_VV');
lm_b_w1 = fitlm(T_wave_1_cleaned, 'D_vis_sm_sm2000_P ~ Com_ND + G_ATC_VV');
lm_c_w1 = fitlm(T_wave_1_cleaned, 'D_vis_sm_sm2000_P ~ G_ATC_VV');
lm_d_w1 = fitlm(T_wave_1_cleaned, 'D_vis_sm_sm2000_P ~ Com_ND');

%% Bootstrapping for the wave 1 only:
% Filter to Wave 1 only
% Extract variables

GABA = T_wave_1_cleaned.G_ATC_VV;
Com_ND = T_wave_1_cleaned.Com_ND;
Dediff = T_wave_1_cleaned.D_vis_sm_sm2000_P;

% Number of bootstrap samples
nBoot = 5000;

indirect_effects = zeros(nBoot, 1);

n = height(T_wave_1_cleaned);

rng(42); % for reproducibility
for i = 1:nBoot
    idx = randsample(n, n, true);
    
    % Resample data
    G_sample = GABA(idx);
    C_sample = Com_ND(idx);
    D_sample = Dediff(idx);
    
    % Path a
    mdl_a = fitlm(G_sample, C_sample);
    a_coef = mdl_a.Coefficients.Estimate(2); % slope for path a

    % Path b
    mdl_b = fitlm([C_sample, G_sample], D_sample);
    b_coef = mdl_b.Coefficients.Estimate(2); % slope for path b

    % Indirect effect
    indirect_effects(i) = a_coef * b_coef;
end

% Compute stats
indirect_mean = mean(indirect_effects);
ci = prctile(indirect_effects, [2.5 97.5]);

% Output results
fprintf('\nBOOTSTRAPPED INDIRECT EFFECT\n');
fprintf('Mean indirect effect: %.4f\n', indirect_mean);
fprintf('95%% CI: [%.4f, %.4f]\n', ci(1), ci(2));

% Optional: Plot distribution
figure;
histogram(indirect_effects, 50, 'FaceColor', [0.4 0.6 1]);
xlabel('Indirect Effect (a × b)');
ylabel('Frequency');
title('Bootstrap Distribution of Indirect Effect');

%% Plot the distribution of the indirect effect
% Compute CI and mean
% Bootstrap distribution plot with 95% CI
figure;
histogram(indirect_effects, 50, 'FaceColor', [0.4 0.6 1]);
xlabel('Indirect Effect (a × b)');
ylabel('Frequency');
title('Bootstrap Distribution of Indirect Effect');

% Compute CI
ci = prctile(indirect_effects, [2.5, 97.5]);
mean_indirect = mean(indirect_effects);

% Plot CI lines
hold on;
yLimits = ylim;
plot([ci(1) ci(1)], yLimits, 'r--', 'LineWidth', 2, 'DisplayName', '95% CI Lower');
plot([ci(2) ci(2)], yLimits, 'r--', 'LineWidth', 2, 'DisplayName', '95% CI Upper');

% Optional: plot mean
plot([mean_indirect mean_indirect], yLimits, 'k-', 'LineWidth', 2, 'DisplayName', 'Mean');
% Add text labels next to each CI line with 3 decimal places
text(ci(1) + 0.01, yLimits(2)*0.9, sprintf('%.3f', ci(1)), ...
    'Color', 'red', 'FontSize', 12, 'HorizontalAlignment', 'left');

text(ci(2) + 0.01, yLimits(2)*0.9, sprintf('%.3f', ci(2)), ...
    'Color', 'red', 'FontSize', 12, 'HorizontalAlignment', 'left');

text(mean_indirect + 0.01, yLimits(2)*0.85, sprintf('Mean = %.3f', mean_indirect), ...
    'Color', 'black', 'FontSize', 12, 'HorizontalAlignment', 'left');
legend('show', 'Location', 'northwest');
box on;

%% filter the only wave one data
% T_filtered_clean = T_filtered;

% T_filtered_clean = T_filtered(T_filtered.Wave_Num == 1, :);

G_VV = T_filtered_clean.G_ATC_VV;
ND = T_filtered_clean.D_vis_sm_sm2000_P;

% using the 3 std to filter the Neural distinctiveness:
ND(ND < 0) = NaN;

VA = T_filtered_clean.NIH_VA_USS;
% VA(VA <= 30.0) = NaN;
FIN_Gaus = T_filtered_clean.FIN_GausLvl;

Com_ND = T_filtered_clean.Com_ND;


% change the timepoints.`
% Com_ND = mean(ND_all{:, 1700:2000}, 2);
% Com_ND(Com_ND == 0) = NaN;
% Com_ND_wave_1 = Com_ND(idx_wave_1, :);

% Define independent (x) and dependent (y) variables
x = Com_ND;
y = ND; 

% Fit linear model
lme_va = fitlme(T_filtered_clean, 'NIH_VA_USS ~ G_ATC_VV + (1 + G_ATC_VV|SubNum)');
disp(lme_va);
p = lme_va.Coefficients.pValue(2);

figure;

% 158/255 201/255 226/255;  % #9EC9E2
blueColor_1 = [158, 201, 226] / 255;  % #3C93C2
blueColor_2 = [60, 147, 194] / 255;  % #3C93C2

purpleColor_1 = [242, 172, 202] / 255;
purpleColor_2 = [233, 86, 148] / 255;
purpleColor_3 = [143, 0, 59] / 255;

greenColor_1 = [64, 173, 90] / 255;
greenColor_2 = [6, 89, 42] / 255;
greenColor_3 = [156, 206, 167] / 255;



% Plot data points
scatter(T_filtered_clean.G_ATC_VV, T_filtered_clean.NIH_VA_USS , 50, greenColor_3, 'filled'); % Blue dots for data
hold on;

% Plot regression line
% Step 1: Create x-values to predict over
x_vals = linspace(min(T_filtered_clean.G_ATC_VV), max(T_filtered_clean.G_ATC_VV), 100)';

% Step 2: Choose an existing subject from your data (check that this ID exists!)
subject_id = T_filtered_clean.SubNum(1); % or use a specific subject like 101

% Step 3: Create prediction table
X_predict = table(x_vals, repmat(subject_id, size(x_vals)), ...
    'VariableNames', {'G_ATC_VV', 'SubNum'});

% Step 4: Predict
y_fit = predict(lme_va, X_predict);

plot(x_vals, y_fit, 'LineWidth', 2);  % Fit line
xlabel('G\_ATC\_VV');
ylabel('Predicted NIH\_VA\_USS');
title('Model Prediction for Subject');
grid on;

% p value for the plot
if p < 0.001
    p_text = 'p < 0.001';
else
    p_text = sprintf('p = %.3f', p);
end

x_pos = min(x) + 0.05 * range(x);
y_pos = max(y) - 0.05 * range(y);

text(x_pos, y_pos, p_text, 'FontSize', 12, 'Color', greenColor_1); % 'k' is black


% Labels and title
xlabel('Visual Acuity');
ylabel('GABA+ level');
title('GABA level x Visual Acuity');
legend('Data', 'Linear Fit', 'Location', 'best');

grid on;
hold off;

%% Next part is the mediation analysis: if com_ND mediate the relationship between G_ATC_VV and ND
% Step 1: Path a: Does GABA predict com_mod distinctiveness
mdl_a = fitlm(G_VV, Com_ND);
a = mdl_a.Coefficients.Estimate(2);
p_a = mdl_a.Coefficients.pValue(2);

% Step 2: Path b & c’: Does com_mod distinctiveness predict cortical distinctiveness
mdl_b = fitlm([G_VV Com_ND], ND);
b = mdl_b.Coefficients.Estimate(3);
c_prime = mdl_b.Coefficients.Estimate(2);
p_b = mdl_b.Coefficients.pValue(3);
p_c_prime = mdl_b.Coefficients.pValue(2);

%% Compute the Indirect Effect (Bootstrapping)
G_VV = T_filtered_clean.G_ATC_VV;
ND = T_filtered_clean.D_vis_sm_sm2000_P;
Com_ND = T_filtered_clean.Com_ND;

num_boot = 5000;
boot_indirect = zeros(num_boot,1);

for i = 1:num_boot
    idx = randi(height(G_VV), height(G_VV), 1); % Resample with replacement
    mdl_a_boot = fitlm(G_VV(idx), Com_ND(idx));
    mdl_b_boot = fitlm([G_VV(idx) Com_ND(idx)], ND(idx));
    
    boot_indirect(i) = mdl_a_boot.Coefficients.Estimate(2) * mdl_b_boot.Coefficients.Estimate(3);
end

% Confidence Interval
ci = prctile(boot_indirect, [2.5, 97.5]);
disp(['Indirect Effect: ', num2str(mean(boot_indirect)), ' (95% CI: ', num2str(ci(1)), ', ', num2str(ci(2)), ')']);

%% demographic: 
% Assuming T is your main table and already loaded in the workspace

% Gender count (1 = Male, 0 = Female)
num_male = sum(T_filtered.Sex_M1 == 1);
num_female = sum(T_filtered.Sex_M1 == 2);

% Average Education
age_data = T_filtered.("Age_DuringParticipation_");
avg_age = mean(age_data, 'omitnan');
std_age = std(age_data, 'omitnan');
median_age = median(age_data, 'omitnan');
min_age = min(age_data, [], 'omitnan');
max_age = max(age_data, [], 'omitnan');
avg_education = mean(T_filtered.Education, 'omitnan');
std_education = std(T_filtered.Education, 'omitnan');

edu_data = T_filtered.Education;
avg_edu = mean(edu_data, 'omitnan');
std_edu = std(edu_data, 'omitnan');
median_edu = median(edu_data, 'omitnan');
min_edu = min(edu_data, [], 'omitnan');
max_edu = max(edu_data, [], 'omitnan');

% Race summary: 1 = White, others = Other
race_white = sum(T_filtered.Race == 1);
race_other = sum(T_filtered.Race ~= 1);

%count wave 1 and wave 2:
num_wave1 = sum(T_filtered.Wave_Num == 1);
num_wave2 = sum(T_filtered.Wave_Num == 2);

%% 
% Get Neural Distinctiveness values
ND_nofilter = T_filtered.G_ATC_VV;

% Remove NaNs
ND_clean = ND_nofilter(~isnan(ND_nofilter));

% Compute mean and standard deviation
mu = mean(ND_clean);
sigma = std(ND_clean);

% Plot histogram
figure; 
histogram(ND_clean, 20);  % Adjust bins as needed
xlabel('Neural Distinctiveness');
ylabel('Frequency');
title('Distribution of Neural Distinctiveness');

% Overlay vertical lines
hold on;
xline(mu, '--k', 'Mean', 'LabelOrientation', 'horizontal', 'FontSize', 10);
xline(mu + 3*sigma, '--r', '+3 SD', 'LabelOrientation', 'horizontal', 'FontSize', 10);
xline(mu - 3*sigma, '--r', '-3 SD', 'LabelOrientation', 'horizontal', 'FontSize', 10);


%% % Assume T_filtered_clean has SubNum, G_ATC_VV, Com_ND
% Step 1: Plot individual lines
figure; hold on;
subjects = unique(T_filtered_clean.SubNum);
for i = 1:length(subjects)
    subData = T_filtered_clean(T_filtered_clean.SubNum == subjects(i), :);
    plot(subData.G_ATC_VV, subData.Com_ND, '-o', 'Color', [0.7 0.7 0.7]); % Light gray
end

% Step 2: Plot fixed effect line
x_vals = linspace(min(T_filtered_clean.G_ATC_VV), max(T_filtered_clean.G_ATC_VV), 100);
fixed_intercept = lme_a.Coefficients.Estimate(1);
fixed_slope     = lme_a.Coefficients.Estimate(2);
y_vals = fixed_intercept + fixed_slope * x_vals;

plot(x_vals, y_vals, 'r-', 'LineWidth', 2); % Red fixed effect line

xlabel('GABA (G_ATC_VV)');
ylabel('Common ND');
title('Individual Slopes + Fixed Effect');


%% Bootstrapping clusters:
nBoot = 1000;
indirect_effects = nan(nBoot, 1);
fail_count = 0;

subjects = unique(T_filtered_clean.SubNum);
nSubj = length(subjects);

rng(42);  % Reproducibility

for b = 1:nBoot
    sampled_subs = datasample(subjects, nSubj, 'Replace', true);
    T_boot = [];

    for s = 1:nSubj
        subj_data = T_filtered_clean(T_filtered_clean.SubNum == sampled_subs(s), :);
        subj_data.SubNum(:) = s;  % Reset subject index
        T_boot = [T_boot; subj_data];
    end

    success = false;

    % First try LME
    try
        lme_a = fitlme(T_boot, 'Com_ND ~ G_ATC_VV + (1 | SubNum)', 'FitMethod', 'ML');
        lme_b = fitlme(T_boot, 'D_vis_sm_sm2000_P ~ Com_ND + G_ATC_VV + (1 | SubNum)', 'FitMethod', 'ML');

        a = lme_a.Coefficients.Estimate('G_ATC_VV');
        b_coef = lme_b.Coefficients.Estimate('Com_ND');

        indirect_effects(b) = a * b_coef;
        success = true;

    catch
        warning("Bootstrap %d: LME failed. Trying fitlm fallback...", b);
    end

    % If LME failed, try fallback fitlm
    if ~success
        try
            mdl_a = fitlm(T_boot, 'Com_ND ~ G_ATC_VV');
            mdl_b = fitlm(T_boot, 'D_vis_sm_sm2000_P ~ Com_ND + G_ATC_VV');

            a = mdl_a.Coefficients.Estimate(2);  % G_ATC_VV
            b_coef = mdl_b.Coefficients.Estimate(2);  % Com_ND

            indirect_effects(b) = a * b_coef;
        catch
            fail_count = fail_count + 1;
            warning("Bootstrap %d: fallback also failed.", b);
        end
    end
end

% Clean and summarize
indirect_effects = indirect_effects(~isnan(indirect_effects));
ci_low = prctile(indirect_effects, 2.5);
ci_high = prctile(indirect_effects, 97.5);
pval_indirect = mean(indirect_effects < 0) * 2;  % two-tailed

fprintf('\nIndirect effect a*b: %.4f\n', mean(indirect_effects));
fprintf('95%% CI: [%.4f, %.4f]\n', ci_low, ci_high);
fprintf('Bootstrap p-value: %.4f\n', pval_indirect);
fprintf('%d out of %d bootstraps failed completely.\n', fail_count, nBoot);

%% % Assuming your table is named T_filtered_wave_1

% Find rows where visual distinctiveness is negative
idx_negative_vis = T_filtered_wave_1.D_vis_sm_sm2000_P < 0;

% Extract SubNum values for those rows
subnums_with_negative_vis = unique(T_filtered_wave_1.SubNum(idx_negative_vis));

% Display the list
disp('Subjects with visual distinctiveness < 0:');
disp(subnums_with_negative_vis);



%% Look at the corr between 2 waves:
% Step 1: Identify subjects who completed both W1 and W2
subs_w1 = unique(T_filtered_clean(T_filtered_clean.Wave_Num == 1, :).SubNum);
subs_w2 = unique(T_filtered_clean(T_filtered_clean.Wave_Num == 2, :).SubNum);
subs_both = intersect(subs_w1, subs_w2);

% Step 2: Keep W2 data from those subjects
T_wave2_both = T_filtered_clean(T_filtered_clean.Wave_Num == 2 & ismember(T_filtered_clean.SubNum, subs_both), :);

% Step 3: Keep W1 data from subjects who did NOT complete W2
subs_w1_only = setdiff(subs_w1, subs_both);
T_wave1_only = T_filtered_clean(T_filtered_clean.Wave_Num == 1 & ismember(T_filtered_clean.SubNum, subs_w1_only), :);

% Step 4: Combine them
T_combined = [T_wave2_both; T_wave1_only];

%% Grab the W1 and W2 data:

% Count the number of occurrences of each SubNum
[unique_ids, ~, id_idx] = unique(T_filtered.SubNum);
id_counts = histcounts(id_idx, 1:(numel(unique_ids)+1));

% Find SubNums that appear more than once
repeated_ids = unique_ids(id_counts > 1);

% Grab all rows with those SubNums
rows_with_duplicates = ismember(T_filtered.SubNum, repeated_ids);
T_bothwaves = T_filtered(rows_with_duplicates, :);

% Extract Wave 1 and Wave 2 data
wave1 = T_bothwaves(T_bothwaves.Wave_Num == 1, :);
wave2 = T_bothwaves(T_bothwaves.Wave_Num == 2, :);

% Find common subjects
common_ids = intersect(wave1.SubNum, wave2.SubNum);
wave1_common = wave1(ismember(wave1.SubNum, common_ids), :);
wave2_common = wave2(ismember(wave2.SubNum, common_ids), :);

% Rename variables for clarity
wave1_common = renamevars(wave1_common, ...
    {'G_ATC_VV', 'Com_ND', 'D_vis_sm_sm2000_P'}, ...
    {'G_ATC_VV_W1', 'Com_ND_W1', 'D_vis_W1'});

wave2_common = renamevars(wave2_common, ...
    {'G_ATC_VV', 'Com_ND', 'D_vis_sm_sm2000_P'}, ...
    {'G_ATC_VV_W2', 'Com_ND_W2', 'D_vis_W2'});

% Join the tables on SubNum
T_joined = innerjoin(wave1_common(:, {'SubNum', 'G_ATC_VV_W1', 'Com_ND_W1', 'D_vis_W1'}), ...
                     wave2_common(:, {'SubNum', 'G_ATC_VV_W2', 'Com_ND_W2', 'D_vis_W2'}), ...
                     'Keys', 'SubNum');

% Correlations with NaNs removed
[r1, p1] = corr(T_joined.G_ATC_VV_W1, T_joined.G_ATC_VV_W2, 'Rows', 'complete');
[r2, p2] = corr(T_joined.Com_ND_W1, T_joined.Com_ND_W2, 'Rows', 'complete');
[r3, p3] = corr(T_joined.D_vis_W1, T_joined.D_vis_W2, 'Rows', 'complete');

% Display results
fprintf('Correlation of G_ATC_VV: r = %.3f, p = %.4f\n', r1, p1);
fprintf('Correlation of Com_ND:   r = %.3f, p = %.4f\n', r2, p2);
fprintf('Correlation of D_vis:    r = %.3f, p = %.4f\n', r3, p3);


%% 
T_filtered_clean = T_filtered;

T_filtered_clean = rmmissing(T_filtered_clean, 'DataVariables', {'FAold_Cog1_Speed', 'PHQ_9', 'SubNum'});

T_filtered_clean.FAold_Cog1_Speed = str2double(T_filtered_clean.FAold_Cog1_Speed);
T_filtered_clean.PHQ_9 = str2double(T_filtered_clean.PHQ_9);

% Step 2: Convert the cleaned cell array of strings to double
T_filtered_clean.FAold_Cog1_Speed = str2double(cleaned_vals);

lme_a_1 = fitlme(T_filtered_clean, 'FAold_Cog1_Speed ~ PHQ_9 + (1 + PHQ_9|SubNum)');









