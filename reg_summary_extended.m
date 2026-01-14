%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an updated version of the original codes created by the authors above, modified for this paper:
% Why is GABA related to neural distinctiveness? A computational account of age-related neural dedifferentiation.
% Quan Zhou [1] and Thad A. Polk [1]
% Department of Psychology, University of Michigan, Ann Arbor, MI, USA
% Correspondence:
% Quan Zhou (Violet) <violetz@umich.edu>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ========================================================================
%  (WIDE) — average last 600 iterations (ND_4401..ND_5000)
%  Modify pipeline to work with: case_a_final_WIDE_data_format.csv
%
%  Assumptions:
%   - CSV has columns: subj_id, GABA_level, Dvis_value, ND_2000 ... ND_5000
%   - Com_ND = mean of the last 600 iterations (4401–5000 inclusive)
%   - Then run the same regressions + optional bootstrap mediation
%   - Keep your outlier filtering logic (3-SD) but updated column names
% ========================================================================

%% ----------------------------
% 0) Read WIDE CSV
% ----------------------------
T = readtable('C:\Users\violetz\Documents\UM\Projects\AD_Model\Paper_materials\Neubiol.Aging_submit\Revision\Results\Extended_Model\Summary\case_b_final_WIDE_data_format.csv');
T.subj_id = string(T.subj_id);   % ensure string type

% delete the subjects who diagnosed as mci = mindm; old = mindo; 
% old wave 2 = mindb;

idx_mindm = startsWith(T.subj_id, "mindm", 'IgnoreCase', true);
fprintf('Removing %d subjects starting with "mindm"\n', sum(idx_mindm));

% delete the subject who has no goggle fit into the fMRI scanner mindo212

% Remove specific bad subject(s)
bad_ids = ["mindo212"];   % add more IDs if needed: ["mindo212","mindoXXX"]
idx_bad = ismember(lower(T.subj_id), lower(bad_ids));
fprintf('Removing %d subjects in bad_ids\n', sum(idx_bad));
T(idx_bad, :) = [];

%% ----------------------------
% 1) Define "last 600" iteration range and compute Com_ND
%    4401..5000 = 600 values
% ----------------------------
iter_start = 4401;
iter_end   = 5000;

nd_cols = strcat("ND_", string(iter_start:iter_end));

% Safety check: ensure all required columns exist
missing_cols = nd_cols(~ismember(nd_cols, string(T.Properties.VariableNames)));
if ~isempty(missing_cols)
    error('Missing ND columns in the CSV (showing up to 10): %s', strjoin(missing_cols(1:min(10,end)), ", "));
end

T.Com_ND = mean(T{:, nd_cols}, 2, 'omitnan');

disp('Com_ND computed as mean across ND_4401–ND_5000 (last 600 iterations).');

%% ----------------------------
% 2) Outlier checks (update variable names)
% ----------------------------
% fMRI neural distinctiveness (Y)
D_mu  = mean(T.Dvis_value, 'omitnan');
D_sig = std(T.Dvis_value,  'omitnan');
D_lower_bound = D_mu - 3*D_sig;
D_upper_bound = D_mu + 3*D_sig;

% GABA (X)
G_mu  = mean(T.GABA_level, 'omitnan');
G_sig = std(T.GABA_level,  'omitnan');
G_lower_bound = G_mu - 3*G_sig;
G_upper_bound = G_mu + 3*G_sig;

%% ----------------------------
% 3) Apply filters
%   - Keep Dvis_value > 0 and < mean+3SD (match your previous rule)
%   - Keep GABA within +/-3SD
% ----------------------------
idx_keep = (T.Dvis_value > D_lower_bound) & (T.Dvis_value < D_upper_bound);
T = T(idx_keep, :);

idx_keep = (T.GABA_level > G_lower_bound) & (T.GABA_level < G_upper_bound);
T = T(idx_keep, :);

% Remove rows where fMRI distinctiveness equals 0
idx_nonzero = T.Dvis_value ~= 0;
fprintf('Removing %d rows with Dvis_value == 0\n', sum(~idx_nonzero));
T = T(idx_nonzero, :);

%% ----------------------------
% 4) Regression models (updated names)
%   X = GABA_level
%   M = Com_ND
%   Y = Dvis_value
% ----------------------------
lm_a_w1 = fitlm(T, 'Com_ND ~ GABA_level');                 % a path
lm_b_w1 = fitlm(T, 'Dvis_value ~ Com_ND + GABA_level');    % b + c'
lm_c_w1 = fitlm(T, 'Dvis_value ~ GABA_level');             % total effect c
lm_d_w1 = fitlm(T, 'Dvis_value ~ Com_ND');                 % Y ~ M only
%% 

%% ----------------------------
% 5) Bootstrap mediation (optional)
% ----------------------------
vars = {'GABA_level','Com_ND','Dvis_value'};
Tmed = rmmissing(T(:, vars));

X = Tmed.GABA_level;
M = Tmed.Com_ND;
Y = Tmed.Dvis_value;

n = height(Tmed);
nBoot = 5000;

indirect_effects = zeros(nBoot, 1);
rng(1); % reproducibility

for b = 1:nBoot
    idx = randsample(n, n, true);

    Xb = X(idx);
    Mb = M(idx);
    Yb = Y(idx);

    % a path: M ~ X
    mdl_a = fitlm(Xb, Mb);
    a = mdl_a.Coefficients.Estimate(2);

    % b path: Y ~ M + X
    mdl_b = fitlm([Mb Xb], Yb, 'VarNames', {'M','X','Y'});
    b_path = mdl_b.Coefficients.Estimate(2); % coefficient for M

    indirect_effects(b) = a * b_path;
end

ci = prctile(indirect_effects, [2.5 97.5]);
indirect_mean = mean(indirect_effects);

fprintf('\nBootstrap indirect effect (a*b): %.6f\n', indirect_mean);
fprintf('95%% CI: [%.6f, %.6f]\n', ci(1), ci(2));

if (ci(1) > 0) || (ci(2) < 0)
    fprintf('Mediation supported (CI excludes zero).\n');
else
    fprintf('Mediation NOT supported (CI includes zero).\n');
end

%% 
% ----------------------------
% % 6) (Optional) Save a cleaned table for later regression
% % ----------------------------
% out_dir = 'C:\Users\violetz\Documents\UM\Projects\AD_Model\Paper_materials\Neubiol.Aging_submit\Revision\Results\Extended_Model\Case_A\Final';
% if ~exist(out_dir, 'dir'); mkdir(out_dir); end
% 
% out_csv = fullfile(out_dir, 'case_a_cleaned_with_ComND_last600.csv');
% writetable(T, out_csv);
% fprintf('\nCleaned table written to:\n%s\n', out_csv);

