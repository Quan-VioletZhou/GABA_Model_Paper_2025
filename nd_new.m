%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is to calculate the cosine similarity depending on individual GABA level for this paper:
% Why is GABA related to neural distinctiveness? A computational account of age-related neural dedifferentiation.
% Quan Zhou [1] and Thad A. Polk [1]
% Department of Psychology, University of Michigan, Ann Arbor, MI, USA
% Correspondence:
% Quan Zhou (Violet) <violetz@umich.edu>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0, 'DefaultFigureVisible', 'off'); 

% Get the current date and time
timestamp = datestr(now, 'mm_dd_HH.MM'); % Format: month_day_hour.minute

base_dir = 'D:\GABA-AD\Paper_results\Final';

% Define the directory name
save_dir = fullfile(base_dir, ['results_case_b_stimulation_final' timestamp]);

% Check if the directory exists, if not, create it
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% Load empirical GABA data
% T = readtable('mindo_gaba_vis.csv'); 
T = readtable('CombinedMasterGLX_new.csv');

T_age = T(T.AgeCategory_Young1 == 2 & T.Subgroup == 1, :);
T_filtered = T_age(:, {'Subject', 'G_ATC_VV', 'G_ATC_RVV', 'D_vis_sm_sm2000_P'});

GABA_values = T_filtered.G_ATC_RVV;  % Extract GABA values for all individuals

% num_individuals = 3;
num_individuals = length(T_filtered.Subject);
iter = 2000;  % Set number of iterations
num_repeats = 3;  % Run each subject 3 times and average results

% Normalize GABA level between 0 (low GABA) and 1 (high GABA)
GABA_min = min(GABA_values, [], 'omitnan');
GABA_max = max(GABA_values, [], 'omitnan');

% GABA_values_sim = [GABA_min:0.1:GABA_max];

% Initialize storage for cosine similarity results
cos_avg_last600 = zeros(num_individuals, 1);  
cos_avg_1400_1450 = zeros(num_individuals, 1);
cos_avg_1450_1500 = zeros(num_individuals, 1);
cos_avg_1400_1600 = zeros(num_individuals, 1);  
cos_avg_1800_2000 = zeros(num_individuals, 1);

% Initialize the full cos_over_time matrix
cos_over_time = zeros(num_individuals, iter); 
% Intialize the full matrix to save each subjects for each repeats
cos_sim_all_repeats = zeros(num_individuals, num_repeats, iter);

for subj = 1:num_individuals
    if isnan(GABA_values(subj))
        disp(['Skipping Subject ' num2str(subj) ' due to NaN GABA level.']);
        continue;  % Skip this subject if GABA level is NaN
    end

    % Initialize arrays for storing results across repeats
    temp_last600 = zeros(num_repeats, 1);
    temp_1400_1500 = zeros(num_repeats, 1);
    temp_1400_1450 = zeros(num_repeats, 1);
    temp_1800_2000 = zeros(num_repeats, 1);
    
    % Initialize temporary storage for cos_over_time
    temp_cos_over_time = zeros(num_repeats, iter);

    for repeat = 1:num_repeats
        
        
        retrain = true;  % Initialize retraining flag
        while retrain  % Keep retrying until cosine similarity is acceptable
            disp(['Processing Subject ' num2str(subj) ', Run ' num2str(repeat)]);
       
        % Set individual's GABA value --< CHange THIS!!
        GABA_level = GABA_values(subj);
        % GABA_level = GABA_max;
        % GABA_scaled = 0.3;

        GABA_scaled = (GABA_level - GABA_min) / (GABA_max - GABA_min); 
        
        disp(['Processing Subject ' num2str(subj) ', Run ' num2str(repeat) ' with GABA: ' num2str(GABA_level) 'with GABA_scaled' num2str(GABA_scaled)]);
        % Run main script for this subject
        main_d1_gaba;  

        % Reset per run
        cos_sim_over_time = zeros(1, iter);  

        for t = 1:iter
            num_patterns = size(ss_activity_per_number, 2);
            cos_sim_matrix = zeros(num_patterns, num_patterns);

            for i = 1:num_patterns
                for j = 1:num_patterns
                    if norm(ss_activity_per_number(:, i, t)) > 0 && norm(ss_activity_per_number(:, j, t)) > 0
                        cos_sim_matrix(i, j) = dot(ss_activity_per_number(:, i, t), ss_activity_per_number(:, j, t)) / ...
                                               (norm(ss_activity_per_number(:, i, t)) * norm(ss_activity_per_number(:, j, t)));
                    else
                        cos_sim_matrix(i, j) = 0;
                    end
                end
            end

            cos_sim_over_time(t) = mean(cos_sim_matrix(:)); % Store cosine similarity per iteration
        end
        
        cos_trial_700_990 = mean(cos_sim_over_time(700:990)); 
        % Check if retraining is needed

        if cos_trial_700_990 > 0.11
            disp(['Retraining Subject ' num2str(subj) ' (cos_trial_700_990 too high: ' num2str(cos_trial_700_990) ')']);
            retrain = true;  % Repeat the simulation
        else
            retrain = false;  % Accept the trial and store results
        
        cos_trial_1400_1500 = mean(cos_sim_over_time(1400:1500));

        % **Print Cosine Similarity for 1400-1600**
        disp(['   cos_sim_1400_1500 = ' num2str(cos_trial_1400_1500)]);

        % Store cosine averages for different time windows
        temp_last600(repeat) = mean(cos_sim_over_time(iter-599:iter));  
        temp_1400_1450(repeat) = mean(cos_sim_over_time(1400:1450));  
        temp_1450_1500(repeat) = mean(cos_sim_over_time(1400:1500));  
        temp_1400_1600(repeat) = mean(cos_sim_over_time(1400:1600));
        temp_1800_2000(repeat) = mean(cos_sim_over_time(1800:2000));  

        % Store full cosine similarity over time for this repeat
        temp_cos_over_time(repeat, :) = cos_sim_over_time;

        cos_sim_all_repeats(subj, repeat, :) = cos_sim_over_time;
        end
        end
end



    % Average across 5 runs
    cos_avg_last600(subj) = mean(temp_last600);
    cos_avg_1400_1450(subj) = mean(temp_1400_1450);
    cos_avg_1450_1500(subj) = mean(temp_1450_1500);
    cos_avg_1400_1600(subj) = mean(temp_1400_1600);
    cos_avg_1800_2000(subj) = mean(temp_1800_2000);
    
    % Average cos_over_time across 5 runs
    cos_over_time(subj, :) = mean(temp_cos_over_time, 1);

end

% Create table with subject IDs and computed cosine averages
T_results = table((1:num_individuals)', cos_avg_last600, cos_avg_1400_1450, cos_avg_1450_1500, cos_avg_1400_1600, cos_avg_1800_2000, ...
    'VariableNames', {'Subject', 'cos_avg_last600', 'cos_avg_1400_1450', 'cos_avg_1450_1500', 'cos_avg_1400_1600','cos_avg_1800_2000'});

%Save cosine similarity table
writetable(T_results, fullfile(save_dir, 'cosine_similarity_results.csv'));

%Save cos_over_time as a CSV file
writematrix(cos_over_time, fullfile(save_dir, 'cosine_similarity_matrix.csv'));

% Display message
disp(['Processing Completed! Results saved in ', save_dir]);


% To save the cosine similarity matrix for each subject and each repeat:
% Reshape: subjects × (repeats * time)
all_rows = [];

for subj = 1:num_individuals
    for rep = 1:num_repeats
        % Extract 1x2000 vector
        cos_row = squeeze(cos_sim_all_repeats(subj, rep, :))';
        
        % Create a row with subject, repeat, and the 2000 values
        row = [{sprintf('sub_%d', subj)}, rep, num2cell(cos_row)];
        
        % Append to overall storage
        all_rows = [all_rows; row];
    end
end

% Create variable names
var_names = [{'Subject', 'Repeat'}, arrayfun(@(x) sprintf('Iter%d', x), 1:iter, 'UniformOutput', false)];

% Convert to table
T_cos_long = cell2table(all_rows, 'VariableNames', var_names);

% Save to CSV
writetable(T_cos_long, fullfile(save_dir, 'cosine_similarity_all_repeats_longformat.csv'));

