%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Special Issue: Cortical Circuitry and Synaptic Dysfunctions in Alzheimer's Disease and Other Dementias
% Alzheimer's disease as a result of stimulus reduction in a GABA-A deficient brain: a neurocomputational model.
% Mariana A. Aguiar-Furucho [1,3] Francisco J. R. Peláez [2,3] 
% [1] Engineering, Neuroscience and Bio-Inspired Systems Study Group (GENeSis), Department of Electrotechnics (DAELT), 
% Universidade Tecnológica Federal do Paraná (UTFPR), Paraná, 80230-901, Brazil. 
% [2] Center of Mathematics, Computation, and Cognition (CMCC). Universidade Federal do ABC.
% [3] Center for Neuroscience and Behavior, Institute of Psychology, University of São Paulo, São Paulo, Brazil. 
% Correspondence: 
% Mariana A. Aguiar-Furucho <marianafurucho@utfpr.edu.br>
% Francisco J. R. Peláez <francisco.pelaez@ufabc.edu.br>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an updated version of the original codes created by the authors above, modified for this paper:
% Why is GABA related to neural distinctiveness? A computational account of age-related neural dedifferentiation.
% Quan Zhou [1] and Thad A. Polk [1]
% Department of Psychology, University of Michigan, Ann Arbor, MI, USA
% Correspondence:
% Quan Zhou (Violet) <violetz@umich.edu>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% clear all;
% close all;
% clc;
iter       = 5000;
pause(0.1);

% [P,col,lin]=inputnumbers_old;
% [P,col,lin]=inputnumbers_positive;
[P,col,lin]=inputnumbers_noise_20(0.1);
% Using the permutation input within-category for each model_id
% [P, col, lin] = inputnumbers_balanced;
% [P, col, lin] = inputnumbers_clean_cat(0.2);
[n_entradas, padroes] = size(P);

% Define the network structure
sensory=n_entradas;
tc=n_entradas;
% spiny=padroes;
spiny =20;
% inhibit=padroes;
inhibit = spiny;
numero_total_de_neuronas=sensory+tc+spiny+inhibit;

% initialize network properties
steepness = 40.*ones(numero_total_de_neuronas,1);
steepness(sensory+1:sensory+tc) = 40;

% generate the connections
camadas = 4;
% We add the random seed to ensure the reproducibility
[w, mascara] = gennet_con_4_capas(sensory,tc,spiny,inhibit);
n_neuronios = size(w, 1);

% Initialize learning parameter, weights, and shift
incw = zeros(size(w)); % store weights update
shift = 0.061.* ones(n_neuronios, 1); % initialize the shift parameter
shift_inicial = shift; 
output = zeros(n_neuronios, 1);
output_antes = output;

% define learning rules
% fator_aprendiz = 0.0019; % original learning rate for weights
fator_aprendiz = 0.0015;
% fator_aprendiz = 0.0038;
velocidade_deslocamento = 0.0199.*(ones(size(shift))); % how fast neurons shift their activation thresh
velocidade_deslocamento(sensory+1:sensory+tc) = 0.0199; % fixed thresh shift speed for tc

% set iterations and graphing parameters
% iter = 2000;
iter_graph = 1;
inter_totais = 1;
graf_shift = zeros(1,iter);
graf_output = zeros(1,iter);

% Parameters for GABA function:
alpha = 5;
threshold = 0.8;

% inputs for STM and pruning
input_reduction = 2; % strength of input is reduced
PX=[1 0 1  1 0 1  0 1 0  1 0 1  1 0 1]; 
PX = PX/sum(PX);

P_STM=P;
P_STM(:,randi(randi([1 10])))=PX;
P_STM_SR = P_STM/input_reduction;

P0=P;
P1=P;
P_SR=P/input_reduction;
[conexoes75, conexoes100, decaimentoConexoes] = pruningParametros(w,sensory,spiny,tc);
time_init=datetime;

% Initiate an input record:
recorded_inputs = zeros(size(P,1), size(P,2), iter);

disp('Processing Koniocortex .....')

%% VZ added %% Initialize a matrix to store SS neuron activations at each iteration
%  Initialize SS activation storage (SS neurons x Numbers x Epochs)
ss_activity_per_number = zeros(spiny, padroes, iter);
graf_incw = zeros(1, iter);

for i = 1:iter
% 
%     Determine if we are at a visualization checkpoint
    % if ismember(i, ceil([5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90] * iter / 100))
    if ismember(i, ceil([5 50 100] * iter / 100))
        % Display network memory visualization
        hfig1 = figure('Menubar','none','Toolbar','none','NumberTitle','off', ...
        'Name',['Training at ', num2str((i/iter)*100), '%'], ...
        'Position',[100 100 800 400]);
        memories1_thresh(w, sensory, tc, spiny, col, lin); % Visualize the memories
        set(hfig1, 'Color', 'w');
        disp(['Wait - Processing at ', num2str((i/iter)*100), '% .']);

        % Display input P_temp at this stage
        hfig2 = figure('Menubar','none','Toolbar','none','NumberTitle','off',...
            'Name',['Input at ', num2str((i/iter)*100), '%']);
        imagesc(P_temp, [0 0.16]); % Visualize current input
        colormap(gray); % Use grayscale for better contrast
        colorbar;
        xlabel('Patterns');
        ylabel('Features');
        title(['Input Matrix P at ', num2str((i/iter)*100), '% Training']);
        set(hfig2, 'Color', 'w');

    end


   P_temp = P;

   % % Apply the STM effect at 70% without SR
   %  if i >= ceil(70 * iter / 100)
   %      P_temp = P_STM;  % STM effect from 70% onward
   %  else
   %      P_temp = P;
   %  end

   % Apply SR at 60% 
      % if i >= ceil(60 * iter / 100)
      %     P_temp = P_SR;
      % else
      %     P_temp = P;
      % end

   % % Apply input reduction at 60% and STM effect at 70%
   %  if i >= ceil(70 * iter / 100)
   %      P_temp = P_STM_SR;  % STM effect from 70% onward
   %  elseif (i>=ceil(60*iter/100))
   %      P_temp = P_SR;  % Weaker input effect from 60% to 69%
   %  else
   %      P_temp = P;  % Normal input before 60%
   %  end
    
    % P_temp = P;  % Normal input before 60%

    % if i == ceil(60 * iter / 100) || i == ceil(70 * iter / 100)
    % disp(['Iteration ' num2str(i) ': Checking P_temp'])
    % disp(P_temp(1:5, 1:5))
    % end

    % S_R at 60% and STM at 85% 
    % if ((i>=ceil(60*iter/100)) && (i<ceil(85*iter/100)))
    %      P_temp = P_SR; % change the input into vague ones
    % elseif (i>=ceil(85*iter/100))
    %     P_temp = P_STM_SR; % add ramdom X
    % end
    
    recorded_inputs(:,:,i) = P_temp;
    
    %% PRUNING
    % aplica_Pruning_Huttenlocher;
    
    for j = 1:1:padroes
        output = zeros(n_neuronios, 1);
        output_antes = output;
        output(1:n_entradas,1) = P_temp(:,j);
        
        for k = 1:1:(2*camadas) % w updates twice per layer 
            w_old = w;
            w = w + incw;
            
            %if mod(i, 500) == 0 % Print every 500 iterations
            % disp(['Iteration ', num2str(i), ' | Sum(abs(w)): ', num2str(sum(abs(w(:))))]);
            %end

            % if (i>=ceil(75*iter/100))
            %     w = w.* mask_pruning;
            % end
            
            % Starting of the Divisive Normalization manipulation
            % Normal DN function: 
            norm_1=output(sensory+1:sensory+tc)'*ones(size(output(sensory+1:sensory+tc)));
            if  norm_1==0
                norm_1=1;        
            end
            % % fprintf('Computed norm_1 for Subject %d: %.6f\n', subj, norm_1);
            % 
            % % %% GABA-A -------------
            % % no GABA impairment:(update tc activation using the
            % % norm-division method:
            output(sensory+1:sensory+tc)=output(sensory+1:sensory+tc)./norm_1;
            % 
            % % % % GABA reduced at 50% 
            % if (i<ceil(50*iter/100))
            %    output(sensory+1:sensory+tc)=output(sensory+1:sensory+tc)./norm_1;
            % else       

            % GABA_scaled = 0.45; 
            % inhibtion factor using linear function:      
            % inhibition_factor = GABA_scaled * norm_1 + (1 - GABA_scaled) * 1;
            % 
            % % Apply updated inhibition factor
            % output(sensory+1:sensory+tc) = output(sensory+1:sensory+tc) ./ inhibition_factor;
            % end

            % % GABA reduced at 50% and recovered at 75%  
            % if (i<ceil(50*iter/100))
            %     output(sensory+1:sensory+tc)=output(sensory+1:sensory+tc)./norm_1;
            % else
            %     if((i>=ceil(50*iter/100)) &&  (i<=ceil(75*iter/100)))
            %        output(sensory+1:sensory+tc)=output(sensory+1:sensory+tc);
            %     else
            %         output(sensory+1:sensory+tc)=output(sensory+1:sensory+tc)./norm_1;
            %     end
            % end
            % END GABA-A ----------
            
            %% a [50x1] is the net input for each neuron
            a = w * output; % a --> the net input
            
            % Output [50x1] of the neurons using sigmoid activation model
            output = 1./(1+exp(-steepness.*(a-5.*shift)));
               
            % update and assign the external inputs to sensory neurons
            output(1:n_entradas,1) = P_temp(:,j);
            
            % last layer inhibition neurons update using the linear
            % function

            %% Manipulting the lateral inhibiton starts here:
            % Normal lateral inhibtion w/o manipulation:

            % output((n_neuronios-inhibit+1:n_neuronios),1) = 1.3.*a((n_neuronios-inhibit+1:n_neuronios),1);
            
            
            % output((n_neuronios-inhibit+1:n_neuronios),1) = 1.4.*a((n_neuronios-inhibit+1:n_neuronios),1);

            % Maniputation of the lateral inhibiton starting from 50% 
            if (i<ceil(50*iter/100))
               output((n_neuronios-inhibit+1:n_neuronios),1) = 1.3.*a((n_neuronios-inhibit+1:n_neuronios),1);
            else  
            
            GABA_scaled = 0.5;
            lateral_inhib_gain  = 1.3 * GABA_scaled;
            output((n_neuronios-inhibit+1:n_neuronios),1) = lateral_inhib_gain.*a((n_neuronios-inhibit+1:n_neuronios),1);
            end

            % incw: update the weights = Hebbian learning term - small
            % decay factor
            % Hebbian lerning term: delta W = Output(t) x output(t-1)
            % Small decay term: ((1+0.05).ones)*output(t-1) 
            incw = (fator_aprendiz .* (output * output_antes' - ((1+0.05).*ones(size(output)) * output_antes') .* w)) .* mascara;
            graf_incw(i) = sum(abs(incw(:)));  
            
            % % MEMANTINE -----------------
            %             if (i<=ceil(55*iter/100))
            %                  incw = (fator_aprendiz .* (output * output_antes' - ((1+0.05).*ones(size(output)) * output_antes') .* w)) .* mascara;
            %             else
            %                 incw = zeros(size(incw));
            %             end
            % % END MEMANTINA --------------
            % 
            % shift update (A.7)
            shift = (velocidade_deslocamento .* output + shift) ./ (1 + velocidade_deslocamento);
            
            %% ACh --------------
            %             if (i == ceil(50*iter/100))
            %                 shift_freezing = shift;
            %             end
            %             if(i>=ceil(50*iter/100))
            %                 shift(sensory+1:sensory+tc)= shift_freezing(sensory+1:sensory+tc);
            %             end
            %% END ACh ----------
            ss_activity_per_number(:, j, i) = output(sensory+ tc +1:sensory+tc+spiny);
            output_antes = output;
            % %% make a record of spiny neurons
            % 
            if i == ceil(50 * iter / 100) || i == ceil(60 * iter / 100) || ...
            i == ceil(70 * iter / 100) || i == ceil(80 * iter / 100)
            spiny_output = output(sensory + tc + 1:sensory + tc + spiny);
            save(['spiny_output_number_' num2str(j-1) '_epoch_' num2str((i/iter)*100) 'percent.mat'], 'spiny_output');
            end

            %% Record the 
            % if i == iter
            % % Each column j = population response (all spiny neurons) for pattern j
            % ss_population(:, j) = output(sensory + tc + 1 : sensory + tc + spiny);
            % end

        end
    end

    graf_shift(i)  = (ones(size(shift(sensory+1:(n_neuronios-inhibit))'))*shift(sensory+1:(n_neuronios-inhibit)))/(n_neuronios-sensory-inhibit);
    graf_output(i) = (ones(size(output(sensory+1:(n_neuronios-inhibit))'))*output(sensory+1:(n_neuronios-inhibit)))/(n_neuronios-sensory-inhibit);
    
end


hfig = figure('Menubar','none','Toolbar','none','NumberTitle','off','Name','100%','Position',[100 100 800 400]);
memories1_thresh(w,sensory,tc,spiny,col,lin);
set(hfig,'Color','w');

x_graf = linspace(1,100,iter);
hfig = figure('Color','w','Menubar','none','Toolbar','none','NumberTitle','off','Name','100%','Position',[100 100 800 400]);
plot(x_graf,graf_output,'r');
hold on;
plot(x_graf,graf_shift,'b');
hold off;
legend('Output','Shift');
xlabel('Training (%)')
nomeFig = 'Iter_';
disp('Processing Completed !!! ');

figure('Color','w', ...
       'Units','pixels', ...
       'Position',[100 100 800 400]);
plot(1:iter, graf_incw);
xlabel('Iteration');
ylabel('Total Absolute Weight Change');
title('Evolution of Weight Updates Over Training');

% drawnow;                            % let MATLAB finish layout
% set(hfig3, 'Position', [100 1000 800 400]);   % force 2:1 at the end
