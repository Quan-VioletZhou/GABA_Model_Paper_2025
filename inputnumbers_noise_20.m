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

function [P, col, lin] = inputnumbers_noise_20(noise_level)
% [P, col, lin] = inputnumbers_noise(noise_level)
%
% Generates 10 Face and 10 House patterns (20 total) on a 5x3 grid,
% with rectified Gaussian noise added to every single cell.
%
% noise_level : standard deviation of the Gaussian noise (before rectification)
%
% Example:
%   [P,~,~] = inputnumbers_noise();        % default noise σ = 0.1
%   [P,~,~] = inputnumbers_noise(0.15);    % lighter noise
if nargin < 1
    noise_level = 0.1;   % nice realistic default
end
 
rng(42);               % always perfectly reproducible

mu = 2.0;               % strength of category signal & individual tag
n_per_cat  = 5;        % <--- CHANGED: 10 exemplars per category
n_total = 2 * n_per_cat; % 20 total patterns

P = zeros(15, n_total); % <--- CHANGED: Matrix size is 15 rows x 20 columns (patterns)

%% ====================================================================
%  Category A — Faces (left column + one unique middle tag) + noise
%% ====================================================================
for i = 1:n_per_cat
    p = zeros(15,1);
    
    % Category signal (Shared Feature for Faces)
    p(1:3:13) = mu;                  % left column: positions 1,4,7,10,13
    
    % Individual "grandmother" tag (Item-wise Variance)
    % Cycles through the 5 middle rows twice to cover 10 items
    p(2 + 3*mod(i-1, 5)) = 0.5 * mu;   % uses mod(i-1, 5) to cycle through rows 2,5,8,11,14
    
    % Add Gaussian noise and rectify → strictly non-negative
    p = p + noise_level * randn(15,1);
    p = max(p, 0);                   % biological rectification
    
    P(:,i) = p;
end
%% ====================================================================
%  Category B — Houses (right column + one unique middle tag) + noise
%% ====================================================================
for i = 1:n_per_cat
    p = zeros(15,1);
    
    % Category signal (Shared Feature for Houses)
    p(3:3:15) = mu;                  % right column: 3,6,9,12,15
    
    % Individual "grandmother" tag (Item-wise Variance)
    % Cycles through the 5 middle rows twice to cover 10 items
    p(2 + 3*mod(i-1, 5)) = 0.5 * mu;       % uses mod(i-1, 5) to cycle through rows 2,5,8,11,14
    
    % Add Gaussian noise and rectify → strictly non-negative
    p = p + noise_level * randn(15,1);
    p = max(p, 0);
    
    % Place in the second half of the P matrix (columns 11 through 20)
    P(:,n_per_cat + i) = p;
end
%% Unit-norm normalization (for cosine similarity)
P = P ./ sqrt(sum(P.^2, 1));

col = 3;
lin = 5;

%% figure plot 

figure('Color','w','Position',[100 100 800 500]);

t = tiledlayout(2,6,'TileSpacing','compact','Padding','compact');

% ---- Category 1: top row tiles 1..5 ----
for i = 1:5
    ax = nexttile(t, i);  % force tile index
    imagesc(ax, reshape(P(:,i),3,5)', [0 0.65]);
    axis(ax,'image','off');
    colormap(ax, parula);
    title(ax, sprintf('%d', i), 'FontWeight','bold');
    ax.LooseInset = [0 0 0 0];
end

% ---- Category 2: bottom row tiles 7..11 ----
for j = 1:5
    i = 5 + j;            % patterns 6..10
    ax = nexttile(t, 6 + j);   % tiles 7..11
    imagesc(ax, reshape(P(:,i),3,5)', [0 0.65]);
    axis(ax,'image','off');
    colormap(ax, parula);
    title(ax, sprintf('%d', j), 'FontWeight','bold'); % label 1..5 again
    ax.LooseInset = [0 0 0 0];
end

% ---- Row category labels ----
% ---- Vertical row category labels ----
annotation('textbox', [0.01 0.62 0.05 0.2], ...
    'String', 'Category 1', ...
    'FontSize', 14, ...
    'FontWeight', 'bold', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'Rotation', 90);

annotation('textbox', [0.01 0.18 0.05 0.2], ...
    'String', 'Category 2', ...
    'FontSize', 14, ...
    'FontWeight', 'bold', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'Rotation', 90);

% ---- Colorbar in the rightmost column (tile 6 or 12 works) ----
% Create it once, then dock it to the east of the layout:
drawnow;  % ensure graphics tree finished updating

cb = colorbar;
cb.Layout.Tile = 'east';
cb.TickDirection = 'out';
cb.Box = 'off';

ylabel(cb, 'Pixel intensity (normalized/ 2)', ...
    'FontSize', 12, 'FontWeight', 'bold');



