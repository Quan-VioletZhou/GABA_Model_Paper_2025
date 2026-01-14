function [] = memories1_thresh(w, sensory, tc, spiny, col, lin, varargin)
% memories1_thresh
% Visualize TC->SS weight "memories" with an adjustable threshold.
%
% USAGE (keeps old behavior by default):
%   memories1_thresh(w, sensory, tc, spiny, col, lin)
%
% OPTIONAL:
%   memories1_thresh(w, sensory, tc, spiny, col, lin, 'thresh', 0.1)
%   memories1_thresh(w, sensory, tc, spiny, col, lin, 'thresh', 0)
%   memories1_thresh(w, sensory, tc, spiny, col, lin, 'scale', 5)
%   memories1_thresh(w, sensory, tc, spiny, col, lin, 'applyThresh', false)

% ----------------------------
% Defaults (match original behavior)
% ----------------------------
p = inputParser;
p.addParameter('thresh', 0.4, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('scale', 10, @(x) isnumeric(x) && isscalar(x));
p.addParameter('applyThresh', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});

thresh      = p.Results.thresh;
scaleFactor = p.Results.scale;
applyThresh = p.Results.applyThresh;

% ----------------------------
% Core computation (unchanged)
% ----------------------------
for k = 1:spiny

    % TC -> SS weights
    mem  = w(1:tc, sensory+tc+1:sensory+tc+spiny);
    mem1 = reshape(mem(:,k), col, lin);

    % SS -> TC weights (kept zeroed, same as original)
    pc  = w(sensory+tc+1:sensory+tc+spiny, sensory+1:sensory+tc);
    pc1 = reshape(pc(k,:), col, lin);
    pc1 = zeros(size(pc1));

    subplot(1, spiny, k);

    % subplot(2, ceil(spiny/2), k);

    % ----------------------------
    % Visualization transform
    % ----------------------------
    mem_plot = mem1' * scaleFactor;

    if applyThresh && thresh > 0
        mem_plot = mem_plot .* (mem_plot > thresh);
    end

    numbergraph1(mem_plot, -pc1', 1, 0);
    axis off;
end
end
