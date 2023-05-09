function PlotXZ(v,L,z,basC,t)
% Name: PlotXZ
% Author: Daniel O'Hara
% Date: 10/05/2019
%
% Description: Plot model landscape.

% Input:
%     v:        Model input structure.
%     L:        Variable-lithology structure.
%     z:        Model elevations.
%     basC:     Model structured basin variable.
%     t:        Current model time.

% Reference:
% O'Hara, D., Karlstrom, L., & Roering, J. J. (2019). Distributed landscape 
%    response to localized uplift and the fragility of steady states. Earth
%    and Planetary Science Letters, 506, 243-254.

if v.allPlot == 0
    % Plot all landscape with no distinction.
    plot(v.x,z)
elseif v.allPlot == 1
    % Plot only channels with distinction.
    plot(v.x(basC.singleBCNodesI{1}),z(basC.singleBCNodesI{1}))
    hold on
    for i = 2:length(basC.singleBCNodesI)
        plot(v.x(basC.singleBCNodesI{i}),z(basC.singleBCNodesI{i}))
    end
    hold off
elseif v.allPlot ==2
    % Plot entire landscape with channel distinction.
    plot(v.x,z)
    hold on
    for i = 1:length(basC.singleBCNodesI)
        plot(v.x(basC.singleBCNodesI{i}),z(basC.singleBCNodesI{i}))
    end
    hold off
end

% Plot layer (if applicable)
if L.runLithologyModel && L.plotLithology
    yl = ylim;
    xl = xlim;
    ti1 = find(~isnan(L.zLp),1);
    ti2 = find(~isnan(L.zLp),1,'last');
    ti3 = find(~isnan(L.zLn),1,'last');
    ti4 = find(~isnan(L.zLn),1);
    lX = [v.x(ti1:ti2),fliplr(v.x(ti4:ti3)),v.x(ti1)];
    lZ = [L.zLp(ti1:ti2),fliplr(L.zLn(ti4:ti3)),L.zLp(ti1)];
    hold on
    plot(lX,lZ,'--k')
    hold off
    xlim(xl)
    ylim(yl)
end

% Add labels and title.
xlabel('X (m)')
ylabel('Z (m)')
title({sprintf('Timestep = %.0f; d_x = %.0f',t,v.dx);sprintf('X_r = %.4f; Z_r = %.4f',basC.singleBRidgeX(1),basC.singleBRidgeZ(1))})
% axis image
drawnow

end

