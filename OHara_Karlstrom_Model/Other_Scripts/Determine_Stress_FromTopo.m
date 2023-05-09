function [sigmas,inputs] = Determine_Stress_FromTopo(z,dz,x,dx,minDepth,maxDepth,sL)

%% Create input structure
% From user
inputs.dz = dz;
inputs.x = x;
inputs.dx = dx;
inputs.Topo = z;

% Others
inputs.G = 10e9;
inputs.lambda = 20e9;
inputs.rho = 2700;

inputs.maxDepth = maxDepth;
inputs.minDepth = minDepth;

inputs.sigmaP_scale = 2000;
inputs.sigmaP_cutByX = 200;%400;
inputs.sigmaP_cutByZ = 25;%50;
inputs.magScale1 = [-.01,.1];
inputs.magScale2 = [0,.5];
inputs.plotTitle = 'tmp';

inputs.kernel_radius = 2e4;
inputs.Analysis = 1;

inputs.VectorWidth = 2;
inputs.sourceLocs = sL;

inputs.topoMaxHeight = max(z)+200;
inputs.makePlots = 0;

scriptPath = '.\..\..\';
addpath(genpath(scriptPath))
sigmas = Calculate_Crustal_Stress_2D(inputs);
end
