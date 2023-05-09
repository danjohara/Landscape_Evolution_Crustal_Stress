function [x,l,z,basC,xcp] = CreateInitialProfile(v,er)
% Name: CreateInitialProfile
% Author: Daniel O'Hara
% Date: 10/05/2019
%
% Description: Creates analytical solution for initial landscape given
% inputs.

% Input:
%     v:        Model input structure.
%     er:       Model erosion structure.

% Output:
%     x:      Model grid.
%     l:      Model grid distances from ridge.
%     z:      Model elevations.
%     basC:   Model basin variables.
%     xcp:    Critical length for fluvial regime.

% Reference:
% O'Hara, D., Karlstrom, L., & Roering, J. J. (2019). Distributed landscape 
%    response to localized uplift and the fragility of steady states. Earth
%    and Planetary Science Letters, 506, 243-254.

% Error Check
if v.xr0<v.xRange(1) || v.xr0 > v.xRange(2)
    error('Ridge outside of domain')
end

if isfield(v,'prevFv')
    x = v.prevFv.v{end}.x;
    l = v.prevFv.lt(end,:);
    z = v.prevFv.zt(end,:);
    basC = v.prevFv.basCs{end};
    xcp = v.prevFv.v{end}.lc;
    return;
end



% Establish variables
x = v.xRange(1):v.dx:v.xRange(2);
l = x - v.xr0;
zln = v.sinkHeight; % Left basin base level.
zlp = v.sinkHeight; % Rigt basin base level.
un = er.u0; % Background uplift that creates left basin height.
h = er.h;
m = er.m;
k = er.k;
ka = er.ka;
D = er.D;
xcn = -(D/(k*ka^m))^(1/(1+h*m)); % Left basin critical length.
xln = -abs(v.xr0-v.xRange(1)); % Left basin length.
xlp = abs(v.xr0-v.xRange(2)); % Right basin length.
xcp = (D/(k*ka^m))^(1/(1+h*m)); % Right basin critical length.
up = (zln - zlp + un*(abs(xln)^(1-h*m) - abs(xcn)^(1-h*m))/(k*ka^m*(1-h*m)) ...
    + un*xcn^2/(2*D)) / ((xlp^(1-h*m) - xcp^(1-h*m))/(k*ka^m*(1-h*m)) ...
    + xcp^2/(2*D)); % Uplift rate required to set right basin to same analytic
                    % height as left basin while keeping length set by
                    % ridge location.

z1 = zln + un*(abs(xln)^(1-h*m) - abs(l).^(1-h*m))./(k*ka^m*(1-h*m));
    % Left basin analytic fluvial landscape solution (valid only on
    % [xln,xcn].
z2 = zlp + up*(abs(xlp)^(1-h*m) - abs(l).^(1-h*m))./(k*ka^m*(1-h*m));
    % Right basin analytic fluvial landscape solution (valid only on
    % [xcp,xlp].
z3 = zln + un*(abs(xln)^(1-h*m) - abs(xcn)^(1-h*m))/(k*ka^m*(1-h*m)) + un*xcn^2/(2*D) - un*l.^2/(2*D);
    % Left basin analytic hillslope solution (valid only on [xcn,0].
z4 = zlp + up*(abs(xlp)^(1-h*m) - abs(xcp)^(1-h*m))/(k*ka^m*(1-h*m)) + up*xcp^2/(2*D) - up*l.^2/(2*D);
    % Right basin analytic hillslope solution (valid only on [0,xcp].
z = [z1(l<=xcn),z3(((l>xcn).*(l<=0))==1),z4(((l>0).*(l<xcp))==1),z2(l>xcp)];
    % Combine analytic solution based on their domains.

% Analytic values for ridge and critical length heights.
topZ = zln + un*(abs(xln)^(1-h*m) - abs(xcn)^(1-h*m))/(k*ka^m*(1-h*m)) + un*xcn^2/(2*D);
zcn = zln + un*(abs(xln)^(1-h*m) - abs(xcn)^(1-h*m))/(k*ka^m*(1-h*m));
zcp = zlp + up*(abs(xlp)^(1-h*m) - abs(xcp)^(1-h*m))/(k*ka^m*(1-h*m));
ind = 1:length(x);

% Save basin statistics
basC.singleBRidgeX = [v.xr0;v.xr0];           
    % Basin ridge X value.
basC.singleBRidgeZ = [topZ;topZ];       
    % Basin ridge Z value.
basC.singleBCHeadX = [v.xr0+xcn;v.xr0+xcp];   
    % Basin channel Head X value.
basC.singleBCHeadZ = [zcn;zcp];                
    % Basin channel Head Z value.
basC.singleBCHeadI = [max(ind(l<=xcn));min(ind(l>=xcp))];           
    % Basin highest Channel Node index.
basC.singleBCHeadE = [un;up];            
    % Basin highest Channel node erosion.
basC.singleBBLX = [v.xRange(1);v.xRange(end)];    
    % Basin base Level X value.
basC.singleBBLZ = [v.sinkHeight;v.sinkHeight];  
    % Basin base Level Z value.
basC.singleBBLI = [ind(1);ind(end)];     
    % Basin lowest Channen Node index.
basC.singleBHI = [max(ind(l<0));min(ind(l>0))]; 
    % Basin highest node.
basC.singleBCNodesI = [{ind(l<=xcn)};{ind(l>=xcp)}];      
    % Basin channel node indexes.
basC.singleBHSNodesI = [{ind(((l>xcn).*(l<0))==1)};{ind(((l>0).*(l<xcp))==1)}];   
    % Basin hillslope node indexes.
basC.singleBRNodesI = [{ind(l==0)};{ind(l==0)}];    
    % Basin ridgetop node index.
basC.singleBLengthsI = [{l(l<0)};{l(l>0)}];     
    % Basin node lengths.
basC.erosionRates = [un.*(x<=v.xr0),up.*(x>=v.xr0)];
basC.erosionRates(basC.erosionRates==0) = [];
    % Basin erosion rates.
basC.totalErosion = basC.erosionRates;
    % Timestep total erosion (= rate at initial).
end