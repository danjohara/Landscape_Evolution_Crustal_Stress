function basC = UpdateBasins(v,z,basC,eT,eR)
% Name: UpdateBasins
% Author: Daniel O'Hara
% Date: 10/05/2019
%
% Description: Update basin information to decide hillslope and fluvial
% regimes. Basin configurations found by iteratively starting from the
% highest node and going down the basin until the slope is either 0 or a
% different sign.

% Input:
%     v:    Model input structure.
%     z:    Model elevation.
%     basC: Current model basin structure.
%     eT:   Total timestep erosion.
%     eR:   Timestep erosion rate.

% Output:
%     basC:   New model basin structure.

% Reference:
% O'Hara, D., Karlstrom, L., & Roering, J. J. (2019). Distributed landscape 
%    response to localized uplift and the fragility of steady states. Earth
%    and Planetary Science Letters, 506, 243-254.

%% Initialize basC

oldBasC = basC;
basC = [];
basC.singleBRidgeX  = [];
basC.singleBRidgeZ = [];
basC.singleBCHeadX = [];
basC.singleBCHeadZ = [];
basC.singleBCHeadI = [];
basC.singleBCHeadE = [];
basC.singleBBLX = [];
basC.singleBBLZ = [];
basC.singleBBLI = [];
basC.singleBHI = [];
basC.singleBCNodesI = {};
basC.singleBHSNodesI = {};
basC.singleBRNodesI = {};
basC.singleBLengthsI = {};
basC.erosionRates = eR;
basC.totalErosion = eT;

slopeThresh = 1e-4; % was 1e-6.

%% Locate divides and minimums.
% Create collector and flag variables.
allDivides = [];
allMins = [];

foundAll = zeros(size(z));
tmpZ = z;

% Start at highest node, and continue down basin in both sides until a
% different slope sign is found. Afterwards, add highest node to divide
% variable and lowerest node to minimum variable.
while sum(foundAll==0)>0
    [~,tmpD] = find(tmpZ==max(tmpZ),1);
    stopLeft = 0;
    stopRight = 0;
    counterLeft = tmpD;
    counterNoPlatLeft = tmpD;
    ignorePlatLeft = 0;
    hitPlatLeft = 0;
    counterRight = tmpD;
    counterNoPlatRight = tmpD;
    ignorePlatRight = 0;
    hitPlatRight = 0;
    
    while stopLeft == 0 || stopRight ==0
        foundAll(counterLeft) = 1;
        tmpZ(counterLeft) = 0;
        if stopLeft == 0 && counterLeft ~= 1 && (z(counterLeft) >= z(counterLeft-1) || abs(z(counterLeft)-z(counterLeft-1)) < 1e-3)
            if hitPlatLeft == 0 && abs(z(counterLeft)-z(counterLeft-1))/v.dx > slopeThresh
                counterNoPlatLeft = counterNoPlatLeft - 1;
            else
                hitPlatLeft = 1;
            end
            counterLeft = counterLeft-1;
        elseif stopLeft == 0 && counterLeft == 1
            ignorePlatLeft = 1;
            stopLeft = 1;
        else
            stopLeft = 1;
        end
        
        foundAll(counterRight) = 1;
        tmpZ(counterRight) = 0;
        if stopRight == 0 && counterRight ~= length(z) && (z(counterRight) >= z(counterRight+1) || abs(z(counterRight) - z(counterRight+1)) < 1e-3)
            if hitPlatRight == 0 && abs(z(counterRight)-z(counterRight+1))/v.dx > slopeThresh
                counterNoPlatRight = counterNoPlatRight + 1;
            else
                hitPlatRight = 1;
            end
            counterRight = counterRight+1;
        elseif stopRight == 0 && counterRight == length(z)
            ignorePlatRight = 1;
            stopRight = 1;
        else
            stopRight = 1;
        end
    end
    
    if ignorePlatLeft == 0
        counterLeft = counterNoPlatLeft;
    end
    
    if ignorePlatRight == 0
        counterRight = counterNoPlatRight;
    end
    
    if counterLeft ~= tmpD
        allMins = [allMins,counterLeft];
        allDivides = [allDivides,tmpD];
    end
    
    if counterRight ~= tmpD
        allMins = [allMins,counterRight];
        allDivides = [allDivides,tmpD];
    end
end

% Sort divides and mins.
allDivides = [sortrows(allDivides')]';
allMins = [sortrows(allMins')]';

%% Create basC
index = 1:length(z);
for i = 1:length(allDivides)
    basC.singleBRidgeX = [basC.singleBRidgeX;v.x(allDivides(i))];
    basC.singleBRidgeZ = [basC.singleBRidgeZ;z(allDivides(i))];
  
    tmpLength = abs(v.x-v.x(allDivides(i)));
    if allMins(i) < allDivides(i)
        tmpLength((v.x>=v.x(allMins(i))).*(v.x<=v.x(allDivides(i))) == 0) = NaN;
    elseif allMins(i) > allDivides(i)
        tmpLength((v.x>=v.x(allDivides(i))).*(v.x<=v.x(allMins(i))) == 0) = NaN;
    end
    
    difLength = tmpLength - v.lc;
    difLength(difLength < 0) = NaN;
    [~,CHI] = find(difLength == min(difLength));
    
    if ~isempty(CHI)
        basC.singleBCHeadX = [basC.singleBCHeadX; v.x(CHI)];
        basC.singleBCHeadZ = [basC.singleBCHeadZ; z(CHI)];
        basC.singleBCHeadI = [basC.singleBCHeadI; index(CHI)];
        basC.singleBCHeadE = [basC.singleBCHeadE; eR(CHI)];
    else
        basC.singleBCHeadX = [basC.singleBCHeadX; NaN];
        basC.singleBCHeadZ = [basC.singleBCHeadZ; NaN];
        basC.singleBCHeadI = [basC.singleBCHeadI; NaN];
        basC.singleBCHeadE = [basC.singleBCHeadE; NaN];
    end
    
    basC.singleBBLX = [basC.singleBBLX, v.x(allMins(i))];
    basC.singleBBLZ = [basC.singleBBLZ, z(allMins(i))];
    basC.singleBBLI = [basC.singleBBLI, index(allMins(i))];
    
    difLength = tmpLength;
    difLength(difLength == 0) = NaN;
    [~,CHI] = find(difLength == min(difLength));
    
    if ~isempty(CHI)
        basC.singleBHI = [basC.singleBHI, index(CHI)];
    else
        basC.singleBHI = [basC.singleBHI, NaN];
    end
    
    difLength = tmpLength - v.lc;
    [~,CI] = find(difLength >= 0);
    [~,HI] = find(difLength < 0);
    
    if ~isempty(HI)
        basC.singleBCNodesI = [basC.singleBCNodesI;{CI}];
    else
        basC.singleBCNodesI = [basC.singleBCNodesI;{[]}];
    end
    
    if ~isempty(HI)
        basC.singleBHSNodesI = [basC.singleBHSNodesI;{HI}];
    else
        basC.singleBHSNodesI = [basC.singleBHSNodesI;{[]}];
    end
    
    basC.singleBRNodesI = [basC.singleBRNodesI;{index(allDivides(i))}];
    basC.singleBLengthsI = [basC.singleBLengthsI;{tmpLength(~isnan(tmpLength))}];
end
end