function v = createUpliftField(v,P,er,U)
% Name: createUpliftField
% Author: Daniel O'Hara
% Date: 10/05/2019
%
% Description: Generate the model uplift field, using either the uplift
% perturbation equation (from O'Hara et al., 2019), or just the provided
% uplift variable.

% Input:
%     v:    Model input structure.
%     P:    model uplift perturbation structure.
%     er:   Model erosion structure.

% Output:
%     v:   Model input structure with updated uplift field.

    v.u0 = er.u0;

    if U.runUpliftGradientModel
        v.u = linspace(U.leftU,U.rightU,length(v.x));
        if P.runPerturbationModel
            v.u = v.u + + P.umax*(1 - ((v.x-P.xP)/P.R).^P.b);
        end
    elseif P.runPerturbationModel
        v.u = ones(size(v.x)).*v.u0 + P.umax*(1 - ((v.x-P.xP)/P.R).^P.b); % Uplift function
        v.u(abs(v.x-P.xP) > P.R) = v.u0; % Uplift function bounds
    else
        v.u = v.u0;
    end
end