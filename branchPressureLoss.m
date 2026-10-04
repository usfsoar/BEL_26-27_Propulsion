
function result = branchPressureLoss( ...
    mdot, rho, mu, L, D, roughness, ...
    K, Cv, Cd, Ainj, dz, options)
% branchPressureLoss
% Calculates total pressure loss through one complete feed branch.
%
% Includes:
%   - straight pipe
%   - fittings
%   - valve
%   - injector
%   - gravity
%
% Optional:
%   fFixed = assumed Darcy friction factor
%
% If fFixed is supplied, pipePressureLoss does not require viscosity.

arguments
    mdot (1,1) double {mustBeNonnegative, mustBeFinite}
    rho (1,1) double {mustBePositive, mustBeFinite}

    % Can be NaN when fixed friction factor is used
    mu (1,1) double

    L (1,1) double {mustBeNonnegative, mustBeFinite}
    D (1,1) double {mustBePositive, mustBeFinite}
    roughness (1,1) double {mustBeNonnegative, mustBeFinite}

    K (1,1) double {mustBeNonnegative, mustBeFinite}
    Cv (1,1) double {mustBePositive, mustBeFinite}
    Cd (1,1) double {mustBePositive, mustBeFinite}
    Ainj (1,1) double {mustBePositive, mustBeFinite}

    dz (1,1) double {mustBeFinite}

    options.fFixed (1,1) double = NaN
end


%% Pipe pressure loss

pipe = pipePressureLoss( ...
    mdot, rho, mu, L, D, roughness, ...
    fFixed=options.fFixed);


%% Fitting pressure loss

fitting = fittingPressureLoss( ...
    K, rho, pipe.v);


%% Valve pressure loss

valve = valvePressureLoss( ...
    mdot, rho, Cv);


%% Injector pressure loss

injector = injectorPressureLoss( ...
    mdot, rho, Cd, Ainj);


%% Gravity pressure change

g = 9.80665;      % [m/s^2]

dP_gravity = rho * g * dz;


%% Total branch pressure loss

dP_total = ...
    pipe.dP + ...
    fitting.dP + ...
    valve.dP + ...
    injector.dP + ...
    dP_gravity;


%% Outputs

result.dP_pipe = pipe.dP;
result.dP_fitting = fitting.dP;
result.dP_valve = valve.dP;
result.dP_injector = injector.dP;
result.dP_gravity = dP_gravity;

result.dP_total = dP_total;

result.v = pipe.v;
result.Re = pipe.Re;
result.f = pipe.f;
result.regime = pipe.regime;

end