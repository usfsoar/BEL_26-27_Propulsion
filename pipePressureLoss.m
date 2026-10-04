function result = pipePressureLoss(mdot, rho, mu, L, D, roughness, options)
% pipePressureLoss
% Calculates pressure loss through a straight pipe using Darcy-Weisbach.
%
% Normal mode:
%   Uses viscosity -> Reynolds number -> friction factor
%
% Fixed-f mode:
%   If fFixed is supplied, Reynolds number and viscosity are not needed.
%
% Inputs:
%   mdot      = mass flow rate [kg/s]
%   rho       = density [kg/m^3]
%   mu        = dynamic viscosity [Pa*s]
%               Use NaN if fFixed is supplied
%   L         = pipe length [m]
%   D         = pipe internal diameter [m]
%   roughness = pipe roughness [m]
%
% Optional name-value input:
%   fFixed    = assumed Darcy friction factor
%
% Example:
%   pipePressureLoss(..., fFixed=0.04)

arguments
    mdot (1,1) double {mustBeNonnegative, mustBeFinite}
    rho (1,1) double {mustBePositive, mustBeFinite}
    mu (1,1) double
    L (1,1) double {mustBeNonnegative, mustBeFinite}
    D (1,1) double {mustBePositive, mustBeFinite}
    roughness (1,1) double {mustBeNonnegative, mustBeFinite}

    options.fFixed (1,1) double = NaN
end


%% Pipe velocity

A = pi * D^2 / 4;

v = mdot / (rho * A);


%% Friction factor

if ~isnan(options.fFixed)

    % ----- Fixed friction-factor mode -----

    if ~isfinite(options.fFixed) || options.fFixed <= 0
        error("pipePressureLoss:InvalidFixedF", ...
            "fFixed must be a positive, finite Darcy friction factor.");
    end

    f = options.fFixed;

    % Reynolds number is not calculated in this mode
    Re = NaN;

    regime = "fixed friction factor";

else

    % ----- Reynolds-number mode -----

    if ~isfinite(mu) || mu <= 0
        error("pipePressureLoss:InvalidViscosity", ...
            "mu must be positive and finite when fFixed is not supplied.");
    end

    Re = rho * v * D / mu;


    if Re == 0

        f = 0;
        regime = "no flow";


    elseif Re < 2300

        % Laminar flow
        f = 64 / Re;
        regime = "laminar";


    elseif Re < 4000

        % Transitional flow
        % Linear interpolation between laminar and turbulent estimates

        relativeRoughness = roughness / D;

        fLaminar = 64 / Re;

        fTurbulent = 1 / ...
            (-1.8 * log10( ...
            (relativeRoughness / 3.7)^1.11 + 6.9 / Re))^2;

        weight = (Re - 2300) / (4000 - 2300);

        f = (1 - weight) * fLaminar + weight * fTurbulent;

        regime = "transitional";


    else

        % Turbulent flow - Haaland equation

        relativeRoughness = roughness / D;

        f = 1 / ...
            (-1.8 * log10( ...
            (relativeRoughness / 3.7)^1.11 + 6.9 / Re))^2;

        regime = "turbulent";

    end
end


%% Darcy-Weisbach pressure loss

dP = f * (L / D) * (rho * v^2 / 2);


%% Outputs

result.dP = dP;
result.v = v;
result.Re = Re;
result.f = f;
result.regime = regime;

end