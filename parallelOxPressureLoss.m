function result = parallelOxPressureLoss( ...
    mdotOx, rho, mu, L, D, roughness, ...
    K, Cv, Cd, AinjTotal, dz, options)
% parallelOxPressureLoss
%
% Calculates pressure loss through two identical parallel
% oxidizer feed branches.
%
% Assumptions:
%   - identical left and right branches
%   - total oxidizer mass flow splits evenly
%   - total injector area splits evenly
%   - pressure drop is the same through each branch
%
% Optional:
%   fFixed = assumed Darcy friction factor
%
% For the current N2O model, use fFixed = 0.04.

arguments
    mdotOx (1,1) double {mustBeNonnegative, mustBeFinite}
    rho (1,1) double {mustBePositive, mustBeFinite}

    % Can be NaN when fixed friction factor is used
    mu (1,1) double

    L (1,1) double {mustBeNonnegative, mustBeFinite}
    D (1,1) double {mustBePositive, mustBeFinite}
    roughness (1,1) double {mustBeNonnegative, mustBeFinite}

    K (1,1) double {mustBeNonnegative, mustBeFinite}
    Cv (1,1) double {mustBePositive, mustBeFinite}
    Cd (1,1) double {mustBePositive, mustBeFinite}
    AinjTotal (1,1) double {mustBePositive, mustBeFinite}

    dz (1,1) double {mustBeFinite}

    options.fFixed (1,1) double = NaN
end


%% Split flow between identical branches

mdotBranch = mdotOx / 2;

% Split total oxidizer injector flow area between branches
AinjBranch = AinjTotal / 2;


%% Calculate one branch

branch = branchPressureLoss( ...
    mdotBranch, ...
    rho, ...
    mu, ...
    L, ...
    D, ...
    roughness, ...
    K, ...
    Cv, ...
    Cd, ...
    AinjBranch, ...
    dz, ...
    fFixed=options.fFixed);


%% Outputs

result.mdot_left = mdotBranch;
result.mdot_right = mdotBranch;

result.left = branch;
result.right = branch;

% Parallel branches see the SAME pressure drop.
% Do not add left and right pressure losses.
result.dP_total = branch.dP_total;

end