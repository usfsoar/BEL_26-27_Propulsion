%% vehicleGeometry.m
function geom = vehicleGeometry()

inch = .0254; %inch to meter 

% Tube internal diameters 
geom.oxLeft.D = .430 * inch;
geom.oxRight.D = .430 * inch;
geom.fuel.D = .305 * inch;

%Unknown until CAD is available
geom.oxLeft.L = .762; %assumptions of 30" pipe lengths
geom.oxRight.L = .762;
geom.fuel.L = .762;

end

