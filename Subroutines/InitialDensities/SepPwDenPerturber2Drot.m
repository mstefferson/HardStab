% PwDenPerturber2Drot.m
% Perturbs input rhp
% Perturbations of the form
% \sum epsilon_k exp( i k_x x ) + \sum epsilon_k exp( i k_y y ) ...
%       + \sum epsilon_k exp( i k_m phi )

function [rho] = SepPwDenPerturber2Drot(rho,ParamObj,GridObj)

%Perturbation coeff
Coeff = ParamObj.WeightPos;
%Change in x
for i = -ParamObj.NumModesX:ParamObj.NumModesX
    
    
    rho = rho + (Coeff + sqrt(-1) .* Coeff) ...
        .* exp( sqrt(-1) .* GridObj.kx(ParamObj.Nx/2+1+i) .* GridObj.x3D ) ; ...
end

%Change in y
for i = -ParamObj.NumModesY:ParamObj.NumModesY
     
    rho = rho + (Coeff + sqrt(-1) .* Coeff) ...
        .* exp( sqrt(-1) .* GridObj.ky(ParamObj.Ny/2+1+i) .* GridObj.y3D ) ; ...
end

Coeff = ParamObj.WeightAng;
%Change in phi
for i = -ParamObj.NumModesM:ParamObj.NumModesM
     
    rho = rho + (Coeff + sqrt(-1) .* Coeff) ...
        .* exp( sqrt(-1) .* GridObj.km(ParamObj.Nm/2+1+i) .* GridObj.phi3D ) ; ...
end

% Take real part
rho = real(rho);

% Fix negative rho if that happened.
[rho] = FixNegDenFnc(rho);

% Integrate first along the depth of matrix w.r.t theta, then across the
% columns w.r.t x, then down the rows w.r.t. y
CurrentNorm = trapz_periodic(GridObj.y,trapz_periodic(GridObj.x,trapz_periodic(GridObj.phi,rho,3),2),1);
rho = real(rho .* ParamObj.Norm ./ CurrentNorm);

end