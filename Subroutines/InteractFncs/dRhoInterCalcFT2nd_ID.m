% Function: dRhoInteract_FT_calc
% Description: Uses the 3nd virial coefficient to calcuate the change in
% density due to hard rod interactions in k-space.
%
% Heavy approximations
% Unlike v1, takes a divergence of the entire flux. Takes the derivative of
% the product of functions. Doesn't distribute the derivative.
% Takes in parameter and grid object as inputs
function [NegDivFluxExcess_FT] = dRhoInterCalcFT2nd_ID(rho,ParamObj,GridObj,DiffMobObj)
%%%%%%%%%%%%%%%%%%%Hard rod interactions%%%%%%%%%%%%%%%%%%%%%%%%%%

% A drastic approximation of the excess chemical potential from the
% 3rd virial coefficient
MuEx_FT = - 1 / 2 .* (pi * 5 / 2 + 3) .* ParamObj.L_r.* fftshiftfft(rho .* rho);
%     MuEx    = real(ifftn(ifftshift(MuEx_FT)));

%Takes its derivative in k-space
dMuEx_dx_FT   =     sqrt(-1) .* GridObj.kx3D .*  MuEx_FT;
dMuEx_dy_FT   =     sqrt(-1) .* GridObj.ky3D .*  MuEx_FT;
dMuEx_dphi_FT =     sqrt(-1) .* GridObj.km3D .*  MuEx_FT;

%Excess chemical potential derivative in real space
%Mayer function derivative in real-space
dMuEx_dx   =  real(ifftn(ifftshift(dMuEx_dx_FT)));
dMuEx_dy   =  real(ifftn(ifftshift(dMuEx_dy_FT)));
dMuEx_dphi =  real(ifftn(ifftshift(dMuEx_dphi_FT)));

%Do the hard disk interaction portion of the PDE in real space then FT it
% Isolate the seperate parts and call them some arbitrary function. We
% will Fourier transform these functions to solve this in Fourier space
%
% Take the divergence of the product of functions. Call these products
% random variables

jx = - DiffMobObj.Mob_pos .* rho .* dMuEx_dx;    %Flux in the x direction with isostropic diffusion
jy = - DiffMobObj.Mob_pos .* rho .* dMuEx_dy;    %Flux in the y direction with isostropic diffusion
jm = - DiffMobObj.Mob_rot .* rho .* dMuEx_dphi;  %Flux in the angular direction with isostropic diffusion

%Fourier transform these
Jx_FT = fftshift(fftn(jx));
Jy_FT = fftshift(fftn(jy));
Jm_FT = fftshift(fftn(jm));

% Calculate the - divergence of the interaction flux
NegDivFluxExcess_FT = - sqrt(-1) .* ( GridObj.kx3D .* Jx_FT + ...
    GridObj.ky3D .* Jy_FT + GridObj.km3D .* Jm_FT );

