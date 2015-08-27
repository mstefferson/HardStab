% Function: dRhoInteract_FT_calc
% Description: Uses the 2nd virial coefficient to calcuate the change in
% density due to hard rod interactions in k-space.
% 
% Unlike v1, takes a divergence of the entire flux. Takes the derivative of
% the product of functions. Doesn't distribute the derivative.
% Takes in parameter and grid object as inputs
function [NegDivFluxExcess_FT] = dRhoInterCalcFT_AI(rho_FT, ParamObj,GridObj,DiffMobObj)
    Nm = ParamObj.Nm;
%%%%%%%%%%%%%%%%%%%Hard rod interactions%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Excess chemical potential in position space is a convolution. In k-space, it is a
    %product. Given by the function derivative of the excess free energy w.r.t.
    %the density profile
     % Calculate rho
     rho = real( ifftn(ifftshift( rho_FT ) ) );
     
     % Just integrating over angle. Already integrated over space to get kernal.
     % Just Fourier transform w.r.t. to angle
     Kern_FT  = fftshift( fft( abs( sin( GridObj.phi3D) ),[], 3) );
     
     %Now includes the correct scale
     MuEx_FT = (2 * pi ) / ParamObj.Nm .* ParamObj.Tmp .* ParamObj.L_rod ^ 2 ...
         .* Kern_FT .* rho_FT;
     
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
    
    jx = - rho .* dMuEx_dx;    %Flux in the x direction with isostropic diffusion
    jy = - rho .* dMuEx_dy;    %Flux in the y direction with isostropic diffusion
    jm = - rho .* dMuEx_dphi;  %Flux in the angular direction with isostropic diffusion
    
    %Fourier transform these
    Jx_FT = fftshift(fftn(jx));
    Jy_FT = fftshift(fftn(jy));
    Jm_FT = fftshift(fftn(jm));
    
   %Replicate the k matrices to do a bunch of terms at once
   kx2Drep = repmat(GridObj.kx2D,[1 1 Nm-4]);
   ky2Drep = repmat(GridObj.ky2D,[1 1 Nm-4]);
    
    %Calculate the -diverance of the flux in Fourier space. ;
    %Lower Boundary terms 
% i = 1
NegDivFluxExcess_FT(:,:,1) = ...
    - sqrt(-1) .* GridObj.kx2D .* DiffMobObj.Mob_par  / 4 .* ( Jx_FT(:,:,Nm-1) + Jx_FT(:,:,3) + 2 .* Jx_FT(:,:,1) ) ;
NegDivFluxExcess_FT(:,:,1) =  NegDivFluxExcess_FT(:,:,1) ...
    - sqrt(-1) .* GridObj.kx2D .* DiffMobObj.Mob_perp / 4 .* ( Jx_FT(:,:,Nm-1) + Jx_FT(:,:,3) - 2 .* Jx_FT(:,:,1) ) ;
NegDivFluxExcess_FT(:,:,1) =  NegDivFluxExcess_FT(:,:,1) ...
    - GridObj.kx2D .* (DiffMobObj.Mob_par - DiffMobObj.Mob_perp) / 4 .* ( Jy_FT(:,:,Nm-1) - Jy_FT(:,:,3) );

NegDivFluxExcess_FT(:,:,1) =  NegDivFluxExcess_FT(:,:,1) ...
    - sqrt(-1) .* GridObj.ky2D .* DiffMobObj.Mob_par  / 4 .* ( Jy_FT(:,:,Nm-1) + Jy_FT(:,:,3) - 2 .* Jy_FT(:,:,1) ) ;
NegDivFluxExcess_FT(:,:,1) =  NegDivFluxExcess_FT(:,:,1) ...
    - sqrt(-1) .* GridObj.ky2D .* DiffMobObj.Mob_perp / 4 .* ( Jy_FT(:,:,Nm-1) + Jy_FT(:,:,3) + 2 .* Jy_FT(:,:,1) ) ;
NegDivFluxExcess_FT(:,:,1) =  NegDivFluxExcess_FT(:,:,1) ...
    - GridObj.ky2D .* (DiffMobObj.Mob_par - DiffMobObj.Mob_perp) / 4 .* ( Jx_FT(:,:,Nm-1) - Jx_FT(:,:,3) );
    % i = 2
NegDivFluxExcess_FT(:,:,2) = ...
    - sqrt(-1) .* GridObj.kx2D .* DiffMobObj.Mob_par  / 4 .* ( Jx_FT(:,:,Nm) + Jx_FT(:,:,4) + 2 .* Jx_FT(:,:,2) ) ;
NegDivFluxExcess_FT(:,:,2) =  NegDivFluxExcess_FT(:,:,2) ...
    - sqrt(-1) .* GridObj.kx2D .* DiffMobObj.Mob_perp / 4 .* ( Jx_FT(:,:,Nm) + Jx_FT(:,:,4) - 2 .* Jx_FT(:,:,2) ) ;
NegDivFluxExcess_FT(:,:,2) =  NegDivFluxExcess_FT(:,:,2) ...
    - GridObj.kx2D .* (DiffMobObj.Mob_par - DiffMobObj.Mob_perp) / 4 .* ( Jy_FT(:,:,Nm) - Jy_FT(:,:,4) );

NegDivFluxExcess_FT(:,:,2) =  NegDivFluxExcess_FT(:,:,2) ...
    - sqrt(-1) .* GridObj.ky2D .* DiffMobObj.Mob_par  / 4 .* ( Jy_FT(:,:,Nm) + Jy_FT(:,:,4) - 2 .* Jy_FT(:,:,2) ) ;
NegDivFluxExcess_FT(:,:,2) =  NegDivFluxExcess_FT(:,:,2) ...
    - sqrt(-1) .* GridObj.ky2D .* DiffMobObj.Mob_perp / 4 .* ( Jy_FT(:,:,Nm) + Jy_FT(:,:,4) + 2 .* Jy_FT(:,:,2) ) ;
NegDivFluxExcess_FT(:,:,2) =  NegDivFluxExcess_FT(:,:,2) ...
    - GridObj.ky2D .* (DiffMobObj.Mob_par - DiffMobObj.Mob_perp) / 4 .* ( Jx_FT(:,:,Nm) - Jx_FT(:,:,4) );
    
%     keyboard
    %No boundary terms
    if Nm >4
        
    NegDivFluxExcess_FT(:,:,3:Nm-2) = ...
        - sqrt(-1) .* kx2Drep .* DiffMobObj.Mob_par  / 4 .* ( Jx_FT(:,:,1:Nm-4) + Jx_FT(:,:,5:Nm) + 2 .* Jx_FT(:,:,3:Nm-2) ) ;
    NegDivFluxExcess_FT(:,:,3:Nm-2) =  NegDivFluxExcess_FT(:,:,3:Nm-2) ...
        - sqrt(-1) .* kx2Drep .* DiffMobObj.Mob_perp / 4 .* ( Jx_FT(:,:,1:Nm-4) + Jx_FT(:,:,5:Nm) - 2 .* Jx_FT(:,:,3:Nm-2) ) ;
    NegDivFluxExcess_FT(:,:,3:Nm-2) =  NegDivFluxExcess_FT(:,:,3:Nm-2) ...
        - kx2Drep .* (DiffMobObj.Mob_par - DiffMobObj.Mob_perp) / 4 .* ( Jy_FT(:,:,1:Nm-4) - Jy_FT(:,:,5:Nm) ); 
    
    NegDivFluxExcess_FT(:,:,3:Nm-2) =  NegDivFluxExcess_FT(:,:,3:Nm-2) ...
        - sqrt(-1) .* ky2Drep .* DiffMobObj.Mob_par  / 4 .* ( Jy_FT(:,:,1:Nm-4) + Jy_FT(:,:,5:Nm) - 2 .* Jy_FT(:,:,3:Nm-2) ) ;
    NegDivFluxExcess_FT(:,:,3:Nm-2) =  NegDivFluxExcess_FT(:,:,3:Nm-2) ...
        - sqrt(-1) .* ky2Drep .* DiffMobObj.Mob_perp / 4 .* ( Jy_FT(:,:,1:Nm-4) + Jy_FT(:,:,5:Nm) + 2 .* Jy_FT(:,:,3:Nm-2) ) ;
    NegDivFluxExcess_FT(:,:,3:Nm-2) =  NegDivFluxExcess_FT(:,:,3:Nm-2) ...
        - ky2Drep .* (DiffMobObj.Mob_par - DiffMobObj.Mob_perp) / 4 .* ( Jx_FT(:,:,1:Nm-4) - Jx_FT(:,:,5:Nm) );
    
    end
    
    
    %Upper Boundary terms
    % i = N-1
    NegDivFluxExcess_FT(:,:,Nm-1) = ...
        - sqrt(-1) .* GridObj.kx2D .* DiffMobObj.Mob_par  / 4 .* ( Jx_FT(:,:,Nm-3) + Jx_FT(:,:,1) + 2 .* Jx_FT(:,:,Nm-1) ) ;
    NegDivFluxExcess_FT(:,:,Nm-1) =  NegDivFluxExcess_FT(:,:,Nm-1) ...
        - sqrt(-1) .* GridObj.kx2D .* DiffMobObj.Mob_perp / 4 .* ( Jx_FT(:,:,Nm-3) + Jx_FT(:,:,1) - 2 .* Jx_FT(:,:,Nm-1) ) ;
    NegDivFluxExcess_FT(:,:,Nm-1) =  NegDivFluxExcess_FT(:,:,Nm-1) ...
        - GridObj.kx2D .* (DiffMobObj.Mob_par - DiffMobObj.Mob_perp) / 4 .* ( Jy_FT(:,:,Nm-3) - Jy_FT(:,:,1) ); 
    
    NegDivFluxExcess_FT(:,:,Nm-1) =  NegDivFluxExcess_FT(:,:,Nm-1) ...
        - sqrt(-1) .* GridObj.ky2D .* DiffMobObj.Mob_par  / 4 .* ( Jy_FT(:,:,Nm-3) + Jy_FT(:,:,1) - 2 .* Jy_FT(:,:,Nm-1) ) ;
        NegDivFluxExcess_FT(:,:,Nm-1) =  NegDivFluxExcess_FT(:,:,Nm-1) ...
        - sqrt(-1) .* GridObj.ky2D .* DiffMobObj.Mob_perp / 4 .* ( Jy_FT(:,:,Nm-3) + Jy_FT(:,:,1) + 2 .* Jy_FT(:,:,Nm-1) ) ;
    NegDivFluxExcess_FT(:,:,Nm-1) =  NegDivFluxExcess_FT(:,:,Nm-1) ...
        - GridObj.ky2D .* (DiffMobObj.Mob_par - DiffMobObj.Mob_perp) / 4 .* ( Jx_FT(:,:,Nm-3) - Jx_FT(:,:,1) );
       
    % i = N
    NegDivFluxExcess_FT(:,:,Nm) = ...
        - sqrt(-1) .* GridObj.kx2D .* DiffMobObj.Mob_par  / 4 .* ( Jx_FT(:,:,Nm-2) + Jx_FT(:,:,2) + 2 .* Jx_FT(:,:,Nm) ) ;
    NegDivFluxExcess_FT(:,:,Nm) =  NegDivFluxExcess_FT(:,:,Nm) ...
        - sqrt(-1) .* GridObj.kx2D .* DiffMobObj.Mob_perp / 4 .* ( Jx_FT(:,:,Nm-2) + Jx_FT(:,:,2) - 2 .* Jx_FT(:,:,Nm) ) ;
    NegDivFluxExcess_FT(:,:,Nm) =  NegDivFluxExcess_FT(:,:,Nm) ...
        - GridObj.kx2D .* (DiffMobObj.Mob_par - DiffMobObj.Mob_perp) / 4 .* ( Jy_FT(:,:,Nm-2) - Jy_FT(:,:,2) ); 
    
    
    NegDivFluxExcess_FT(:,:,Nm) =  NegDivFluxExcess_FT(:,:,Nm) ...
        - sqrt(-1) .* GridObj.ky2D .* DiffMobObj.Mob_par  / 4 .* ( Jy_FT(:,:,Nm-2) + Jy_FT(:,:,2) - 2 .* Jy_FT(:,:,Nm) ) ;
    NegDivFluxExcess_FT(:,:,Nm) =  NegDivFluxExcess_FT(:,:,Nm) ...
        - sqrt(-1) .* GridObj.ky2D .* DiffMobObj.Mob_perp / 4 .* ( Jy_FT(:,:,Nm-2) + Jy_FT(:,:,2) + 2 .* Jy_FT(:,:,Nm) ) ;
    NegDivFluxExcess_FT(:,:,Nm) =  NegDivFluxExcess_FT(:,:,Nm) ...
        - GridObj.ky2D .* (DiffMobObj.Mob_par - DiffMobObj.Mob_perp) / 4 .* ( Jx_FT(:,:,Nm-2) - Jx_FT(:,:,2) );
    
    %Add the C(k) term last
     NegDivFluxExcess_FT = NegDivFluxExcess_FT - DiffMobObj.Mob_rot .* sqrt(-1) .* GridObj.km3D .* Jm_FT;
    
%      dRhoHDInteract_FT   =  NegDivFluxExcess_FT .* delta_t;
    
%      keyboard