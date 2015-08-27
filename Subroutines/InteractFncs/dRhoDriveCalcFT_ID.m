function [dRhoDrive_FT] = dRhoDriveCalcFT_ID(rho,v,phi3D,kx3D,ky3D)

vx = v .* cos(phi3D);
vy = v .* sin(phi3D);
dRhoDrive_FT = - sqrt(-1) .* kx3D .* fftshift(fftn( rho .* vx ) ) + ... ;
    - sqrt(-1) .* ky3D .* fftshift(fftn( rho .* vy ) );

end

