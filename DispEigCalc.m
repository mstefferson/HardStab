function [eigVals] = ...
    DispEigCalc(DiffMobObj,GridObj,ParamObj,Interactions,Diffusion,...
    SparseMat,kxHolder,kyHolder)


Nm = ParamObj.Nm;

% Prefactors
DiagPreFac = -( DiffMobObj.D_par + DiffMobObj.D_perp) / 2 .* ...
    ( GridObj.kx(kxHolder) .^ 2 + GridObj.ky(kyHolder) .^ 2 )...
    - DiffMobObj.D_rot * GridObj.km .^2 ;
Mp1PreFac = -ParamObj.vD ./ 2 .* ( ...
    sqrt(-1) .* GridObj.kx(kxHolder) + GridObj.ky(kyHolder) );
Mm1PreFac = -ParamObj.vD ./ 2 .* ( ...
    sqrt(-1) .* GridObj.kx(kxHolder) - GridObj.ky(kyHolder) );

Mp2PreFac =  ( DiffMobObj.D_par - DiffMobObj.D_perp ) / 4 .* ...
    ( GridObj.ky(kyHolder) - sqrt(-1) .* GridObj.kx(kxHolder)  ) .^2;
Mm2PreFac = ( DiffMobObj.D_par - DiffMobObj.D_perp ) / 4 .* ...
    ( GridObj.ky(kyHolder) + sqrt(-1) .* GridObj.kx(kxHolder)  ) .^2;

% keyboard
% Build diagonal
if Diffusion
    DiffContrib = DiagPreFac;
else
    DiffContrib = zeros(1,ParamObj.Nm);
end

% rho_0 = c / (2*pi)
if Interactions
    %     tic
    Fm = MayerFncDiffBtwPntsCalc(...
        ParamObj.Nx, ParamObj.Ny, ParamObj.Nm, ...
        ParamObj.Lx, ParamObj.Ly, GridObj.dx,...
        GridObj.dy, GridObj.dphi, ParamObj.L_rod);
    Fm_FT = fftshift(fftn( Fm ) );
    %     toc
    %     keyboard
    IntPreFac = pi * ParamObj.Lx * ParamObj.Ly / ...
        ( ParamObj.L_rod * ParamObj.Nx * ParamObj.Ny * ParamObj.Nm)...
        * ParamObj.bc;
    IntDiagContrib = -IntPreFac .* DiagPreFac .*  ...
        reshape(Fm_FT(kxHolder,kyHolder,:),[1,ParamObj.Nm]);
    % Anisotropic diffusion
    Mp2vec = Mp2PreFac .* ( ones(1,ParamObj.Nm) ...
        -IntPreFac .* reshape(Fm_FT(kxHolder,kyHolder,:),[1,ParamObj.Nm])  );
    Mm2vec = Mm2PreFac .* ( ones(1,ParamObj.Nm) ...
        -IntPreFac .* reshape(Fm_FT(kxHolder,kyHolder,:),[1,ParamObj.Nm])  );
    
else
    Fm_FT = zeros(ParamObj.Nx,ParamObj.Ny,ParamObj.Nm);
    IntDiagContrib  = zeros(1,ParamObj.Nm);
    % Anisotropic diffusion
    Mp2vec = Mp2PreFac .* ones(1,ParamObj.Nm);
    Mm2vec = Mm2PreFac .* ones(1,ParamObj.Nm);
    
end %Interactions

MdiagVec = DiffContrib + IntDiagContrib;
% keyboard

% Driving
Mp1vec = Mp1PreFac .* ones(1,ParamObj.Nm);
Mm1vec = Mm1PreFac .* ones(1,ParamObj.Nm);

if SparseMat == 1
    M = spdiags( [Mm2vec' Mm1vec'  MdiagVec' Mp1vec' Mp2vec'],...
        -2:2,Nm,Nm);
    % PBC driving
    M(1,ParamObj.Nm) = Mm1PreFac;
    M(ParamObj.Nm,1) = Mp1PreFac;
    % PBC aniso
    M(1,ParamObj.Nm-1) = Mp2PreFac .* ...
        (1 -IntPreFac .* Fm_FT(kxHolder,kyHolder,1) );
    M(ParamObj.Nm-1,1) = Mp2PreFac .* ...
        (1 -IntPreFac .* Fm_FT(kxHolder,kyHolder,Nm-1)  );
    M(2,ParamObj.Nm) = Mp2PreFac .* ...
        (1 -IntPreFac .* Fm_FT(kxHolder,kyHolder,2) );
    M(ParamObj.Nm,2) = Mp2PreFac .* ...
        (1 -IntPreFac .* Fm_FT(kxHolder,kyHolder,Nm)  );
    
    [eigVals] = eigs(M,Nm-1,'lm');
      
else
    
    M =  diag( MdiagVec') +...
        diag( Mp1vec(1:Nm-1)',1) + diag( Mm1vec(1:Nm-1)',-1) +...
        diag( Mp2vec(1:Nm-2)',2) + diag( Mp2vec(1:Nm-2)',-2) ;
    % PBC
    M(1,ParamObj.Nm) = Mm1PreFac;
    M(ParamObj.Nm,1) = Mp1PreFac;
    % PBC aniso
    M(1,ParamObj.Nm-1) = Mp2PreFac .* ...
        (1 -IntPreFac .* Fm_FT(kxHolder,kyHolder,1) );
    M(ParamObj.Nm-1,1) = Mp2PreFac .* ...
        (1 -IntPreFac .* Fm_FT(kxHolder,kyHolder,Nm-1)  );
    M(2,ParamObj.Nm) = Mp2PreFac .* ...
        (1 -IntPreFac .* Fm_FT(kxHolder,kyHolder,2) );
    M(ParamObj.Nm,2) = Mp2PreFac .* ...
        (1 -IntPreFac .* Fm_FT(kxHolder,kyHolder,Nm)  );
    
    [eigVals] = eig(M);
end



end
% toc