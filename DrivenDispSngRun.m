%% General dispersion relation maker proof of concept

function DrivenDispSngRun
CurrentDir = pwd;
addpath( genpath( CurrentDir) );

Interactions = 1;
Diffusion    = 1;
SparseMat    = 0;
AnisoDiff    = 0;
PerturbGen   = 0;
SaveMe       = 0;

ModeX    = 1;
ModeY    = 1;

vD    = 0;
bcE   = 1.4; % Equilbirum bc
bcP   = 1.6; % Perturbation bc

Nx    = 2^6;
Ny    = Nx;
Nm    = 2^6;

kxHolder = Nx/2+1 + ModeX;
kyHolder = Ny/2+1 + ModeY;

L_rod   = 1;
Lx      = 10;
Ly      = Lx;
Mob_0   = 1;

Mob_par  = 2 * Mob_0;
Mob_perp = Mob_0;
Mob_rot = 6 *Mob_par/(L_rod^2);

D_par = Mob_par;

if AnisoDiff;
    D_perp = Mob_perp;
    D_rot = Mob_rot;
else
    D_par  = 1;
    D_perp = D_par;
    D_rot  = D_par;
end


DiffMobObj = struct('Mob_par', Mob_par,'D_par',D_par,'D_perp',D_perp,...
    'Mob_rot', Mob_rot,'D_rot',D_rot);

ParamObj = struct('Nx',Nx,'Ny',Ny,'Nm',Nm,'Lx',Lx,'Ly',Ly,'L_rod',L_rod,...
    'bcE',bcE,'bcP',bcP,'vD',vD);

GridObj = DispGridMaker(...
    ParamObj.Nx,ParamObj.Ny,ParamObj.Nm,ParamObj.Lx,ParamObj.Ly);

if PerturbGen
[eigVecs,eigVals] = ...
    DispEigCalcGen(DiffMobObj,GridObj,ParamObj,...
    kxHolder,kyHolder,Interactions);
else
    [eigVecs,eigVals] = ...
    DispEigCalcIsoSS(DiffMobObj,GridObj,ParamObj,Interactions,Diffusion,...
    SparseMat,kxHolder,kyHolder);
end

max( real( diag( eigVals ) ) )
% keyboard
eigValsGenR = sort( real( diag(eigVals) ) );
eigValsGenI = sort( imag( diag(eigVals) ) );

% plot( real( ( diag(eigVals) ) ) )
% keyboard
% plotyy(1:Nm, [ eigValsGenR ],1:Nm, [ eigValsGenI ])
legend('real','imag')

[FmDistBtwnPts] = ...
    MayerFncDiffBtwPntsCalc(Nx, Ny, Nm, Lx, Ly, GridObj.dx, GridObj.dy, GridObj.dphi, L_rod);
    
    Fm_FT = fftshift( fftn( FmDistBtwnPts ) );
 
 %%  
 kxHolder = Nx/2  + 1 - 1;
 kyHolder = Ny/2  + 1;
 kmHolder = Nm/2  + 1 - 2;
 
 EigContrFixed = -( GridObj.kx(kxHolder )^2 + GridObj.ky( kyHolder )^2 + GridObj.km(kmHolder)^2 ) .* ...
 (1 -  ( ParamObj.bcP * pi * Lx * Ly) / (Nx * Ny * Nm * L_rod ^2) * Fm_FT(kxHolder, kyHolder ,kmHolder) )   

%%
 kxHolder = Nx/2  + 1;
 kyHolder = Ny/2  + 1;
 
Range = 11;
 LowB = Nm/2 + 1 - (Range - 1)/2;
 UpB = Nm/2 + 1 + (Range - 1)/2;

 EigContr =  ...
     -(1 -  ( ParamObj.bcP * pi * Lx * Ly) / (Nx * Ny * Nm * L_rod ^2) *...
     Fm_FT(kxHolder, kyHolder ,LowB:UpB) );   

 
plot( reshape( real( EigContr), [1, Range ] ) )
%%


Range    = 11;
LowB     = Nx/2 + 1 - (Range-1) / 2 ;
UpperB   = Nx/2 + 1 + (Range-1) / 2;
kmHolder = Nm/2 +3 ;

xSubSet  = LowB:UpperB;

EigKgrid = zeros(Range,Range);

for i = 1:Range
    for j = 1:Range
        
        EigKgrid(i,j) = -(1 -  ( ParamObj.bcP * pi * Lx * Ly) / (Nx * Ny * Nm * L_rod ^2) *...
            real( Fm_FT( xSubSet(i) , xSubSet(j) , kmHolder ) ) );
        
    end
end

% real( EigKgrid) > 0


pcolor(EigKgrid)
colorbar

keyboard

%%
function GridObj = DispGridMaker(Nx,Ny,Nm,Lx,Ly)

dx   = Lx/Nx;
dy   = Ly/Ny;
dphi = 2*pi/Nm;
% Make vectors and grids
% x                = ( 0: dx : Lx - dx );
x                = ( -Lx/2 : dx: Lx/2 - dx);
% y                = ( 0: dy : Ly - dy );
y                = ( -Ly/2 : dy: Ly/2 - dy);
phi              = ( 0: dphi: (2*pi - dphi) );
% Make k-space spacings
dkx          = 2*pi/Lx;
dky          = 2*pi/Ly;
dkm          = 1;
% Make k vectors and grids
kx_max           = pi/dx;            %Maximum spatial k-vector allowed by grid
ky_max           = pi/dy;            %Maximum spatial k-vector allowed by grid
km_max           = pi / dphi;        %Maximum angular k-vector
kx               = ( -kx_max: dkx: (kx_max - dkx) );
ky               = ( -ky_max: dky: (ky_max - dky) );
km               = ( -km_max: dkm: (km_max - dkm) );

%Put it all in an object
GridObj = struct('x',x,'y',y,'phi',phi,'dx',dx,'dy',dy,'dphi',dphi,...
    'kx',kx,'ky',ky,'km',km);
