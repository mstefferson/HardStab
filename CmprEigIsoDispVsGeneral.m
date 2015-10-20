%% General dispersion relation maker proof of concept

function CmprEigIsoDispVsGeneral
CurrentDir = pwd;
addpath( genpath( CurrentDir) );

Interactions = 1;
Diffusion    = 1;
AnisoDiff    = 1;
SparseMat    = 0;
SaveMe       = 0;

ModeX    = 0;
ModeY    = 1;

vD    = 100;
bcE   = 1.6; % Equilbirum bc
bcP   = 1.6; % Perturbation bc

Nx    = 2^8;
Ny    = Nx;
Nm    = 2^8;

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

% [M,DiffContrib,IntDiagContrib,eigVals] = DispEigCalcTemp(DiffMobObj,GridObj,ParamObj,Interactions,Diffusion,...
%     SparseMat ,kxHolder,kyHolder);
% disp('Old meth');disp( max( real( eigVals ) ) )

[Miso,IntDiagContrib,eigValsIso] = ...
    DispEigCalcIsoSS(DiffMobObj,GridObj,ParamObj,Interactions,Diffusion,...
    SparseMat,kxHolder,kyHolder);

[Mgen,Mint,eigValsGen] = ...
    DispEigCalcGen(DiffMobObj,GridObj,ParamObj,...
    kxHolder,kyHolder,Interactions);

% keyboard
eigValsIsoR = sort( real( eigValsIso) );
eigValsGenR = sort( real( eigValsGen) );

eigValsIsoI = sort( imag( eigValsIso) );
eigValsGenI = sort( imag( eigValsGen) );


deltaEigR = (eigValsGenR - eigValsIsoR);
deltaEigI = (eigValsGenI - eigValsIsoI);

 'Max error real scaled ='
    [MaxErReal, MaxIndR ] = max( abs( real( Miso(:) ) - real(Mgen(:) ) ) );
    MaxErReal / real( Mgen(MaxIndR) )
'Max error imag scaled= '
    [MaxErImag, MaxIndI ] = max( abs( imag( Miso(:) ) - imag(Mgen(:) ) ) );
    MaxErImag / imag( Mgen(MaxIndI) )
    
    
subplot(2,1,1)
% keyboard
plotyy(1:Nm, [ eigValsGenR, eigValsIsoR ],...
    1:Nm, [ eigValsGenI, eigValsIsoI ])
subplot(2,1,2)
plot(1:Nm,deltaEigR,1:Nm,deltaEigI)
legend('real','imag')

keyboard

if 0
    
    figure()
    pcolor( real(Miso) - real(Mgen) )
    'error real ='
    [MaxErReal, MaxIndR ] = max( abs( real( Miso(:) ) - real(Mgen(:) ) ) );
    MaxErReal / real( Mgen(MaxIndR))
   
    figure()
    pcolor( imag(Miso) - imag(Mgen) )
    'error imag= '
    [MaxErImag, MaxIndI ] = max( abs( imag( Miso(:) ) - imag(Mgen(:) ) ) );
    MaxErImag / imag( Mgen(MaxIndI))
    
    keyboard
    pcolor( real(Mint) )
    pcolor( real(Mint) )
    IndTemp = 15;
    IndTempR = 15;
    IndTempC = 15;
    Miso(IndTempR,IndTempC) - Mgen(IndTempR,IndTempC);
    Miso(IndTempC,IndTempR) - Mgen(IndTempC,IndTempR);
    Mgen(IndTempR,IndTempR);
    Miso(IndTempR,IndTempC);
    Mint(IndTempR,IndTempC) - IntDiagContrib(IndTempC);

end
% keyboard
disp('New meth');disp(max( real( eigValsGen ) ) );


if 0

Nc = 10;
[Coeff_best, ~] = CoeffCalcExpCos2D(Nc,GridObj.phi,ParamObj.bcE); % Calculate coeff
fE = DistBuilderExpCos2Dsing(Nc,GridObj.phi,Coeff_best);        % Build equil distribution
fE_FT = fftshift( fft( fE ) );
[Coeff_best, ~] = CoeffCalcExpCos2D(Nc,GridObj.phi,ParamObj.bcP); % Calculate coeff
fP = DistBuilderExpCos2Dsing(Nc,GridObj.phi,Coeff_best);        % Build equil distr
fP_FT = fftshift( fft( fP ) );    
figure()

subplot(2,1,1)
plot(GridObj.phi,fE, GridObj.phi,fP )
legend('Equilbrium', 'Perturbation','location','best')
subplot(2,1,2)
plot(1:Nm, real(fE_FT), 1:Nm, real(fP_FT) )
legend('Equilbrium', 'Perturbation','location','best')


% disp( max(real(eigVals) ) )
% keyboard
figure()
subplot(2,1,1)
% plot(1:Nm, sort( real( eigVals ) ), 1:Nm,sort(real(eigValsGen)) )
plot(1:Nm,sort(real(eigValsGen)) )
% legend('Using constant', 'Using Eq. Dist','location','best')

title('Real')
subplot(2,1,2)
% plot(1:Nm, sort( imag( eigVals ) ), 1:Nm,sort(imag(eigValsGen)) )
plot(1:Nm,sort(imag(eigValsGen)) )
% legend('Using constant', 'Using Eq. Dist','location','best')
title('Imag')
end

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
