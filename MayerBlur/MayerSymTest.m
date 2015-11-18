
ModeX    = 1;
ModeY    = 0;

vD    = 0;
bcE   = 1.4; % Equilbirum bc
bcP   = 1.4; % Perturbation bc

Nx    = 2^7;
Ny    = 2^7;
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
D_perp = Mob_perp;
D_rot = Mob_rot;


DiffMobObj = struct('Mob_par', Mob_par,'D_par',D_par,'D_perp',D_perp,...
    'Mob_rot', Mob_rot,'D_rot',D_rot);

ParamObj = struct('Nx',Nx,'Ny',Ny,'Nm',Nm,'Lx',Lx,'Ly',Ly,'L_rod',L_rod,...
    'bcE',bcE,'bcP',bcP,'vD',vD);


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



Fm = MayerFncDiffBtwPntsCalc(...
        ParamObj.Nx, ParamObj.Ny, ParamObj.Nm, ...
        ParamObj.Lx, ParamObj.Ly, GridObj.dx,...
        GridObj.dy, GridObj.dphi, ParamObj.L_rod);
    
%     %%
%     Mind  = 19;
%     Mind2 = Mind + Nm/2;
%     pcolor( reshape( Fm(:,:,Mind2), Nx, Ny )' )
    
    %%
    Fm_FT = fftshift(fftn( Fm ) );

    mInd = Nm/2;
    Sum = 0;
    for i = 1:Nx
        for  j = 1:Ny
       Sum = Sum + Fm_FT(i,j,mInd);     
        end
    end     
   Sum;
    
     xMode = 3;
     yMode = 6;  
     mMode = 1;
     
     mInd  = Nm/2+1 + mMode;
     kxInd = Nx/2+1 + xMode;
     kyInd = Ny/2+1 + yMode;
     
     mInd2  = Nm/2+1 - mMode;
     kxInd2 = kxInd;
     kyInd2 = kyInd;

     
     GridObj.km(mInd)
     
     Fm_FT(kxInd,kyInd,mInd)
     Fm_FT(kxInd2,kyInd2,mInd2)
     %%
     xMode = 0;
     yMode = 0;  
     
     kxInd = Nx/2+1 + xMode;
     kyInd = Ny/2+1 + yMode;
     
     figure     
     
     plot( 1:Nm, real( reshape( Fm_FT(kxInd,kyInd,:) , Nm, 1 ) ),...
           1:Nm, imag( reshape( Fm_FT(kxInd,kyInd,:) , Nm, 1 ) ))
    legend('real','imag')
       plot( 1:Nm, imag( reshape( Fm_FT(kxInd,kyInd,:) , Nm, 1 ) ))
       
%%
xInd = 1;
yInd = 2;
FmFixedPos =  reshape( Fm(xInd,yInd,:), Nm,1 );

plot(phi,FmFixedPos)

FmFixedPos_FT = fftshift(fft(FmFixedPos));

FmFixedPos_FT(:);

%%
Fm_FTtest = Fm_FT;

for  i =  1:Nx
    for j = 1:Ny
        for k = 1:Nm
            
            if mod(k,2) ~= 0
                FM_FTtest = 0;
            end
        end
    end
end

FmTest = ifftn(ifftshift( Fm_FTtest ) );

deltaFM = Fm - FmTest;

 %%

 
 Sum = 0;
 Int = 1; 
 for i = 0:N-1
     Sum =  Sum + exp( -2*pi*sqrt(-1)*i/N * Int);
 end
 
 Sum
                
    
    
