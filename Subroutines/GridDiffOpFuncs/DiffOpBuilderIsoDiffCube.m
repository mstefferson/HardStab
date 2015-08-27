% Build the diffusion operator

function [Lop_kcube] = DiffOpBuilderIsoDiffCube(DiffMobObj,GridObj)

%%%%%%%%%%%%%%%%%%Diagonal operator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lop_kcube = -( DiffMobObj.D_pos .* ( GridObj.kx3D.^2 + GridObj.ky3D.^2 )...            
             + DiffMobObj.D_rot .* GridObj.km3D.^2 );

end
