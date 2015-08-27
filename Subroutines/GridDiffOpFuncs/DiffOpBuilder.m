% Build the diffusion operator

function [Lop] = DiffOpBuilder(DiffMobObj,GridObj,Nm,N2,N3)

%%%%%%%%%%%%%%%%%%Diagonal operator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lop_kcube = -( ( (DiffMobObj.D_par + DiffMobObj.D_perp)./ 2 ) .* (GridObj.kx3D.^2 + GridObj.ky3D.^2) + DiffMobObj.D_rot .* GridObj.km3D.^2 );

%Diagonal matrix part of the operator (no interactions)
Lop_DiagMtx = spdiags( reshape( Lop_kcube, N3, 1 ), 0, N3, N3 );

%%%%%%%%%%%%%%%%%Off diagonal mode coupling%%%%%%%%%%%%%%%%%%%%%%%%%%
%Handle the cross terms of the PDE
CpMplus2Vec  = reshape( repmat( DiffMobObj.CfMplus2,  [1 1 Nm]), 1, N3 );
CpMminus2Vec = reshape( repmat( DiffMobObj.CfMminus2, [1 1 Nm]), 1, N3 );

% Put them in the matrix
Mplus2Mtx  = spdiags( [ zeros(1,2*N2) CpMplus2Vec( 1:(N2*(Nm-2)) ) ]', 2*N2, N3, N3 );
Mminus2Mtx = spdiags( [ CpMminus2Vec(1:(N2*(Nm-2))) zeros(1,2*N2) ]' , -2*N2, N3, N3 );

%Periodic terms
Mplus2Mtx =  Mplus2Mtx  + spdiags( [ CpMplus2Vec( N2*(Nm-2)+1:N3 ) zeros( 1, N2*(Nm-2) ) ]' , -(N2*(Nm-2)), N3, N3 ); %m = N coupling to m = 2
Mminus2Mtx = Mminus2Mtx + spdiags( [ zeros( 1, N2*(Nm-2) ) CpMminus2Vec( 1:2*N2 ) ]', N2*(Nm-2), N3, N3 ) ;                  %m = 1 coupling to m = N-1

%%%%%%%%%%%%%%%%%%%%%%Non-interacting Propagator%%%%%%%%%%%%%%%%%%%%%%%%%%
Lop = Lop_DiagMtx + Mplus2Mtx + Mminus2Mtx;

end
