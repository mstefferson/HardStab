%%
kmTemp = 2;
Perturb = ParamObj.vD^2/(2*(1-IntPreFac*Fm_FT(kxHolder,kyHolder,Nm/2+3)))* ...
    (GridObj.kx(kxHolder)^2 + GridObj.ky(kyHolder)^2) / (1-4*kmTemp^2)

(1-IntPreFac*Fm_FT(kxHolder,kyHolder,Nm/2+3))
IntPreFac*Fm_FT(kxHolder,kyHolder,Nm/2+3)
PreFac = -( DiffMobObj.D_par + DiffMobObj.D_perp) / 2 .* ...
    ( GridObj.kx(kxHolder) .^ 2 + GridObj.ky(kyHolder) .^ 2 )...
    - DiffMobObj.D_rot * kmTemp .^2
EigNoDrive = PreFac * (1-IntPreFac*Fm_FT(kxHolder,kyHolder,Nm/2+3));
EigNoDrive + Perturb
fprintf('Max eig w/ drive = %.7e \n', max( real( eig(M) )  ) );
fprintf('Max eig no drive = %.7e \n', max( real( MdiagVec )  ) );
fprintf('Drive - No Drive = %.7e \n', ....
    max( real( eig(M) )  ) -max( real( MdiagVec )  )  )