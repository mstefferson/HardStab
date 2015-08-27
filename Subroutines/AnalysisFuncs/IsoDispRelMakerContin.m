function [lambda] = IsoDispRelMakerContin(kx,ky,km,D_pos,D_rot,bc)

lambda = - ( D_pos .* ( kx.^2 + ky.^2 )  ...
           + D_rot .* km.^2 ) .* ...
             ( 1 -  ( 1 + (-1).^( km ) ) .* ...
               bc ./ ( km.^2 - 1 + 0.0000001 ) );

end