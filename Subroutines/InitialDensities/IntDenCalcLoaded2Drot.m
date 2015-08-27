function [ rho ] = IntDenCalcLoaded2Drot( VarNmString )
%INTDEN2DROTCALCLOADED Initial density is loaded from a saved rho
rho = 0; 
% keyboard
Name2Load = ...
    sprintf('C:/Users/MWS/Documents/MATLAB/research/Bg/diffFT/hRddft/Subroutines/InitialDensities/SavedRhos/%s',...
             VarNmString);
load(Name2Load)

end

