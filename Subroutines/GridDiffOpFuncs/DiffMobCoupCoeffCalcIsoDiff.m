% Calculates the mobility and diffusion tensor based on dimensions of rod
% and fluid

function [DiffMobObj] = DiffMobCoupCoeffCalcIsoDiff(T,Mob_pos,Mob_rot)


% Use Einstein diffusion relations
D_pos = Mob_pos * T; 
D_rot = Mob_rot * T; 

% Make Diff object
DiffMobObj = struct('Mob_pos', Mob_pos,'D_pos',D_pos,'Mob_rot', Mob_rot,'D_rot',D_rot);
end