%%
% Grid
Nx = 1000;
Lbox = 10;
x  = linspace(-Lbox/2,Lbox/2,Nx);

% LJ stuff
epsilon = 0.1;
sigma = 1;
rm = 2^(1/6) * sigma;

LJv = epsilon .* ( (rm ./ abs(x) ) .^ (12) - 2 .* (rm ./ abs(x) ) .^ (6) );
HSv = zeros(1,Nx);
HSv( Nx/2 - floor(Nx * sigma / Lbox):Nx/2 + floor(Nx * sigma / Lbox) ) = 10000;

% may fncs
FmLJ = exp(-LJv) - 1;
FmHS = exp(-HSv) - 1;
% Blurred Mayer fnc
lambda     = 0.1;
FmBlurred =  zeros(1,Nx);
dx   = x(2) - x(1);
Gauss = 1/( lambda * sqrt(2 * pi) ) .* exp( -(x-Lbox/2).^2 ./ (2* lambda^2) ) ...
    +1/( lambda * sqrt(2 * pi) ) .* exp( -(x+Lbox/2).^2 ./ (2* lambda^2) );

for i = 1:Nx
    for ik = 1:Nx
        delta = i-ik+1;
        if delta <= 0
            delta = delta + Nx; 
        end
        FmBlurred(i) = FmBlurred(i) + FmHS(ik).* Gauss(delta) * dx;
    end
end

figure()
subplot(2,1,1)
plot(x,HSv,x,LJv)
Ax = gca;
Ax.YLim = [-1 2];


subplot(2,1,2)
plot(x,FmHS,x,FmLJ,x,FmBlurred)
Ax = gca;
Ax.YLim = [-1 2];
