% Gaussian Blurrer
%% 1D

N     = 10000;
Nint  = 3000;
FmBlurred = zeros(1,N);
Fm = zeros(1,N);
Fm( 1:Nint ) = - 1;
Fm( 2*Nint:end ) = - 1;

Lbox   = 10;
lambda = 0.1;

x = linspace(0,Lbox,N);

dx = x(2)-x(1);
Gauss = 1/( lambda * sqrt(2 * pi) ) .* exp( -x.^2 ./ (2* lambda^2) ) ...
    +1/( lambda * sqrt(2 * pi) ) .* exp( -(x-Lbox).^2 ./ (2* lambda^2) );
trapz(x,Gauss)
trapz_periodic(x,Gauss)

for i = 1:N
    for ik = 1:N
        delta = i-ik+1;
        if delta <= 0
            delta = delta + N; 
        end
        FmBlurred(i) = FmBlurred(i) + Fm(ik).* Gauss(delta) * dx;
    end
end
figure()
plot(x,Fm,x,FmBlurred)
legend('step','blurred step')
Ax = gca;
Ax.YLim = [-1 1];
Ax.LineWidth = 0.5;

%%
%% 2D

N     = 100;
FmBlurred = zeros(N,N);
Fm = zeros(N,N);
Fm(1:ceil(N/4),1:ceil(N/4) ) = -1;

Lbox   = 10;
lambda = 0.1;

x = linspace(0,10,N);
[x2d,y2d] = meshgrid(x);
dx = x(2)-x(1);

Gauss = ( 1/ ( lambda * sqrt(2 * pi) ) ) .^2 .* (...
    exp( -( (x2d).^2 + (y2d).^2  ) ./ (2* lambda^2) ) + ...
    exp( -( (x2d).^2 + (y2d-Lbox).^2  ) ./ (2* lambda^2) ) +...
    exp( -( (x2d-Lbox).^2 + (y2d).^2  ) ./ (2* lambda^2) ) +...
    exp( -( (x2d-Lbox).^2 + (y2d-Lbox).^2  ) ./ (2* lambda^2) ) ) ;

trapz(trapz(x,Gauss,1),x,2)

for i = 1:N
    for j = 1:N
        for ik = 1:N
            deltaI = i-ik+1;
            
            if deltaI <= 0
                deltaI = deltaI + N;
            end
            
            for jk = 1:N                
                deltaJ = j-jk+1;
  
                if deltaJ <= 0
                    deltaJ = deltaJ + N;
                end
                
                FmBlurred(i,j) = FmBlurred(i,j) + ...
                    Fm(ik,jk).* Gauss(deltaI,deltaJ) * dx;
            end
        end
    end
end
 
figure()
subplot(2,1,1)
pcolor(Fm)
subplot(2,1,2)
pcolor(FmBlurred)

%% 3D

N         = 32;
FmBlurred = zeros(N,N,N);
Fm = zeros(N,N,N);
Fm(1:ceil(N/4),1:ceil(N/4),1:ceil(N/4)  ) = -1;

Lbox   = 10;
lambda = 0.1;

x = linspace(0,10,N);
y = x;
phi = linspace(0,2*pi,N);

[x2d,y2d,phi3d] = meshgrid(x,y,phi);
dx = x(2)-x(1);

Gauss = 1/( lambda * sqrt(2 * pi) ) .*...
    exp( -( (x3d).^2 + (y3d).^2  ) ./ (2* lambda^2) ) + ...
    1/( lambda * sqrt(2 * pi) ) .*...
    exp( -( (x2d-Lbox).^2 + (y2d-Lbox).^2  ) ./ (2* lambda^2) );


for i = 1:N
    for j = 1:N
        for ik = 1:N
            deltaI = i-ik+1;
            
            if deltaI <= 0
                deltaI = deltaI + N;
            end
            
            for jk = 1:N                
                deltaJ = j-jk+1;
  
                if deltaJ <= 0
                    deltaJ = deltaJ + N;
                end
                
                FmBlurred(i,j) = FmBlurred(i,j) + ...
                    Fm(ik,jk).* Gauss(deltaI,deltaJ) * dx;
            end
        end
    end
end
 
figure()
subplot(2,1,1)
pcolor(Fm)
subplot(2,1,2)
pcolor(FmBlurred)

