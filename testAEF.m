% Commputes table 2 in the paper E Graf and A Townsend:
% A Hidden Variable Resultant Method for the Polynomial Multiparameter
% Eigenvalue Problem, 2025.
% 
% Test for multipolyeig on AEF Problem from
% A. Pons, S. Gutschmidt, Multiparameter spectral analysis for aeroelastic instability problems, 
% J. Appl. Mech. 85 (6) (2018) 061011. doi:10.1115/ 1.4039671.
%
% Uses quad_twopareig from B. Plestenjak, Multipareig, MATLAB Central File Exchange.
% https://www.mathworks.com/matlabcentral/fileexchange/ 47844-multipareig.
%
% E Graf and A Townsend, 2025.

clearvars;
rng(1);

%%%Matrix dimensions
nf = 2;
ng = 2;

%%%Polynomial degrees in each variable (+1)
mx = 2;
my = 3;

%%%Raw data
mu = 20;
r = .4899;
omegah = .5642;
omegatheta = 1.4105;
zetah = 1.4105;
zetatheta = 2.3508;
rtheta = -.1;
a = -.2;
G0 = 1/mu.*[1 a; a 1/8+a^2];
G1 = 1/mu.*[-2*i 2*i*(1-a); -i*(1+2*a) i*a*(1-2*a)];
G2 = 1/mu.*[0 2; 0 1+2*a];
M0 = [1 -rtheta; -rtheta r^2];
D0 = [2*i*zetah*omegah 0; 0 2*i*r^2*zetatheta*omegatheta];
K0 = [omegah^2 0; 0 r^2*omegatheta^2];

%%%Set up input for multipolyeig
F = zeros(2,2,2,3);
F(:,:,2,1) = K0;
F(:,:,1,1) = M0+G0+1/2*G2;
F(:,:,1,3) = 1/2*G2;
F(:,:,1,2) = G1;

%Solve
%%%Random transformation
[Q,~] = qr(randn(2));
FF = F;
FF(:,:,1,3) = Q(1,1)^2*F(:,:,1,3);
FF(:,:,3,1) = Q(2,1)^2*F(:,:,1,3);
FF(:,:,1,1) = F(:,:,1,1) + (Q(2,1)^2+Q(1,1)^2-1)*F(:,:,1,3);
FF(:,:,2,2) = 4*Q(1,1)*Q(2,1)*F(:,:,1,3);
FF(:,:,1,2) = Q(1,1)*F(:,:,1,2) + Q(1,2)*F(:,:,2,1);
FF(:,:,2,1) = Q(2,1)*F(:,:,1,2) + Q(2,2)*F(:,:,2,1);

GG = conj(FF);

roots = multipolyeig({FF;GG});

%%%Transform back
roots = [-roots(:,1) roots(:,2)]*Q;

x = roots(:,1);
y = roots(:,2);


%%Now solve with multipareig
G = conj(F);
A1 = F(:,:,1,1) - F(:,:,1,3);
B1 = F(:,:,2,1);
C1 = F(:,:,1,2);
D1 = zeros(2,2);
E1 = F(:,:,2,2);
F1 = 2*F(:,:,1,3);
A2 = G(:,:,1,1) - G(:,:,1,3);
B2 = G(:,:,2,1);
C2 = G(:,:,1,2);
D2 = zeros(2,2);
E2 = G(:,:,2,2);
F2 = 2*G(:,:,1,3);


%%%
opts.refine = 0;
opts.fast=1;
opts.singular = 1;
[x2,y2] = quad_twopareig(A1,B1,C1,D1,E1,F1,A2,B2,C2,D2,E2,F2,opts);

%%%Check
disp("multipolyeig:");
for i = 1:size(y)
    disp(strcat("Root ", string(i), ": (", num2str(x(i),15), ", ", num2str(y(i),15), ")"));
end
res = 0;
for i = 1:size(y)
    %%%Set Up Vandermonde
    vx = zeros(mx,1);
    vy = zeros(my,1);
    for j = 1:mx
        vx(j) = cos((j-1)*acos(x(i)));
    end
    for j = 1:my
        vy(j) = cos((j-1)*acos(y(i)));
    end
    
    Ft = zeros(nf,nf);
    Gt = zeros(ng,ng);
    for k = 1:mx
        for l = 1:my
            Ft = Ft + vx(k)*vy(l)*F(:,:,k,l);
            Gt = Gt + vx(k)*vy(l)*G(:,:,k,l);
        end
    end
    disp(strcat("Residual ",string(i),"; min(svd(F)): ",string(min(svd(Ft))),"; min(svd((G))): ",string(min(svd((Gt))))));
    res = res + min(svd(Ft)) + min(svd(Gt));
end
disp(strcat("Average residual: ", string(res/(2*size(y,1)))));

disp("multipareig:");
for i = 1:size(y2)
    disp(strcat("Root ", string(i), ": (", num2str(x2(i),15), ", ", num2str(y2(i),15), ")"));
end
res = 0;
for i = 1:size(y2)
    %%%Set Up Vandermonde
    vx = zeros(mx,1);
    vy = zeros(my,1);
    for j = 1:mx
        vx(j) = cos((j-1)*acos(x2(i)));
    end
    for j = 1:my
        vy(j) = cos((j-1)*acos(y2(i)));
    end
    
    Ft = zeros(nf,nf);
    Gt = zeros(ng,ng);
    for k = 1:mx
        for l = 1:my
            Ft = Ft + vx(k)*vy(l)*F(:,:,k,l);
            Gt = Gt + vx(k)*vy(l)*G(:,:,k,l);
        end
    end
    disp(strcat("Residual ",string(i),"; min(svd(F)): ",string(min(svd(Ft))),"; min(svd((G))): ",string(min(svd((Gt))))));
    res = res + min(svd(Ft)) + min(svd(Gt));
end
disp(strcat("Average residual: ", string(res/(2*size(y2,1)))));

