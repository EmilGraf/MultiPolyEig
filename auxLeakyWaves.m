function sols = auxLeakyWaves(L2,L1,L0,M,R1,R2,k1,k2,mu)

% Solves an instance of the 3PEP in eq. (25) in the paper E Graf and A Townsend:
% A Hidden Variable Resultant Method for the Polynomial Multiparameter
% Eigenvalue Problem, 2025.
%
% Based on eig_LeakySolid from RandomJointEig, https://github.com/borplestenjak/RandomJointEig,
% Bor Plestenjak, Hauke Gravenkamp, Daniel Kiefer 2024.
%
% Output
%  - sols: [k beta eta]

% E Graf and A Townsend, 2025

% Set up in monomial basis

% A{4,1} = -(L0 + mu^2*M - k1^2*R3 -k2^2*R4);
A{4,1} = -(L0 + mu^2*M);
A{4,2} =  L1;
A{4,3} =  R1;
A{4,4} =  R2;
A{4,5} =  L2;

% xi1^2 = -xi3*(k1^2+xi3), xi3 = lambda^2
A{2,1} = -[0 -k1^2 ;0 0];
A{2,2} = zeros(2);
A{2,3} = eye(2);
A{2,4} = zeros(2);
A{2,5} = [0 -1;1 0];

% xi2^2 = -xi3*(k2^2+xi3), xi3 = lambda^2
A{3,1} = -[0 -k2^2;0 0];
A{3,2} = zeros(2);
A{3,3} = zeros(2);
A{3,4} = eye(2);
A{3,5} = [0 -1;1 0];

% Convert to Chebyshev

F1 = zeros([size(L1) 2 2 3]);
F2 = zeros(2,2,2,2,3);
F3 = zeros(2,2,2,2,3);

F1(:,:,1,1,1) = -A{4,1}+.5*A{4,5};
F1(:,:,2,1,1) = A{4,3};
F1(:,:,1,2,1) = A{4,4};
F1(:,:,1,1,2) = A{4,2};
F1(:,:,1,1,3) = .5*A{4,5};

F2(:,:,1,1,1) = -A{3,1}+.5*A{3,5};
F2(:,:,2,1,1) = A{3,3};
F2(:,:,1,2,1) = A{3,4};
F2(:,:,1,1,2) = A{3,2};
F2(:,:,1,1,3) = .5*A{3,5};

F3(:,:,1,1,1) = -A{2,1}+.5*A{2,5};
F3(:,:,2,1,1) = A{2,3};
F3(:,:,1,2,1) = A{2,4};
F3(:,:,1,1,2) = A{2,2};
F3(:,:,1,1,3) = .5*A{2,5};

evblockv = [0 40 20];
options.res = 100;

F1t = F1;
F2t = F2;
F3t = F3;

% Solve

for ii = 1:3
options.evblock = evblockv(ii);

if ii > 1
    %%%Random shift
    [Q,~] = qr(randn(3));
    F1t(:,:,2,1,1) = Q(1,1)*F1(:,:,2,1,1) + Q(1,2)*F1(:,:,1,2,1) + Q(1,3)*F1(:,:,1,1,2);
    F1t(:,:,1,2,1) = Q(2,1)*F1(:,:,2,1,1) + Q(2,2)*F1(:,:,1,2,1) + Q(2,3)*F1(:,:,1,1,2);
    F1t(:,:,1,1,2) = Q(3,1)*F1(:,:,2,1,1) + Q(3,2)*F1(:,:,1,2,1) + Q(3,3)*F1(:,:,1,1,2);
    F1t(:,:,3,1,1) = Q(1,3)^2*F1(:,:,1,1,3);
    F1t(:,:,1,3,1) = Q(2,3)^2*F1(:,:,1,1,3);
    F1t(:,:,1,1,3) = Q(3,3)^2*F1(:,:,1,1,3);
    F1t(:,:,1,2,2) = 4*Q(2,3)*Q(3,3)*F1(:,:,1,1,3);
    F1t(:,:,2,1,2) = 4*Q(1,3)*Q(3,3)*F1(:,:,1,1,3);
    F1t(:,:,2,2,1) = 4*Q(1,3)*Q(2,3)*F1(:,:,1,1,3);
    
    F2t(:,:,2,1,1) = Q(1,1)*F2(:,:,2,1,1) + Q(1,2)*F2(:,:,1,2,1) + Q(1,3)*F2(:,:,1,1,2);
    F2t(:,:,1,2,1) = Q(2,1)*F2(:,:,2,1,1) + Q(2,2)*F2(:,:,1,2,1) + Q(2,3)*F2(:,:,1,1,2);
    F2t(:,:,1,1,2) = Q(3,1)*F2(:,:,2,1,1) + Q(3,2)*F2(:,:,1,2,1) + Q(3,3)*F2(:,:,1,1,2);
    F2t(:,:,3,1,1) = Q(1,3)^2*F2(:,:,1,1,3);
    F2t(:,:,1,3,1) = Q(2,3)^2*F2(:,:,1,1,3);
    F2t(:,:,1,1,3) = Q(3,3)^2*F2(:,:,1,1,3);
    F2t(:,:,1,2,2) = 4*Q(2,3)*Q(3,3)*F2(:,:,1,1,3);
    F2t(:,:,2,1,2) = 4*Q(1,3)*Q(3,3)*F2(:,:,1,1,3);
    F2t(:,:,2,2,1) = 4*Q(1,3)*Q(2,3)*F2(:,:,1,1,3);
    
    F3t(:,:,2,1,1) = Q(1,1)*F3(:,:,2,1,1) + Q(1,2)*F3(:,:,1,2,1) + Q(1,3)*F3(:,:,1,1,2);
    F3t(:,:,1,2,1) = Q(2,1)*F3(:,:,2,1,1) + Q(2,2)*F3(:,:,1,2,1) + Q(2,3)*F3(:,:,1,1,2);
    F3t(:,:,1,1,2) = Q(3,1)*F3(:,:,2,1,1) + Q(3,2)*F3(:,:,1,2,1) + Q(3,3)*F3(:,:,1,1,2);
    F3t(:,:,3,1,1) = Q(1,3)^2*F3(:,:,1,1,3);
    F3t(:,:,1,3,1) = Q(2,3)^2*F3(:,:,1,1,3);
    F3t(:,:,1,1,3) = Q(3,3)^2*F3(:,:,1,1,3);
    F3t(:,:,1,2,2) = 4*Q(2,3)*Q(3,3)*F3(:,:,1,1,3);
    F3t(:,:,2,1,2) = 4*Q(1,3)*Q(3,3)*F3(:,:,1,1,3);
    F3t(:,:,2,2,1) = 4*Q(1,3)*Q(2,3)*F3(:,:,1,1,3);
end

Ft = {F1t;F2t;F3t};

% Get candidate roots
r = multipolyeig(Ft,options);

if ii > 1
    r = r*Q;
end

if ii == 1
    lambda = r;
else
    lambda = [lambda; r];
end
end

k = lambda(:,3)/1i;
sols = [k lambda(:,1)./k lambda(:,2)./k];


% take sols with beta and eta close to sqrt(k1^2-k^2) and sqrt(k1^2-k^2)
    test_sr = abs(log(abs(sqrt(k1^2-sols(:,1).^2)./sols(:,2)))) + abs(log(abs(sqrt(k2^2-sols(:,1).^2)./sols(:,3))));
    ind = find(test_sr<1e-2);

    if ~isempty(ind)
        sols = sols(ind,:);
    else
        sols = sols(1,:)*nan;
    end
