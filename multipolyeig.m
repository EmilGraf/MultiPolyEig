function roots = multipolyeig(F,options)
%%% Global Solver for Polynomial Multiparameter Eigenvalue Problems
%%% For d=1,2,3, the algorithm returns the common eigenvalues of F = cell(F1,...,Fd),
%%% where Fi are d-parameter polynomial matrices, given as d+2 dimensional
%%% tensors where the last d dimensions correspond to variables
%%% and the first two to matrix dimensions.
%%% The algorithm is based on the hidden variable Dixon resultant method from
%%% E. Graf and A. Townsend, "A Hidden Variable Resultant Method for the
%%% Polynomial Multiparameter Eigenvalue Problem," 2025.

%%% Options:
%%% evblock: 1/evblock indicates the percentage of eigenvector entry ratios to consider
%%% when calculating cordinates x_1,...,x_{d-1}
%%% res: tolerance in final residual check

if nargin < 2
    options = [];
end

if isfield (options,'evblock'), evblock = options.evblock; else, evblock = 20; end
if isfield (options,'res'), res = options.res; else, res = 1e-8; end


d = size(F,1);

%%%Extract dimensions, last d columns give degree in each variable, first 2
%%%give matrix dimensions (which must match)
dims = zeros(d,d+2);
for i=1:d
    dims(i,:) = size(F{i});
end

%%%Get candidate roots
roots = runPolyEig(F,evblock);

%%%Residual Check
r = NaN*zeros(size(roots,1),1);
v = cell(d,1);
for j = 1:size(roots,1)
    for i = 1:d
        v{i} = flip(vandercheb(roots(j,i),dims(1,2+i)));
    end
    for i = 1:d
        Fi = F{i};
        for k = 1:d
            Fi = reshape(reshape(Fi,[],dims(1,d+3-k))*v{d+1-k},dims(i,1:d+2-k));
        end
        r(j) = max(r(j),min(svd(Fi)));
    end
end
ind = r < res;
roots = roots(ind,:);
end

function roots = runPolyEig(F,evblock)
%%%Run the eigenvalue algorithm

d = size(F,1);
%%%Extract dimensions, last d columns give degree in each variable, first 2
%%%give matrix dimensions (which must match)
dims = zeros(d,d+2);
for i=1:d
    dims(i,:) = size(F{i});
end

if min(dims(1,3:end)) < 3
    roots = runPolyEiglin(F);
    return;
end

if d==1
    roots = singChebPEP(flip(F{1},3));
else
%%%Construct Hidden Variable Resultant
R = tensorDixon(F,d,dims);

%%%Solve via colleague linearization
[lambda,V] = singChebPEP(flip(R,3));

%%%Extract roots
roots = zeros(size(lambda,1),d);
roots(:,d) = lambda;

%%%Detect generic null space
    v = flip(vandercheb(randn,size(R,3)));
    %%%Evaluate
    a = size(R,1);
    T = reshape(reshape(R,[],size(R,3))*v,a,a);
    [~,S,W] = svd(T);
    idx = find(abs(diag(S)/S(1,1))>1e-13,1,'last');
    N = W(:,idx+1:end);

%%%d=2
if d==2
if dims(1,3) < 3
    for i = 1:size(lambda)
        if (isfinite(lambda(i)))
            v = flip(vandercheb(lambda(i),dims(1,d+2)));
        
            %%%Evaluate
            F1 = reshape(reshape(F{1},[],dims(1,d+2))*v,dims(1,1:d+1));
            F2 = reshape(reshape(F{2},[],dims(1,d+2))*v,dims(2,1:d+1));
        
            %%%Candidate roots
            r1 = singChebPEP(flip(F1,3));
            r2 = singChebPEP(flip(F2,3));
        
            %%%Closest
            if size(r1,1) > 0 && size(r2,1) > 0
                T = r1-r2.';
                [~,I] = min(abs(T),[],'all');
                [a,b] = ind2sub(size(T),I);
                roots(i,1) = mean([r1(a,:);r2(b,:)]);
            end
        end
    end
else
    j = prod(dims(:,1),'all');
    Vt = V(end-j+1:end,:);
    Vt(abs(Vt) < 1e-13) = nan;
    tmp = V(end-2*j+1:end-j,:)./Vt;
    VV1 = max(abs(N(end-j+1:end,:)),[],2);
    VV2 = max(abs(N(end-2*j+1:end-j,:)),[],2);
    if size(N,2) > 0
        ind = VV1 < 1e-13 & VV2 < 1e-13;
    else
        ind = ones(size(tmp,1));
    end

    tmp = tmp(ind,:);
    Vt = Vt(ind,:);
    [~,ind2] = max(abs(Vt),[],1);
    if size(ind2,1) > 0
        ind2 = sub2ind(size(tmp),ind2,1:size(Vt,2));
        roots(:,1) = tmp(ind2);
    end
end
end

if d==3
    if dims(1,3) == 2
        %%%Linear case, we can not extract from eigenvector
        for i = 1:size(lambda,1)
            %%%Evaluate
            v = flip(vandercheb(lambda(i),dims(1,5)));
            F1 = permute(reshape(reshape(F{1},[],dims(1,5))*v,dims(1,1:4)),[1 2 4 3]);
            F2 = permute(reshape(reshape(F{2},[],dims(2,5))*v,dims(2,1:4)),[1 2 4 3]);
            F3 = permute(reshape(reshape(F{3},[],dims(3,5))*v,dims(3,1:4)),[1 2 4 3]);
        
            %%%Candidiate roots
            R = tensorDixon({F1,F2},2,dims(1:2,[1 2 4 3]));
            r1 = singChebPEP(flip(R,3));
            R = tensorDixon({F2,F3},2,dims(2:3,[1 2 4 3]));
            r2 = singChebPEP(flip(R,3)).';
            R = tensorDixon({F1,F3},2,dims([1 3],[1 2 4 3]));
            tmp = singChebPEP(flip(R,3));
            r3 = zeros(1,1,size(tmp,1));
            r3(1,1,:) = tmp;
        
            %%%Select closest
            if size(r1,1) > 0 && size(r2,1) > 0 && size(r3,1) > 0
                T = abs(r1-r2)+abs(r2-r3)+abs(r1-r3);
                [~,I] = min(T,[],'all');
                [a,b,c] = ind2sub(size(T),I);
                roots(i,1) = median([r1(a,:);r2(:,b);r3(:,:,c)]);
            end
        end
    else
        j = prod(dims(:,1),'all');
        Vt = V(end-j+1:end,:);
        %Vt(abs(Vt) < 1e-6) = nan;
        tmp = V(end-2*j+1:end-j,:)./Vt;
        VV1 = max(abs(N(end-j+1:end,:)),[],2);
        VV2 = max(abs(N(end-2*j+1:end-j,:)),[],2);
        if size(N,2) > 0
            ind = VV1 < 1e-11 & VV2 < 1e-11;
        else
            ind = ones(size(tmp,1));
        end
        tmp = tmp(ind,:);
        Vt = Vt(ind,:);
        k = ceil(size(Vt,1)/evblock);
        [~,ind2] = maxk(abs(Vt),k,1);
        if size(ind2,1) > 0
            ind2 = sub2ind(size(tmp),ind2,repmat(1:size(Vt,2),k,1));
            tmp = tmp(ind2);
            roots(:,1) = mean(tmp,1,'omitnan');
        end
    end

    if dims(1,4) == 2
        %%%Linear case, we can not extract from eigenvector
        for i = 1:size(lambda,1)
            %%%Evaluate
            v = flip(vandercheb(lambda(i),dims(1,5)));
            F1 = reshape(reshape(F{1},[],dims(1,5))*v,dims(1,1:4));
            F2 = reshape(reshape(F{2},[],dims(2,5))*v,dims(2,1:4));
            F3 = reshape(reshape(F{3},[],dims(3,5))*v,dims(3,1:4));
        
            %%%Candidiate roots
            R = tensorDixon({F1,F2},2,dims(1:2,1:4));
            r1 = singChebPEP(flip(R,3));
            R = tensorDixon({F2,F3},2,dims(2:3,1:4));
            r2 = singChebPEP(flip(R,3)).';
            R = tensorDixon({F1,F3},2,dims([1 3],1:4));
            tmp = singChebPEP(flip(R,3));
            r3 = zeros(1,1,size(tmp,1));
            r3(1,1,:) = tmp;
        
            %%%Select closest
            if size(r1,1) > 0 && size(r2,1) > 0 && size(r3,1) > 0
                T = abs(r1-r2)+abs(r2-r3)+abs(r1-r3);
                [~,I] = min(T,[],'all');
                [a,b,c] = ind2sub(size(T),I);
                roots(i,2) = median([r1(a,:);r2(:,b);r3(:,:,c)]);
            end
        end
    else 
        j = prod(dims(:,1),'all')*2*(dims(1,3)-1);
        Vt = V(end-j+1:end,:);
        %Vt(abs(Vt) < 1e-6) = nan;
        tmp = V(end-2*j+1:end-j,:)./Vt;
        VV1 = max(abs(N(end-j+1:end,:)),[],2);
        VV2 = max(abs(N(end-2*j+1:end-j,:)),[],2);
        if size(N,2) > 0
            ind = VV1 < 1e-11 & VV2 < 1e-11;
        else
            ind = ones(size(tmp,1));
        end
        tmp = tmp(ind,:);
        Vt = Vt(ind,:);
        k = ceil(size(Vt,1)/evblock);
        [~,ind2] = maxk(abs(Vt),k,1);
        if size(ind2,1) > 0
            ind2 = sub2ind(size(tmp),ind2,repmat(1:size(Vt,2),k,1));
            tmp = tmp(ind2);
            roots(:,2) = mean(tmp,1,'omitnan');
        end
    end
end
end
end



function roots = runPolyEiglin(F)
%%%Run the eigenvalue algorithm, fast method for partially linear problem

d = size(F,1);
%%%Extract dimensions, last d columns give degree in each variable, first 2
%%%give matrix dimensions (which must match)

dims = zeros(d,d+2);
for i=1:d
    dims(i,:) = size(F{i});
end

if d==1
    roots = singChebPEP(flip(F{1},3));
elseif d==2
%%%d=2
%%%Construct Hidden Variable Resultant
if dims(1,4)==2
    rev = 0;
    R = tensorDixon(F,d,dims);
else
    F{1} = permute(F{1},[1 2 4 3]);
    F{2} = permute(F{2},[1 2 4 3]);
    rev = 1;
    R = tensorDixon(F,d,dims(:,[1 2 4 3]));
end

%%%Solve via colleague linearization
[lambda,V] = singChebPEP(flip(R,3));

%%%Extract roots
roots = zeros(size(lambda,1),d);
roots(:,d) = lambda;

%%%project eigenvectors
    v = flip(vandercheb(randn,size(R,3)));
    %%%Evaluate
    a = size(R,1);
    T = reshape(reshape(R,[],size(R,3))*v,a,a);
    [~,S,W] = svd(T);
    idx = find(abs(diag(S)/S(1,1))>1e-13,1,'last');
    N = W(:,idx+1:end);
    for i = 1:size(V,2)
       for j = 1:size(N,2)
           V(:,i) = V(:,i) - (V(:,i).'*N(:,j))*N(:,j);
       end
    end
    j = prod(dims(:,1),'all');
    tmp = V(end-2*j+1:end-j,:)./V(end-j+1:end,:);
    roots(:,1) = median(tmp,1,"omitnan");


    Vt = V(end-j+1:end,:);
    Vt(abs(Vt) < 1e-13) = nan;
    tmp = V(end-2*j+1:end-j,:)./Vt;
    VV1 = max(abs(N(end-j+1:end,:)),[],2);
    VV2 = max(abs(N(end-2*j+1:end-j,:)),[],2);
    if size(N,2) > 0
        ind = VV1 < 1e-13 & VV2 < 1e-13;
    else
        ind = ones(size(tmp,1));
    end

    tmp = tmp(ind,:);
    Vt = Vt(ind,:);
    [~,ind2] = max(abs(Vt),[],1);
    if size(ind2,1) > 0
        ind2 = sub2ind(size(tmp),ind2,1:size(Vt,2));
        roots(:,1) = tmp(ind2);
    end

if rev
    roots = flip(roots,2);
end

elseif d==3
%%%d==3;
%%%Construct Hidden Variable Resultant
rev = dims(1,3:5) == 2;
numlin = sum(rev);
if isequal(rev,[0 0 1])
    R = tensorDixon(F,d,dims);
elseif isequal(rev,[0 1 0])
    for i = 1:3
        F{i} = permute(F{i},[1 2 3 5 4]);
    end
    R = tensorDixon(F,d,dims(:,[1 2 3 5 4]));
elseif isequal(rev,[1 0 0])
    for i = 1:3
        F{i} = permute(F{i},[1 2 4 5 3]);
    end
    dims = dims(:,[1 2 4 5 3]);
    R = tensorDixon(F,d,dims);
elseif isequal(rev,[0 1 1])
    R1 = tensorDixon(F,d,dims);
    for i = 1:3
        F{i} = permute(F{i},[1 2 3 5 4]);
    end
    R2 = tensorDixon(F,d,dims(:,[1 2 3 5 4]));
elseif isequal(rev,[1 1 0])
    for i = 1:3
        F{i} = permute(F{i},[1 2 5 4 3]);
    end
    R1 = tensorDixon(F,d,dims(:,[1 2 5 4 3]));
    for i = 1:3
        F{i} = permute(F{i},[1 2 3 5 4]);
    end
    R2 = tensorDixon(F,d,dims(:,[1 2 5 3 4]));
elseif isequal(rev,[1 0 1])
    for i = 1:3
        F{i} = permute(F{i},[1 2 4 3 5]);
    end
    R1 = tensorDixon(F,d,dims(:,[1 2 4 3 5]));
    for i = 1:3
        F{i} = permute(F{i},[1 2 3 5 4]);
    end
    R2 = tensorDixon(F,d,dims(:,[1 2 4 5 3]));
end

%%%Solve via colleague linearization
if numlin == 1
    %one linear var
    [lambda,V] = singChebPEP(flip(R,3));

%%%Extract roots
roots = zeros(size(lambda,1),d);
roots(:,d) = lambda;

%%%project eigenvectors
    v = flip(vandercheb(randn,size(R,3)));
    %%%Evaluate
    a = size(R,1);
    T = reshape(reshape(R,[],size(R,3))*v,a,a);
    [~,S,W] = svd(T);
    idx = find(abs(diag(S)/S(1,1))>1e-13,1,'last');
    N = W(:,idx+1:end);
    for i = 1:size(V,2)
       for j = 1:size(N,2)
           V(:,i) = V(:,i) - (V(:,i).'*N(:,j))*N(:,j);
       end
    end
        j = prod(dims(:,1),'all');
        Vt = V(end-j+1:end,:);
        Vt(abs(Vt) < 1e-13) = nan;
        tmp = V(end-2*j+1:end-j,:)./Vt;
        VV1 = max(abs(N(end-j+1:end,:)),[],2);
        VV2 = max(abs(N(end-2*j+1:end-j,:)),[],2);
        if size(N,2) > 0
            ind = VV1 < 1e-13 & VV2 < 1e-13;
        else
            ind = ones(size(tmp,1));
        end

        tmp = tmp(ind,:);
        Vt = Vt(ind,:);
        [~,ind2] = max(abs(Vt));
        if size(ind2,1) > 0
            ind2 = sub2ind(size(tmp),ind2,1:size(Vt,2));
            %roots(:,1) = median(tmp,1,'omitnan');
            roots(:,1) = tmp(ind2);
        end

    
        j = prod(dims(:,1),'all')*2*(dims(1,3)-1);
        tmp = V(end-2*j+1:end-j,:)./V(end-j+1:end,:);

        VV1 = max(abs(N(end-j+1:end,:)),[],2);
        VV2 = max(abs(N(end-2*j+1:end-j,:)),[],2);
        if size(N,2) > 0
            ind = VV1 < 1e-13 & VV2 < 1e-13;
        else
            ind = ones(size(tmp,1));
        end
        %tmp = tmp(ind,:);

        tmp = tmp(ind,:);
        Vt = Vt(ind,:);
        [~,ind2] = max(abs(Vt));
        if size(ind2,1) > 0
            ind2 = sub2ind(size(tmp),ind2,1:size(Vt,2));
            %roots(:,1) = median(tmp,1,'omitnan');
            roots(:,2) = tmp(ind2);
        end

        %roots(:,2) = median(tmp,1);
else
    %Two linear vars
    [lambda1,V1] = singChebPEP(flip(R1,3));
    [lambda2,V2] = singChebPEP(flip(R2,3));

%%%Extract roots
roots = zeros(size(lambda1,1),d);
roots(:,3) = lambda1;

%%%First coordinate
%%%project eigenvectors
    v = flip(vandercheb(randn,size(R1,3)));
    %%%Evaluate
    a = size(R1,1);
    T = reshape(reshape(R1,[],size(R1,3))*v,a,a);
    [~,S,W] = svd(T);
    idx = find(abs(diag(S)/S(1,1))>1e-13,1,'last');
    N = W(:,idx+1:end);
    for i = 1:size(V1,2)
       for j = 1:size(N,2)
           V1(:,i) = V1(:,i) - (V1(:,i).'*N(:,j))*N(:,j);
       end
    end
        j = prod(dims(:,1),'all');
        Vt = V1(end-j+1:end,:);
        %Vt(abs(Vt) < 1e-6) = nan;
        tmp = V1(end-2*j+1:end-j,:)./Vt;
        VV1 = max(abs(N(end-j+1:end,:)),[],2);
        VV2 = max(abs(N(end-2*j+1:end-j,:)),[],2);
        if size(N,2) > 0
            ind = VV1 < 1e-13 & VV2 < 1e-13;
        else
            ind = ones(size(tmp,1));
        end
        tmp = tmp(ind,:);
        Vt = Vt(ind,:);
        k = ceil(size(Vt,1)/20);
        [~,ind2] = maxk(abs(Vt),k,1);
        if size(ind2,1) > 0
            ind2 = sub2ind(size(tmp),ind2,repmat(1:size(Vt,2),k,1));
            tmp = tmp(ind2);
            roots(:,1) = median(tmp,1,'omitnan');
            %roots(:,1) = tmp(ind2);
        end



%%%Second coordinate matching
v = flip(vandercheb(randn,size(R2,3)));
    %%%Evaluate
    a = size(R2,1);
    T = reshape(reshape(R2,[],size(R2,3))*v,a,a);
    [~,S,W] = svd(T);
    idx = find(abs(diag(S)/S(1,1))>1e-13,1,'last');
    N = W(:,idx+1:end);
    for i = 1:size(V2,2)
       for j = 1:size(N,2)
           V2(:,i) = V2(:,i) - (V2(:,i).'*N(:,j))*N(:,j);
       end
    end
        j = prod(dims(:,1),'all');
        Vt = V2(end-j+1:end,:);
        %Vt(abs(Vt) < 1e-6) = nan;
        tmp = V2(end-2*j+1:end-j,:)./Vt;
        VV1 = max(abs(N(end-j+1:end,:)),[],2);
        VV2 = max(abs(N(end-2*j+1:end-j,:)),[],2);
        if size(N,2) > 0
            ind = VV1 < 1e-13 & VV2 < 1e-13;
        else
            ind = ones(size(tmp,1));
        end
        % tmp = tmp(ind,:);
        % match = median(tmp,1,'omitnan');
        % M = match.' - roots(:,1).';
        % [~,ind] = min(abs(M));
        % roots(:,2) = lambda2(ind);

        tmp = tmp(ind,:);
        Vt = Vt(ind,:);
        k = ceil(size(Vt,1)/20);
        [~,ind2] = maxk(abs(Vt),k,1);
        %rtmp = roots;
        if size(ind2,1) > 0
            ind2 = sub2ind(size(tmp),ind2,repmat(1:size(Vt,2),k,1));
            tmp = tmp(ind2);
            match = median(tmp,1,'omitnan');
            M = match.' - roots(:,1).';
            [~,ind] = min(abs(M));
            roots(:,1) = mean([roots(:,1) match(ind).'],2);
            roots(:,2) = lambda2(ind);
            % rtmp(:,2) = lambda2(ind(2,:));
            % roots = [roots;rtmp];
        end

end

if isequal(rev,[0 1 0])
    roots = roots(:,[1 3 2]);
elseif isequal(rev,[1 0 0])
    roots = roots(:,[3 1 2]);
elseif isequal(rev,[0 1 1])
    roots = roots(:,[1 2 3]);
elseif isequal(rev,[1 1 0])
    roots = roots(:,[3 2 1]);
elseif isequal(rev,[1 0 1])
    roots = roots(:,[2 1 3]);
end

end

end


function R = tensorDixon(F,d,dims)
%%%Constructs the tensor Dixon matrix of F, hiding the variable xd

%%%Find max degrees
m = max(dims(:,3:d+2),[],1);

%%%Block size
b = prod(dims(:,1),'all');

%%%Pretend we have a generic system of d-degree m
k = prod(m(1:d-1)-ones(1,d-1),'all')*b*factorial(d-1);

%%%Form Dixon resultant at Cheb points
dy = sum(dims(:,d+2))-d+1;
R = zeros(k,k,dy);

xd = chebpts(dy);
for i = 1:dy
    R(:,:,i) = scalarTensorDixon(F,dims,m,d,xd(i));
end

%%%Convert back to coefficients
R = mvals2coeffs(R,dy);

%%%Cutoff Negligible R
nrm = norm(R(:,:,1),'fro');
for ii=size(R,3):-1:1
    if norm(R(:,:,ii),'fro')/nrm > 10*eps,    break;    end
end
R = R(:,:,1:ii);

end

function R = scalarTensorDixon(F,dims,m,d,xd)
%%%Constructs the tensor Dixon matrix of F at a specific value of xd (d matrix polynomials in d-1
%%%variables)

%%%Evaluate F at xd
for i = 1:d
    %%%Get Chebyshev vandermonde
    v = flip(vandercheb(xd,dims(i,d+2)));

    %%%Evaluate
    F{i} = reshape(reshape(F{i},[],dims(i,d+2))*v,dims(i,1:d+1));
end

%%%Form matrix at specific value of xd

if d==2
    %R = kron2(F,d,m);
    R = kron2(F,m);
elseif d==3
    R = kron3(F);
end

end

function v = vandercheb(xd,d)
% Chebyshev vandermonde vector
v = zeros(d,1);
for i = 1:d 
    v(end-i+1) = cos((i-1)*acos(xd));
end
end

function A = mvals2coeffs(B,d)
%Matrix vals2coeffs, univariate
A = zeros([size(B,1:2),d]);
    for i = 1:size(B,1)
        for j = 1:size(B,2)
            values = reshape(B(i,j,:),[],1);
            n = size(values, 1);
            tmp = [values(n:-1:2,:) ; values(1:n-1,:)];
            coeffs = ifft(tmp);
            coeffs = coeffs(1:n,:);
            coeffs(2:n-1,:) = 2*coeffs(2:n-1,:);
            A(i,j,:) = coeffs.';
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%Kronecker Methods%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Construct Kronecker products that generate the Dixon function
%%%Given F, returns the coefficient matrix of the Dixon function

function D = kron2(F,m)
    %%%Construct numerator:
    F1 = F{1};
    F2 = F{2};

    [a1, b1, c1] = size(F1);
    [a2, b2, c2] = size(F2);

    F1 = reshape(permute(F1, [3 2 1]), [], a1*b1);
    F2 = reshape(permute(F2, [3 2 1]), [], a2*b2);
    T1 = reshape(permute(kron(F1,F2), [2 1]),a1*b1,a2*b2,c1,c2);
    T1 = flip(flip(permute(reshape(permute(reshape(permute(reshape(T1,a2,a2,[],c1,c2),[1 3 2 4 5]),b1,b2,b1,[],c1,c2),[1 2 4 3 5 6]),a1*a2,[],c1,c2),[2 1 4 3]),3),4);
    R = T1 - permute(T1,[1 2 4 3]);

    %%%"Divide" by s-t

    %%%Size of block matrix
    k = m(1)-1;

    %%%Block size
    b = size(R,1);

    D = zeros(b,b,k,k);
    R = 2*R;

    if k == 1
        D = R(:,:,1,2); 
    else
        D(:,:,1,:) = R(:,:,1,2:end);
        D(:,:,2,:) = R(:,:,2,2:end)+cat(4,D(:,:,1,2:end-1),2*D(:,:,1,end),zeros(b,b,1,1))+cat(4,zeros(b,b,1,1),D(:,:,1,1:end-1));
        for i = 3:k                                    % backwards substitution
            D(:,:,i,:) = R(:,:,i,2:end)-D(:,:,i-2,1:end)+cat(4,D(:,:,i-1,2:end-1),2*D(:,:,i-1,end),zeros(b,b,1,1))+cat(4,zeros(b,b,1,1),D(:,:,i-1,1:end-1));
        end
        D(:,:,k,:) = D(:,:,k,:)/2;
    end


    %%%Reshape into matrix
    D = reshape(permute(D,[1 3 2 4]),size(D,1)*size(D,3),size(D,2)*size(D,4));

end


function D = kron3(F)
    %%%Construct numerator:
    F1 = F{1};
    F2 = F{2};
    F3 = F{3};
    [n1, n2] = size(F1,[3 4]);
    R = zeros([size(F1, [1 2]).*size(F2,[1 2]).*size(F3,[1 2]) n1 2*n2-1 2*n1-1 n2]);


    for i1 = 1:n1
        for i2 = 1:n2
            for j1 = 1:n1
                for j2 = 1:n2
                    for k1 = 1:n1
                        for k2 = 1:n2
                            %%%Take product from given indices
                            A = kron(kron(F1(:,:,i1,i2),F2(:,:,j1,j2)),F3(:,:,k1,k2))/4;

                            %%%F1(x1,x2) F2({x}1,x2) F3({x}1,{x}2)
                            R(:,:,n1-i1+1,2*n2-i2-j2+1,2*n1-j1-k1+1,n2-k2+1) = R(:,:,n1-i1+1,2*n2-i2-j2+1,2*n1-j1-k1+1,n2-k2+1) + A;
                            R(:,:,n1-i1+1,2*n2-abs(i2-j2)-1,2*n1-j1-k1+1,n2-k2+1) = R(:,:,n1-i1+1,2*n2-abs(i2-j2)-1,2*n1-j1-k1+1,n2-k2+1) + A;
                            R(:,:,n1-i1+1,2*n2-i2-j2+1,2*n1-abs(j1-k1)-1,n2-k2+1) = R(:,:,n1-i1+1,2*n2-i2-j2+1,2*n1-abs(j1-k1)-1,n2-k2+1) + A;
                            R(:,:,n1-i1+1,2*n2-abs(i2-j2)-1,2*n1-abs(j1-k1)-1,n2-k2+1) = R(:,:,n1-i1+1,2*n2-abs(i2-j2)-1,2*n1-abs(j1-k1)-1,n2-k2+1) + A;

                            %%%F1({x}1,x2) F2({x}1,{x}2) F3(x1,x2)
                            R(:,:,n1-k1+1,2*n2-i2-k2+1,2*n1-i1-j1+1,n2-j2+1) = R(:,:,n1-k1+1,2*n2-i2-k2+1,2*n1-i1-j1+1,n2-j2+1) + A;
                            R(:,:,n1-k1+1,2*n2-abs(i2-k2)-1,2*n1-i1-j1+1,n2-j2+1) = R(:,:,n1-k1+1,2*n2-abs(i2-k2)-1,2*n1-i1-j1+1,n2-j2+1) + A;
                            R(:,:,n1-k1+1,2*n2-i2-k2+1,2*n1-abs(i1-j1)-1,n2-j2+1) = R(:,:,n1-k1+1,2*n2-i2-k2+1,2*n1-abs(i1-j1)-1,n2-j2+1) + A;
                            R(:,:,n1-k1+1,2*n2-abs(i2-k2)-1,2*n1-abs(i1-j1)-1,n2-j2+1) = R(:,:,n1-k1+1,2*n2-abs(i2-k2)-1,2*n1-abs(i1-j1)-1,n2-j2+1) + A;

                            %%%F1({x}1,{x}2) F2(x1,x2) F3({x}1,x2)
                            R(:,:,n1-j1+1,2*n2-j2-k2+1,2*n1-i1-k1+1,n2-i2+1) = R(:,:,n1-j1+1,2*n2-j2-k2+1,2*n1-i1-k1+1,n2-i2+1) + A;
                            R(:,:,n1-j1+1,2*n2-abs(j2-k2)-1,2*n1-i1-k1+1,n2-i2+1) = R(:,:,n1-j1+1,2*n2-abs(j2-k2)-1,2*n1-i1-k1+1,n2-i2+1) + A;
                            R(:,:,n1-j1+1,2*n2-j2-k2+1,2*n1-abs(i1-k1)-1,n2-i2+1) = R(:,:,n1-j1+1,2*n2-j2-k2+1,2*n1-abs(i1-k1)-1,n2-i2+1) + A;
                            R(:,:,n1-j1+1,2*n2-abs(j2-k2)-1,2*n1-abs(i1-k1)-1,n2-i2+1) = R(:,:,n1-j1+1,2*n2-abs(j2-k2)-1,2*n1-abs(i1-k1)-1,n2-i2+1) + A;

                            %%%F1(x1,x2) F2({x}1,{x}2) F3({x}1,x2)
                            R(:,:,n1-i1+1,2*n2-i2-k2+1,2*n1-j1-k1+1,n2-j2+1) = R(:,:,n1-i1+1,2*n2-i2-k2+1,2*n1-j1-k1+1,n2-j2+1) - A;
                            R(:,:,n1-i1+1,2*n2-abs(i2-k2)-1,2*n1-j1-k1+1,n2-j2+1) = R(:,:,n1-i1+1,2*n2-abs(i2-k2)-1,2*n1-j1-k1+1,n2-j2+1) - A;
                            R(:,:,n1-i1+1,2*n2-i2-k2+1,2*n1-abs(j1-k1)-1,n2-j2+1) = R(:,:,n1-i1+1,2*n2-i2-k2+1,2*n1-abs(j1-k1)-1,n2-j2+1) - A;
                            R(:,:,n1-i1+1,2*n2-abs(i2-k2)-1,2*n1-abs(j1-k1)-1,n2-j2+1) = R(:,:,n1-i1+1,2*n2-abs(i2-k2)-1,2*n1-abs(j1-k1)-1,n2-j2+1) - A;

                            %%%F1({x}1,x2) F2(x1,x2) F3({x}1,{x}2)
                            R(:,:,n1-j1+1,2*n2-i2-j2+1,2*n1-i1-k1+1,n2-k2+1) = R(:,:,n1-j1+1,2*n2-i2-j2+1,2*n1-i1-k1+1,n2-k2+1) - A;
                            R(:,:,n1-j1+1,2*n2-abs(i2-j2)-1,2*n1-i1-k1+1,n2-k2+1) = R(:,:,n1-j1+1,2*n2-abs(i2-j2)-1,2*n1-i1-k1+1,n2-k2+1) - A;
                            R(:,:,n1-j1+1,2*n2-i2-j2+1,2*n1-abs(i1-k1)-1,n2-k2+1) = R(:,:,n1-j1+1,2*n2-i2-j2+1,2*n1-abs(i1-k1)-1,n2-k2+1) - A;
                            R(:,:,n1-j1+1,2*n2-abs(i2-j2)-1,2*n1-abs(i1-k1)-1,n2-k2+1) = R(:,:,n1-j1+1,2*n2-abs(i2-j2)-1,2*n1-abs(i1-k1)-1,n2-k2+1) - A;

                            %%%F1({x}1,{x}2) F2({x}1,x2) F3(x1,x2)
                            R(:,:,n1-k1+1,2*n2-j2-k2+1,2*n1-i1-j1+1,n2-i2+1) = R(:,:,n1-k1+1,2*n2-j2-k2+1,2*n1-i1-j1+1,n2-i2+1)  - A;
                            R(:,:,n1-k1+1,2*n2-abs(j2-k2)-1,2*n1-i1-j1+1,n2-i2+1) = R(:,:,n1-k1+1,2*n2-abs(j2-k2)-1,2*n1-i1-j1+1,n2-i2+1) - A;
                            R(:,:,n1-k1+1,2*n2-j2-k2+1,2*n1-abs(i1-j1)-1,n2-i2+1) = R(:,:,n1-k1+1,2*n2-j2-k2+1,2*n1-abs(i1-j1)-1,n2-i2+1) - A;
                            R(:,:,n1-k1+1,2*n2-abs(j2-k2)-1,2*n1-abs(i1-j1)-1,n2-i2+1) = R(:,:,n1-k1+1,2*n2-abs(j2-k2)-1,2*n1-abs(i1-j1)-1,n2-i2+1) - A;
                        end
                    end
                end
            end
        end
    end

    %%%"Divide" by \prod s_i-t_i

    % %%%Size of block matrix
    % k1 = m(2)-1;
    % k2 = 2*k1;

    %%%Block size
    b = size(R,1);

    D = zeros(b,b,n1-1,2*n1-2,2*n2-1,n2);
    R = 2*permute(R,[1 2 3 5 4 6]);

    %%%First dimension
    if n1 == 2
        D = R(:,:,1,2:end,:,:); 
    else
        D(:,:,1,:,:,:) = R(:,:,1,2:end,:,:);
        D(:,:,2,:,:,:) = R(:,:,2,2:end,:,:)+cat(4,D(:,:,1,2:end-1,:,:),2*D(:,:,1,end,:,:),zeros(b,b,1,1,2*n2-1,n2))+cat(4,zeros(b,b,1,1,2*n2-1,n2),D(:,:,1,1:end-1,:,:));
        for i = 3:n1-1                                    % backwards substitution
            D(:,:,i,:,:,:) = R(:,:,i,2:end,:,:)-D(:,:,i-2,1:end,:,:)+cat(4,D(:,:,i-1,2:end-1,:,:),2*D(:,:,i-1,end,:,:),zeros(b,b,1,1,2*n2-1,n2))+cat(4,zeros(b,b,1,1,2*n2-1,n2),D(:,:,i-1,1:end-1,:,:));
        end
        D(:,:,n1-1,:,:,:) = D(:,:,n1-1,:,:,:)/2;
    end

    %%%Second Dimension
    R = 2*permute(D,[1 2 3 4 6 5]);
    D = zeros(b,b,n1-1,2*n1-2,n2-1,2*n2-2);
    if n2 == 2
        D = R(:,:,:,:,1,2:end); 
    else
        D(:,:,:,:,1,:) = R(:,:,:,:,1,2:end);
        D(:,:,:,:,2,:) = R(:,:,:,:,2,2:end)+cat(6,D(:,:,:,:,1,2:end-1),2*D(:,:,:,:,1,end),zeros(b,b,n1-1,2*n1-2,1,1))+cat(6,zeros(b,b,n1-1,2*n1-2,1,1),D(:,:,:,:,1,1:end-1));
        for i = 3:n2-1                                   % backwards substitution
            D(:,:,:,:,i,:) = R(:,:,:,:,i,2:end)-D(:,:,:,:,i-2,1:end)+cat(6,D(:,:,:,:,i-1,2:end-1),2*D(:,:,:,:,i-1,end),zeros(b,b,n1-1,2*n1-2,1,1))+cat(6,zeros(b,b,n1-1,2*n1-2,1,1),D(:,:,:,:,i-1,1:end-1));
        end
        D(:,:,:,:,n2-1,:) = D(:,:,:,:,n2-1,:)/2;
    end
    D = permute(D,[1 2 3 4 6 5]);


    %%%Reshape into matrix
    %D = permute(D,[1 2 3 5 4 6]);
    D = reshape(permute(D,[1 3 5 2 4 6]),size(D,1)*size(D,3)*size(D,5),size(D,2)*size(D,4)*size(D,6));
end


%%%%%%%%%%%%%%%%%%%%%%%%%Eigenvalue Method%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lambda,X] = singChebPEP(c)
% Finds roots of a singular Cheb eigenvalue problem via projection. Coefficients are ordered  highest degree down.
%
% c is a nxnxk array.  k-1 is the degree, n is the matrix size.

if size(c,3) ==2  % linear case
    opts.method = 'project';
    [~,~,Z,lambda,X,~,~,V,~,~,~,~] = singgep(c(:,:,2),-c(:,:,1),opts);
    ind = Z(:,8)==1;
    lambda = lambda(ind);
    X = X(:,ind);
    X = [zeros(size(V,1)-size(X,1),size(X,2)); X];
    X = V*X; 
else
    %%%Project pencil
    k=length(c(1,1,:)); n=length(c(:,:,1));

    %%%Normal rank
    nrank = 0;
    for ii = 1:6
        v = vandercheb(randn,k);
        nrank = max(nrank,rank(reshape(reshape(c,[],k)*v,n,n)));
    end
    
    [U,~] = qr(randn(n));
    [V,~] = qr(randn(n));
    
    c = permute(reshape(permute(reshape(permute(reshape(U'*reshape(c,n,[]),n,n,k),[1 3 2]),[],n)*V,[2 1]),n,n,k),[2 1 3]);
    
    ca = c(1:end-nrank,end-nrank+1:end,:);
    cb = c(end-nrank+1:end,1:end-nrank,:);
    c = c(end-nrank+1:end,end-nrank+1:end,:);
    
    ns = n;
    n = nrank;
    
    for ii=2:k
        c(:,:,ii)=c(:,:,ii)*(-.5); % coefficients
    end
    c(:,:,3) = c(:,:,3)+.5*c(:,:,1);
    
    oh = .5*ones(n*(k-2),1);
    % form colleague matrix A,B:
    A = diag(oh,n)+diag(oh,-n);
    A(end-n+1:end,end-2*n+1:end-n) = eye(n);
    
    for ii=1:k-1
        A(1:n,(ii-1)*n+1:ii*n) = c(:,:,ii+1);
    end
    B=eye(size(A)); B(1:n,1:n)=c(:,:,1);
    
    [X,lambda,Y] = eig(A,B);
    lambda = diag(lambda);
    
    %%%Check for finite regular eigenvalues
    X = X(end-n+1:end,:);
    Y = Y(1:n,:);
    ind = ~isinf(lambda);
    lambda = lambda(ind);
    X = X(:,ind);
    Y = Y(:,ind);
    ind = zeros(size(lambda));
    for j = 1:size(lambda)
        v = vandercheb(lambda(j),k);
        c1 = reshape(reshape(ca,[],k)*v,ns-n,n)*X(:,j);
        c2 = Y(:,j)'*reshape(reshape(cb,[],k)*v,n,ns-n);
        ind(j) = (norm(c1) < 1e-8) && (norm(c2) < 1e-8);
    end
    ind = (ind == 1);
    
    %lambda = lambda(ind);
    %X = X(:,ind);
    
    %%%augment eigenvectors
    X = [zeros(size(V,1)-size(X,1),size(X,2)); X];
    X = V*X;
end


end


