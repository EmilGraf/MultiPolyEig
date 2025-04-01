% Plots Figure 1 in the paper E Graf and A Townsend:
% A Hidden Variable Resultant Method for the Polynomial Multiparameter
% Eigenvalue Problem, 2025.
%
% Based on ExampleBrassTeflon from RandomJointEig, https://github.com/borplestenjak/RandomJointEig,
% Bor Plestenjak, Hauke Gravenkamp, Daniel Kiefer 2024.
%
% We use data from H. Gravenkamp, B. Plestenjak, D. A. Kiefer, and J. Elias: Computation 
% of leaky waves in layered structures coupled to unbounded media by 
% exploiting multiparameter eigenvalue problems. arXiv:2404.15277, 2024.
% Requires RandomJointEig, including data from example_brassTeflon_Data.
%
% E Graf and A Townsend, 2025.

clearvars;

% load data
load example_brassTeflon_Data.mat
R1 = R{1};
R2 = R{2};

w = 2*pi*linspace(.01, 3, 150).'; % frequency

rng(2)

n = size(E0,1);
NN = nan(length(w), 2^3*n);

% multipolyeig
k1t = cell(length(w));
kyL1t = k1t;
kyS1t = k1t;
% this loop is set up to be replaced with parfor if desired
for i = 1:length(w)
    disp(strcat("run: ", string(i)));
    kappa = w(i)./c;
    kappal = kappa(1);
    kappat = kappa(2);
    sols = auxLeakyWaves(E0,E1,-E2,M,R1,R2,kappal,kappat,w(i));
    k1t{i}  = sols(:,1);
    kyL1t{i}  = sols(:,2);
    kyS1t{i}  = sols(:,3);
end

k1 = NN + 1i*NN;
kyL1 = k1;
kyS1 = k1;

for i = 1:length(w)
    k1(i,1:numel(k1t{i}))  = k1t{i};
    kyL1(i,1:numel(kyL1t{i}))  = kyL1t{i};
    kyS1(i,1:numel(kyS1t{i}))  = kyS1t{i};
end

% MultiParEig
if contains(path, ['OldMultiParEig', pathsep])
    rmpath OldMultiParEig
end
addpath NewMultiParEig

% options for new multipareig as in ExampleBrassTeflon
% this section exactly as in ExampleBrassTeflon
opts = [];
opts.refine = 0;
opts.rand_orth = 1;
opts.twosidedRQ = 1;
opts.solver = 'eig';

k2 = NN + 1i*NN;
kyL2 = k2;
kyS2 = k2;
for i = 1:length(w)
    kappa = w(i)./c;
    kappal = kappa(1);
    kappat = kappa(2);
    [lambda2, tmp_lambda2] = eig_LeakySolid(E0,E1,-E2,M,R1,R2,kappal,kappat,w(i),opts);
    k2(i,1:numel(lambda2))  = lambda2;
    kyL2(i,1:numel(lambda2))  = tmp_lambda2(:,2);
    kyS2(i,1:numel(lambda2))  = tmp_lambda2(:,3);
end

fh1 = w/2/pi.*ones(size(k1));
fh2 = w/2/pi.*ones(size(k2));

% filter modes for positive real parts
% keep original data k1,k2
kt1 = k1;
kt2 = k2;
rm1 = (real(kyL1)>-1e-2)|(real(kyS1)>-1e-2);
kt1(rm1) = nan + 1i*nan;
rm2 = (real(kyL2)>-1e-2)|(real(kyS2)>-1e-2);
kt2(rm2) = nan + 1i*nan;

% plot with settings from ExampleBrassTeflon
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultLineLineWidth', 1.75);
set(groot, 'defaultFigurePosition', [100, 50, 1000, 800]);
set(groot, 'defaultAxesLabelFontSizeMultiplier', 2);  % Sets the multiplier for axis label font size
set(groot, 'defaultAxesTitleFontSizeMultiplier', 1.25);
figure 
hold on
plot(fh1(:), real(kt1(:)), 'k.', 'MarkerSize',18,'DisplayName','multipolyeig')
plot(fh2(:), real(kt2(:)), 'ro', 'MarkerSize',6,'DisplayName','multipareig');
ylabel('wavenumber in rad/mm')
xlabel('frequency in Mhz')
legend('MultiPolyEig','multipareig','Location','southwest')
ylim([0.8, 2.9]) 
xlim([0.5, 3])