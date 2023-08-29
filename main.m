%% Preliminaries
clc
close all
clear

addpath('utils')
forig = double(rgb2gray(imread('ncat.png')));

% padding
forig = [forig; zeros(1, size(forig, 1))];
forig = [forig, zeros(size(forig, 1) ,1)];
forig = forig/max(forig(:));

% create grid
N = size(forig, 1)-1;
x = linspace(-N/2, N/2-1, N);
[X, Y] = meshgrid(x, x);

% limited angular range
Phi = deg2rad(65);

% create Curvelet tiling
nscales = 4;
[W, COORD] = create_tiling(N, nscales, Phi, 'meyer');

% define analysis and synthesis operators
Psi  = @(f) transform(f, W, COORD, N, false, 'analysis_loc');
PsiT = @(c) transform(c, W, COORD, N, false, 'synthesis_loc');

% define Radon transform and adjoint
thetas_rad = linspace(-pi/2, pi/2*(1-1/180), 180);
thetas_rad(abs(thetas_rad)>Phi) = 0;
thetas_rad = thetas_rad(find(thetas_rad));
thetas = rad2deg(thetas_rad);
A  = @(f) radon(f, thetas);
AT = @(g) iradon(g, thetas, N+1, 'None');
L = 25; % constant for Chambolle-Pock algorithm, should be <= ||A||_2

% create data
rng(123)
p = 1e3; % number of incident photons
mu = 0.02; % mass attenuation coefficient for water
data = A(forig);
data = exp(-data .* mu);
data = -log(poissrnd(data * p) /p)*1/mu;

%% sparse ell1 regularization

% initialization
Nsyn = 200; % number of iterations
c0 = Psi(zeros(N+1));
h  = zeros(N+1);

% regularization parameter
alpha = 0.01;

[crec_sparse_ell1, error] = sparse_ell1(c0, data, h, A, AT, Psi, PsiT, N+1, L, alpha, Nsyn, forig, 'None');
frec_sparse_ell1 = PsiT(crec_sparse_ell1);

figure, imagesc(frec_sparse_ell1), colormap bone, axis off, axis square, title('sparse ell1 recon')

%% total variation regularization

% initialization
beta = 0.01; % regularization parameter
Ntv = 500; % number of iterations
x0 = zeros(N+1);
frec_tv = tv(x0, data, A, AT, N, L, beta, Ntv, forig, 1);

figure, imagesc(frec_tv), colormap bone, axis off, axis square, title('TV recon')

%% hybrid sparse ell1-TV regularization

alpha = 0.1; % regularization parameter for ell1
beta = 0.1; % regularization parameter for TV
Ncomb = 500; % number of iterations
crec_hybrid = hybrid(c0, data, A, AT, Psi, PsiT, N+1, alpha, beta, L, Ncomb, forig);
frec_hybrid = PsiT(crec_hybrid);

figure, imagesc(frec_hybrid), colormap bone, axis off, axis square, title('hybrid recon')

%% complementary sparse ell1-TV regularization
Nsyn = 200; % number of iterations for ell1 part
Ntv  = 200; % number of iterations for TV part
Nouter = 2; % number of outer iterations
c0 = Psi(zeros(N+1));
f0 = zeros(N+1);
alpha = 0.009;
beta = 0.008;
flag = 0;
for k = 1:Nouter
    [crec, err_syn] = sparse_ell1(c0, data, h, A, AT, Psi, PsiT, N+1, L, alpha, Nsyn, forig, flag);
    [frec, err_tv] = tv(f0, A(PsiT(crec)), A, AT, N, L, beta, Ntv, forig, flag);
    h = A(frec);
    c0 = crec;
    f0 = frec;
    flag = 1;
    beta = beta/2;
end

figure, imagesc(frec), colormap bone, axis off, axis square, title('comp. ell1-TV recon')
