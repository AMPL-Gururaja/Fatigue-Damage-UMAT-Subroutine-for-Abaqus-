% Full Eshelby-Mori-Tanaka model for a short-fibre thermoplastic
% (CF-ABS) with scalar isotropic damage on the matrix. Given a measured
% composite damage curve D_c(N), recovers the RVE-averaged matrix damage
% D_m(N) under the assumption that only the matrix damages.
%
% Matrix (ABS)  : isotropic
% Fiber (carbon): transversely isotropic (axis of symmetry = fiber axis),
%                 assumed elastic and undamaged throughout (D_f = 0)
%
% Uses:
%   - full fibre length spectrum (32 bins, volume-weighted)
%   - full orientation distribution from phi and theta marginals
%
% Requires:  fiber_distribution.xlsx, caf_damage_data_65uts.xlsx
%            in the current MATLAB working directory.

clear; clc; close all;

%% =====================================================================
%% 1. Material properties 
%% =====================================================================
E_m0 = 2.30e9;        % ABS matrix Young's modulus, Pa
nu_m = 0.35;          % ABS Poisson's ratio (preserved under scalar damage)

% Carbon fibre: transversely isotropic, symmetry axis = fibre axis.
EL_f   = 230e9;       % axial (fibre-direction) modulus, Pa
ET_f   = 15e9;        % transverse modulus, Pa
GLT_f  = 15e9;        % axial shear modulus, Pa
nuLT_f = 0.20;        % major Poisson ratio (axial load -> transverse strain)
nuT_f  = 0.25;        % transverse Poisson ratio (in the 1-2 plane)
% inside transversely_isotropic_stiffness as G12 = ET/(2*(1+nuT)).

E_c0 = 8.3e9;           % Measured pristine composite modulus, Pa
d_f  = 8.0;           % Fibre diameter, um

% 2a. Eshelby tensor large-alpha (continuous fibre) limit -- depends only
%     on matrix properties and inclusion shape, unaffected by the fibre
%     stiffness model.
S_inf = eshelby_prolate(1e5, nu_m);

% 2b. Transversely isotropic fibre stiffness
%     Symmetry axis is local x3, so expect:
%       C1111 == C2222                 (isotropy in the 1-2 plane)
%       C1212 == (C1111 - C1122)/2     (in-plane shear tied to in-plane E,nu)
Cf_check = transversely_isotropic_stiffness(EL_f, ET_f, GLT_f, nuLT_f, nuT_f);

%% =====================================================================
%% 3. Load fibre length & orientation distributions
%% =====================================================================
raw = readmatrix('fiber_distribution.xlsx', 'Range', 'A3:K34');

L_lo = raw(:,1);  L_hi = raw(:,2);  L_cnt   = raw(:,3);
phi_lo = raw(:,5); phi_hi = raw(:,6); phi_cnt = raw(:,7);
th_lo = raw(:,9); th_hi  = raw(:,10); th_cnt = raw(:,11);

L_mid   = 0.5*(L_lo + L_hi);           % um
phi_mid = deg2rad(0.5*(phi_lo + phi_hi));
th_mid  = deg2rad(0.5*(th_lo  + th_hi));

% Volume-weighted aspect ratio distribution (fibre volume ~ count*length)
ar   = L_mid / d_f;
ar_w = L_cnt .* L_mid;   ar_w = ar_w / sum(ar_w);

% Normalised marginal ODFs; joint = product 
p_phi = phi_cnt / sum(phi_cnt);
p_th  = th_cnt  / sum(th_cnt);

fprintf('\nFibre data: %d fibres, mean L = %.1f um, vol-weighted mean AR = %.1f\n', ...
    sum(L_cnt), sum(L_mid.*L_cnt)/sum(L_cnt), sum(ar.*ar_w));

%% =====================================================================
%% 4. Calibrate V_f against the measured pristine composite modulus
%% =====================================================================
Vf_cal = bisect_Vf(E_c0, 0.03, 0.35, E_m0, nu_m, ...
                   EL_f, ET_f, GLT_f, nuLT_f, nuT_f, ...
                   ar, ar_w, phi_mid, p_phi, th_mid, p_th);
Cbar0  = MT_effective_stiffness(E_m0, nu_m, EL_f, ET_f, GLT_f, nuLT_f, nuT_f, Vf_cal, ...
                                ar, ar_w, phi_mid, p_phi, th_mid, p_th);
Ec_check = young_from_C(Cbar0, 1);
fprintf('\nCalibrated Vf = %.4f  =>  E_c0 (model) = %.3f GPa  (target %.3f)\n', ...
    Vf_cal, Ec_check/1e9, E_c0/1e9);

%% =====================================================================
%% 5. Sweep D_m and build the D_m -> D_c mapping
%% =====================================================================
Dm_grid = linspace(0, 0.6, 61);
Ec_grid = zeros(size(Dm_grid));
fprintf('\nSweeping D_m ...');
for i = 1:length(Dm_grid)
    Cbar = MT_effective_stiffness(E_m0*(1-Dm_grid(i)), nu_m, ...
                                  EL_f, ET_f, GLT_f, nuLT_f, nuT_f, Vf_cal, ...
                                  ar, ar_w, phi_mid, p_phi, th_mid, p_th);
    Ec_grid(i) = young_from_C(Cbar, 1);
end
Dc_grid = 1 - Ec_grid / Ec_grid(1);
fprintf(' done.\n');

%% =====================================================================
%% 6. Load composite damage & invert pointwise
%% =====================================================================
T = readmatrix('caf_damage_data_65uts.xlsx', 'Range', 'A3:B1108');
T = T(all(~isnan(T),2), :);
N_meas  = T(:,1);
Dc_meas = T(:,2);

Dm_meas = interp1(Dc_grid, Dm_grid, Dc_meas, 'pchip');
ratio   = Dm_meas ./ max(Dc_meas, 1e-9);

writematrix([N_meas, Dc_meas, Dm_meas, ratio], ...
    'matrix_damage_recovered_matlab.csv');

%% =====================================================================
%% 7. Plots
%% =====================================================================
figure('Position', [80 80 1500 430], 'Color', 'w');

subplot(1,3,1);
plot(Dm_grid, Dc_grid, 'k-', 'LineWidth', 2); hold on;
plot([0 0.6], [0 0.6], 'r--');
xlabel('Matrix damage D_m'); ylabel('Composite damage D_c');
title(sprintf('MT-predicted D_m \\rightarrow D_c\n(CF-ABS, V_f = %.2f, mean AR = %.1f)', ...
    Vf_cal, sum(ar.*ar_w)));
legend({'MT model', 'D_c = D_m'}, 'Location', 'best'); grid on;

subplot(1,3,2);
plot(N_meas, Dc_meas, 'b-', 'LineWidth', 1); hold on;
plot(N_meas, Dm_meas, 'r-', 'LineWidth', 1);
xlabel('Cycles N'); ylabel('Damage');
title('Composite vs matrix damage over life');
legend({'D_c (measured)', 'D_m (MT-recovered)'}, 'Location', 'best'); grid on;

subplot(1,3,3);
mask = Dc_meas > 0.005;
plot(N_meas(mask), ratio(mask), 'g-', 'LineWidth', 1); hold on;
yline(1, ':k');
xlabel('Cycles N'); ylabel('D_m / D_c');
title('Damage amplification factor'); grid on;
ylim([0, max(3, max(ratio(mask))*1.05)]);

saveas(gcf, 'matrix_damage_recovery_matlab.png');

fprintf('\nAmplification D_m/D_c (mask Dc > 0.005):\n');
fprintf('  mean = %.3f, min = %.3f, max = %.3f\n', ...
    mean(ratio(mask)), min(ratio(mask)), max(ratio(mask)));
fprintf('  final N = %d: D_c = %.4f,  D_m = %.4f,  ratio = %.2f\n', ...
    round(N_meas(end)), Dc_meas(end), Dm_meas(end), ratio(end));

%% =====================================================================
%% ==================  LOCAL FUNCTIONS  ================================
%% =====================================================================

function S = eshelby_prolate(alpha, nu)
% Prolate spheroid (alpha > 1) Eshelby tensor, aligned along x3,
% embedded in an isotropic matrix of Poisson ratio nu.
% Depends only on matrix properties and inclusion aspect ratio -- NOT
% on the inhomogeneity (fibre) stiffness
% Reference: Mura (1987). Verified against alpha->inf fibre limit.
    if alpha < 1.0001, alpha = 1.0001; end
    a2 = alpha^2;
    g  = alpha/(a2-1)^1.5 * (alpha*sqrt(a2-1) - acosh(alpha));

    S1111 = 3/(8*(1-nu))*a2/(a2-1) + 1/(4*(1-nu))*(1-2*nu - 9/(4*(a2-1)))*g;
    S3333 = 1/(2*(1-nu))*(1-2*nu + (3*a2-1)/(a2-1) - (1-2*nu + 3*a2/(a2-1))*g);
    S1122 = 1/(4*(1-nu))*(a2/(2*(a2-1)) - (1-2*nu + 3/(4*(a2-1)))*g);
    S1133 = 1/(2*(1-nu))*(-a2/(a2-1) + 0.5*(3*a2/(a2-1) - (1-2*nu))*g);
    S3311 = 1/(2*(1-nu))*((2*nu-1) - 1/(a2-1) + (1-2*nu + 3/(2*(a2-1)))*g);
    S1212 = 1/(4*(1-nu))*(a2/(2*(a2-1)) + (1-2*nu - 3/(4*(a2-1)))*g);
    S1313 = 1/(4*(1-nu))*(1-2*nu - (a2+1)/(a2-1) - 0.5*(1-2*nu - 3*(a2+1)/(a2-1))*g);

    S = zeros(3,3,3,3);
    S = fill_sym(S,1,1,1,1,S1111);
    S = fill_sym(S,2,2,2,2,S1111);
    S = fill_sym(S,3,3,3,3,S3333);
    S = fill_sym(S,1,1,2,2,S1122);   S = fill_sym(S,2,2,1,1,S1122);
    S = fill_sym(S,1,1,3,3,S1133);   S = fill_sym(S,2,2,3,3,S1133);
    S = fill_sym(S,3,3,1,1,S3311);   S = fill_sym(S,3,3,2,2,S3311);
    S = fill_sym(S,1,2,1,2,S1212);
    S = fill_sym(S,1,3,1,3,S1313);   S = fill_sym(S,2,3,2,3,S1313);
end

function T = fill_sym(T,i,j,k,l,v)
% Fill minor-symmetric components T_ijkl = T_jikl = T_ijlk
    T(i,j,k,l) = v;  T(j,i,k,l) = v;
    T(i,j,l,k) = v;  T(j,i,l,k) = v;
end

function C = isotropic_stiffness(E, nu)
% Isotropic 4th-order stiffness C_ijkl. Retained for the matrix, and
% for anyone who wants to compare fibre models by re-enabling it.
    lam = E*nu / ((1+nu)*(1-2*nu));
    mu  = E / (2*(1+nu));
    C   = zeros(3,3,3,3);
    d   = eye(3);
    for i = 1:3
      for j = 1:3
        for k = 1:3
          for l = 1:3
            C(i,j,k,l) = lam*d(i,j)*d(k,l) + mu*(d(i,k)*d(j,l) + d(i,l)*d(j,k));
          end
        end
      end
    end
end

function C = transversely_isotropic_stiffness(EL, ET, GLT, nuLT, nuT)
% Transversely isotropic 4th-order stiffness, symmetry axis = local x3
% (fibre axis). This is evaluated in the LOCAL fibre frame, i.e. before
% rotate4() maps it into the global (measured phi/theta) frame -- exactly
% the frame in which a transversely isotropic tensor is naturally defined.
%
% Inputs (5 independent constants):
%   EL   - axial (fibre-direction) Young's modulus
%   ET   - transverse Young's modulus
%   GLT  - axial (longitudinal-transverse) shear modulus
%   nuLT - major Poisson ratio (axial load -> transverse strain)
%   nuT  - transverse Poisson ratio (in the 1-2 plane)
%
% Derived (not independent):
%   G12  = ET / (2*(1+nuT))            transverse-plane shear modulus
%   nu31 = nuLT * ET/EL                 minor Poisson ratio, by reciprocity

    G12 = ET / (2*(1+nuT));

    % Compliance in Voigt order [11 22 33 23 13 12], axis 3 = fibre axis
    S = zeros(6,6);
    S(1,1) = 1/ET;      S(2,2) = 1/ET;      S(3,3) = 1/EL;
    S(1,2) = -nuT/ET;   S(2,1) = S(1,2);
    S(1,3) = -nuLT/EL;  S(3,1) = S(1,3);    % reciprocity: nu13/E1 = nu31/E3
    S(2,3) = -nuLT/EL;  S(3,2) = S(2,3);
    S(4,4) = 1/GLT;     % 23 shear
    S(5,5) = 1/GLT;     % 13 shear
    S(6,6) = 1/G12;     % 12 shear (transverse plane)

    Cv = inv(S);

    % Voigt -> 4th-order tensor
    idx = [1 1; 2 2; 3 3; 2 3; 1 3; 1 2];
    C = zeros(3,3,3,3);
    for I = 1:6
        for J = 1:6
            a=idx(I,1); b=idx(I,2); c=idx(J,1); d=idx(J,2);
            v = Cv(I,J);
            C(a,b,c,d)=v; C(b,a,c,d)=v; C(a,b,d,c)=v; C(b,a,d,c)=v;
        end
    end
end

function I = identity4()
% Symmetric 4th-order identity I_ijkl = 0.5*(d_ik d_jl + d_il d_jk).
    I = zeros(3,3,3,3);
    d = eye(3);
    for i = 1:3
      for j = 1:3
        for k = 1:3
          for l = 1:3
            I(i,j,k,l) = 0.5*(d(i,k)*d(j,l) + d(i,l)*d(j,k));
          end
        end
      end
    end
end

function R = dcon4(A, B)
% Double contraction (A:B)_ijkl = A_ijmn B_mnkl (4th-order tensors).
% MATLAB column-major reshape puts (i,j) on rows with i fast, (k,l) on
% cols with k fast, matching the tensor-contraction index convention.
    R = reshape(reshape(A,9,9) * reshape(B,9,9), 3,3,3,3);
end

function Tinv = inv4_mandel(T)
% Invert a 4th-order tensor that maps symmetric->symmetric, via the
% Mandel 6x6 representation (isotropic C is singular in the naive
% 9x9 view because minor-symmetric rows are redundant; Mandel cures it).
    idx = [1 1; 2 2; 3 3; 2 3; 1 3; 1 2];
    w   = [1 1 1 sqrt(2) sqrt(2) sqrt(2)];
    M = zeros(6,6);
    for I = 1:6
        for J = 1:6
            a=idx(I,1); b=idx(I,2); c=idx(J,1); d=idx(J,2);
            M(I,J) = w(I)*w(J) * T(a,b,c,d);
        end
    end
    Minv = inv(M);
    Tinv = zeros(3,3,3,3);
    for I = 1:6
        for J = 1:6
            a=idx(I,1); b=idx(I,2); c=idx(J,1); d=idx(J,2);
            v = Minv(I,J) / (w(I)*w(J));
            Tinv(a,b,c,d) = v;  Tinv(b,a,c,d) = v;
            Tinv(a,b,d,c) = v;  Tinv(b,a,d,c) = v;
        end
    end
end

function Trot = rotate4(T, R)
% Rotate 4th-order tensor: T'_ijkl = R_ip R_jq R_kr R_ls T_pqrs
% via T_rot = A * T_mat * A' with A = kron(R,R) (9x9). Verified in
% Python against einsum.
    A = kron(R, R);
    Trot = reshape(A * reshape(T,9,9) * A.', 3,3,3,3);
end

function R = fiber_rotation(phi, theta)
% Rotation matrix with column 3 = fibre-axis unit vector n.
% Convention: theta = azimuth in xy-plane, phi = elevation from xy-plane.
% Transverse orientation of columns 1,2 is arbitrary (fibre has trans. iso.,
% so the choice of e1,e2 within the transverse plane does not matter).
    n = [cos(phi)*cos(theta); cos(phi)*sin(theta); sin(phi)];
    if abs(n(3)) < 0.9, v = [0;0;1]; else, v = [1;0;0]; end
    e1 = cross(v, n);  e1 = e1 / norm(e1);
    e2 = cross(n, e1);
    R  = [e1, e2, n];
end

function Cbar = MT_effective_stiffness(Em, num, EL, ET, GLT, nuLT, nuT, Vf, ...
                                       ar, ar_w, phi, p_phi, th, p_th)
% Full 4th-order Mori-Tanaka with orientation & aspect-ratio averaging.
% Matrix: isotropic. Fibre: transversely isotropic (built in local
% fibre frame, then rotated per (phi,theta) bin).
% Outer loop: aspect ratio (Eshelby tensor changes with alpha).
% Inner loops: phi, theta bins (each gets a rotation of A_dil_local).
% Returns effective stiffness Cbar (3x3x3x3), symmetrised (major sym).
    Cm  = isotropic_stiffness(Em, num);
    Cf  = transversely_isotropic_stiffness(EL, ET, GLT, nuLT, nuT);
    dC  = Cf - Cm;
    Cmi = inv4_mandel(Cm);
    I4  = identity4();

    Adil_avg = zeros(3,3,3,3);
    idx_p = find(p_phi > 0);
    idx_t = find(p_th  > 0);

    for ia = 1:length(ar)
        if ar_w(ia) <= 0, continue; end
        S     = eshelby_prolate(ar(ia), num);
        term  = dcon4(dcon4(S, Cmi), dC);
        Adil  = inv4_mandel(I4 + term);     % dilute in fibre frame
        for ip = idx_p'
            for it = idx_t'
                Rrot = fiber_rotation(phi(ip), th(it));
                Ag   = rotate4(Adil, Rrot);
                Adil_avg = Adil_avg + ar_w(ia) * p_phi(ip) * p_th(it) * Ag;
            end
        end
    end

    denom = (1-Vf) * I4 + Vf * Adil_avg;
    AMT   = dcon4(Adil_avg, inv4_mandel(denom));
    Cbar  = Cm + Vf * dcon4(dC, AMT);

    % Symmetrise (major symmetry)
    Cbar_sym = zeros(3,3,3,3);
    for i = 1:3
      for j = 1:3
        for k = 1:3
          for l = 1:3
            Cbar_sym(i,j,k,l) = 0.5 * (Cbar(i,j,k,l) + Cbar(k,l,i,j));
          end
        end
      end
    end
    Cbar = Cbar_sym;
end

function E = young_from_C(Cbar, dir)
% Young's modulus in direction 'dir' (=1,2,3) from Cbar.
    Sc = inv4_mandel(Cbar);
    E  = 1 / Sc(dir, dir, dir, dir);
end

function Vf = bisect_Vf(Etarget, lo, hi, Em, num, EL, ET, GLT, nuLT, nuT, ...
                        ar, arw, phi, pp, th, pt)
% Bisection for Vf such that the pristine model E_x equals Etarget.
    for it = 1:30
        mid  = 0.5*(lo + hi);
        Cbar = MT_effective_stiffness(Em, num, EL, ET, GLT, nuLT, nuT, mid, ...
                                      ar, arw, phi, pp, th, pt);
        Emid = young_from_C(Cbar, 1);
        if Emid < Etarget, lo = mid; else, hi = mid; end
        if abs(Emid - Etarget)/Etarget < 1e-3, Vf = mid; return; end
    end
    Vf = mid;
end