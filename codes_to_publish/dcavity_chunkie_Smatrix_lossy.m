%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  dcavity_chunkie_Smatrix_lossy.m
%
%  Scattering matrix of a LOSSY (impedance-wall) D-shaped reverberation cavity,
%  computed with a dense boundary-integral-equation (BIE) solve using chunkie.
%  Identical to dcavity_chunkie_Smatrix_lossless.m EXCEPT the outer enclosure
%  (arc + chord segments) carries a finite-conductivity Leontovich IMPEDANCE
%  boundary condition instead of PEC -- so the S-matrix is SUB-unitary (the walls
%  absorb).
%
%  PHYSICS
%    2D scalar Helmholtz, TM (scalar H_z), time convention exp(-i*omega*t).
%    PEC wall               -> Neumann:  dH_z/dn = 0.
%    Leontovich impedance   -> Robin:    dH_z/dn + gamma*H_z = 0 ,
%       gamma = i*omega_phys*eps0*Zs ,
%       Zs    = sqrt(i*omega_phys*mu0/(sigma + i*omega_phys*eps0))   (surface impedance).
%    Impedance on the arc + the chord segments; PEC on the lead walls + ellipses.
%
%  METHOD  (same as the lossless script)
%    Single-layer BIE, mixed Robin/Neumann operator scaled to (I + K) form and
%    assembled as a per-edge chunkie kernel matrix (RCIP-correct at corners).
%    Drive each lead; project onto the TEM mode on a cut-line; S = a_out / a_in.
%    Validation: manufactured-solution gate; absorption per channel > 0.
%
%  OUTPUT (into dcavity_chunkie_Smatrix_lossy_outputs/):
%    S_lossy_chunkie.mat / .csv   (the 4x4 S-matrix + metadata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

%% ---------------------------------------------------------------- chunkie path
% Locate chunkie: prefer a local copy/submodule next to these scripts (the
% published repo ships chunkie as a git submodule -- run `git submodule update
% --init`), else fall back to a development copy under $HOME.
here_ = fileparts(mfilename('fullpath')); if isempty(here_), here_ = pwd; end
cands_ = {fullfile(here_,'chunkie'), fullfile(here_,'..','chunkie'), ...
          fullfile(here_,'..','..','chunkie'), fullfile(getenv('HOME'),'main_projects','Numerics','chunkie')};
chunkie_root = '';
for ic_ = 1:numel(cands_), if exist(fullfile(cands_{ic_},'startup.m'),'file'), chunkie_root = cands_{ic_}; break; end, end
assert(~isempty(chunkie_root), 'chunkie not found -- run: git submodule update --init --recursive  (see README).');
olddir = pwd; cd(chunkie_root); startup; cd(olddir);

thisdir = fileparts(mfilename('fullpath'));  if isempty(thisdir), thisdir = pwd; end
outdir  = fullfile(thisdir,'dcavity_chunkie_Smatrix_lossy_outputs');
if ~exist(outdir,'dir'), mkdir(outdir); end

%% ------------------------------------------------------------------ parameters
c_light = 299792458;
eps0 = 8.8541878128e-12;  mu0 = 1.25663706212e-6;

omega   = 0.97*(2*pi);
lambda  = 2*pi/omega;
f_hz    = c_light/lambda;
omega_phys = c_light*omega;

sigma_cond = 50;                    % LOSSY: wall conductivity [S/m]  (50% absorption regime)
Zs    = sqrt( 1i*omega_phys*mu0 / (sigma_cond + 1i*omega_phys*eps0) );
gamma = 1i*omega_phys*eps0*Zs;
fprintf('LOSSY D-cavity: sigma=%.3g  Zs=%.4g%+.4gi Ohm  |Zs|/Z0=%.3e  gamma=%.4g%+.4gi\n',...
        sigma_cond, real(Zs),imag(Zs), abs(Zs)/sqrt(mu0/eps0), real(gamma),imag(gamma));

% --- D-cavity + leads + ellipse geometry (same as the lossless script) ---
g_.R = 5.75; g_.d = g_.R/2; g_.L = 4.0; g_.W = 0.5;
g_.lead_ycen = [3.7, -2.9, 1.7, -1.6];        % 4 leads => port order 1..4
g_.ell_semi  = [0.75, 0.5];  g_.ell_rot = 110;
g_.ell_cens  = [ 1.2, 0.0, -1.3 ; 2.4, 4.0, 1.5 ];            % 3 ellipses
g_.numell    = size(g_.ell_cens,2);
g_.maxchunklen = lambda/4;  g_.legk = 16;
numch = numel(g_.lead_ycen);

%% --------------------------------------------------------- build the geometry
t0 = tic;
cgrph = build_cavity_chunkgraph(g_);
t_geom = toc(t0);
npt   = cgrph.npt;  nedge = numel(cgrph.echnks);

% Impedance edges = chord flats [1 5 9 13 17] + arc [18] (= 2+4*numch):
imp_edges = [1 + 4*(0:numch), 2 + 4*numch];
imp_pt = edge_point_mask(cgrph, imp_edges);
fprintf('Geometry: %d edges (%d impedance: [%s]), npt=%d\n', ...
        nedge, numel(imp_edges), num2str(imp_edges), npt);

%% ------------------------------- per-edge kernel matrix, scaled to (I + K) form
sp = kernel('h','sp',omega);  sk = kernel('h','s',omega);
kern_imp = 2*sp + (2*gamma)*sk;   % Robin (impedance) wall
kern_pec = 2*sp;                  % PEC (Neumann) wall
Kmat(nedge,nedge) = kernel();
for i = 1:nedge
    ki = kern_pec;  if ismember(i,imp_edges), ki = kern_imp; end
    for j = 1:nedge, Kmat(i,j) = ki; end
end

%% ------------------------------------------------- assemble system matrix (timed)
t0 = tic;
sysmat = chunkermat(cgrph, Kmat) + eye(npt);
t_assemble = toc(t0);
fprintf('System matrix assembled in %.2f s ; cond = %.2e\n', t_assemble, cond(sysmat));

%% --------------------------------------- VALIDATION GATE: manufactured solution
mms_relerr = robin_mms_check(cgrph, sysmat, omega, gamma, imp_pt);
fprintf('Manufactured-solution (Robin/Neumann) check: rel err = %.3e\n', mms_relerr);
assert(mms_relerr < 1e-6, 'Impedance BIE failed MMS gate (rel err %.2e).', mms_relerr);

%% ---------------------------------------------------- solve for the S-matrix (timed)
t0 = tic;
[Lf,Uf,Pf] = lu(sysmat);
ch = build_channels(cgrph, g_, omega);
field_coefs = zeros(numch);  field_der_coefs = zeros(numch);
for j = 1:numch
    rhs = zeros(npt,1);  rhs(ch(j).leadend_inds) = 2.0;
    sigma = Uf\(Lf\(Pf*rhs));
    for i = 1:numch
        u_cut  = chunkerkerneval(cgrph, sk, sigma, ch(i).XY).';
        ux_cut = chunkerkerneval(cgrph, sp, sigma, ch(i).XY_ptinfo).';
        field_coefs(i,j)     = ch(i).proj * u_cut.';
        field_der_coefs(i,j) = ch(i).proj * ux_cut.';
    end
end
beta  = arrayfun(@(s) s.beta, ch);
a_in  = (field_coefs + (1./(1i*beta(:))).*field_der_coefs)/2;
a_out = (field_coefs - (1./(1i*beta(:))).*field_der_coefs)/2;
S = a_out / a_in;
t_solve = toc(t0);
t_total = t_geom + t_assemble + t_solve;

%% ----------------------------------------------------------------- report
unit_defect = norm(S'*S - eye(numch));
absorb = 1 - sum(abs(S).^2,1);
fprintf('\n=== chunkie S-matrix (LOSSY impedance, TM, %d ports) ===\n', numch);  disp(S);
fprintf('|S| =\n');  disp(abs(S));
fprintf('||S^H S - I|| = %.3e   (lossy => >0)\n', unit_defect);
fprintf('absorption per input channel = [%s]\n', num2str(absorb,'%.4f '));
fprintf('symmetry ||S - S^T|| = %.3e  (reciprocal cavity)\n', norm(S - S.'));
if any(absorb < -1e-6)
    warning('Negative absorption => gamma sign wrong for this convention.');
end
fprintf('\n--- timing [s] ---  geom %.3f | assemble %.3f | solve %.3f | TOTAL %.3f (npt=%d)\n',...
        t_geom, t_assemble, t_solve, t_total, npt);

%% ----------------------------------------------------------------- save
writematrix(S, fullfile(outdir,'S_lossy_chunkie.csv'));
meta = struct('solver','chunkie-BIE','pol','TM','bc','impedance+PEC','nports',numch, ...
              'omega',omega,'lambda',lambda,'f_hz',f_hz,'sigma_cond',sigma_cond, ...
              'Zs',Zs,'gamma',gamma,'imp_edges',imp_edges, ...
              'lead_ycen',g_.lead_ycen,'ell_cens',g_.ell_cens,'ell_semi',g_.ell_semi, ...
              'ell_rot',g_.ell_rot,'unit_defect',unit_defect,'absorption',absorb,'npt',npt, ...
              't_geom',t_geom,'t_assemble',t_assemble,'t_solve',t_solve,'t_total',t_total);
save(fullfile(outdir,'S_lossy_chunkie.mat'),'S','meta');
fprintf('Saved S_lossy_chunkie.{mat,csv} in %s\n', outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            LOCAL  FUNCTIONS  (identical to the lossless script)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cgrph = build_cavity_chunkgraph(g_)
R=g_.R; d=g_.d; L=g_.L; W=g_.W; yc=g_.lead_ycen;
a60=deg2rad(60); a300=deg2rad(300);
Atop=R*[cos(a60);sin(a60)]; Abot=R*[cos(a300);sin(a300)];
xin=d; xout=d+L;  ys = sort(yc);
verts = [xin; Abot(2)];
for m = 1:numel(ys)
    yb=ys(m)-W/2; yt=ys(m)+W/2;
    verts = [verts, [xin;yb], [xout;yb], [xout;yt], [xin;yt]]; %#ok<AGROW>
end
verts = [verts, [xin; Atop(2)]];   nv = size(verts,2);
edgeinc = [1:nv ; circshift(1:nv,-1)];
edgefuns = cell(1,nv);
for i=1:nv-1, edgefuns{i} = @(t)(1-t(:).').*verts(:,i)+t(:).'.*verts(:,i+1); end
edgefuns{nv} = @(t) R*[cos(a60+t(:).'*(a300-a60)); sin(a60+t(:).'*(a300-a60))];
Rr=[cosd(g_.ell_rot),-sind(g_.ell_rot); sind(g_.ell_rot),cosd(g_.ell_rot)];
for q=1:g_.numell
    cq=g_.ell_cens(:,q);  edgeinc=[edgeinc, nan(2,1)]; %#ok<AGROW>
    edgefuns{end+1}=@(t) cq + Rr*[ g_.ell_semi(1)*cos(2*pi*t(:).'); ...
                                  -g_.ell_semi(2)*sin(2*pi*t(:).') ]; %#ok<AGROW>
end
cp=[]; cp.eps=1e-10; cp.nover=0; cp.maxchunklen=g_.maxchunklen;
pr=[]; pr.k=g_.legk;
cgrph = chunkgraph(verts, edgeinc, edgefuns, cp, pr);
end

%-----------------------------------------------------------------------------
function mask = edge_point_mask(cgrph, edge_list)
np = arrayfun(@(e) e.npt, cgrph.echnks);
mask = false(cgrph.npt,1);  off = 0;
for e = 1:numel(np)
    idx = off+(1:np(e));
    if ismember(e, edge_list), mask(idx) = true; end
    off = off + np(e);
end
end

%-----------------------------------------------------------------------------
function relerr = robin_mms_check(cgrph, sysmat, omega, gamma, imp_pt)
xy = cgrph.r(:,:); nn = cgrph.n(:,:);
z0 = [ -14.0 ; 9.0 ];
r0 = @(X) vecnorm(X - z0);
uman  = @(X) (1i/4)*besselh(0,1,omega*r0(X));
gradu = @(X) (1i/4)*(-omega)*besselh(1,1,omega*r0(X)).*((X-z0)./r0(X));
un = sum(nn.*gradu(xy),1).';
ub = uman(xy).';
data = un + gamma*(imp_pt.*ub);
sigma = sysmat\(2*data);
targ = [ -3.0, -1.0,  0.0,  1.0, -2.0 ;
         -3.0, -2.0, -3.0, -2.0, -1.0 ];
sk = kernel('h','s',omega);
ucalc = chunkerkerneval(cgrph, sk, sigma, targ);
utrue = uman(targ).';
relerr = norm(ucalc-utrue)/norm(utrue);
end

%-----------------------------------------------------------------------------
function ch = build_channels(cgrph, g_, omega)
boundary_xy = cgrph.r(:,:);
xout = g_.d + g_.L;  W = g_.W;  yc = g_.lead_ycen;
xcut = g_.d + 0.8*g_.L;  Ny = 21;
for j = 1:numel(yc)
    cen = yc(j);
    ch(j).leadend_inds = ( abs(boundary_xy(1,:)-xout)<1e-5 & ...
                           abs(boundary_xy(2,:)-cen)<=W/2 + 1e-9 ).';
    ylin = cen + linspace(-0.4*W, 0.4*W, Ny);
    ch(j).XY = [ xcut*ones(1,Ny) ; ylin ];
    pi_=[]; pi_.r=ch(j).XY; pi_.d=repmat([0;1],1,Ny);
    pi_.d2=zeros(2,Ny); pi_.n=repmat([1;0],1,Ny);
    ch(j).XY_ptinfo = pi_;
    phi = ones(Ny,1);  ch(j).proj = (phi.')/(phi.'*phi);
    ch(j).beta = omega;
end
end
