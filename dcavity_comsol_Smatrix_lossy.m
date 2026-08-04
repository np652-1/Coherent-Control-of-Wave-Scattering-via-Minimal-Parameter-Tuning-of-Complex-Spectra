%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  dcavity_comsol_Smatrix_lossy.m
%
%  Scattering matrix of the LOSSY (impedance-wall) D-shaped reverberation cavity,
%  computed with COMSOL Multiphysics via LiveLink for MATLAB.  This is the
%  finite-element counterpart of dcavity_chunkie_Smatrix_lossy.m and is used to
%  cross-check the boundary-integral solver.
%
%  PHYSICS / GEOMETRY  (identical to the chunkie script)
%    2D, TM polarization: the Electromagnetic Waves, Frequency Domain (ewfd)
%    interface with an in-plane E-field (=> the scalar wave field is H_z).
%    D-shaped cavity: circular arc (R, 60->300 deg) + flat chord wall at x=R/2
%    carrying 4 TEM waveguide leads/ports; NELL PEC ellipse scatterers inside.
%    Lossy walls: a Leontovich IMPEDANCE boundary (user-defined conductivity
%    sigma) on the arc + the chord segments; PEC on the lead walls + ellipses.
%    => a sub-unitary 4x4 S-matrix (the walls absorb).
%
%  HOW TO RUN
%    1) Start a COMSOL server in a terminal:   comsol mphserver -port 2036
%       (or launch MATLAB itself via:          comsol mphserver matlab )
%    2) Set COMSOL_MLI below to your COMSOL "mli" directory, then run this file.
%
%  OUTPUT (into dcavity_comsol_Smatrix_lossy_outputs/):
%    S_lossy_comsol.mat / .csv   (the 4x4 S-matrix + metadata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

%% -------------------------------------------------- connect to COMSOL LiveLink
COMSOL_MLI  = '/usr/local/comsol61/multiphysics/mli';   % <-- your COMSOL LiveLink "mli" dir
COMSOL_PORT = 2036;                                     % <-- port of your `comsol mphserver`
addpath(COMSOL_MLI);
import com.comsol.model.*
import com.comsol.model.util.*
% connect to the running COMSOL server, unless this MATLAB is already linked
% (e.g. it was launched via `comsol mphserver matlab`):
try, ModelUtil.getModelNames(); linked = true; catch, linked = false; end
if ~linked, mphstart(COMSOL_PORT); end

thisdir = fileparts(mfilename('fullpath')); if isempty(thisdir), thisdir = pwd; end
outdir  = fullfile(thisdir,'dcavity_comsol_Smatrix_lossy_outputs');
if ~exist(outdir,'dir'), mkdir(outdir); end

%% ---------------------------------------------------------------- parameters
P = struct();
c_light  = 299792458;
P.omega  = 0.97*(2*pi);
P.lambda = 2*pi/P.omega;
P.f_hz   = c_light/P.lambda;
P.R   = 5.75;  P.d = P.R/2;  P.L = 4.0;  P.W = 0.5;
P.lead_ycen = [3.7, -2.9, 1.7, -1.6];            % 4 leads (port order 1..4)
P.ell_semi  = [0.75, 0.5];
P.ell_rot_deg = 110;
P.ell_cens  = [ 1.2, 0.0, -1.3 ; 2.4, 4.0, 1.5 ];   % 3 ellipses
P.minppw = 20;                                    % mesh: elements per wavelength
P.sigma_cond = 50;                                % LOSSY: impedance-wall conductivity [S/m]
nP = numel(P.lead_ycen);

fprintf('COMSOL LOSSY D-cavity: f=%.6e Hz, lambda=%.5f, sigma=%.3g, %d ports\n', ...
        P.f_hz, P.lambda, P.sigma_cond, nP);

%% --------------------------------------------------- build the model (timed)
t0 = tic;
model = build_model(P);
t_build = toc(t0);
fprintf('Model build (geom+mesh+setup): %.3f s\n', t_build);

% diagnostics: expected surface impedance
eps0 = 8.8541878128e-12; mu0 = 1.25663706212e-6;  omega_phys = c_light*P.omega;
Zs_expect = sqrt(1i*omega_phys*mu0/(P.sigma_cond + 1i*omega_phys*eps0));
fprintf('Expected Zs = %.4g%+.4gi Ohm  (|Zs|/Z0=%.3e)\n', ...
        real(Zs_expect),imag(Zs_expect), abs(Zs_expect)/sqrt(mu0/eps0));

%% ----------------------------------------- full S-matrix (excite each lead, timed)
t0 = tic;
S = compute_Smatrix(model, nP);
t_solve = toc(t0);

%% ----------------------------------------------------------------- report
unit_defect = norm(S'*S - eye(nP));
absorb = 1 - sum(abs(S).^2,1);
fprintf('\n=== COMSOL S-matrix (LOSSY impedance, TM, %d ports) ===\n', nP);  disp(S);
fprintf('|S| =\n'); disp(abs(S));
fprintf('||S^H S - I|| = %.3e   (lossy => >0)\n', unit_defect);
fprintf('absorption per input channel = [%s]\n', num2str(absorb,'%.4f '));
fprintf('\n--- timing [s] ---  build %.3f | solve(%d excitations) %.3f | TOTAL %.3f\n',...
        t_build, nP, t_solve, t_build+t_solve);

%% ----------------------------------------------------------------- save
writematrix(S, fullfile(outdir,'S_lossy_comsol.csv'));
meta = struct('solver','COMSOL-ewfd','pol','TM-inplane','bc','impedance+PEC','nports',nP, ...
              'omega',P.omega,'lambda',P.lambda,'f_hz',P.f_hz,'sigma_cond',P.sigma_cond, ...
              'lead_ycen',P.lead_ycen,'ell_cens',P.ell_cens,'ell_semi',P.ell_semi, ...
              'ell_rot_deg',P.ell_rot_deg,'minppw',P.minppw, ...
              'unit_defect',unit_defect,'absorption',absorb, ...
              't_build',t_build,'t_solve',t_solve,'t_total',t_build+t_solve);
save(fullfile(outdir,'S_lossy_comsol.mat'),'S','meta');
fprintf('Saved S_lossy_comsol.{mat,csv} in %s\n', outdir);
% model left in workspace; free with ModelUtil.remove('Model') when done.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            LOCAL  FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = compute_Smatrix(model, nP)
% Excite each lead in turn (that lead = port 1 / active input, the rest passive)
% and read the S-parameters ewfd.S(k)1 from the global-evaluation table.
S = zeros(nP,nP);  order = 1:nP;
for j = 1:nP
    others = order(order~=j);  sel = [j, others];   % make lead j the driven port
    for k = 1:nP
        model.component('comp1').physics('ewfd').feature(sprintf('port%d',k)).selection.named(sprintf('leadbox%d',sel(k)));
    end
    model.sol('sol1').runAll;
    model.result.numerical('gev1').setResult;
    tbl = mphtable(model,'tbl1');
    row = tbl.data(end,:);                 % [freq, S11, S21, ..., S(nP)1]
    for k = 1:nP, S(sel(k), j) = row(1+k); end
    fprintf('  excited lead %d done\n', j);
end
end

%-----------------------------------------------------------------------------
function model = build_model(P)
import com.comsol.model.*
import com.comsol.model.util.*
model = ModelUtil.create('Model');
model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 2);
nL = numel(P.lead_ycen);

model.param.set('R', num2str(P.R)); model.param.set('d','R/2');
model.param.set('Lw', num2str(P.L)); model.param.set('Ww', num2str(P.W));

% ----- geometry: arc + chord -> solid; subtract-in lead rectangles; cut ellipses
g = model.component('comp1').geom('geom1');
g.create('ca1','CircularArc');
g.feature('ca1').set('r','R'); g.feature('ca1').set('angle1',60); g.feature('ca1').set('angle2',300);
g.create('ls1','LineSegment');
g.feature('ls1').set('specify1','coord'); g.feature('ls1').set('coord1',{'R/2' 'sqrt(3)*R/2'});
g.feature('ls1').set('specify2','coord'); g.feature('ls1').set('coord2',{'R/2' '-sqrt(3)*R/2'});
g.create('csol1','ConvertToSolid'); g.feature('csol1').selection('input').set({'ca1' 'ls1'});
ycen = P.lead_ycen;
leadnames = cell(1,nL);
for k=1:nL
    rn = sprintf('r%d',k); leadnames{k}=rn; g.create(rn,'Rectangle');
    g.feature(rn).set('size',{'Lw' 'Ww'}); g.feature(rn).set('base','center');
    g.feature(rn).set('pos',{'d + Lw/2', num2str(ycen(k))});
end
ellnames = cell(1,size(P.ell_cens,2));
for q=1:size(P.ell_cens,2)
    en = sprintf('e%d',q); ellnames{q}=en; g.create(en,'Ellipse');
    g.feature(en).set('semiaxes',[P.ell_semi(1) P.ell_semi(2)]);
    g.feature(en).set('pos',[P.ell_cens(1,q) P.ell_cens(2,q)]);
    g.feature(en).set('rot', P.ell_rot_deg);
end
g.create('dif1','Difference');
g.feature('dif1').selection('input').set([{'csol1'} leadnames]);
g.feature('dif1').selection('input2').set(ellnames);
g.run;

xin = P.d;  xout = P.d + P.L;  arcY = P.R*sqrt(3)/2;
% ----- lead-end (port) edge selections
for k=1:nL
    sn = sprintf('leadbox%d',k);
    model.component('comp1').selection.create(sn,'Box');
    model.component('comp1').selection(sn).set('entitydim',1);
    model.component('comp1').selection(sn).set('xmin', xout-0.05); model.component('comp1').selection(sn).set('xmax', xout+0.05);
    model.component('comp1').selection(sn).set('ymin', ycen(k)-0.30); model.component('comp1').selection(sn).set('ymax', ycen(k)+0.30);
    model.component('comp1').selection(sn).set('condition','allvertices');
end

% ----- IMPEDANCE-wall edge selection = arc + (nL+1) flat chord segments
model.component('comp1').selection.create('imp_arc','Box');
model.component('comp1').selection('imp_arc').set('entitydim',1);
model.component('comp1').selection('imp_arc').set('xmin',-P.R-0.2); model.component('comp1').selection('imp_arc').set('xmax',-4.5);
model.component('comp1').selection('imp_arc').set('ymin',-3.5); model.component('comp1').selection('imp_arc').set('ymax',3.5);
model.component('comp1').selection('imp_arc').set('condition','intersects');
ys = sort(ycen);
fy_lo = [-arcY, ys + P.W/2];        % lower edges of the (nL+1) flat segments
fy_hi = [ys - P.W/2,  arcY];        % upper edges
flatnames = cell(1,nL+1);
for q=1:nL+1
    sn = sprintf('imp_f%d',q); flatnames{q}=sn;
    model.component('comp1').selection.create(sn,'Box');
    model.component('comp1').selection(sn).set('entitydim',1);
    model.component('comp1').selection(sn).set('xmin',xin-0.03); model.component('comp1').selection(sn).set('xmax',xin+0.03);
    model.component('comp1').selection(sn).set('ymin',fy_lo(q)-0.1); model.component('comp1').selection(sn).set('ymax',fy_hi(q)+0.1);
    model.component('comp1').selection(sn).set('condition','inside');
end
model.component('comp1').selection.create('impsel','Union');
model.component('comp1').selection('impsel').set('entitydim',1);
model.component('comp1').selection('impsel').set('input',[{'imp_arc'} flatnames]);

% ----- material (vacuum domain; the wall loss is a boundary property below)
model.component('comp1').material.create('mat1','Common');
model.component('comp1').material('mat1').propertyGroup('def').set('relpermittivity',{'1'});
model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability',{'1'});
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity',{'0'});

% ----- mesh
hmax = P.lambda/P.minppw;
mesh1 = model.component('comp1').mesh.create('mesh1');
ftri1 = mesh1.create('ftri1','FreeTri'); ftri1.create('size1','Size');
ftri1.feature('size1').set('custom','on');
ftri1.feature('size1').set('hmax', num2str(hmax)); ftri1.feature('size1').set('hmaxactive', true);
mesh1.run;

% ----- physics: ewfd (TM), 4 TEM ports, impedance BC on the outer enclosure
phys = model.component('comp1').physics.create('ewfd','ElectromagneticWavesFrequencyDomain','geom1');
phys.prop('components').set('components','inplane');     % TM: in-plane E => scalar H_z
for k=1:nL
    pn=sprintf('port%d',k); phys.create(pn,'Port',1);
    phys.feature(pn).selection.named(sprintf('leadbox%d',k));
    phys.feature(pn).set('PortType','Rectangular'); phys.feature(pn).set('PortModeType','TEM');
end
phys.create('imp1','Impedance',1);
phys.feature('imp1').selection.named('impsel');
phys.feature('imp1').set('DisplacementFieldModel','RelativePermittivity');
phys.feature('imp1').set('epsilonr_mat','userdef'); phys.feature('imp1').set('epsilonr',1);
phys.feature('imp1').set('murbnd_mat','userdef');   phys.feature('imp1').set('murbnd',1);
phys.feature('imp1').set('sigmabnd_mat','userdef'); phys.feature('imp1').set('sigmabnd', P.sigma_cond);

% ----- study + solver (direct, complex-valued frequency-domain)
model.study.create('std1'); model.study('std1').create('freq','Frequency');
model.study('std1').feature('freq').set('plist', num2str(P.f_hz));
model.study('std1').feature('freq').set('punit','Hz');
model.study('std1').setGenConv(false); model.study('std1').setGenPlots(false);
model.sol.create('sol1'); model.sol('sol1').study('std1');
model.sol('sol1').create('st1','StudyStep'); model.sol('sol1').feature('st1').set('study','std1');
model.sol('sol1').feature('st1').set('studystep','freq');
model.sol('sol1').create('v1','Variables'); model.sol('sol1').feature('v1').set('control','freq');
model.sol('sol1').create('s1','Stationary'); model.sol('sol1').feature('s1').set('stol','1e-6');
model.sol('sol1').feature('s1').create('p1','Parametric');
model.sol('sol1').feature('s1').feature.remove('pDef');
model.sol('sol1').feature('s1').feature('p1').set('pname',{'freq'});
model.sol('sol1').feature('s1').feature('p1').set('plistarr', num2str(P.f_hz));
model.sol('sol1').feature('s1').feature('p1').set('punit',{'Hz'});
model.sol('sol1').feature('s1').feature('aDef').set('complexfun', true);
model.sol('sol1').feature('s1').create('fc1','FullyCoupled');
model.sol('sol1').feature('s1').create('d1','Direct');
model.sol('sol1').feature('s1').feature('d1').set('linsolver','pardiso');
model.sol('sol1').feature('s1').feature('fc1').set('linsolver','d1');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').attach('std1');

% ----- S-parameter global evaluation:  ewfd.S11, S21, ..., S(nL)1
exprs = arrayfun(@(k) sprintf('ewfd.S%d1',k), 1:nL, 'uni',0);
model.result.numerical.create('gev1','EvalGlobal'); model.result.numerical('gev1').set('data','dset1');
model.result.numerical('gev1').set('expr', exprs);
model.result.table.create('tbl1','Table'); model.result.numerical('gev1').set('table','tbl1');
end
