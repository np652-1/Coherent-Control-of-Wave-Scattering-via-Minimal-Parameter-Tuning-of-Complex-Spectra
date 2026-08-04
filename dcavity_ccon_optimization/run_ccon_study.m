%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  run_ccon_study.m   (6 CCON determinant-FOM studies, ONE unified job array)
%
%  ** LOSSLESS re-run: identical to statistical_study_2 / statistical_study_3,
%     EXCEPT ALL cavity walls are PEC -- no impedance boundary condition anywhere
%     (gamma=0; see ccon_config C.lossless).  Same 7 cases (RTTT,NDTT,NNDD,RRRR,
%     RDDT,RDTT,NDDD), same random N-subset positions / random freq / port
%     permutation / det(C'C) FOM / v6 backend, so results drop straight into the
%     analysis alongside the lossy studies for a loss vs no-loss comparison. **
%
%  Generalizes statistical_study_1_NDDD/to_upload/statistical_study_v2.m to a
%  family of figures of merit
%             F(theta) = det(C'*C),   C = S[rows,cols]  (a sub-matrix of S),
%  one per CCON letter word (RTTT, NDTT, NNDD, RDDT, RDTT, NDDD; RRRR dropped; see
%  ccon_cases.m and ali_paper.pdf).  Same un-squared functional form as the
%  completed NDDD study.  The VALIDATED v6 fast backend (Schur + precomputed field-eval +
%  multipole M_CE/M_EC) and its matching fast adjoint gradient are reused
%  VERBATIM; only the FOM layer (fom_value / fom_cotangents / fom_Sbar) is new.
%
%  Per realization (reproducible from row.seed): random N_scat-subset of the 12
%  sites, random frequency in [0.96,0.98] f_c, and a RANDOM port-role permutation
%  (which physical lead plays each CCON letter).  The N_scat ellipse angles are
%  optimized: multi-start fmincon (box [0,2pi], analytic adjoint gradient,
%  5*N_scat restarts) + fminsearch polish (200*N_scat) -- identical to results_4.
%
%  USAGE
%    run_ccon_study()                 % production: array task -> manifest rows
%    run_ccon_study('verify')         % adjoint-vs-FD gradient gate, all 6 cases
%    run_ccon_study('verify','RTTT')  % ... just one case
%    run_ccon_study('smoke')          % one tiny optimization per case (quick)
%    run_ccon_study('smoke','NNDD')   % ... just one case
%
%  OUTPUTS (production)
%    results/<CASE>/summary_k####.mat       -- compact res struct (for analysis)
%    results/<CASE>/intermediate_k####.mat  -- full per-eval log (angles,F,S)
%  Both carry the case name + permutation, so outputs are trivially disentangled.
%
%  chunkie expected at ../../chunkie (dev fallback ~/main_projects/Numerics/chunkie).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = run_ccon_study(mode, varargin)
if nargin<1 || isempty(mode), mode='run'; end
ensure_chunkie();
switch lower(mode)
    case 'run',     ccon_run_production();
    case 'verify',  ok=ccon_verify_gradient(varargin{:}); if nargout, varargout{1}=ok; end
    case 'smoke',   ccon_smoke(varargin{:});
    case 'example', ccon_example(varargin{:});      % restricted 0/50/90% loss demo (~10 min)
    otherwise, error('run_ccon_study: unknown mode "%s"',mode);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        shared config
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = ccon_config()
C.c=299792458; C.eps0=8.8541878128e-12; C.mu0=1.25663706212e-6;
C.R=5.75; C.d=C.R/2; C.L=4.0; C.W=0.5; C.lead_ycen=[3.7,-2.9,1.7,-1.6];
C.ell_semi=[0.75,0.5]; C.numch=4;
% LOSSLESS study: ALL cavity walls are PEC -- no impedance BC anywhere.
% This forces gamma=0 in precompute_realization (zero surface impedance).
% Identical to statistical_study_2 / _3 in every other respect; sigma is left
% only as an unused informational value (it plays no role when lossless=true).
C.lossless=true; C.sigma=Inf;   % PEC walls, gamma=0
C.nSites=12;
C.Positions=[1.2,2.4; 0,4; -1.3,1.5; -3.3,2.2; 0.2,-1.6; -3.1,-2.6; ...
             -4.3,-0.1; 1.2,-4.2; 1.3,0.4; -1.4,-3.9; -1.8,-0.4; -1.8,3.7].';
C.f_c=C.c/(2*C.W); C.f_min=0.96*C.f_c; C.f_max=0.98*C.f_c;
C.mcl_factor=4; C.mp_tol=1e-12;
end

function opt = ccon_opt(quick)
if quick
    opt.restart_per_dim=3; opt.restart_cap=6; opt.polish_fev_per_dim=20;
    opt.fmc=optimoptions('fmincon','Algorithm','interior-point','SpecifyObjectiveGradient',true,...
              'Display','off','MaxIterations',8,'MaxFunctionEvaluations',60,'OptimalityTolerance',1e-6);
    opt.fms=optimset('Display','off','TolFun',1e-8,'TolX',1e-4);
else
    % IDENTICAL to Ali / results_4
    opt.restart_per_dim=5; opt.restart_cap=inf; opt.polish_fev_per_dim=200;
    opt.fmc=optimoptions('fmincon','Algorithm','interior-point','SpecifyObjectiveGradient',true,...
              'Display','off','MaxIterations',500,'MaxFunctionEvaluations',10000,'OptimalityTolerance',1e-10);
    opt.fms=optimset('Display','off','TolFun',1e-12,'TolX',1e-6);
end
end

function dirs = ccon_dirs(sub)
if nargin<1, sub='results'; end
dirs.main = fileparts(mfilename('fullpath')); if isempty(dirs.main), dirs.main=pwd; end
dirs.results = fullfile(dirs.main, sub);
ensure_dir(dirs.results);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        PRODUCTION runner (one unified job array over the global manifest)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ccon_run_production()
base_seed=2026; R_per=10;
C = ccon_config(); opt = ccon_opt(false); cases = ccon_cases();
dirs = ccon_dirs('results');
man = job_manifest(R_per, base_seed); R = numel(man);

setup_pool();
[myrows, tag] = array_partition(R); nrow = numel(myrows);
fprintf('\n=== CCON study%s: %d of R=%d realizations (6 cases (RRRR dropped), LOSSLESS/PEC walls) ===\n', tag, nrow, R);
t0=tic;
for ii=1:nrow
    row = man(myrows(ii));
    r_ = process_realization(row, C, opt, dirs, cases);
    fprintf('  [%s k=%d] N=%d perm=%s C=S[%s;%s] best-restart=%.4g polished=%.4g (%.0fs) %s\n',...
        r_.name, r_.k, r_.N_scat, mat2str(r_.perm), mat2str(r_.rows), mat2str(r_.cols),...
        r_.F_best_restart, r_.F_polished, r_.walltime, r_.err);
end
fprintf('task wall time: %.1f s.  DONE.\n', toc(t0));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        ONE REALIZATION (serial inside; parfor over restarts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = process_realization(row, C, opt, dirs, cases)
ensure_chunkie(); tR=tic;
rng(row.seed,'twister');
N = row.N; cdef = cases(row.ci);
cfg = draw_realization(C, N, cdef);
n_restarts = min(opt.restart_per_dim*N, opt.restart_cap);
lb=zeros(1,N); ub=2*pi*ones(1,N);
fms = optimset(opt.fms,'MaxIter',opt.polish_fev_per_dim*N,'MaxFunEvals',opt.polish_fev_per_dim*N);
caseDir = fullfile(dirs.results, row.name); ensure_dir(caseDir);
try
    P = precompute_realization(C, cfg);
catch ME
    res = base_res(row,cfg,N); res.err = ME.message; res.walltime = toc(tR);
    save(fullfile(caseDir,sprintf('summary_k%04d.mat',row.k)),'res','-v7.3');
    return;
end
x0s = 2*pi*rand(n_restarts, N);
xr=zeros(n_restarts,N); fr=zeros(n_restarts,1); rlog=cell(1,n_restarts);
fom=cfg.fom; fmc=opt.fmc;
parfor r=1:n_restarts
    [xr(r,:),fr(r),rlog{r}] = run_restart(P, fom, x0s(r,:), lb, ub, fmc);
end
[~,ib]=min(fr); xbest=xr(ib,:);
[xpol,fpol,plog] = run_polish(P, fom, xbest, fms); xpol=mod(xpol,2*pi);
[~,S_best] = fom_only(P, xpol, fom);

% per-evaluation LOG: restart evals tagged rid=r, polish rid=-1
LOG=struct('rid',[],'x',{{}},'F',[],'S',{{}});
for r=1:n_restarts
    m=numel(rlog{r}.F);
    LOG.rid=[LOG.rid, r*ones(1,m)]; LOG.x=[LOG.x, rlog{r}.x]; LOG.F=[LOG.F, rlog{r}.F]; LOG.S=[LOG.S, rlog{r}.S];
end
mp=numel(plog.F);
LOG.rid=[LOG.rid, -ones(1,mp)]; LOG.x=[LOG.x, plog.x]; LOG.F=[LOG.F, plog.F]; LOG.S=[LOG.S, plog.S];

res = base_res(row,cfg,N);
res.omega=P.omega; res.gamma=P.gamma; res.mcl=P.mcl;
res.x_best=xbest; res.x_polished=xpol; res.F_polished=fpol; res.F_best_restart=min(fr);
res.S_best=S_best; res.walltime=toc(tR); res.err='';
save_intermediate(caseDir,row,C,cfg,P,LOG,x0s,xr,fr,xbest,xpol,fpol);
save(fullfile(caseDir,sprintf('summary_k%04d.mat',row.k)),'res','-v7.3');
end

function res = base_res(row,cfg,N)
res=struct('name',row.name,'ci',row.ci,'k',row.k,'gi',row.gi,'N_scat',N,'np',[],'Nmax',row.Nmax,...
    'word',cfg.fom.word,'rows',cfg.fom.rows,'cols',cfg.fom.cols,'perm',cfg.fom.perm,...
    'positions',cfg.positions,'sel',cfg.sel,'f',cfg.f,'seed',row.seed,...
    'omega',NaN,'gamma',NaN,'mcl',NaN,'x_best',nan(1,N),'x_polished',nan(1,N),...
    'F_polished',NaN,'F_best_restart',NaN,'S_best',nan(4,4),'walltime',NaN,'err','');
end

% ----- one multi-start restart (parfor worker; returns its eval log) -----
function [xopt,fopt,lg] = run_restart(P, fom, x0, lb, ub, fmc)
lg=struct('x',{{}},'F',[],'S',{{}});
[xopt,fopt]=fmincon(@obj,x0,[],[],[],[],lb,ub,[],fmc);
    function [F,g]=obj(x), [F,g,S]=fom_and_grad(P,x,fom); lg.x{end+1}=x(:); lg.F(end+1)=F; lg.S{end+1}=S; end
end
% ----- fminsearch polish (serial; returns its eval log) -----
function [xopt,fopt,lg] = run_polish(P, fom, x0, fms)
lg=struct('x',{{}},'F',[],'S',{{}});
[xopt,fopt]=fminsearch(@obj,x0,fms);
    function F=obj(x), [F,S]=fom_only(P,x,fom); lg.x{end+1}=x(:); lg.F(end+1)=F; lg.S{end+1}=S; end
end

function save_intermediate(caseDir,row,C,cfg,P,LOG,x0s,xr,fr,xbest,xpol,fpol)
n=numel(LOG.F); N=cfg.N_scat; nch=C.numch;
eval_x=zeros(N,n); for i=1:n, eval_x(:,i)=LOG.x{i}; end
eval_S=zeros(nch,nch,n); for i=1:n, eval_S(:,:,i)=LOG.S{i}; end
realiz=struct(); realiz.name=row.name; realiz.ci=row.ci; realiz.k=row.k; realiz.gi=row.gi;
realiz.seed=row.seed; realiz.cfg=cfg;
realiz.omega=P.omega; realiz.gamma=P.gamma; realiz.mcl=P.mcl; realiz.lambda=P.lambda;
realiz.eval_rid=LOG.rid; realiz.eval_x=eval_x; realiz.eval_F=LOG.F; realiz.eval_S=eval_S;
realiz.restart_inits=x0s; realiz.restart_x=xr; realiz.restart_F=fr;
realiz.x_best=xbest; realiz.x_polished=xpol; realiz.F_polished=fpol;
C_global=C; %#ok<NASGU>
save(fullfile(caseDir,sprintf('intermediate_k%04d.mat',row.k)),'realiz','C_global','-v7.3');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        VERIFY: adjoint gradient vs central finite difference (the gate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = ccon_verify_gradient(only_name)   % gate over all 7 cases (incl. NDDD)
C = ccon_config(); cases = ccon_cases(); base_seed=2026;
if nargin>=1 && ~isempty(only_name)
    cases = cases(strcmpi({cases.name},only_name));
    assert(~isempty(cases),'no such case "%s"',only_name);
end
Ncheck=3; dth=1e-6; tol=1e-5; buf='';
prn('CCON adjoint-gradient vs central finite difference (det(C''C) FOM)\n');
prn('================================================================\n');
maxerr_all=0;
for ci=1:numel(cases)
    cdef=cases(ci); N=min(Ncheck,cdef.Nmax);
    rng(base_seed+777+ci,'twister');
    cfg=draw_realization(C,N,cdef);
    P=precompute_realization(C,cfg);
    ang=2*pi*rand(1,N);
    [F0,gA]=fom_and_grad(P,ang,cfg.fom);
    gF=zeros(1,N);
    for q=1:N
        ap=ang; ap(q)=ap(q)+dth; am=ang; am(q)=am(q)-dth;
        gF(q)=(fom_only(P,ap,cfg.fom)-fom_only(P,am,cfg.fom))/(2*dth);
    end
    relerr=abs(gA-gF)./max(abs(gF),eps); me=max(relerr); maxerr_all=max(maxerr_all,me);
    prn('\n%-5s  N=%d  C=S[%s;%s]  F=%.6g\n',cdef.name,N,mat2str(cfg.fom.rows),mat2str(cfg.fom.cols),F0);
    for q=1:N
        prn('  angle %d: adjoint=%+.8g  FD=%+.8g  rel.err=%.2e\n',q,gA(q),gF(q),relerr(q));
    end
    prn('  max rel.err (%s) = %.2e\n',cdef.name,me);
end
prn('\nOVERALL max rel.err = %.2e   (tol %.0e)  ->  %s\n',maxerr_all,tol,tern(maxerr_all<tol,'PASS','FAIL'));
% write the accumulated report robustly (-batch fopen to mfilename dir can fail)
dirs = ccon_dirs('results'); outp=fullfile(dirs.main,'gradient_check.txt');
fid=fopen(outp,'w'); if fid<0, outp=fullfile(pwd,'gradient_check.txt'); fid=fopen(outp,'w'); end
if fid>0, fprintf(fid,'%s',buf); fclose(fid); fprintf('[wrote %s]\n',outp);
else, warning('could not open gradient_check.txt for writing'); end
ok = maxerr_all<tol;
if ~ok, warning('GRADIENT GATE FAILED: max rel.err %.2e >= tol %.0e',maxerr_all,tol); end
    function prn(varargin), s=sprintf(varargin{:}); fprintf('%s',s); buf=[buf s]; end %#ok<AGROW>
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        SMOKE: one tiny optimization per case (writes files, checks F drops)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ccon_smoke(only_name)
C = ccon_config(); opt = ccon_opt(true); cases = ccon_cases(); base_seed=2026;
sel = 1:numel(cases);
if nargin>=1 && ~isempty(only_name)
    sel = find(strcmpi({cases.name},only_name)); assert(~isempty(sel),'no such case "%s"',only_name);
end
dirs = ccon_dirs('results_smoke');
setup_pool();
fprintf('\n=== CCON smoke (quick) ===\n');
for ci = sel
    N=2;
    row=struct('gi',ci,'ci',ci,'name',cases(ci).name,'N',N,'k',N,...
               'seed',base_seed+5000+ci,'Nmax',cases(ci).Nmax);
    r_ = process_realization(row, C, opt, dirs, cases);
    fprintf('  %-5s N=%d  polished F=%.4g  best-restart F=%.4g (%.0fs) %s\n',...
        r_.name, r_.N_scat, r_.F_polished, r_.F_best_restart, r_.walltime, r_.err);
end
fprintf('smoke done; files under %s\n', dirs.results);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        EXAMPLE: restricted 0% / 50% / 90%-loss CCON optimization (~10 min)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ccon_example(wordname, N_scat)
%CCON_EXAMPLE  A small, self-contained demonstration of the CCON optimization,
%  sized to run on a common laptop in ~10 minutes (NO cluster, NO Parallel
%  Computing Toolbox needed).
%
%  For each of three wall-loss levels (lossless / 50% / 90%) it draws ONE fixed
%  realization -- a random N-ellipse arrangement + random frequency + random
%  port-role permutation, all pinned by a fixed seed -- then minimizes the CCON
%  figure of merit
%        F(theta) = det(C'C),   C = S[rows,cols]
%  over the N ellipse rotation angles, using the SAME fast v6 forward solver and
%  analytic adjoint gradient as the production study, in "quick" optimizer mode
%  (a few restarts + a short polish).  It prints, and bar-plots, the FOM before
%  vs after optimization at each loss level.
%
%  Cost is dominated by the ONE-TIME boundary-matrix precompute (~30-90 s per
%  loss level depending on the machine); the per-evaluation solves are ~0.5 s.
%
%  Usage:  run_ccon_study('example')            % default: word NNDD, N=2 ellipses
%          run_ccon_study('example','NDDD',3)   % a different word / ellipse count
if nargin<1 || isempty(wordname), wordname='NNDD'; end   % headline CCON word (n_p=2)
if nargin<2 || isempty(N_scat),   N_scat  =2;      end    % 2 ellipses = 2 tunable angles
seed = 2026;  n_rand_before = 8;  disp_floor = 1e-12;     % FOM display floor (dB)

ensure_chunkie();
cases = ccon_cases(); opt = ccon_opt(true);              % QUICK optimizer settings
ci = find(strcmpi({cases.name}, wordname), 1);
assert(~isempty(ci), 'unknown CCON word "%s" (have: %s)', wordname, strjoin({cases.name},','));
cdef = cases(ci);
dirs = ccon_dirs('results_example');

levels = { 'lossless (0%)', true,  Inf     ; ...
           '50% loss',      false, 50      ; ...
           '90% loss',      false, 0.51527 };
nlev = size(levels,1);
res = struct('label',{},'sigma',{},'F_before',{},'F_after',{},'ratio',{},'S_after',{},'angles',{});

fprintf('\n=== CCON optimization example: word %s, N=%d ellipses, quick optimizer ===\n', wordname, N_scat);
fprintf('(3 loss levels; each = one fixed realization, optimized over %d ellipse angles)\n', N_scat);
t_all = tic;
for L = 1:nlev
    C = ccon_config();                                   % base D-cavity config
    C.lossless = levels{L,2};  C.sigma = levels{L,3};    % override the wall-loss level
    rng(seed,'twister');                                 % SAME realization at every loss level
    cfg = draw_realization(C, N_scat, cdef);
    tL = tic;
    P  = precompute_realization(C, cfg);                 % one-time boundary precompute
    % --- FOM BEFORE optimization: median over a few random angle configs ---
    Fr = zeros(1,n_rand_before);
    for t = 1:n_rand_before, Fr(t) = fom_only(P, 2*pi*rand(1,N_scat), cfg.fom); end
    F_before = median(Fr);
    % --- OPTIMIZE: multi-start fmincon (adjoint gradient) + fminsearch polish ---
    lb = zeros(1,N_scat); ub = 2*pi*ones(1,N_scat);
    nrs = min(opt.restart_per_dim*N_scat, opt.restart_cap);
    x0s = 2*pi*rand(nrs, N_scat); xr = zeros(nrs,N_scat); fr = zeros(nrs,1);
    for r = 1:nrs, [xr(r,:),fr(r)] = run_restart(P, cfg.fom, x0s(r,:), lb, ub, opt.fmc); end
    [~,ib] = min(fr);
    fms = optimset(opt.fms,'MaxIter',opt.polish_fev_per_dim*N_scat,'MaxFunEvals',opt.polish_fev_per_dim*N_scat);
    [xpol,F_after] = run_polish(P, cfg.fom, xr(ib,:), fms); xpol = mod(xpol,2*pi);
    [~,S_after] = fom_only(P, xpol, cfg.fom);
    ratio = F_before / max(F_after, realmin);
    fprintf('  [%-13s sigma=%-8.4g] f=%.4f GHz  F_before=%.4e  F_after=%.4e  drop=%.3gx  (%.0fs)\n',...
            levels{L,1}, C.sigma, cfg.f/1e9, F_before, F_after, ratio, toc(tL));
    res(L) = struct('label',levels{L,1},'sigma',C.sigma,'F_before',F_before,...
                    'F_after',F_after,'ratio',ratio,'S_after',S_after,'angles',xpol);
end
fprintf('example total wall time: %.0f s\n', toc(t_all));

% ---- save summary + before/after bar chart ----
save(fullfile(dirs.results,'ccon_example_summary.mat'),'res','wordname','N_scat','seed');
try
    Fb = max([res.F_before], disp_floor);  Fa = max([res.F_after], disp_floor);
    fig = figure('Color','w','Position',[100 100 760 520],'Visible','off');
    bar(10*log10([Fb(:) Fa(:)])); grid on;
    set(gca,'XTick',1:nlev,'XTickLabel',{res.label},'FontSize',12);
    ylabel('CCON FOM   10\cdotlog_{10} det(C^HC)   [dB]');
    legend({'random start (median)','optimized'},'Location','best');
    title(sprintf('CCON "%s" optimization (N=%d ellipses) vs wall loss', wordname, N_scat));
    exportgraphics(fig, fullfile(dirs.results,'ccon_example_before_after.png'),'Resolution',150);
    close(fig);
    fprintf('wrote ccon_example_before_after.png\n');
catch ME, warning('plot failed: %s', ME.message); end
fprintf('example done; summary + figure under %s\n', dirs.results);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        realization config + FOM bookkeeping (NEW: random port-role perm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cfg = draw_realization(C, N_scat, cdef)
sel=randperm(C.nSites,N_scat); cfg.positions=C.Positions(:,sel); cfg.N_scat=N_scat; cfg.sel=sel;
cfg.f=C.f_min+(C.f_max-C.f_min)*rand;
perm=randperm(C.numch);                 % random assignment of leads to CCON letters
cfg.fom=make_fom(cdef,perm);
end
function fom = make_fom(cdef, perm)
fom.word=cdef.word; fom.perm=perm;
fom.rows=perm(cdef.rows); fom.cols=perm(cdef.cols);
fom.name=sprintf('%s: min det(C''C), C=S[%s;%s]',cdef.word,mat2str(fom.rows),mat2str(fom.cols));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        FOM + cotangent (NEW generic det(C'C)) ; adjoint reused verbatim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = fom_value(S, fom)
Cm=S(fom.rows,fom.cols); F=real(det(Cm'*Cm));
end
function [F,S] = fom_only(P, angles, fom)
out=solve_v6(P,angles); S=out.S; F=fom_value(S,fom);
end
function G = fom_Sbar(S, fom)
% G = dF/dS (Wirtinger) for F = det(C'C), C=S[rows,cols].
% dF/dC = conj(C * adj(C'C)),  adj(M)=det(M) inv(M)  (finite at the optimum).
Cm=S(fom.rows,fom.cols); Gm=Cm'*Cm; ad=ccon_adjugate(Gm);
G=zeros(size(S)); G(fom.rows,fom.cols)=conj(Cm*ad);
end
function A=ccon_adjugate(M)
n=size(M,1);
if n==1, A=1;
elseif n==2, A=[M(2,2) -M(1,2); -M(2,1) M(1,1)];
else, A=det(M)*inv(M); end %#ok<MINV>
end
function [F,Wfc,Wfdc] = fom_cotangents(S,A,fom,omega)
% Propagate G=dF/dS back to the field-coefficient cotangents (FC,FDC), the only
% place the FOM enters the adjoint.  Derived from, and matched to, the original
% per-entry leakage cotangent loop:
%   Wfc  =  (1/2)   (I - S.') G (ain.')^{-1}
%   Wfdc = -(1/2iw) (I + S.') G (ain.')^{-1}
nch=size(S,1); F=fom_value(S,fom); G=fom_Sbar(S,fom);
iw=1/(1i*omega); AinvT=(A.')\eye(nch); I=eye(nch);
Wfc  =  0.5    *(I-S.')*G*AinvT;
Wfdc = -0.5*iw *(I+S.')*G*AinvT;
end
function lam = adjoint_v6(P, out, Wfc, Wfdc)
QvalE=P.Proj*out.MtarE_val; QderE=P.Proj*out.MtarE_sp;
GC=P.QvalC.'*Wfc+P.QderC.'*Wfdc; GE=QvalE.'*Wfc+QderE.'*Wfdc;
rhsE=GE-out.Mmpe.'*(P.MCmp.'*(P.iMCC.'*GC));
lamE=out.PE.'*(out.LE.'\(out.UE.'\rhsE));
lamC=P.iMCC.'*(GC-P.MmpC.'*(out.Memp.'*lamE));
lam=zeros(P.npt,P.numch); lam(P.Cidx,:)=lamC; lam(P.Eidx,:)=lamE;
end
function [F, grad, S] = fom_and_grad(P, angles, fom)
out=solve_v6(P,angles); S=out.S;
[F,Wfc,Wfdc]=fom_cotangents(S,out.ain,fom,P.omega);
lam=adjoint_v6(P,out,Wfc,Wfdc);
grad=grad_multipole(P,out,lam,Wfc,Wfdc);
end
function grad = grad_multipole(P, out, lam, Wfc, Wfdc)
nch=P.numch; k=P.omega; grad=zeros(1,P.g_.numell);
sig=out.sig; sigC=out.sig(P.Cidx,:); Xe=out.Xe; Ne=out.Ne; we=out.we; lamC=lam(P.Cidx,:);
for q=1:P.g_.numell
    ml=P.mlocs{q}; rg=P.Eloc{q}; cq=P.g_.ell_cens(:,q);
    MCmp_q=P.MCmp(:,ml); MmpC_q=P.MmpC(ml,:); Mmpe_q=out.Mmpe(ml,rg); Memp_q=out.Memp(rg,ml); imv=1i*P.mvec(:);
    Eqg=P.Eglob{q}; lamEq=lam(Eqg,:); sigEq=sig(Eqg,:);
    vq=[-(Xe(2,rg)-cq(2));(Xe(1,rg)-cq(1))]; nqe=Ne(:,rg); dnqe=[-nqe(2,:);nqe(1,:)];
    acc=0;
    for j=1:nch
        a=MCmp_q.'*lamC(:,j); b=Mmpe_q*sigEq(:,j); t8=sum(a.*(-imv).*b);
        cc=Memp_q.'*lamEq(:,j); d=MmpC_q*sigC(:,j); t7=sum(cc.*(imv).*d);
        acc=acc-(t7+t8);
    end
    for qp=1:P.g_.numell
        if qp==q, continue; end
        rgp=P.Eloc{qp}; Eqpg=P.Eglob{qp};
        DM7=2*dSpmat_cross(Xe(:,rg),Ne(:,rg),Xe(:,rgp),we(rgp),k,vq,dnqe,zeros(2,numel(rgp)));
        DM8=2*dSpmat_cross(Xe(:,rgp),Ne(:,rgp),Xe(:,rg),we(rg),k,zeros(2,numel(rgp)),zeros(2,numel(rgp)),vq);
        for j=1:nch, acc=acc-( lam(Eqg,j).'*(DM7*sig(Eqpg,j)) + lam(Eqpg,j).'*(DM8*sig(Eqg,j)) ); end
    end
    for j=1:nch
        dUv=dS_eval(P.tinfo.r,Xe(:,rg),we(rg),sigEq(:,j),k,vq);
        dUs=dSp_eval(P.tinfo.r,Xe(:,rg),we(rg),sigEq(:,j),k,vq);
        acc=acc+Wfc(:,j).'*(P.Proj*dUv)+Wfdc(:,j).'*(P.Proj*dUs);
    end
    grad(q)=2*real(acc);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        cluster / pool plumbing  (verbatim from statistical_study_v2.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ensure_chunkie()
if ~isempty(which('chunkermat')), return; end
sd=fileparts(mfilename('fullpath'));
cands={fullfile(sd,'..','..','chunkie'), fullfile(sd,'..','chunkie'), ...
       fullfile(getenv('HOME'),'main_projects','Numerics','chunkie')};
for i=1:numel(cands)
    if exist(cands{i},'dir')
        addpath(genpath(cands{i}));
        if ~isempty(which('chunkermat')), return; end
    end
end
error('chunkie not found (looked in ../../chunkie, ../chunkie, ~/main_projects/Numerics/chunkie)');
end
function [myk, tag] = array_partition(R)
ts=getenv('SLURM_ARRAY_TASK_ID');
if isempty(ts), myk=1:R; tag=''; return; end
tid=str2double(ts);
ntasks=str2double(getenv('SLURM_ARRAY_TASK_COUNT'));
mn=str2double(getenv('SLURM_ARRAY_TASK_MIN')); mx=str2double(getenv('SLURM_ARRAY_TASK_MAX'));
if isnan(mn), mn=0; end
if isnan(ntasks)||ntasks<1, if ~isnan(mx), ntasks=mx-mn+1; else, ntasks=1; end, end
myk=(tid-mn+1):ntasks:R; tag=sprintf('_task%03d',tid);
end
function setup_pool()
if ~isempty(gcp('nocreate')), return; end
if ~license('test','Distrib_Computing_Toolbox'), return; end
n=str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(n)||n<2, n=str2double(getenv('SLURM_CPUS_ON_NODE')); end
if isnan(n)||n<2, n=feature('numcores'); end
try
    pc=parcluster('local'); jid=getenv('SLURM_JOB_ID');
    if ~isempty(jid), jsl=fullfile(tempdir,['matlab_jsl_' jid]); if ~exist(jsl,'dir'),mkdir(jsl);end, pc.JobStorageLocation=jsl; end
    pc.NumWorkers=n; parpool(pc,n); fprintf('parpool: %d workers\n',n);
catch ME
    fprintf('parpool setup failed (%s); parfor runs serially.\n',ME.message);
end
end
function ensure_dir(d), if ~exist(d,'dir'), try mkdir(d); catch, end, end, end
function fprintf_both(fid,varargin), fprintf(varargin{:}); if fid>0, fprintf(fid,varargin{:}); end, end
function s=tern(c,a,b), if c, s=a; else, s=b; end, end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        PRECOMPUTE (v6)  -- verbatim from statistical_study_v2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = precompute_realization(C, cfg)
omega=2*pi*cfg.f/C.c; oph=C.c*omega;
if isfield(C,'lossless') && C.lossless
    gamma=0;                                   % PEC walls: zero surface impedance (lossless)
else
    Zs=sqrt(1i*oph*C.mu0/(C.sigma+1i*oph*C.eps0)); gamma=1i*oph*C.eps0*Zs;
end
g_=struct('R',C.R,'d',C.d,'L',C.L,'W',C.W,'lead_ycen',C.lead_ycen,'ell_semi',C.ell_semi,...
          'ell_cens',cfg.positions,'numell',cfg.N_scat,'maxchunklen',(2*pi/omega)/C.mcl_factor,'legk',16);
P=precompute_v6(g_, omega, gamma, C.numch, zeros(1,cfg.N_scat), C.mp_tol);
P.C=C; P.cfg=cfg; P.lambda=2*pi/omega; P.mcl=g_.maxchunklen;
end

function P = precompute_v6(g_, omega, gamma, numch, base_ang, mp_tol)
sp=kernel('h','sp',omega); sk=kernel('h','s',omega);
cg0=build_cav(g_,base_ang); nedge=numel(cg0.echnks); npt=cg0.npt;
imp_edges=[1+4*(0:numch),2+4*numch]; ncav_edge=2+4*numch;
np=arrayfun(@(e)e.npt,cg0.echnks); estart=[0 cumsum(np(:).')]; eidx=@(e)(estart(e)+1):estart(e+1);
Cidx=1:estart(ncav_edge+1); Eidx=(estart(ncav_edge+1)+1):npt; nC=numel(Cidx); nE=numel(Eidx);
Eglob=cell(1,g_.numell); Eloc=cell(1,g_.numell); off=0;
for q=1:g_.numell, rg=eidx(ncav_edge+q); Eglob{q}=rg; Eloc{q}=off+(1:numel(rg)); off=off+numel(rg); end
impc=false(nC,1); for e=imp_edges, impc(eidx(e))=true; end
wall=wts_all(cg0); Xc=cg0.r(:,Cidx); Nc=cg0.n(:,Cidx); wc=wall(Cidx);
Kmat=make_Kmat(nedge,imp_edges,sp,sk,gamma);
Mfull0=chunkermat(cg0,Kmat)+eye(npt); MCC=Mfull0(Cidx,Cidx); iMCC=inv(MCC); %#ok<MINV>
MEiEi=cell(1,g_.numell); for q=1:g_.numell, MEiEi{q}=Mfull0(Eglob{q},Eglob{q}); end
ch=build_channels(cg0,g_,omega); Ny=size(ch(1).XY,2); Ntar=numch*Ny;
tinfo.r=[ch(:).XY]; tinfo.d=repmat([0;1],1,Ntar); tinfo.d2=zeros(2,Ntar); tinfo.n=repmat([1;0],1,Ntar);
MtarC_val=chunkerkernevalmat(cg0,sk,tinfo); MtarC_val=MtarC_val(:,Cidx);
MtarC_sp =chunkerkernevalmat(cg0,sp,tinfo); MtarC_sp =MtarC_sp (:,Cidx);
RC=zeros(nC,numch); for j=1:numch, RC(idxmap(ch(j).leadend_inds,Cidx),j)=2.0; end
Proj=zeros(numch,Ntar); for i=1:numch, Proj(i,(i-1)*Ny+(1:Ny))=ch(i).proj; end
beta=arrayfun(@(s)s.beta,ch).';
r_enc=max(g_.ell_semi); pmode=0;
for q=1:g_.numell
    dmin=min(sqrt(sum((Xc-g_.ell_cens(:,q)).^2,1))); ratio=r_enc/dmin;
    assert(ratio<1,'ellipse %d not separated from cavity (ratio=%.3f)',q,ratio);
    pmode=max(pmode,max(ceil(log(mp_tol)/log(ratio)),ceil(omega*r_enc)+5));
end
mvec=-pmode:pmode; nm=numel(mvec); Nmp=g_.numell*nm; mlocs=cell(1,g_.numell);
for q=1:g_.numell, mlocs{q}=(q-1)*nm+(1:nm); end
MCmp=zeros(nC,Nmp); MmpC=zeros(Nmp,nC);
for q=1:g_.numell, cq=g_.ell_cens(:,q);
    MCmp(:,mlocs{q})=build_Acmp(Xc,Nc,impc,cq,mvec,omega,gamma);
    MmpC(mlocs{q},:)=build_Dmpc(Xc,wc,cq,mvec,omega);
end
P=struct(); P.omega=omega; P.gamma=gamma; P.sp=sp; P.sk=sk; P.g_=g_; P.base_ang=base_ang;
P.Cidx=Cidx; P.Eidx=Eidx; P.Eglob=Eglob; P.Eloc=Eloc; P.nC=nC; P.nE=nE; P.npt=npt; P.numch=numch;
P.Xc=Xc; P.Nc=Nc; P.wc=wc; P.impc=impc; P.iMCC=iMCC; P.MEiEi=MEiEi; P.Kmat=Kmat; P.ncav_edge=ncav_edge;
P.ch=ch; P.tinfo=tinfo; P.Proj=Proj; P.beta=beta; P.RC=RC; P.Ntar=Ntar; P.Ny=Ny;
P.MtarC_val=MtarC_val; P.MtarC_sp=MtarC_sp; P.mvec=mvec; P.mlocs=mlocs; P.Nmp=Nmp; P.MCmp=MCmp; P.MmpC=MmpC;
P.iMCC_RC=iMCC*RC; P.iMCC_MCmp=iMCC*MCmp; P.QvalC=Proj*MtarC_val; P.QderC=Proj*MtarC_sp;
P.EC_val=MtarC_val*iMCC; P.EC_sp=MtarC_sp*iMCC; P.U0_val=P.EC_val*RC; P.U0_sp=P.EC_sp*RC;
P.Psi_val=P.EC_val*MCmp; P.Psi_sp=P.EC_sp*MCmp; P.r_mp=MmpC*(iMCC*RC); P.W=MmpC*(iMCC*MCmp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        SOLVE (v6)  -- verbatim from statistical_study_v2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = solve_v6(P, ang)
cg=build_cav(P.g_,ang); wall=wts_all(cg); k=P.omega;
Xe=cg.r(:,P.Eidx); Ne=cg.n(:,P.Eidx); we=wall(P.Eidx);
Mmpe=zeros(P.Nmp,P.nE); Memp=zeros(P.nE,P.Nmp);
for q=1:P.g_.numell, cq=P.g_.ell_cens(:,q); rg=P.Eloc{q}; ml=P.mlocs{q};
    Mmpe(ml,rg)=build_Bmpe(Xe(:,rg),we(rg),cq,P.mvec,k);
    Memp(rg,ml)=build_Cemp(Xe(:,rg),Ne(:,rg),cq,P.mvec,k);
end
MEE=assemble_MEE(P,Xe,Ne,we,k); SE=MEE-Memp*(P.W*Mmpe); [LE,UE,PE]=lu(SE);
SigE=UE\(LE\(PE*(-(Memp*P.r_mp)))); sigC=P.iMCC_RC-P.iMCC_MCmp*(Mmpe*SigE);
ell_chnkr=merge(cg.echnks(P.ncav_edge+1:end));
MtarE_val=chunkerkernevalmat(ell_chnkr,P.sk,P.tinfo); MtarE_sp=chunkerkernevalmat(ell_chnkr,P.sp,P.tinfo);
gmp=Mmpe*SigE; U_=P.U0_val-P.Psi_val*gmp+MtarE_val*SigE; UX=P.U0_sp-P.Psi_sp*gmp+MtarE_sp*SigE;
FC=P.Proj*U_; FDC=P.Proj*UX;
ain=(FC+(1./(1i*P.beta)).*FDC)/2; aout=(FC-(1./(1i*P.beta)).*FDC)/2; S=aout/ain;
sig=zeros(P.npt,P.numch); sig(P.Cidx,:)=sigC; sig(P.Eidx,:)=SigE;
out=struct('S',S,'ain',ain,'sig',sig,'SigE',SigE,'Mmpe',Mmpe,'Memp',Memp,'LE',LE,'UE',UE,'PE',PE,...
           'MtarE_val',MtarE_val,'MtarE_sp',MtarE_sp,'Xe',Xe,'Ne',Ne,'we',we);
end
function MEE = assemble_MEE(P,Xe,Ne,we,k)
MEE=zeros(P.nE);
for i=1:P.g_.numell, MEE(P.Eloc{i},P.Eloc{i})=P.MEiEi{i};
    for j=1:P.g_.numell, if i==j, continue; end
        MEE(P.Eloc{i},P.Eloc{j})=block_EC(Xe(:,P.Eloc{i}),Ne(:,P.Eloc{i}),Xe(:,P.Eloc{j}),we(P.Eloc{j}),k); end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        derivative kernels  -- verbatim from statistical_study_v2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=dS_eval(Xt,Xs,ws,sig,k,vs)
nt=size(Xt,2);out=zeros(nt,1);
for i=1:nt,d=Xs-Xt(:,i);r=sqrt(sum(d.^2,1));
 ker=(-1i*k/4)*besselh(1,1,k*r).*(sum(d.*vs,1)./r);out(i)=sum(ker(:).*sig(:).*ws(:));end
end
function out=dSp_eval(Xt,Xs,ws,sig,k,vs)
out=dSp_cross(Xt,repmat([1;0],1,size(Xt,2)),Xs,ws,sig,k,zeros(2,size(Xt,2)),zeros(2,size(Xt,2)),vs);
end
function out=dSp_cross(XT,nT,XS,wS,sigS,k,vT,dnT,vS)
NT=size(XT,2);out=zeros(NT,1);c1=-1i*k/4;
for i=1:NT,d=XT(:,i)-XS;r=sqrt(sum(d.^2,1));H0=besselh(0,1,k*r);H1=besselh(1,1,k*r);
 nd=nT(:,i).'*d;dd=vT(:,i)-vS;ddel=sum(d.*dd,1);ndel=(dnT(:,i).'*d)+sum(nT(:,i).*dd,1);
 term=c1*(nd.*ddel.*(k*H0./r.^2-2*H1./r.^3)+H1.*ndel./r);out(i)=sum(term(:).*sigS(:).*wS(:));end
end
function M=dSpmat_cross(XT,nT,XS,wS,k,vT,dnT,vS)
NT=size(XT,2);NS=size(XS,2);M=zeros(NT,NS);c1=-1i*k/4;
for i=1:NT,d=XT(:,i)-XS;r=sqrt(sum(d.^2,1));H0=besselh(0,1,k*r);H1=besselh(1,1,k*r);
 nd=nT(:,i).'*d;dd=vT(:,i)-vS;ddel=sum(d.*dd,1);ndel=(dnT(:,i).'*d)+sum(nT(:,i).*dd,1);
 M(i,:)=c1*(nd.*ddel.*(k*H0./r.^2-2*H1./r.^3)+H1.*ndel./r).*wS(:).';end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% multipole factor builders -- verbatim %%%
function A=build_Acmp(Xc,Nc,impc,cq,mvec,k,gamma)
dx=Xc(1,:)-cq(1);dy=Xc(2,:)-cq(2);rho=sqrt(dx.^2+dy.^2);al=atan2(dy,dx);
rhat_n=(dx.*Nc(1,:)+dy.*Nc(2,:))./rho;ahat_n=(-sin(al).*Nc(1,:)+cos(al).*Nc(2,:));
gam=(2*gamma)*double(impc(:).');A=zeros(numel(rho),numel(mvec));
for a=1:numel(mvec),m=mvec(a);
    Hm=besselh(m,1,k*rho);Hmp=(besselh(m-1,1,k*rho)-besselh(m+1,1,k*rho))/2;
    e=exp(1i*m*al);O=(1i/4)*Hm.*e;dnO=(1i/4)*e.*(k*Hmp.*rhat_n+(1i*m./rho).*Hm.*ahat_n);
    A(:,a)=(2*dnO+gam.*O).';
end
end
function B=build_Bmpe(Xq,wq,cq,mvec,k)
dx=Xq(1,:)-cq(1);dy=Xq(2,:)-cq(2);rho=sqrt(dx.^2+dy.^2);al=atan2(dy,dx);
B=zeros(numel(mvec),numel(rho));
for a=1:numel(mvec),m=mvec(a);B(a,:)=besselj(m,k*rho).*exp(-1i*m*al).*wq(:).';end
end
function C=build_Cemp(Xq,Nq,cq,mvec,k)
dx=Xq(1,:)-cq(1);dy=Xq(2,:)-cq(2);rho=sqrt(dx.^2+dy.^2);al=atan2(dy,dx);
rhat_n=(dx.*Nq(1,:)+dy.*Nq(2,:))./rho;ahat_n=(-sin(al).*Nq(1,:)+cos(al).*Nq(2,:));
C=zeros(numel(rho),numel(mvec));
for a=1:numel(mvec),m=mvec(a);
    Jm=besselj(m,k*rho);Jmp=(besselj(m-1,k*rho)-besselj(m+1,k*rho))/2;
    e=exp(1i*m*al);dnpsi=e.*(k*Jmp.*rhat_n+(1i*m./rho).*Jm.*ahat_n);C(:,a)=(2*dnpsi).';
end
end
function D=build_Dmpc(Xc,wc,cq,mvec,k)
dx=Xc(1,:)-cq(1);dy=Xc(2,:)-cq(2);rho=sqrt(dx.^2+dy.^2);al=atan2(dy,dx);
D=zeros(numel(mvec),numel(rho));
for a=1:numel(mvec),m=mvec(a);D(a,:)=(1i/4)*besselh(m,1,k*rho).*exp(-1i*m*al).*wc(:).';end
end
function M=block_EC(Xt,Nt,Xs,ws,k)
nt=size(Xt,2);M=zeros(nt,size(Xs,2));
for i=1:nt,d=Xt(:,i)-Xs;r=sqrt(sum(d.^2,1));
    M(i,:)=2*((-1i*k/4)*besselh(1,1,k*r).*((Nt(:,i).'*d)./r)).*ws(:).';end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% shared helpers -- verbatim %%%%%%
function K=make_Kmat(nedge,imp_edges,sp,sk,gamma)
ki=2*sp+(2*gamma)*sk; kp=2*sp; K(nedge,nedge)=kernel();
for i=1:nedge,kk=kp;if ismember(i,imp_edges),kk=ki;end,for j=1:nedge,K(i,j)=kk;end,end
end
function w=wts_all(cg), w=[]; for e=1:numel(cg.echnks), we=cg.echnks(e).wts; w=[w;we(:)]; end, end %#ok<AGROW>
function loc=idxmap(globlogical,Cidx), g=find(globlogical); [~,loc]=ismember(g,Cidx); loc=loc(loc>0); end
function ch=build_channels(cg,g_,omega)
bxy=cg.r(:,:); xout=g_.d+g_.L; W=g_.W; yc=g_.lead_ycen; xcut=g_.d+0.8*g_.L; Ny=21;
for j=1:numel(yc), cen=yc(j);
    ch(j).leadend_inds=(abs(bxy(1,:)-xout)<1e-5 & abs(bxy(2,:)-cen)<=W/2+1e-9).'; %#ok<AGROW>
    ylin=cen+linspace(-0.4*W,0.4*W,Ny); ch(j).XY=[xcut*ones(1,Ny);ylin]; %#ok<AGROW>
    pp=[]; pp.r=ch(j).XY; pp.d=repmat([0;1],1,Ny); pp.d2=zeros(2,Ny); pp.n=repmat([1;0],1,Ny);
    ch(j).XY_ptinfo=pp; phi=ones(Ny,1); ch(j).proj=(phi.')/(phi.'*phi); ch(j).beta=omega; %#ok<AGROW>
end
end
function cg=build_cav(g_,ell_ang)
R=g_.R;d=g_.d;L=g_.L;W=g_.W;yc=g_.lead_ycen;
a60=deg2rad(60);a300=deg2rad(300);Atop=R*[cos(a60);sin(a60)];Abot=R*[cos(a300);sin(a300)];
xin=d;xout=d+L;ys=sort(yc);verts=[xin;Abot(2)];
for m=1:numel(ys),yb=ys(m)-W/2;yt=ys(m)+W/2;verts=[verts,[xin;yb],[xout;yb],[xout;yt],[xin;yt]];end %#ok<AGROW>
verts=[verts,[xin;Atop(2)]];nv=size(verts,2);edgeinc=[1:nv;circshift(1:nv,-1)];edgefuns=cell(1,nv);
for i=1:nv-1,edgefuns{i}=@(t)(1-t(:).').*verts(:,i)+t(:).'.*verts(:,i+1);end
edgefuns{nv}=@(t)R*[cos(a60+t(:).'*(a300-a60));sin(a60+t(:).'*(a300-a60))];
for q=1:g_.numell,cq=g_.ell_cens(:,q);th=ell_ang(q);Rr=[cos(th) -sin(th);sin(th) cos(th)];
 edgeinc=[edgeinc,nan(2,1)]; %#ok<AGROW>
 edgefuns{end+1}=@(t)cq+Rr*[g_.ell_semi(1)*cos(2*pi*t(:).');-g_.ell_semi(2)*sin(2*pi*t(:).')];end %#ok<AGROW>
cp=[];cp.eps=1e-10;cp.nover=0;cp.maxchunklen=g_.maxchunklen;pr=[];pr.k=g_.legk;
cg=chunkgraph(verts,edgeinc,edgefuns,cp,pr);
end
