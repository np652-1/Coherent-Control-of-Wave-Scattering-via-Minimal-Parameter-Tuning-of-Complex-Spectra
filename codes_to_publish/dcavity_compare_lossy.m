%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  dcavity_compare_lossy.m
%
%  Cross-check that the boundary-integral (chunkie) solver and COMSOL give the
%  SAME physics for the LOSSY (impedance-wall) D-shaped cavity, at the level of
%  BOTH the scattering matrix AND the interior field.
%
%  The two solvers are implemented here as self-contained local subroutines:
%     [S, U] = chunkie_solve(P, XY)   % dense BIE  (see dcavity_chunkie_Smatrix_lossy.m)
%     [S, U] = comsol_solve (P, XY)   % COMSOL FEM (see dcavity_comsol_Smatrix_lossy.m)
%  Each returns the 4x4 S-matrix S and the total field U at a set of interior
%  evaluation points XY, one column per physical incoming channel m=1..4.
%
%  WHAT IS COMPARED
%    * |S|  and per-channel absorption 1 - sum|S(:,m)|^2.  We compare the
%      gauge-invariant magnitude |S| (and absorption): the two solvers use a
%      different per-port phase reference and the opposite time convention
%      (e^{-iwt} vs e^{+jwt}), under which |S| is invariant.  Expected ~1e-3
%      (COMSOL mesh-limited).
%    * interior field per channel.  "Unit incoming at port m" is normalized
%      differently by the two codes, so the interior fields agree UP TO ONE
%      COMPLEX SCALAR per channel; we fit it by complex least squares and
%      compare.  Expected ~1e-3.
%
%  HOW TO RUN
%    Start a COMSOL server:  comsol mphserver -port 2036   (set COMSOL_MLI below),
%    put chunkie on the path (auto), then run this file.
%
%  OUTPUT (into dcavity_compare_lossy_outputs/):
%    S_comparison.png, field_compare_channel_#.png, compare_data.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

thisdir = fileparts(mfilename('fullpath')); if isempty(thisdir), thisdir = pwd; end
outdir  = fullfile(thisdir,'dcavity_compare_lossy_outputs');
if ~exist(outdir,'dir'), mkdir(outdir); end

%% ------------------------------------------------------------------ parameters
P = struct();
P.c_light = 299792458; P.eps0 = 8.8541878128e-12; P.mu0 = 1.25663706212e-6;
P.omega   = 0.97*(2*pi);
P.lambda  = 2*pi/P.omega;
P.f_hz    = P.c_light/P.lambda;
P.R = 5.75; P.d = P.R/2; P.L = 4.0; P.W = 0.5;
P.lead_ycen = [3.7, -2.9, 1.7, -1.6];
P.ell_semi  = [0.75, 0.5];  P.ell_rot_deg = 110;
P.ell_cens  = [ 1.2, 0.0, -1.3 ; 2.4, 4.0, 1.5 ];
P.numell    = size(P.ell_cens,2);
P.maxchunklen = P.lambda/4;  P.legk = 16;
P.minppw    = 20;
P.sigma_cond = 50;                 % LOSSY: impedance-wall conductivity [S/m]
P.lossless   = false;
P.numch = numel(P.lead_ycen);
% chunkie surface admittance for the impedance walls:
oph = P.c_light*P.omega;
P.Zs    = sqrt(1i*oph*P.mu0/(P.sigma_cond+1i*oph*P.eps0));
P.gamma = 1i*oph*P.eps0*P.Zs;
% COMSOL LiveLink:
P.COMSOL_MLI = '/usr/local/comsol61/multiphysics/mli';   % <-- your COMSOL "mli" dir
P.COMSOL_PORT = 2036;
% interior evaluation grid:
P.bdr_delta = 0.30;   % keep eval points this far from walls/ellipses (chunkie near-field)
P.ngrid     = 180;

fprintf('=== D-cavity chunkie-vs-COMSOL compare (LOSSY sigma=%g, %d ports) ===\n', P.sigma_cond, P.numch);

%% --------------------------------------------------- interior evaluation grid
[XY, valid, X, Y, xg, yg] = build_eval_grid(P);
fprintf('interior grid: %d valid points (of %d)\n', numel(valid), numel(X));

%% --------------------------------------------------- run the two solvers
fprintf('\n-- chunkie BIE --\n');   [Sb, Ub] = chunkie_solve(P, XY);
fprintf('\n-- COMSOL FEM --\n');    [Sc, Uc] = comsol_solve (P, XY);

%% --------------------------------------------------- compare S-matrices
fprintf('\n================ |S| COMPARISON ================\n');
dA = abs(abs(Sb)-abs(Sc));
fprintf('|S| (chunkie):\n'); disp(abs(Sb));
fprintf('|S| (COMSOL):\n');  disp(abs(Sc));
fprintf('max | |S_chunkie|-|S_comsol| |   = %.3e\n', max(dA(:)));
fprintf('rel. Frobenius |S| difference    = %.3e\n', norm(dA,'fro')/norm(abs(Sc),'fro'));
aB = 1-sum(abs(Sb).^2,1); aC = 1-sum(abs(Sc).^2,1);
fprintf('absorption  chunkie = [%s]\n', num2str(aB,'%.4f '));
fprintf('absorption  COMSOL  = [%s]\n', num2str(aC,'%.4f '));
fprintf('absorption max|diff| = %.3e\n', max(abs(aB-aC)));
plot_S_comparison(Sb, Sc, outdir, 'LOSSY (sigma=50) impedance cavity');

%% --------------------------------------------------- compare interior fields
fprintf('\n================ FIELD COMPARISON (per channel) ================\n');
[gamma_fit, relerr] = compare_fields(P, Ub, Uc, valid, X, Y, xg, yg, outdir);
fprintf('max field relerr over channels = %.3e\n', max(relerr));

%% --------------------------------------------------- save
save(fullfile(outdir,'compare_data.mat'),'Sb','Sc','Ub','Uc','gamma_fit','relerr', ...
     'X','Y','valid','P');
fprintf('\nDONE.  Figures + data in %s\n', outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     SUBROUTINE:  chunkie BIE solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S, U] = chunkie_solve(P, XY)
% Dense single-layer BIE.  Returns the 4x4 S-matrix and the interior total field
% U (numel(XY) x numch), one column per physical incoming channel.
persistent started
if isempty(started)
    croot = ''; sd_ = fileparts(mfilename('fullpath'));
    cds = {fullfile(sd_,'chunkie'), fullfile(sd_,'..','chunkie'), ...
           fullfile(sd_,'..','..','chunkie'), fullfile(getenv('HOME'),'main_projects','Numerics','chunkie')};
    for ic=1:numel(cds), if exist(fullfile(cds{ic},'startup.m'),'file'), croot=cds{ic}; break; end, end
    assert(~isempty(croot),'chunkie not found -- run: git submodule update --init --recursive (see README).');
    od=pwd; cd(croot); startup; cd(od); started=true;
end
omega = P.omega; numch = P.numch;
sp = kernel('h','sp',omega);  sk = kernel('h','s',omega);
cg = build_cav(P);  nedge = numel(cg.echnks);  npt = cg.npt;
imp_edges = [1+4*(0:numch), 2+4*numch];
gamma = P.gamma; if P.lossless, gamma = 0; end
Kmat = make_Kmat(nedge, imp_edges, sp, sk, gamma);
M = chunkermat(cg, Kmat) + eye(npt);
ch = build_channels(cg, P, omega);
[Lf,Uf,Pf] = lu(M);
sigmas = cell(1,numch);
for j = 1:numch
    rhs = zeros(npt,1); rhs(ch(j).leadend_inds) = 2.0;
    sigmas{j} = Uf\(Lf\(Pf*rhs));
end
% S-matrix from projected in/out amplitudes on the lead cut-lines:
beta = arrayfun(@(s)s.beta,ch).';
fc = zeros(numch); fdc = zeros(numch);
for j = 1:numch
    for i = 1:numch
        u  = chunkerkerneval(cg, sk, sigmas{j}, ch(i).XY).';
        ux = chunkerkerneval(cg, sp, sigmas{j}, ch(i).XY_ptinfo).';
        fc(i,j)  = ch(i).proj*u.';
        fdc(i,j) = ch(i).proj*ux.';
    end
end
a_in  = (fc + (1./(1i*beta)).*fdc)/2;
a_out = (fc - (1./(1i*beta)).*fdc)/2;
S = a_out/a_in;
% interior field for each PURELY-INCOMING channel m (unit a_in at port m):
a_in_inv = inv(a_in);
U = nan(size(XY,2), numch);
for m = 1:numch
    special = zeros(npt,1);
    for j = 1:numch, special = special + a_in_inv(j,m)*sigmas{j}; end
    U(:,m) = chunkerkerneval(cg, sk, special, XY).';
end
fprintf('chunkie: S extracted, absorption=[%s]\n', num2str(1-sum(abs(S).^2,1),'%.4f '));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     SUBROUTINE:  COMSOL FEM solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S, U] = comsol_solve(P, XY)
% COMSOL ewfd (TM) FEM.  Returns the 4x4 S-matrix and the interior field U at XY,
% one column per physical incoming channel.  Both come out of the SAME 4 solves.
addpath(P.COMSOL_MLI);
import com.comsol.model.*
import com.comsol.model.util.*
try, ModelUtil.getModelNames(); linked=true; catch, linked=false; end
if ~linked, mphstart(P.COMSOL_PORT); end
model = build_model_comsol(P);
nP = P.numch; order = 1:nP;
S = zeros(nP,nP);  U = nan(size(XY,2), nP);
for m = 1:nP
    others = order(order~=m);  sel = [m, others];
    for k = 1:nP
        model.component('comp1').physics('ewfd').feature(sprintf('port%d',k)).selection.named(sprintf('leadbox%d',sel(k)));
    end
    model.sol('sol1').runAll;
    model.result.numerical('gev1').setResult;
    tbl = mphtable(model,'tbl1'); row = tbl.data(end,:);   % [freq, S11, S21, ..., S(nP)1]
    for k = 1:nP, S(sel(k), m) = row(1+k); end
    Hz = mphinterp(model,'ewfd.Hz','coord',XY,'complexout','on');
    U(:,m) = Hz(:);
    fprintf('  COMSOL excitation lead %d done\n', m);
end
fprintf('COMSOL: S extracted, absorption=[%s]\n', num2str(1-sum(abs(S).^2,1),'%.4f '));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     interior grid + field/S comparison helpers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [XY, valid, X, Y, xg, yg] = build_eval_grid(P)
xin=P.d; R=P.R; bd=P.bdr_delta;
xg = linspace(-R+bd, xin-bd, P.ngrid);
yg = linspace(-R*sqrt(3)/2+bd, R*sqrt(3)/2-bd, P.ngrid);
[X,Y] = meshgrid(xg,yg);
inside = (X.^2+Y.^2) <= (R-bd)^2 & X <= xin-bd;         % D-shape body, off the walls
th = deg2rad(P.ell_rot_deg); Rr=[cos(th) -sin(th); sin(th) cos(th)];
a=P.ell_semi(1); b=P.ell_semi(2);
for q=1:P.numell                                         % carve out the ellipses (+ margin)
    dx=X-P.ell_cens(1,q); dy=Y-P.ell_cens(2,q);
    xp = Rr(1,1)*dx+Rr(2,1)*dy;  yp = Rr(1,2)*dx+Rr(2,2)*dy;
    marg = bd/min(a,b);
    inside = inside & ~((xp/a).^2+(yp/b).^2 <= (1+marg)^2);
end
valid = find(inside);
XY = [X(valid).'; Y(valid).'];
end

%-----------------------------------------------------------------------------
function plot_S_comparison(Sb, Sc, outdir, ttl)
dA = abs(abs(Sb)-abs(Sc)); n = size(Sb,1);
f = figure('units','normalized','position',[0.1 0.3 0.85 0.4],'Color','w','Visible','off');
tl = {'|S| chunkie','|S| COMSOL','| |S_{chunkie}|-|S_{COMSOL}| |'};
dat = {abs(Sb), abs(Sc), dA};
for p=1:3
    subplot(1,3,p); imagesc(dat{p}); colorbar; axis square;
    set(gca,'xtick',1:n,'ytick',1:n); title(tl{p});
    for ii=1:n, for jj=1:n, text(jj,ii,sprintf('%.3f',dat{p}(ii,jj)), ...
        'horiz','center','color','w','fontweight','bold','fontsize',8); end; end
end
sgtitle(sprintf('%s: chunkie BIE vs COMSOL  |S|', ttl));
exportgraphics(f, fullfile(outdir,'S_comparison.png'),'Resolution',150); close(f);
end

%-----------------------------------------------------------------------------
function [gamma_fit, relerr] = compare_fields(P, Ub, Uc, valid, X, Y, xg, yg, outdir)
addpath(fileparts(mfilename('fullpath')));   % local redwhiteblue.m (vendored beside these scripts)
if isempty(which('redwhiteblue'))            % dev fallback
    addpath(fullfile(getenv('HOME'),'main_projects','Numerics','utilities'));
end
R=P.R; xin=P.d; th=deg2rad(P.ell_rot_deg);
cavth=linspace(deg2rad(60),deg2rad(300),400);
cavx=[R*cos(cavth), xin, R*cos(deg2rad(60))];
cavy=[R*sin(cavth), -R*sqrt(3)/2, R*sqrt(3)/2];
tt=linspace(0,2*pi,200);
numch=P.numch; gamma_fit=zeros(1,numch); relerr=zeros(1,numch);
for m=1:numch
    A = Ub(:,m); B = Uc(:,m);
    g = (B'*A)/(B'*B);           % complex LSQ: A ~ g*B (one gauge scalar per channel)
    gamma_fit(m)=g; Bs=g*B; relerr(m)=norm(A-Bs)/norm(A);
    Ua=nan(size(X)); Ua(valid)=A;
    Ub2=nan(size(X)); Ub2(valid)=Bs;
    D=nan(size(X)); D(valid)=log10(abs(A-Bs)+realmin);
    amp=max(abs([real(A);real(Bs)]));
    fig=figure('Position',[80 80 1500 460],'Color','w','Visible','off');
    titles={sprintf('Re[H_z] chunkie (chan %d)',m), ...
            sprintf('Re[\\gamma H_z] COMSOL (chan %d)',m), ...
            sprintf('log_{10}|chunkie-COMSOL|  relerr=%.2e',relerr(m))};
    dat={real(Ua),real(Ub2),D};
    for p=1:3
        subplot(1,3,p); hold on; pcolor(xg,yg,dat{p}); shading flat;
        if p<3, colormap(gca,redwhiteblue(-amp,amp,1e2)); caxis([-amp amp]);
        else,   colormap(gca,parula); caxis([-8 max(-3,min(0,max(D(:))))]); end
        colorbar; plot(cavx,cavy,'k-','LineWidth',1);
        for q=1:P.numell
            ptx=P.ell_semi(1)*cos(tt); pty=P.ell_semi(2)*sin(tt);
            ex=P.ell_cens(1,q)+cos(th)*ptx-sin(th)*pty;
            ey=P.ell_cens(2,q)+sin(th)*ptx+cos(th)*pty;
            plot(ex,ey,'k-','LineWidth',1);
        end
        hold off; axis equal tight; xlim([-R xin]); ylim([min(yg) max(yg)]);
        title(titles{p},'FontWeight','bold'); set(gca,'Layer','top');
    end
    sgtitle(sprintf('D-cavity interior field: chunkie vs COMSOL, incoming channel %d',m),'FontWeight','bold');
    print(fig, fullfile(outdir,sprintf('field_compare_channel_%d.png',m)), '-dpng','-r120'); close(fig);
    fprintf('  channel %d: |gauge|=%.3e  relerr=%.3e\n', m, abs(g), relerr(m));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 chunkie geometry / kernel / channel helpers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = make_Kmat(nedge, imp_edges, sp, sk, gamma)
kern_imp = 2*sp + (2*gamma)*sk;  kern_pec = 2*sp;
K(nedge,nedge)=kernel();
for i=1:nedge
    ki=kern_pec; if ismember(i,imp_edges), ki=kern_imp; end
    for j=1:nedge, K(i,j)=ki; end
end
end

function ch = build_channels(cg, P, omega)
bxy=cg.r(:,:); xout=P.d+P.L; W=P.W; yc=P.lead_ycen; xcut=P.d+0.8*P.L; Ny=21;
for j=1:numel(yc)
    cen=yc(j);
    ch(j).leadend_inds=( abs(bxy(1,:)-xout)<1e-5 & abs(bxy(2,:)-cen)<=W/2+1e-9 ).';
    ylin=cen+linspace(-0.4*W,0.4*W,Ny);
    ch(j).XY=[xcut*ones(1,Ny); ylin];
    pp=[]; pp.r=ch(j).XY; pp.d=repmat([0;1],1,Ny); pp.d2=zeros(2,Ny); pp.n=repmat([1;0],1,Ny);
    ch(j).XY_ptinfo=pp;
    phi=ones(Ny,1); ch(j).proj=(phi.')/(phi.'*phi); ch(j).beta=omega;
end
end

function cg = build_cav(P)
R=P.R; d=P.d; L=P.L; W=P.W; yc=P.lead_ycen;
a60=deg2rad(60); a300=deg2rad(300);
Atop=R*[cos(a60);sin(a60)]; Abot=R*[cos(a300);sin(a300)];
xin=d; xout=d+L; ys=sort(yc);
verts=[xin;Abot(2)];
for m=1:numel(ys), yb=ys(m)-W/2; yt=ys(m)+W/2;
    verts=[verts,[xin;yb],[xout;yb],[xout;yt],[xin;yt]]; end %#ok<AGROW>
verts=[verts,[xin;Atop(2)]]; nv=size(verts,2);
edgeinc=[1:nv;circshift(1:nv,-1)]; edgefuns=cell(1,nv);
for i=1:nv-1, edgefuns{i}=@(t)(1-t(:).').*verts(:,i)+t(:).'.*verts(:,i+1); end
edgefuns{nv}=@(t)R*[cos(a60+t(:).'*(a300-a60));sin(a60+t(:).'*(a300-a60))];
Rr=[cosd(P.ell_rot_deg),-sind(P.ell_rot_deg); sind(P.ell_rot_deg),cosd(P.ell_rot_deg)];
for q=1:P.numell
    cq=P.ell_cens(:,q); edgeinc=[edgeinc,nan(2,1)]; %#ok<AGROW>
    edgefuns{end+1}=@(t) cq + Rr*[ P.ell_semi(1)*cos(2*pi*t(:).'); ...
                                  -P.ell_semi(2)*sin(2*pi*t(:).') ]; %#ok<AGROW>
end
cp=[]; cp.eps=1e-10; cp.nover=0; cp.maxchunklen=P.maxchunklen; pr=[]; pr.k=P.legk;
cg=chunkgraph(verts,edgeinc,edgefuns,cp,pr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 COMSOL model builder (impedance walls if lossy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = build_model_comsol(P)
import com.comsol.model.*
import com.comsol.model.util.*
model = ModelUtil.create('Model');
model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 2);
nL = numel(P.lead_ycen);
model.param.set('R', num2str(P.R)); model.param.set('d','R/2');
model.param.set('Lw', num2str(P.L)); model.param.set('Ww', num2str(P.W));
g = model.component('comp1').geom('geom1');
g.create('ca1','CircularArc');
g.feature('ca1').set('r','R'); g.feature('ca1').set('angle1',60); g.feature('ca1').set('angle2',300);
g.create('ls1','LineSegment');
g.feature('ls1').set('specify1','coord'); g.feature('ls1').set('coord1',{'R/2' 'sqrt(3)*R/2'});
g.feature('ls1').set('specify2','coord'); g.feature('ls1').set('coord2',{'R/2' '-sqrt(3)*R/2'});
g.create('csol1','ConvertToSolid'); g.feature('csol1').selection('input').set({'ca1' 'ls1'});
ycen = P.lead_ycen; leadnames = cell(1,nL);
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
for k=1:nL
    sn = sprintf('leadbox%d',k);
    model.component('comp1').selection.create(sn,'Box');
    model.component('comp1').selection(sn).set('entitydim',1);
    model.component('comp1').selection(sn).set('xmin', xout-0.05); model.component('comp1').selection(sn).set('xmax', xout+0.05);
    model.component('comp1').selection(sn).set('ymin', ycen(k)-0.30); model.component('comp1').selection(sn).set('ymax', ycen(k)+0.30);
    model.component('comp1').selection(sn).set('condition','allvertices');
end
% impedance-wall selection (arc + chord flats) -- only used if lossy
if ~P.lossless
    model.component('comp1').selection.create('imp_arc','Box');
    model.component('comp1').selection('imp_arc').set('entitydim',1);
    model.component('comp1').selection('imp_arc').set('xmin',-P.R-0.2); model.component('comp1').selection('imp_arc').set('xmax',-4.5);
    model.component('comp1').selection('imp_arc').set('ymin',-3.5); model.component('comp1').selection('imp_arc').set('ymax',3.5);
    model.component('comp1').selection('imp_arc').set('condition','intersects');
    ys = sort(ycen); fy_lo = [-arcY, ys + P.W/2]; fy_hi = [ys - P.W/2, arcY];
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
end
model.component('comp1').material.create('mat1','Common');
model.component('comp1').material('mat1').propertyGroup('def').set('relpermittivity',{'1'});
model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability',{'1'});
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity',{'0'});
hmax = P.lambda/P.minppw;
mesh1 = model.component('comp1').mesh.create('mesh1');
ftri1 = mesh1.create('ftri1','FreeTri'); ftri1.create('size1','Size');
ftri1.feature('size1').set('custom','on');
ftri1.feature('size1').set('hmax', num2str(hmax)); ftri1.feature('size1').set('hmaxactive', true);
mesh1.run;
phys = model.component('comp1').physics.create('ewfd','ElectromagneticWavesFrequencyDomain','geom1');
phys.prop('components').set('components','inplane');
for k=1:nL
    pn=sprintf('port%d',k); phys.create(pn,'Port',1);
    phys.feature(pn).selection.named(sprintf('leadbox%d',k));
    phys.feature(pn).set('PortType','Rectangular'); phys.feature(pn).set('PortModeType','TEM');
end
if ~P.lossless
    phys.create('imp1','Impedance',1);
    phys.feature('imp1').selection.named('impsel');
    phys.feature('imp1').set('DisplacementFieldModel','RelativePermittivity');
    phys.feature('imp1').set('epsilonr_mat','userdef'); phys.feature('imp1').set('epsilonr',1);
    phys.feature('imp1').set('murbnd_mat','userdef');   phys.feature('imp1').set('murbnd',1);
    phys.feature('imp1').set('sigmabnd_mat','userdef'); phys.feature('imp1').set('sigmabnd', P.sigma_cond);
end
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
exprs = arrayfun(@(k) sprintf('ewfd.S%d1',k), 1:nL, 'uni',0);
model.result.numerical.create('gev1','EvalGlobal'); model.result.numerical('gev1').set('data','dset1');
model.result.numerical('gev1').set('expr', exprs);
model.result.table.create('tbl1','Table'); model.result.numerical('gev1').set('table','tbl1');
end
