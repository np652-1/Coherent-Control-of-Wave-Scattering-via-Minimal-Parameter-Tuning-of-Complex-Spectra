function man = job_manifest(R_per, base_seed)
%JOB_MANIFEST  Deterministic global list of all realizations across the 7 cases.
%
%  One unified SLURM job array runs every row of this manifest; the array index
%  selects exactly one row (one realization).  The manifest is CASE-MAJOR, and
%  within a case N_scat-major (N=1..Nmax, R_per reps each), so that
%      total rows = sum_cases(Nmax)*R_per = (4+4+4+4+8+8+8)*R_per = 40*R_per
%  With the default R_per=10 that is 400 realizations (NDDD is the 7th case).
%
%  Each row carries everything needed to (a) reproduce the realization (seed),
%  (b) route its output to the right per-case folder (name), and (c) disentangle
%  the results afterwards.  Fields:
%     gi    global index (1..R)
%     ci    case index into ccon_cases()
%     name  case name (e.g. 'RTTT')
%     N     number of ellipses for this realization
%     k     index within the case (= (N-1)*R_per + rep), N-major like results_4
%     seed  per-realization RNG seed (independent across cases)
%     Nmax  the case's sweep ceiling

if nargin<1 || isempty(R_per),    R_per    = 10;   end
if nargin<2 || isempty(base_seed), base_seed = 2026; end

cases = ccon_cases();
man = struct('gi',{},'ci',{},'name',{},'N',{},'k',{},'seed',{},'Nmax',{});
gi = 0;
for ci = 1:numel(cases)
    Nmax = cases(ci).Nmax;
    for N = 1:Nmax
        for rep = 1:R_per
            gi = gi + 1;
            k  = (N-1)*R_per + rep;
            man(gi) = struct('gi',gi,'ci',ci,'name',cases(ci).name,'N',N,'k',k,...
                             'seed',base_seed + 1000*ci + k,'Nmax',Nmax);
        end
    end
end
end
