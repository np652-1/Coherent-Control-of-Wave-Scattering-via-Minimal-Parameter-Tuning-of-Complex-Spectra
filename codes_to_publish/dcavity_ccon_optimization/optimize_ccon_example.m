%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  optimize_ccon_example.m
%
%  Restricted CCON (coherent-control) optimization demo for 0% / 50% / 90% wall
%  loss, sized to run on a common laptop in ~10 minutes (no cluster, no Parallel
%  Computing Toolbox required).
%
%  For each of the three wall-loss levels it takes one fixed D-shaped-cavity
%  realization (a random arrangement of N=2 ellipse scatterers + a random
%  frequency + a random assignment of the 4 leads to the CCON letter roles) and
%  optimizes the ellipse rotation angles to minimize a coherent-control figure
%  of merit
%        F(theta) = det(C'C),    C = S[rows,cols]  (a sub-block of the 4x4 S),
%  using the fast v6 boundary-integral forward solver together with its analytic
%  adjoint gradient (both reused verbatim from the production study).  It prints,
%  and bar-plots, the FOM before vs after optimization at each loss level.
%
%  This is a thin wrapper around run_ccon_study('example'); the full study engine
%  -- every CCON word, the production job array, and the adjoint-vs-finite-
%  difference gradient gate -- lives in run_ccon_study.m alongside this file.
%
%  Requires: chunkie (auto-found by run_ccon_study), MATLAB Optimization Toolbox
%  (fmincon).  Outputs -> results_example/.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
addpath(here);                       % so run_ccon_study.m + ccon_cases.m are found

run_ccon_study('example');           % default: word NNDD, N=2 ellipses, 3 loss levels

% --- variations you can try -------------------------------------------------
%   run_ccon_study('example','NDDD',3)   % "dark-state" word, 3 ellipses
%   run_ccon_study('example','RTTT',2)   % single-port reflectionless
%
% --- full study utilities (optional) ----------------------------------------
%   run_ccon_study('verify')   % validate the adjoint gradient vs finite differences
%   run_ccon_study('smoke')    % one quick optimization for every CCON word
%   run_ccon_study('run')      % full production job-array realization (needs a manifest row)
