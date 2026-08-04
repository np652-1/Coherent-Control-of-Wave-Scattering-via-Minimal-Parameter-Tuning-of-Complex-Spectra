function cases = ccon_cases()
%CCON_CASES  The seven CCON figure-of-merit cases for the lossy 4-port NDDD cavity.
%
%  Each case minimizes        F(theta) = det(C'*C)      (real, >= 0)
%  where C = S[rows,cols] is a sub-matrix of the 4x4 scattering matrix S and
%  theta are the ellipse rotation angles.  This is the SAME functional form as
%  the completed NDDD study (whose FOM was |S21|^2+|S31|^2+|S41|^2 = det(C'C)
%  for the column C=(S21;S31;S41)), kept un-squared for consistency.
%
%  rows/cols are DERIVED from the CCON letter word (see ali_paper.pdf, Fig. 1):
%     R = R-tilde (reflectionless): driven, zero outgoing  -> contributes a ROW and a COL
%     N = normal               : driven, free outgoing     -> contributes a COL
%     D = dark                 : no drive, zero outgoing    -> contributes a ROW
%     T = transmission         : no drive, free outgoing    -> contributes nothing
%  i.e.  rows = ports in {R,D} (zero-output),  cols = ports in {N,R} (driven).
%
%  N_max is the sweep ceiling for the number of ellipses; it sits a margin above
%  n_p, the predicted #angles needed to drive F to zero (the floor should appear
%  at N_scat = n_p):
%     square C  (det C = 0, one complex condition)        -> n_p = 2
%     column C  (m x 1, whole column -> 0, m complex)      -> n_p = 2*m
%
%  Returns a struct array with fields: name, word, rows, cols, Nmax, np.
%  NDDD is now included as a NATIVE case in this list (word 'NDDD' ->
%  rows=[2 3 4], cols=1, C=(S21;S31;S41), the same column FOM as the original
%  results_4 leakage study, n_p=6).  It was run separately before; here it is
%  swept in the same unified job array as the other six.

% RRRR is DROPPED in the lossless study: a lossless cavity has a unitary S, so
% F=det(S'S)=|det S|^2 == 1 for all angles -- degenerate (nothing to optimize).
words = {'RTTT','NDTT','NNDD','RDDT','RDTT','NDDD'};
Nmax  = [   4 ,   4 ,   4 ,   8 ,   8 ,   8  ];

cases = struct('name',{},'word',{},'rows',{},'cols',{},'Nmax',{},'np',{});
for i = 1:numel(words)
    [r,c,np] = ccon_word_to_idx(words{i});
    cases(i) = struct('name',words{i},'word',words{i},'rows',r,'cols',c,...
                      'Nmax',Nmax(i),'np',np);
end

% self-check: the derived sub-matrices must match the user's specification
assert(isequal(cases(1).rows,1)        && isequal(cases(1).cols,1),        'RTTT != S11');
assert(isequal(cases(2).rows,2)        && isequal(cases(2).cols,1),        'NDTT != S21');
assert(isequal(cases(3).rows,[3 4])    && isequal(cases(3).cols,[1 2]),    'NNDD != S[34,12]');
assert(isequal(cases(4).rows,[1 2 3])  && isequal(cases(4).cols,1),        'RDDT != (S11;S21;S31)');
assert(isequal(cases(5).rows,[1 2])    && isequal(cases(5).cols,1),        'RDTT != (S11;S21)');
assert(isequal(cases(6).rows,[2 3 4])  && isequal(cases(6).cols,1),        'NDDD != (S21;S31;S41)');
assert(cases(6).np==6 && cases(6).Nmax==8,                                 'NDDD np/Nmax mismatch');
end


function [rows,cols,np] = ccon_word_to_idx(word)
%CCON_WORD_TO_IDX  Map a CCON letter word to the constraint sub-matrix indices.
L = upper(word); rows = []; cols = [];
for p = 1:numel(L)
    switch L(p)
        case 'R', rows(end+1) = p; cols(end+1) = p;   %#ok<AGROW> reflectionless: row & col
        case 'D', rows(end+1) = p;                    %#ok<AGROW> dark: row only
        case 'N', cols(end+1) = p;                    %#ok<AGROW> normal: col only
        case 'T'                                       % transmission: neither
        otherwise, error('ccon_word_to_idx:badletter','bad CCON letter "%s"',L(p));
    end
end
m = numel(rows); n = numel(cols);
% predicted #tunable angles needed to null F=det(C'C):
if n <= 1
    np = 2*m;            % column C: whole column must vanish (m complex conditions)
else
    np = 2;              % square/wide C: det C = 0 (one complex condition)
end
end
