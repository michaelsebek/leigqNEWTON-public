function S = newton_doc_matrices(which)
%NEWTON_DOC_MATRICES  Fixed paper matrices for documentation examples.
%  (The function is used by almost all documentation scripts.)
%
%   S = newton_doc_matrices('HuangSo')
%   S = newton_doc_matrices('MVPS')
%   S = newton_doc_matrices('PanNg')
%   S = newton_doc_matrices
%   S = newton_doc_matrices('all')
%
% Output S is a struct with named matrices and (when available) reference lambdas.
%
% This helper exists only to keep documentation pages concise.

if exist('quaternion','class') ~= 8
    error('Requires MATLAB quaternion class.');
end

q0 = quaternion(0,0,0,0);
q1 = quaternion(1,0,0,0);
qi = quaternion(0,1,0,0);
qj = quaternion(0,0,1,0);
qk = quaternion(0,0,0,1);


% Allow calling with no argument (for interactive exploration).
if nargin < 1 || (isstring(which) && strlength(which)==0) || (ischar(which) && isempty(which))
    which = 'all';
end

which = lower(string(which));
S = struct();

if which == "all"
    S.HuangSo = newton_doc_matrices('HuangSo');
    S.MVPS    = newton_doc_matrices('MVPS');
    S.PanNg   = newton_doc_matrices('PanNg');
    return;
end


switch which
    case "huangso"
        % Example 2.5
        S.A25 = [ q0, q1 + qi; q1 - qi, q0 ];
        S.lam25_true = [ quaternion(sqrt(2),0,0,0); quaternion(-sqrt(2),0,0,0) ];
        % Example 2.6
        S.A26 = [ q0, qi; qj, q1 ];
        S.lam26_true = [ (q1 + qi + qj - qk)/2; (q1 - qi - qj - qk)/2 ];
        % Example 2.7 (sphere)
        S.A27 = [ quaternion(2,0,0,0),  qi; -qi, quaternion(2,0,0,0) ];

    case "mvps"
        S.A19 = [  qi,    q0,   q0; qk, qj, q0; -3*qi, 2*qk, qk ];
        S.A38 = [  q0, qi, q1; 3*qi-qk, q0, q1; qk, (-q1+qj+qk), q0 ];
        S.B51 = S.A38 + qi*eye(3);
        S.A52 = [ qj, q1, q0; 2*qi, (-qk), q1; (2*q1-qi-2*qj), (-q1-qj+qk), (-qi-qk) ];
        S.A55 = [ qk, q0, q0; 3*qi-qj, (-qi), qi; (q1-2*qk), qj, (-qj) ];
        S.A56 = [ (-qi-qj), q0, q0; qk, (-qi), qi; (q1-qi), qj, (-qj) ];

    case "panng"
        q = @(w,x,y,z) quaternion(w,x,y,z);
        a = q(-2, 1, 1, 4);
        b = q( 2, 4, 1, 1);
        c = q( 1, 3, 2, 2);
        d = q(-1, 2, 2, 3);
        S.A = [a b c d; d a b c; c d a b; b c d a];

    otherwise
        error('Unknown matrix set: %s', which);
end
end