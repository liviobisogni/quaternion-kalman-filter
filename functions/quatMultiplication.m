function q_times_p = quatMultiplication(q, p)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% [https://en.wikipedia.org/wiki/Quaternion#Hamilton_product]
% QUATMULTIPLICATION Multiplies the quaternions q and p, thus obtaining their Hamilton product.
% Please note: the quaternion product is NOT commutative!
% ((( NB: [Eq. post-5 Trawny] is not used ))
%
% INPUT:
%   * q,                Quaternion q(k) (q = [q0; q1; q2; q3], q0 is the scalar)        (4 x 1) vector      []
%   * p,                Quaternion p(k) (p = [p0; p1; p2; p3], p0 is the scalar)        (4 x 1) vector      []
%
% OUTPUT:
%   * q_times_p,        (not unitary) Quaternion, product of (in this order) q and p	(4 x 1) vector      []
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    % Check number of arguments
    narginchk(2,2);
    
    if (~isequal(size(q), [4 1]) && ~isequal(size(p), [4 1]))
        error('q and p must be (4 x 1) vectors.');
    end
    if (~isequal(size(p), [4 1]))
        error('p must be a (4 x 1) vector.');
    end
    if (~isequal(size(q), [4 1]))
        error('q must be a (4 x 1) vector.');
    end

    q0p0 = q(1) * p(1);
    q0p1 = q(1) * p(2);
    q0p2 = q(1) * p(3);
    q0p3 = q(1) * p(4);
    q1p0 = q(2) * p(1);
    q1p1 = q(2) * p(2);
    q1p2 = q(2) * p(3);
    q1p3 = q(2) * p(4);
    q2p0 = q(3) * p(1);
    q2p1 = q(3) * p(2);
    q2p2 = q(3) * p(3);
    q2p3 = q(3) * p(4);
    q3p0 = q(4) * p(1);
    q3p1 = q(4) * p(2);
    q3p2 = q(4) * p(3);
    q3p3 = q(4) * p(4);

    % https://en.wikipedia.org/wiki/Quaternion#Hamilton_product
    q_times_p = [q0p0 - q1p1 - q2p2 - q3p3;
                 q0p1 + q1p0 + q2p3 - q3p2;
                 q0p2 - q1p3 + q2p0 + q3p1;
                 q0p3 + q1p2 - q2p1 + q3p0];

    
    %%%%%%%% (([Eq. post-5 Trawny]; not used, since it probably swaps q and p))
    %

end