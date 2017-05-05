function [ psi ] = annulusPsiEqn( psi, f, dr, dtheta, r )
%ANNULUSPSIEQN: Solves delsqr(psi) = f by multigrid

tol = 1e-10;
maxres = 1;
while maxres > tol
    psi = psi_MultiGridV( psi, f, dr, dtheta, r);
    res = psi_residual( psi, f, dr, dtheta, r);
    maxres = max(max(abs(res)));
end

end

