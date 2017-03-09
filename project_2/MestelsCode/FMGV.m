function u = FMGV(u0,f,NV)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Full Multigrid and then NV V-cycles  % 
% delsqr u = f on square               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % initialization
  u0 = FullMG(f);

  % Now do Multigrid V-cycles
  for k=1:NV
    u0=MultigridV(u0,f);
  end

  % and tell us the answer
  error=max(max(residual(u0,f)))
  u = u0;
  
