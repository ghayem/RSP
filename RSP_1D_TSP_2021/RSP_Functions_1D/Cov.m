
function C = Cov(roh, sig,x1,x2)

[DeltaX,DeltaX_Prime] = meshgrid(x2,x1);

C = sig.^2 .* exp( -1/(2*roh^2) * abs(DeltaX-DeltaX_Prime).^2);

% K = size(C);
% C = C + ((sig^2)*1e-6)*eye(K);
% C = C + (1e-6)*eye(K);

end
