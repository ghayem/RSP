
function C = My_2D_Cov(rho, sig,x1,x2)

for i = 1: size(x1,1)
    
    for j = 1: size(x2,1)
        
        Delta = vecnorm(x1(i,:)-x2(j,:));
        C(i,j) = sig.^2 .* exp( -1/(2*rho^2) * Delta.^2);
        
    end
    
end


end
