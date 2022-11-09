

function [J_N] = CostFunction_CDF_Complementary(indSelected, indNonSelected, nbPts, x, ma, Ka,...
                 aSelected, maSelected, Kn, a ,gg, g_step  )



nbX = size(indSelected,1);

%% Generate a_N | vec(a)_N-1)
% indNonSelected = setdiff(1:nbPts, indSelected);

nbPtsTest = 100;

for indN = 1:nbPts
    xN(indN) = x(indN);
    
    if ismember(indN, indNonSelected)
        %indNonSelected   %(randomOrder(1));
        
        mu_aN(indN,1) = ma(indN) + Ka(indSelected, indN).' * ((Ka(indSelected,indSelected)+1e-8*eye(nbX)) \ (aSelected - maSelected));
        s2_aN(indN,1) = Ka(indN,indN) - Ka(indSelected, indN).' / (Ka(indSelected,indSelected)+1e-8*eye(nbX)) * Ka(indSelected, indN);
        
        aN(indN, 1:nbPtsTest) = mu_aN(indN) + sqrt(s2_aN(indN)) * randn(1, nbPtsTest);
        
    else
        mu_aN(indN,1) = a(indN);
        s2_aN(indN,1) = 0;
        
        aN(indN, 1:nbPtsTest) = a(indN);
        
    end
    
end


%% Compute SNR mean

for ind = 1:length(indNonSelected)
    indN = indNonSelected(ind);      %indNonSelected(randomOrder(1));
    
    Rn     = eye(nbX+1) / (Kn([indSelected ; indN], [indSelected ; indN]) + 1e-8*eye(nbX+1));
    SNRExpectationTh(ind) = aSelected.' * Rn(1:nbX,1:nbX) * aSelected ...
        + 2 * aSelected.' * Rn(1:nbX,end) * mu_aN(indN) ...
        + Rn(end,end) * (s2_aN(indN) + mu_aN(indN)^2);
end



%% Compute SNR CDF

SNR_Initial     = aSelected.' / (Kn(indSelected,indSelected)+1e-8*eye(nbX)) * aSelected;
[SNR_ExpectedMin, indexMax] = max(SNRExpectationTh);

for ind = 1:length(indNonSelected)
    indN = indNonSelected(ind);      %indNonSelected(randomOrder(1));
    
    for indTest = 1:nbPtsTest
        SNR_N(indTest) = [aSelected ; aN(indN,indTest)].' / (Kn([indSelected ; indN],[indSelected ; indN]) + 1e-8*eye(nbX+1)) * [aSelected ; aN(indN,indTest)];
    end
    
    
    [SNRHist(:,ind), SNRHistAbsc(:,ind)] = hist(SNR_N, 35);
    SNRHist(:,ind) = SNRHist(:,ind) / (sum(SNRHist(:,ind)) * (SNRHistAbsc(2,ind)-SNRHistAbsc(1,ind)));
    
    
    
    Rn     = eye(nbX+1) / (Kn([indSelected ; indN], [indSelected ; indN]) + 1e-8*eye(nbX+1));
    alpha  = aSelected.' * Rn(1:nbX,1:nbX) * aSelected;
    beta   = aSelected.' * Rn(1:nbX,end);
    gamma  = Rn(end,end);
    lambda = (mu_aN(indN) + beta/gamma)^2 / s2_aN(indN);
    
    
    xSNRTh(:,ind) = linspace(min(SNRHistAbsc(:,ind)), max(SNRHistAbsc(:,ind)), 1e3).';
    SNR_pdfTh(:,ind) = 1/(gamma*s2_aN(indN)) * ncx2pdf(1/(gamma*s2_aN(indN))*(xSNRTh(:,ind) - (alpha - beta^2/gamma)), 1, lambda);
    
    
    probaImprovement(ind,1) = 1-ncx2cdf(1/(gamma*s2_aN(indN))*(SNR_ExpectedMin - (alpha - beta^2/gamma)), 1, lambda);
    probaDegradate(ind,1) = ncx2cdf(1/(gamma*s2_aN(indN))*(SNR_Initial - (alpha - beta^2/gamma)), 1, lambda);
    
    %     for g = 1:4
    %         threshold(g) = .25 * g * (SNR_ExpectedMin-SNR_Initial);
    %         probaDegradate_BIS(ind,g) = ncx2cdf(1/(gamma*s2_aN(indN))*(SNR_Initial + threshold(g) - (alpha - beta^2/gamma)), 1, lambda);
    %     end
    
    
    g = gg;
    threshold = g_step * g * (SNR_ExpectedMin-SNR_Initial);
    probaDegradate_BIS(ind) = ncx2cdf(1/(gamma*s2_aN(indN))*(SNR_Initial + threshold - (alpha - beta^2/gamma)), 1, lambda);

        
end

J_N = 1 - probaDegradate_BIS;

end








