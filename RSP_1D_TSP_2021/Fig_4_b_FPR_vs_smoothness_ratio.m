clear all
close all
clc

tic
%%
% l_ratio_set = [.01:.1:0.91,1:1:10];
nbla = 15;
l_ratio_set = sort([logspace(-2, 1, nbla).';1]);

la = 0.1;

ind_rho = size(l_ratio_set,1);
nbMcMc = 2;
sigmaS = 1;  % power of s

for i = 1:ind_rho
    
    for Mont = 1:nbMcMc
        Mont
        i
        %% Parameters
        nbPts = 3e2;
        x = linspace(0, 1, nbPts).';
        
        % Parameters of n
        ratio_rho = l_ratio_set(i);

        sa2 = 1;
        
        % Parameters of the mean function
        lma  = la;
        sma2 = sa2;
        
        %%% SNR %%%
        SNR_dB = 1;
        SNR = (10^SNR_dB) / 10;
        
        % Parameters of noise
        ln  = ratio_rho*la;
        sn2 = sa2 * sigmaS / sqrt(SNR);
        
        % Covariance function
        kSquare = @(x,s2,l)(s2 * exp(-squareform(pdist(x,'squaredeuclidean')) / (2*l^2)));
        k = @(x,xP,s2,l) (s2 * exp(-(x-xP).^2 / (2*l^2)));
        
        
        
        %% Selection X_(N-1) = [x_1, ... , x_{N-1}]
        nbX = 3;
        
        indSelected = fix([1:nbX]/(nbX+1)*nbPts);
        %indSelected = randperm(nbPts);
        indSelected = sort(indSelected(1:nbX));
        xSelected   = x(indSelected);
        
        
        %% Generate noise
        Kn = kSquare(x, sn2, ln);
        
        %% Generate ma
        Kma = kSquare(x, sma2, lma);
        Lma = chol(Kma + 1e-8*eye(nbPts), 'lower');
        ma  = Lma * randn(nbPts, 1);
        
        maSelected = ma(indSelected);
        
        
        %% Generate a
        Ka = kSquare(x, sa2, la);
        La = chol(Ka + 1e-8*eye(nbPts), 'lower');
        a  = ma + La * randn(nbPts, 1);
        
        % aSelected = ma(indSelected);
        aSelected = a(indSelected);
        
        
        %% Generate a_N | vec(a)_N-1)
        indNonSelected = setdiff(1:nbPts, indSelected);
        
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
        
        %% True SNR
        
        aSelected_True = a(indSelected);
        for ind = 1:length(indNonSelected)
            indN = indNonSelected(ind);      %indNonSelected(randomOrder(1));
            a_total = [aSelected_True;a(indN)];
            x_total = [x(indSelected) ; x(indN)];
            Rn1 = pinv( kSquare(x_total, sn2, ln) );
            SNR_true(ind,1) = a_total' * Rn1 * a_total;
        end
        
        
        %% True SNR(f_hat)
        
        zK = ma(indSelected);
        tol = 1e-10;
        flagPrior = 1;
        xK   = x(indSelected);
        ma_K = ma(indSelected);
        
        
        % 1. Estimation of a(x) from z(x)
        if flagPrior == 1;
            aEst = ma;
            
            [xx_, xy_] = meshgrid(x, x);
            Sigma_aEst = k(xx_, xy_, la, sa2);
            
        else
            [xx_kk, xy_kk] = meshgrid(xK, xK);
            SigmaA_kk = k(xx_kk, xy_kk, la, sa2);
            
            [xx_k, xy_k] = meshgrid(xK, x);
            SigmaA_k  = k(xx_k, xy_k, la, sa2);
            
            aEst = ma + SigmaA_k * (SigmaA_kk \ (zK - ma_K));
            
            [xx_, xy_] = meshgrid(x, x);
            SigmaA_k  = k(xx_, xy_, la, sa2);
            
            Sigma_aEst = SigmaA_k - SigmaA_k * (SigmaA_kk \ SigmaA_k.');
        end
        
        
        % Computation of the true SNR with estimated a(x)
        
        SNR_True_Estimated = zeros(nbPts);
        
        for indSensor1 = 1:nbPts
            for indSensor2 = 1:indSensor1
                indexExtendedSensors = unique([indSelected' ; indSensor1 ; indSensor2]);
                
                xExtended     = x(indexExtendedSensors);
                aExtended     = a(indexExtendedSensors);
                aEst_Extended = ma(indexExtendedSensors);
                
                xDist                   = squareform(pdist(xExtended, 'squaredeuclidean'));
                Rn_Extended_computation = sn2^2 * (exp( -xDist / (2*ln^2) ) + tol*eye(length(xExtended)));
                SigmaA_Extended         = Sigma_aEst(indexExtendedSensors, indexExtendedSensors);
                
                
                aTmp = Rn_Extended_computation \ aExtended;
                SNR_True_Estimated(indSensor1, indSensor2) = ...
                    sigmaS^2 * ( (aTmp.' * aEst_Extended)^2 + aTmp.' * SigmaA_Extended * aTmp ) / ...
                    ( (aEst_Extended.' / Rn_Extended_computation) * aEst_Extended + ...
                    trace(Rn_Extended_computation \ SigmaA_Extended) );
            end
        end
        
        tmp = tril(SNR_True_Estimated, -1);   % Complete the upper part of SNR_True (should be symmetric matrix)
        SNR_True_Estimated = SNR_True_Estimated + tmp.';
        SNR_True_Estimated_diag = diag(SNR_True_Estimated);
        
        
        toc
        
        SNR_True_Estimated_int = SNR_True_Estimated_diag(indSelected(1));
        
        SNR_True_Estimated_diag_norm_int = (SNR_True_Estimated_diag - SNR_True_Estimated_int) / (max(SNR_True_Estimated_diag) - SNR_True_Estimated_int);
        
        % plot(SNR_True_Estimated_diag_norm_int, 'LineWidth', 2.5); hold on
        
        
        %% Compute SNR mean
        
        % SNRExpectationTh = zeros(nbPts,1);
        
        Rn_KK     = eye(nbX) / (Kn(indSelected , indSelected) + 1e-8*eye(nbX));
        Ka_KK = Ka(indSelected , indSelected);
        
        % SNRExpectationTh_int = sigmaS^2 * ( trace(Rn_KK * Ka_KK) + ma_K' * Rn_KK * ma_K);
        % SNRExpectationTh_int = sigmaS^2 * ( ma_K' * Rn_KK * ma_K);
        
        % SNRExpectationTh(indSelected,1) = repmat(SNRExpectationTh_int, nbX, 1);
        
        for ind = 1:length(indNonSelected)
            indN = indNonSelected(ind);      %indNonSelected(randomOrder(1));
            
            Rn_KN     = eye(nbX+1) / (Kn([indSelected , indN], [indSelected , indN]) + 1e-8*eye(nbX+1));
            SNRExpectationTh(ind,1) = sigmaS^2 * ( aSelected.' * Rn_KN(1:nbX,1:nbX) * aSelected ...
                + 2 * aSelected.' * Rn_KN(1:nbX,end) * mu_aN(indN) ...
                + Rn_KN(end,end) * (s2_aN(indN) + mu_aN(indN)^2) );
        end
        
        
        %% Compute SNR CDF
        
        SNR_Initial = sigmaS^2 * ( aSelected.' / (Kn(indSelected,indSelected)+1e-8*eye(nbX)) * aSelected);
        
        [SNR_ExpectedMax, indexMax] = max(SNRExpectationTh);
        
        
        for ind = 1:length(indNonSelected)
            indN = indNonSelected(ind);      %indNonSelected(randomOrder(1));
            
            for indTest = 1:nbPtsTest
                SNR_N(indTest) = [aSelected ; aN(indN,indTest)].' / (Kn([indSelected indN],[indSelected indN]) + 1e-8*eye(nbX+1)) * [aSelected ; aN(indN,indTest)];
            end
            
            
            [SNRHist(:,ind), SNRHistAbsc(:,ind)] = hist(SNR_N, 35);
            SNRHist(:,ind) = SNRHist(:,ind) / (sum(SNRHist(:,ind)) * (SNRHistAbsc(2,ind)-SNRHistAbsc(1,ind)));
            
            
            
            Rn_KN     = eye(nbX+1) / (Kn([indSelected , indN], [indSelected , indN]) + 1e-8*eye(nbX+1));
            alpha  = aSelected.' * Rn_KN(1:nbX,1:nbX) * aSelected;
            beta   = aSelected.' * Rn_KN(1:nbX,end);
            gamma  = Rn_KN(end,end);
            lambda = (mu_aN(indN) + beta/gamma)^2 / s2_aN(indN);
            
            
            xSNRTh(:,ind) = linspace(min(SNRHistAbsc(:,ind)), max(SNRHistAbsc(:,ind)), 1e3).';
            SNR_pdfTh(:,ind) = 1/(gamma*s2_aN(indN)) * ncx2pdf(1/(gamma*s2_aN(indN))*(xSNRTh(:,ind) - (alpha - beta^2/gamma)), 1, lambda);
            
            
            probaImprovement(ind,1) = 1-ncx2cdf(1/(gamma*s2_aN(indN))*(SNR_ExpectedMax - (alpha - beta^2/gamma)), 1, lambda);
            probaDegradate(ind,1) = ncx2cdf(1/(gamma*s2_aN(indN))*(SNR_Initial - (alpha - beta^2/gamma)), 1, lambda);
            
            delta_set = [1,5,10,15,20]';
            g_max = size(delta_set,1);
            
            for g = 1:g_max
                threshold(g) = delta_set(g) * (SNR_ExpectedMax-SNR_Initial);
                probaDegradate_BIS(ind,g) = ncx2cdf(1/(gamma*s2_aN(indN))*(SNR_Initial + threshold(g) - (alpha - beta^2/gamma)), 1, lambda);
                
            end
            
            
        end
        
        %%  normalization
        
        % SNRExpectation
        SNRExpectationTh_norm_int = (SNRExpectationTh - min(SNRExpectationTh)) / (max(SNRExpectationTh) - min(SNRExpectationTh));
        
        % CDF_Complementary
        CDF_Complementary_norm= zeros(size(indNonSelected,2),g_max);
        % normalize
        for g = 1:g_max
            CDF_Complementary(:,g) = 1-probaDegradate_BIS(:,g);
            CDF_min(g) = min(CDF_Complementary(:,g));
            CDF_Complementary_norm(:,g) = (CDF_Complementary(:,g)-CDF_min(g)) / (max(CDF_Complementary(:,g)) - CDF_min(g));
        end
        
        
        
        %% failure percentage
        
        index_failure = find(SNR_True_Estimated_diag_norm_int<0);
        
        % SNRExpectation
        index_select_SNRExpectationTh = find(SNRExpectationTh_norm_int>0);
        failure_SNRExpectationTh(Mont) = 100 * numel(intersect(index_select_SNRExpectationTh,index_failure)) / numel(index_failure);
        
        % CDF_Complementary
        
        for gg=1:g_max
            index_select_CDF_Complementary = find(CDF_Complementary_norm(:,gg)>0);
            failure_CDF_Complementary(Mont,gg) = 100 * numel(intersect(index_select_CDF_Complementary,index_failure)) / numel(index_failure);
        end
        
        
    end
    
    failure_SNRExpectation_average(i) = nanmean(failure_SNRExpectationTh)
%     failure_SNRExpectation_average(i) = nanmean(failure_SNRExpectationTh,'all')
    
    failure_CDF_Complementary_average(i,:) = nanmean(failure_CDF_Complementary,1)
    
    
end

%%
semilogx(l_ratio_set,smoothdata(smoothdata(failure_SNRExpectation_average)),'LineWidth',5);
hold on 
semilogx(l_ratio_set,smoothdata(failure_CDF_Complementary_average,1),'LineWidth',5);
% xlim([0,100])
set(gca,'FontSize',55);
% xticks([0:.1:1.1,5:10:100]')

% legend(...
%     ['$J_{E}(X_N|\hat{a}_K)$'],...
%     ['$J_{P}(X_N,\delta=$',num2str(delta_set(1)),'$|\hat{a}_K)$'],...
%     ['$J_{P}(X_N,\delta=$',num2str(delta_set(2)),'$|\hat{a}_K)$'],...
%     ['$J_{P}(X_N,\delta=$',num2str(delta_set(3)),'$|\hat{a}_K)$'],...
%     ['$J_{P}(X_N,\delta=$',num2str(delta_set(4)),'$|\hat{a}_K)$'],...
%     ['$J_{P}(X_N,\delta=$',num2str(delta_set(5)),'$|\hat{a}_K)$'],...
%     'interpreter','latex','FontSize',55);

grid on
ylabel('$\overline{FPR}[\%]$','interpreter','latex','FontSize',55)
xlabel('$\beta=\rho_n/\rho_a$','interpreter','latex','FontSize',70)

set (gcf, 'Position' , [120, 150, 1800, 1300])


