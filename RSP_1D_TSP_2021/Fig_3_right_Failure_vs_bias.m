
clear all
% close all
clc

tic


%% Parameters
la1  = .2;      % length scale of a(x)
sa21 = .5;      % power of a(x)


nbX = 3;        % number of initial sensors

sigmaS = 2;  % power of s


% set of the uncertaunty power on a(x) -> \hat{a}(x)
nbSigmaA2 = 20;
sma_ratio = logspace(-4, 4, nbSigmaA2).';


nbA    = 50;         % nb of MC for a(x)


ratio_rho = .01;    % ratio between l_n and l_a => ratio_rho = l_n / l_a, with l_N: length scale of the spatial noise


nbMcMc = 10;        % nb of MC simulations for each configuratio, i.e. pair (a(x), uncertainty on \hat{a}(x))


nbPts = 3e2;    % nb of points 'x' in the grid (for optimization)



%%
x = linspace(0, 1, nbPts).';

% Covariance functions
kSquare = @(x,s2,l)(s2 * exp(-squareform(pdist(x,'squaredeuclidean')) / (2*l^2)));
k       = @(x,xP,s2,l) (s2 * exp(-(x-xP).^2 / (2*l^2)));




failure_size_mont = zeros(nbA, nbSigmaA2, nbMcMc);

for inda = 1:nbA    % loop MC on a(x)
    % Generate a
    Ka = kSquare(x, sa21, la1);
    La = chol(Ka + 1e-8*eye(nbPts), 'lower');
    a  = La * randn(nbPts, 1);
    
    
    % Selection X_(N-1) = [x_1, ... , x_{N-1}]
    
    indSelected = fix([1:nbX]/(nbX+1)*nbPts);
    %indSelected = randperm(nbPts);
    indSelected = sort(indSelected(1:nbX));
    xSelected   = x(indSelected);
    
    % aSelected = ma(indSelected);
    aSelected = a(indSelected);
    
    
    figure(1); clf
    axs = axes;
    plot(x,a); hold on
    plot(xSelected, aSelected, 'o');
    
    grid on
    lbl = xlabel('$x$');
    lbl.Interpreter = 'latex';
    lbl = ylabel('$a(x)$');
    lbl.Interpreter = 'latex';
    axs.FontSize = 16;
    
    
    % Failure computation
    for indUncertainty = 1:nbSigmaA2     % loop on uncertainty on a(x)
        figure(10); clf
        subplot(211)
        title('SNR')
        hold on
        
        subplot(212)
        title('SNR normalized')
        hold on
        
        
        
        for Mont = 1:nbMcMc     % loop MC on \hat{a}(x)
            [inda indUncertainty Mont]
            
            % Parameters of \hat{a}
            la  = la1;
            %             sa2 = sa2_set(indUncertainty);
            sa2 = 0.01;
            
            % Parameters of the mean function
            lma  = la;
            %             sma2 = sa2_set(indUncertainty);
            
            sma2 = sma_ratio(indUncertainty)*sa2;   % here sa2_set represents the ratio between "sigma_ma" and "sigma_a"
            
            % Parameters of noise
            ln  = ratio_rho*la;
            sn2 = .01;
            
            
            % 1. Generate ma
            Kma = kSquare(x, sma2, lma);
            Lma = chol(Kma + 1e-8*eye(nbPts), 'lower');
            ma  = a + Lma * randn(nbPts, 1);
            
            maSelected = ma(indSelected);
            
            
            
            figure(2); clf
            axs = axes;
            plt = fill([x ; flipud(x)], [ma + 2*sqrt(sa2) ; flipud(ma) - 2*sqrt(sa2)], .85*[1 1 1]); hold on
            plt.HandleVisibility = 'off';
            plt.LineStyle = 'none';
            
            plt = plot(x,a);
            plt.Color = axs.ColorOrder(1,:);
            plt = plot(xSelected, aSelected, 'o');
            plt.MarkerEdgeColor = axs.ColorOrder(1,:);
            plt.HandleVisibility = 'off';
            
            plt = plot(x, ma);
            plt.Color = axs.ColorOrder(2,:);
            plt = plot(xSelected, maSelected, 'o');
            plt.Color = axs.ColorOrder(2,:);
            plt.HandleVisibility = 'off';
            
            grid on
            
            lbl = xlabel('$x$');
            lbl.Interpreter = 'latex';
            lbl = ylabel('$a(x)$');
            lbl.Interpreter = 'latex';
            axs.FontSize = 16;
            lgd = legend({'$a(x)$' ; '$\hat{a}(x)$'});
            lgd.Interpreter = 'latex';
            drawnow
            
            
            
            % 2. Generate noise
            Kn = kSquare(x, sn2, ln);
            
            
            %         % 3. Generate \hat{a}_N | vec(a)_N-1
            indNonSelected = setdiff(1:nbPts, indSelected);
            %
            %         nbPtsTest = 100;
            %
            %         for indN = 1:nbPts
            %             if ismember(indN, indNonSelected)
            %                 %indNonSelected   %(randomOrder(1));
            %
            %                 mu_aN(indN,1) = ma(indN) + Ka(indSelected, indN).' * ((Ka(indSelected,indSelected)+1e-8*eye(nbX)) \ (aSelected - maSelected));
            %                 s2_aN(indN,1) = Ka(indN,indN) - Ka(indSelected, indN).' / (Ka(indSelected,indSelected)+1e-8*eye(nbX)) * Ka(indSelected, indN);
            %
            %                 aN(indN, 1:nbPtsTest) = mu_aN(indN) + sqrt(s2_aN(indN)) * randn(1, nbPtsTest);
            %
            %             else
            %                 mu_aN(indN,1) = a(indN);
            %                 s2_aN(indN,1) = 0;
            %
            %                 aN(indN, 1:nbPtsTest) = a(indN);
            %
            %             end
            %
            %         end
            
            % 4. Optimal SNR: "perfect knowledge on a(x)"
            aSelected_True = a(indSelected);
            for indN = 1:nbPts
                if ismember(indN, indNonSelected)
                    a_total = [aSelected_True ; a(indN)];
                    x_total = [x(indSelected) ; x(indN)];
                    Rn1 = kSquare(x_total, sn2, ln);
                    
                else
                    a_total = [aSelected_True];
                    x_total = [x(indSelected)];
                    Rn1 = kSquare(x_total, sn2, ln);
                end
                SNR_optimal(indN,1) = sigmaS^2 * a_total' * (Rn1 \ a_total);
            end
            
            
            
            % 5. True SNR(f_hat): "actual SNR with uncertain a(x)"
            flagPrior = 1;
            
            zK     = ma(indSelected);
            tol    = 1e-10;
            xK     = x(indSelected);
            ma_K   = ma(indSelected);
            
            
            % a. Estimation of a(x) from z(x)
            if flagPrior == 1
                aEst = ma;
                
                [xx_, xy_] = meshgrid(x, x);
                Sigma_aEst = k(xx_, xy_, sa2, la);
                
            else
                [xx_kk, xy_kk] = meshgrid(xK, xK);
                SigmaA_kk = k(xx_kk, xy_kk, sa2, la);
                
                [xx_k, xy_k] = meshgrid(xK, x);
                SigmaA_k  = k(xx_k, xy_k, sa2, la);
                
                aEst = ma + SigmaA_k * (SigmaA_kk \ (zK - ma_K));
                
                [xx_, xy_] = meshgrid(x, x);
                SigmaA_k  = k(xx_, xy_, sa2, la);
                
                Sigma_aEst = SigmaA_k - SigmaA_k * (SigmaA_kk \ SigmaA_k.');
            end
            
            
            % b. Computation of the true SNR with estimated a(x)
            SNR_True_Estimated = zeros(nbPts,1);
            
            for indSensor1 = 1:nbPts
                indexExtendedSensors = unique([indSelected' ; indSensor1]);
                
                xExtended     = x(indexExtendedSensors);
                aExtended     = a(indexExtendedSensors);
                aEst_Extended = ma(indexExtendedSensors);
                
                xDist                   = squareform(pdist(xExtended, 'squaredeuclidean'));
                Rn_Extended_computation = kSquare(xExtended, sn2, ln) + tol*sn2*eye(length(xExtended));
                SigmaA_Extended         = Sigma_aEst(indexExtendedSensors, indexExtendedSensors);
                
                
                aTmp = Rn_Extended_computation \ aExtended;
                SNR_True_Estimated(indSensor1,1) = ...
                    sigmaS^2 * ( (aTmp.' * aEst_Extended)^2 + aTmp.' * SigmaA_Extended * aTmp ) / ...
                    ( (aEst_Extended.' / Rn_Extended_computation) * aEst_Extended + ...
                    trace(Rn_Extended_computation \ SigmaA_Extended) );
            end
            
            
            
            toc
            
            SNR_True_Estimated_int = SNR_True_Estimated(indSelected(1));
            
            SNR_True_Estimated_diag_norm_int = (SNR_True_Estimated - SNR_True_Estimated_int) / (max(SNR_True_Estimated) - SNR_True_Estimated_int);
            
            
            figure(3); clf
            axs = axes;
            plot(x, 10*log10(SNR_optimal)); hold on
            plt = plot(xSelected, 10*log10(SNR_optimal(indSelected)), 'o');
            plt.MarkerEdgeColor = axs.ColorOrder(1,:);
            plt.HandleVisibility = 'off';
            plt = plot([0 1], 10*log10(SNR_optimal(indSelected(1)))*ones(2,1), '--');
            plt.Color = axs.ColorOrder(1,:);
            plt.HandleVisibility = 'off';
            
            plt = plot(x, 10*log10(SNR_True_Estimated));
            plt.Color = axs.ColorOrder(2,:);
            plt = plot(xSelected, 10*log10(SNR_True_Estimated(indSelected)), 'o');
            plt.Color = axs.ColorOrder(2,:);
            plt.HandleVisibility = 'off';
            plt = plot([0 1], 10*log10(SNR_True_Estimated(indSelected(1)))*ones(2,1), '--');
            plt.Color = axs.ColorOrder(2,:);
            plt.HandleVisibility = 'off';
            grid on
            
            
            lgd = legend({'Optimal', 'With uncertainty'});
            grid on
            
            lbl = xlabel('$x$');
            lbl.Interpreter = 'latex';
            lbl = ylabel('$SNR$ [dB]');
            lbl.Interpreter = 'latex';
            axs.FontSize = 16;
            
            drawnow
            
            
            figure(10); % tmp
            subplot(211)
            plt = plot(x, 10*log10(SNR_True_Estimated)); hold on
            plt.Color = axs.ColorOrder(1,:);
            grid on
            subplot(212)
            plt = plot(x, SNR_True_Estimated_diag_norm_int); hold on
            plt.Color = axs.ColorOrder(1,:);
            grid on
            drawnow
            
            
            %% Compute SNR mean
            %
            %         % SNRExpectationTh = zeros(nbPts,1);
            %
            %         Rn_KK     = eye(nbX) / (Kn(indSelected , indSelected) + 1e-8*eye(nbX));
            %         Ka_KK = Ka(indSelected , indSelected);
            %
            %         % SNRExpectationTh_int = sigmaS^2 * ( trace(Rn_KK * Ka_KK) + ma_K' * Rn_KK * ma_K);
            %         % SNRExpectationTh_int = sigmaS^2 * ( ma_K' * Rn_KK * ma_K);
            %
            %         % SNRExpectationTh(indSelected,1) = repmat(SNRExpectationTh_int, nbX, 1);
            %
            %         for ind = 1:length(indNonSelected)
            %             indN = indNonSelected(ind);      %indNonSelected(randomOrder(1));
            %
            %             Rn_KN     = eye(nbX+1) / (Kn([indSelected , indN], [indSelected , indN]) + 1e-8*eye(nbX+1));
            %             SNRExpectationTh(ind,1) = sigmaS^2 * ( aSelected.' * Rn_KN(1:nbX,1:nbX) * aSelected ...
            %                 + 2 * aSelected.' * Rn_KN(1:nbX,end) * mu_aN(indN) ...
            %                 + Rn_KN(end,end) * (s2_aN(indN) + mu_aN(indN)^2) );
            %         end
            %
            
            %% Compute SNR CDF
            
            %         SNR_Initial = sigmaS^2 * ( aSelected.' / (Kn(indSelected,indSelected)+1e-8*eye(nbX)) * aSelected);
            %
            %         [SNR_ExpectedMax, indexMax] = max(SNRExpectationTh);
            %
            %
            %         for ind = 1:length(indNonSelected)
            %             indN = indNonSelected(ind);      %indNonSelected(randomOrder(1));
            %
            %             for indTest = 1:nbPtsTest
            %                 SNR_N(indTest) = [aSelected ; aN(indN,indTest)].' / (Kn([indSelected indN],[indSelected indN]) + 1e-8*eye(nbX+1)) * [aSelected ; aN(indN,indTest)];
            %             end
            %
            %
            %             [SNRHist(:,ind), SNRHistAbsc(:,ind)] = hist(SNR_N, 35);
            %             SNRHist(:,ind) = SNRHist(:,ind) / (sum(SNRHist(:,ind)) * (SNRHistAbsc(2,ind)-SNRHistAbsc(1,ind)));
            %
            %
            %
            %             Rn_KN     = eye(nbX+1) / (Kn([indSelected , indN], [indSelected , indN]) + 1e-8*eye(nbX+1));
            %             alpha  = aSelected.' * Rn_KN(1:nbX,1:nbX) * aSelected;
            %             beta   = aSelected.' * Rn_KN(1:nbX,end);
            %             gamma  = Rn_KN(end,end);
            %             lambda = (mu_aN(indN) + beta/gamma)^2 / s2_aN(indN);
            %
            %
            %             xSNRTh(:,ind) = linspace(min(SNRHistAbsc(:,ind)), max(SNRHistAbsc(:,ind)), 1e3).';
            %             SNR_pdfTh(:,ind) = 1/(gamma*s2_aN(indN)) * ncx2pdf(1/(gamma*s2_aN(indN))*(xSNRTh(:,ind) - (alpha - beta^2/gamma)), 1, lambda);
            %
            %
            %             probaImprovement(ind,1) = 1-ncx2cdf(1/(gamma*s2_aN(indN))*(SNR_ExpectedMax - (alpha - beta^2/gamma)), 1, lambda);
            %             probaDegradate(ind,1) = ncx2cdf(1/(gamma*s2_aN(indN))*(SNR_Initial - (alpha - beta^2/gamma)), 1, lambda);
            %
            %             g_max = 4;
            %
            %             for g = 1:g_max
            %                 threshold(g) = .25 * g * (SNR_ExpectedMax-SNR_Initial);
            %                 probaDegradate_BIS(ind,g) = ncx2cdf(1/(gamma*s2_aN(indN))*(SNR_Initial + threshold(g) - (alpha - beta^2/gamma)), 1, lambda);
            %
            %             end
            %
            %
            %         end
            
            %%  normalization
            
            % SNRExpectation
            %         SNRExpectationTh_norm_int = (SNRExpectationTh - min(SNRExpectationTh)) / (max(SNRExpectationTh) - min(SNRExpectationTh));
            %
            %         % CDF_Complementary
            %         CDF_Complementary_norm= zeros(size(indNonSelected,2),g_max);
            %         % normalize
            %         for g = 1:g_max
            %             CDF_Complementary(:,g) = 1-probaDegradate_BIS(:,g);
            %             CDF_min(g) = min(CDF_Complementary(:,g));
            %             CDF_Complementary_norm(:,g) = (CDF_Complementary(:,g)-CDF_min(g)) / (max(CDF_Complementary(:,g)) - CDF_min(g));
            %         end
            %
            
            
            %% failure percentage
            
            index_failure = find(SNR_True_Estimated_diag_norm_int<0);
            failure_size(Mont) = size(index_failure,1);
            
            %         % SNRExpectation
            %         index_select_SNRExpectationTh = find(SNRExpectationTh_norm_int>0);
            %         failure_SNRExpectationTh(Mont) = 100 * numel(intersect(index_select_SNRExpectationTh,index_failure)) / numel(index_failure);
            %
            %         % CDF_Complementary
            %
            %         for gg=1:g_max
            %             index_select_CDF_Complementary = find(CDF_Complementary_norm(:,gg)>0);
            %             failure_CDF_Complementary(Mont,gg) = 100 * numel(intersect(index_select_CDF_Complementary,index_failure)) / numel(index_failure);
            %         end
            %
            
        end
        
        %     failure_SNRExpectation_average(i) = nanmean(failure_SNRExpectationTh)
        %
        %     failure_CDF_Complementary_average(i,:) = nanmean(failure_CDF_Complementary,1)
        failure_size_mont(inda, indUncertainty,:) = failure_size;
    end
    
end

%%
nbSubplotRow = fix(sqrt(nbA));
nbSubplotCol = ceil(nbA/nbSubplotRow);
figure(4); clf
for indPlot = 1:nbA
    
    axs = subplot(nbSubplotRow, nbSubplotCol, indPlot);
    
    plt = fill([sma_ratio.^(1/2) ; flipud(sma_ratio.^(1/2))], 100*[quantile(squeeze(failure_size_mont(indPlot,:,:)).', .25).' ; flipud(quantile(squeeze(failure_size_mont(indPlot,:,:)).', .75).')]/nbPts, .75*[1 1 1]); hold on
    plt.HandleVisibility = 'off';
    plt.LineStyle = 'none';
    
    plt = plot(sma_ratio.^(1/2), 100*median(squeeze(failure_size_mont(indPlot,:,:)).')/nbPts); hold on
    plt.Color = axs.ColorOrder(1,:);
    
    
    grid on
    axs.XScale = 'log';
    
    ylabel('failure [%]')
    xlabel('sigma_a (i.e. ''Uncertainty on a(x)'')')
end


% Global performance
failure_size_Global = reshape(permute(failure_size_mont, [2 1 3]), nbSigmaA2, nbA*nbMcMc);

%%

figure(5); clf
axs = axes;
sma2_vect = sa2.*sma_ratio;
plt = fill([sma2_vect.^(1/2) ; flipud(sma2_vect.^(1/2))], 100*[quantile(failure_size_Global.', .25).' ; flipud(quantile(failure_size_Global.', .75).')]/nbPts, .85*[1 1 1]); hold on
plt.HandleVisibility = 'off';
plt.LineStyle = 'none';

plt = plot(sma2_vect.^(1/2), 100*median(failure_size_Global.')/nbPts,'LineWidth',3.5,'Color',[0 0.4470 0.7410]); hold on
plt.Color = axs.ColorOrder(1,:);


grid on
axs.XScale = 'log';
set(gca,'FontSize',30);
ylabel('Failure [$\%$]','fontsize',35,'interpreter','latex')
xlabel('$\sigma_{b}$','fontsize',40,'interpreter','latex')
% ylim([0,40])
set(gcf,'Position',[100 200 1200 1000])

%%k



%
% figure(1)
% plot(sa2_set,smoothdata(failure_CDF_Complementary_average,1),'LineWidth',2.5);
% hold on
% plot(sa2_set,smoothdata(smoothdata(failure_SNRExpectation_average)),'LineWidth',2.5);
% % xlim([0,100])
% set(gca,'FontSize',16);
% % xticks([0:.1:1.1,5:10:100]')
% legend('$G_{cmpl},\delta=0.25$','$G_{cmpl},\delta=0.50$',...
%     '$G_{cmpl},\delta=0.75$','$G_{cmpl},\delta=1.00$','$\hat{g}(X_N)$','interpreter','latex','FontSize',20);
% grid on
% ylabel('Failure','FontSize',20)
% xlabel('$\rho_a$','interpreter','latex','FontSize',24)




% figure(4); clf
% semilogx(sa2_set,smooth(smooth(sa2_set,failure_size_average)),'LineWidth',2.5)
% grid on
% % xlim(1e-4,5)
% ylabel('\textbf{Size of the Failure Region}','interpreter','latex','FontSize',20)
% xlabel('\textbf{Uncertainty on spatial gain} $a$ $(\sigma_a)$','interpreter','latex','FontSize',24)
% set(gca,'FontSize',20);
%
% figure(5); clf
% fill([sa2_set' ; flipud(sa2_set')], [smooth(failure_size_average'+var(failure_size_mont')'.^(1/2)) ; ...
%     flipud(smooth(failure_size_average'-var(failure_size_mont')'.^(1/2)))], .85*ones(1,3), 'edgeColor', 'k','LineWidth',1)
% hold on
% plot(sa2_set',smooth(failure_size_average'),'b','LineWidth',2.5); hold on
% grid on
% ylabel('\textbf{Size of the Failure Region}','interpreter','latex','FontSize',20)
% xlabel('\textbf{Uncertainty on spatial gain} $a$ $(\sigma_a)$','interpreter','latex','FontSize',24)
% set(gca, 'XScale', 'log')
% set(gca,'FontSize',20);