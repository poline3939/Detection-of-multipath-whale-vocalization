clear
load roc_spec_dist_realdata.mat
load roc_yuan_STFT_realdata.mat
load roc_realdata_spectcorr.mat
N=100;
ii=round(linspace(1,numel(KK2_1),N));
figure;
plot(KK2_2(ii),KK2_1(ii),'m^-','LineWidth',2.5)
hold on
plot(KK2_1(ii),KK2_2(ii),'b*-',KK3_2(ii),KK3_1(ii),'g+-')
xlabel('P_F','FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
ylabel('P_D','FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
axis([0 1 0 1])
legend('STFT','Spectrogram distribution','Spectrogram Correlation')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
hold off
grid on
grid off

N=100;
ii=round(linspace(1,numel(KK1_1),N));
figure
plot(KK1_2(ii),KK1_1(ii),'b*-')

N=100;
ii=round(linspace(1,numel(KK3_1),N));
figure
plot(KK3_2(ii),KK3_1(ii),'g+-')