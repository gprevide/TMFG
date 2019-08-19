figure
size = [50 100 150 200 300 400 500 700 1000 1100 1300 1500 1700 2000];
TMFG_elapsed = [0.375986	0.040929	0.061819	0.098221	0.148019	0.209489	0.281514	0.449734	0.757743	0.860448	1.111806	1.389367	1.69763	2.305995];
TMFGT1_elapsed = [0.250084	0.281018	0.433037	0.706785	0.92695	1.300162	1.600658	2.306769	3.491609	3.895897	4.792208	6.207399	6.67014	8.234642];
TMFGT2_K4_elapsed =	[0.207861	0.148106	0.267811	0.458638	0.768263	1.242158	1.846474	3.250649	6.284583	7.577598	10.162726	13.354423	17.101997	23.131266];
PMFG_elapsed = [0.38917	2.402623	9.465224	17.865617	61.67472	145.778166	289.903112	797.503718	2276.6015	3094.509287	5125.679278	7735.270361	11551.099308	18487.253867];
lx = [0:50:2000];
mdl1 = LinearModel.fit(log(size), log(TMFG_elapsed));
ly1 = exp(mdl1.Coefficients.Estimate(1) + mdl1.Coefficients.Estimate(2) * log(lx))
streq1 = sprintf('log(elapsed time) = %4.2f + %4.2f * log(matrix size)', mdl1.Coefficients.Estimate(1), mdl1.Coefficients.Estimate(2));
mdl2 = LinearModel.fit(log(size), log(TMFGT1_elapsed));
ly2 = exp(mdl2.Coefficients.Estimate(1) + mdl2.Coefficients.Estimate(2) * log(lx))
streq2 = sprintf('log(elapsed time) = %4.2f + %4.2f * log(matrix size)', mdl2.Coefficients.Estimate(1), mdl2.Coefficients.Estimate(2));
mdl3 = LinearModel.fit(log(size), log(TMFGT2_K4_elapsed));
ly3 = exp(mdl3.Coefficients.Estimate(1) + mdl3.Coefficients.Estimate(2) * log(lx))
streq3 = sprintf('log(elapsed time) = %4.2f + %4.2f * log(matrix size)', mdl3.Coefficients.Estimate(1), mdl3.Coefficients.Estimate(2));
mdl4 = LinearModel.fit(log(size), log(PMFG_elapsed));
ly4 = exp(mdl4.Coefficients.Estimate(1) + mdl4.Coefficients.Estimate(2) * log(lx))
streq4 = sprintf('log(elapsed time) = %4.2f + %4.2f * log(matrix size)', mdl4.Coefficients.Estimate(1), mdl4.Coefficients.Estimate(2));

hlines = plot(size, TMFG_elapsed, '^b', size, TMFGT1_elapsed,'+g', size, TMFGT2_K4_elapsed,'xr', size,PMFG_elapsed, 'om',...
    lx, ly1, '-b', lx, ly2, '-g', lx, ly3, '-r', lx, ly4, '-m');
axis([0 2000 0 20000])
set(gca,'xscale','log','yscale','log')
set(gca,'xtick',[0, 50, 200, 500, 1000, 2000],'xticklabel',{'0', '50', '200', '500', '1,000', '2,000'})
set(gca,'ytick',[0.1, 20000],'yticklabel',{'0.1','20,000'})
set(gca,'fontsize',15)
xlabel('$Size \ of \ matrix$','fontsize',15,'interpreter','latex')
ylabel('$Elapsed \ time$','fontsize',15,'interpreter','latex')


set(hlines(1), 'DisplayName', 'TMFG')
set(hlines(2), 'DisplayName', 'TMFGT1')
set(hlines(3), 'DisplayName', 'TMFGT2')
set(hlines(4), 'DisplayName', 'PMFG')
set(hlines(5), 'DisplayName', streq1)
set(hlines(6), 'DisplayName', streq2)
set(hlines(7), 'DisplayName', streq3)
set(hlines(8), 'DisplayName', streq4)

legend1 = legend('show');
set(legend1,...
    'Location','NorthWest',...
    'FontSize',12,...
    'FontName','Cambria')



