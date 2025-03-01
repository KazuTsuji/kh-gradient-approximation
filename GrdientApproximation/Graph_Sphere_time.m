
semilogy(eq_nodes,eqweight,'-o','Color','#00F');
hold on

semilogy(ls_nodes,line_search,'-+','Color','#ff8c00');
hold on

semilogy(pmp_nodes,pmp,'-v','Color','#000');
hold on

semilogy(fc_pmp_nodes,fc_pmp,'-o','Color','#0B0');
hold on

semilogy(gcos_nodes,gcos,'-+','Color','#00ffff');
hold on


semilogy(fc_gcos_nodes,fc_gcos,'-v','Color','#cd853f');
hold on

semilogy(FC_nodes,FC,'-o','Color','#F00');
hold on

semilogy(SBQ_nodes,SBQ_err,'-+','Color','#A2142F')
hold on

xlabel('computational time (s)','FontSize',20)
ylabel('MMD','FontSize',20)
 
hold off


legend({'eq-weight','linesearch','PMP','FC-PMP','gcos','FC-gcos','FC','SBQ'},'FontSize',14,'NumColumns',2)
