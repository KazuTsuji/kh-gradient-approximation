loglog(eq_nodes,eqweight,'-o','Color','#00F');%loglog
hold on

loglog(ls_nodes,line_search,'-+','Color','#ff8c00');%loglog
hold on

loglog(pmp_nodes,pmp,'-v','Color','#000');%loglog
hold on

loglog(fc_pmp_nodes,fc_pmp,'-o','Color','#0B0');%loglog
hold on

loglog(gcos_nodes,gcos,'-+','Color','#00ffff');%loglog
hold on

loglog(fc_gcos_nodes,fc_gcos,'-v','Color','#cd853f');
hold on

loglog(FC_nodes,FC,'-o','Color','#F00');
hold on

loglog(SBQ_nodes,SBQ_err,'-+','Color','#A2142F')
hold on

loglog(legendre_nodes,legendre_error,'-v','Color','#7E2F8E')
hold on


loglog(linspace(1,maxnodes,maxnodes),order_23,'--','Color','#ff1493')
hold on


xlabel('the number of nodes','FontSize',20)
ylabel('MMD','FontSize',20)
 
hold on


legend({'eq-weight','linesearch','PMP','FC-PMP','gcos','FC-gcos','FC','SBQ','Legendre','(1/n)^{5/4}'},'FontSize',14,'NumColumns',2)
set(gca,'FontSize',16); 

print('gcos_matern32_nodes_with_Legendre.pdf', '-dpdf')
hold off