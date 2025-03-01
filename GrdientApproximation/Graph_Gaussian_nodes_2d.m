loglog(eq_nodes,eqweight,'-o','Color','#00F');%loglog
hold on

loglog(ls_nodes,line_search,'-+','Color','#ff8c00');%loglog
hold on

loglog(pmp_nodes,pmp,'-v','Color','#000');%loglog
hold on

loglog(gcos_nodes,gcos,'-+','Color','#00ffff');%loglog
hold on


xlabel('the number of nodes','FontSize',20)
ylabel('MMD','FontSize',20)
 
hold on


legend({'eq-weight','linesearch','PMP','gcos'},'FontSize',16,'NumColumns',2)
set(gca,'FontSize',16); 
%}
hold off