%semilogy(monte_carlo_nodes,monte_carlo_err,'-v','Color','#7E2F8E');
%hold on

%semilogy(eq_nodes,eqweight,'-o','Color','#00F');
%hold on

semilogy(ls_nodes,line_search,'-+','Color','#ff8c00');
hold on


%semilogy(pmp_nodes,pmp,'-v','Color','#000');
%hold on

%semilogy(fc_pmp_nodes,fc_pmp,'-o','Color','#0B0');
%hold on

semilogy(gcos_nodes,gcos,'-+','Color','#00ffff');
hold on
%semilogy(gcos_nodes,gcos,'-+','Color','blue');%スライド用


%semilogy(fc_gcos_nodes,fc_gcos,'-v','Color','#cd853f');
semilogy(fc_gcos_nodes,fc_gcos,'-v','Color','green');
hold on


semilogy(FC_nodes,FC,'-o','Color','#F00');
hold on

semilogy(SBQ_nodes,SBQ_err,'-+','Color','#A2142F')
hold on


%semilogy(linspace(1,maxnodes,maxnodes),order_23,'--','Color','#ff1493')
%hold on
%semilogy(linspace(1,maxnodes,maxnodes),order_25,'--','Color','#ff69b4')
%hold on
%semilogy(linspace(1,maxnodes,maxnodes),order_exp,'--','Color','#ff69b4')
%hold on


%xlabel('the number of nodes','FontSize',20)
xlabel('computational time (s)','FontSize',20)
ylabel('MMD','FontSize',20)
 
hold on


%legend({'Monte carlo','eq-weight','linesearch','PMP','FC-PMP','gcos','FC-gcos','FC','SBQ','(1/n)^{5/4}'},'FontSize',14,'NumColumns',2)
%legend({'Monte carlo','eq-weight','linesearch','PMP','FC-PMP','gcos','FC-gcos','FC','SBQ','(1/n)'},'FontSize',16,'NumColumns',2)
%legend({'Monte carlo','eq-weight','linesearch','PMP','FC-PMP','gcos','FC-gcos','FC','SBQ','(1/n)^{4/3}'},'FontSize',14,'NumColumns',2)
%legend({'eq-weight','linesearch','PMP','FC-PMP','gcos','FC-gcos','FC','SBQ','exp(-(n)^{1/2})'},'FontSize',14,'NumColumns',2)
%legend({'eq-weight','linesearch','Away','PCG','SBQ','BPCG','BPCG-lazified'},'FontSize',20,'NumColumns',2)
%legend({'eq-weight','linesearch','PMP','gcos'},'FontSize',20,'NumColumns',2)
%legend({'eq-weight','linesearch','PMP','FC-PMP','gcos','FC-gcos','FC','SBQ'},'FontSize',18,'NumColumns',2)

%スライド用
%legend({'Kernel Herding','cos最大化','cos最大化+係数最適化','(1/n)^{5/4}(最適レート)'},'FontSize',20,'NumColumns',2)
%legend({'linesearch','gcos','FC-gcos','FC','SBQ','(1/n)^{5/4}(最適レート)'},'FontSize',20,'NumColumns',2)
legend({'linesearch','gcos','FC-gcos','FC','SBQ'},'FontSize',22,'NumColumns',2)
lgd = legend;
lgd.NumColumns = 2; % 
%legend({'Monte carlo','eq-weight','linesearch','PMP','FC-PMP','gcos','FC-gcos','FC','SBQ','(1/n)^{5/4}'},'FontSize',18,'NumColumns',2)
%legend({'eq-weight','linesearch','PMP','FC-PMP','gcos','FC-gcos','FC','SBQ'},'FontSize',20,'NumColumns',2)
%legend({'FC-PMP','FC-gcos','FC'},'FontSize',20,'NumColumns',2)
%txt = {' FC-gcos ','  (係数最適化)', '  \downarrow  '};
txt = {'\leftarrow FC-gcos(係数最適化)'};
%text(200,0.001,txt,'FontSize',16,'FontWeight','bold','Color','black')
text(40,0.001,txt,'FontSize',16,'FontWeight','bold','Color','black')
set(gca,'FontSize',16); 
%}
hold off