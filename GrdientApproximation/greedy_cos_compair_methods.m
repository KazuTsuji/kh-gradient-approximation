%compair the four methods on the convergence speed of the worst case errors for the number of nodes when dim=2

% Gaussian dim =2 
%{
epsilon=1;
n=400;
dim=2;
meanarr=[[0,0]];
partion=[1];
var=[1];
delta=1e-4;
K=30;
T=35;
functype=2;
option=4;
maxtime=600;
maxnodes=100;
%}

% Matern3/2 dim=2
%{
epsilon=1;
n=80;
dim=2;
meanarr=[[0,0]];
partion=[1];
var=[1];
delta=1e-5;
K=20;
T=35;
functype=3;
option=4;
maxtime=14400;
maxnodes=200;
%}

% Sovolev Sphere
%{
epsilon=1;
n=80;
dim=2;
meanarr=[[0,0]];
partion=[1];
var=[1];
delta=1e-5;
K=20;
T=35;
functype=6;
option=4;
maxtime=14400;
maxnodes=200;
%}
%Gaussian dim =3 
epsilon=1;
n=20;
dim=3;
meanarr=[[0,0,0]];
partion=[1];
var=[1];
delta=1e-4;
K=30;
T=35;
functype=2;
option=4;
maxtime=3600;
maxnodes=130;

%Matern3/2 dim=3
%{
epsilon=1;
n=20;
dim=3;
meanarr=[[0,0,0]];
 
iteration=1;
 
partion=[1];
 
var=[1];
delta=0.0001;
K=20;
T=35;
functype=3;
option=4;
maxtime=14400;
maxnodes=200;
%}

points =linspace(-1,1,n) ;
 
c_var=zeros(length(partion),dim);
for jj=1:dim
    for ii=1:length(partion)
        c_var(ii,jj)=1/sqrt(var(ii)) *(sqrt(pi)/2)*(erf(sqrt(var(ii))*(1-meanarr(ii,jj)))-  erf(sqrt(var(ii))*(-1-meanarr(ii,jj)))) ;   
    end
end
mu=zeros(1,n^(dim));
 
for ii=0:n^(dim)-1;
    indexes= func_index(ii,n,dim);
    value=1;    
    if functype==1;
        val=zeros(1,length(partion))+1;
        for jj=1:length(partion)
            for kk=1:dim
                int_value1= erf(sqrt(var(jj))* ( points(indexes(kk)) - (2*meanarr(jj,kk)*var(jj)+epsilon)/(2*var(jj) ) ) ) - erf(sqrt(var(jj))* ( -1 - (2*meanarr(jj,kk)*var(jj)+epsilon)/(2*var(jj) ) ) );
                int1= 1/(c_var(jj,kk)) * exp(-epsilon*points(indexes(kk)) -var(jj)* (meanarr(jj,kk))^2 )* exp((2*meanarr(jj,kk)*var(jj)+epsilon)^2 /(4*var(jj) )) * sqrt(pi)/2 *1/sqrt(var(jj))* int_value1;
                
                int_value2= erf(sqrt(var(jj))* ( 1 - (2*meanarr(jj,kk)*var(jj)- epsilon)/(2*var(jj) ) ) ) - erf(sqrt(var(jj))* ( points(indexes(kk)) - (2*meanarr(jj,kk)*var(jj)-epsilon)/(2*var(jj)) ) );
                int2= 1/(c_var(jj,kk)) * exp(epsilon*points(indexes(kk)) -var(jj)* (meanarr(jj,kk))^2 )* exp((2*meanarr(jj,kk)*var(jj)-epsilon)^2 /(4*var(jj) )) * sqrt(pi)/2 *1/sqrt(var(jj))* int_value2;
                
                val(jj)=val(jj)*(int1+int2);
            end
        end
        value= partion*transpose(val);
       
    elseif functype==2;
        val=zeros(1,length(partion))+1;
        for jj=1:length(partion)
            for kk=1:dim
                int_value= erf(sqrt(var(jj)+epsilon)+ (meanarr(jj,kk)*var(jj) +epsilon*points(indexes(kk)))/sqrt(var(jj)+epsilon) )+ erf(sqrt(var(jj)+epsilon) - (meanarr(jj,kk)*var(jj)+epsilon*points(indexes(kk)))/sqrt(var(jj)+epsilon) );
                val(jj)=val(jj)* 1/(c_var(jj,kk))* exp(-var(jj)*(meanarr(jj,kk))^2 - epsilon*(points(indexes(kk)))^2 + ((meanarr(jj,kk)*var(jj)+epsilon*(points(indexes(kk))))^2)/(var(jj)+epsilon) )* 1/sqrt(var(jj)+epsilon) * (sqrt(pi)/2)*int_value ; 
            end
        end
        value= partion*transpose(val);
    end
    mu(ii+1)=value;
    
    if functype==3;
        p_ar=[];
        for l=1:dim;
            p_ar=[p_ar,points(indexes(l))];
        end
        if dim==2
            ff=@(x,y)(1+sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2 ) ) .*exp(-sqrt( (x-p_ar(1)).^2+(y-p_ar(2)).^2 ));
            mu(ii+1)=(0.5)^dim *integral2(ff,-1,1,-1,1,'RelTol',1e-10);
        elseif dim==3
            ff=@(x,y,z)(1+sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2 ) ) .*exp(-sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2));
            mu(ii+1)=(0.5)^dim *integral3(ff,-1,1,-1,1,-1,1,'RelTol',1e-6);
        end
    
    elseif functype==4;
        p_ar=[];
        for l=1:dim;
            p_ar=[p_ar,points(indexes(l))];
        end
        if dim==2
            ff=@(x,y)(1+sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2 )+((x-p_ar(1)).^2+(y-p_ar(2)).^2 )/3 ) .*exp(-sqrt( (x-p_ar(1)).^2+(y-p_ar(2)).^2 ));
            mu(ii+1)=(0.5)^dim *integral2(ff,-1,1,-1,1,'RelTol',1e-10);
        elseif dim==3
            ff=@(x,y,z)(1+sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2 )+((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2 )/3 ) .*exp(-sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2));
            mu(ii+1)=(0.5)^dim *integral3(ff,-1,1,-1,1,-1,1,'RelTol',1e-10);
        end
        
    elseif functype==5;
        p_ar=[];
        for l=1:dim;
            p_ar=[p_ar,points(indexes(l))];
        end
        if dim==2
            ff=@(x,y) exp(-epsilon*sqrt( (x-p_ar(1)).^2+(y-p_ar(2)).^2 ));
            mu(ii+1)=(0.5)^dim *integral2(ff,-1,1,-1,1);
        elseif dim==3
            ff=@(x,y,z) exp(-epsilon*sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2));
            mu(ii+1)=(0.5)^dim *integral3(ff,-1,1,-1,1,-1,1);
        end
    elseif functype==6
        mu(ii+1)=4/3;      
    end
end

order_1=[];
order_2=[];
num_nodes=[];


output_monte_carlo= monte_carlo(maxnodes,epsilon,dim,functype,option);
len_monte_carlo=  int64(length(output_monte_carlo)/2);
monte_carlo_err =[];
monte_carlo_nodes=[];

fact=0;
stride=5;

for i=1: len_monte_carlo
    m= int64(output_monte_carlo(i));
    if fact*stride < m
        monte_carlo_err = [monte_carlo_err, output_monte_carlo(len_monte_carlo +i)];
        monte_carlo_nodes=[monte_carlo_nodes,m];
        order_1=[order_1,1/double(m)];
        order_2=[order_2,1/sqrt(double(m))];
        fact=fact+1;
    end
end


output_linesearch= linesearch(maxtime,maxnodes,epsilon,n,dim,points,c_var,mu,var,meanarr,partion,functype,option);
len_ls=  int64(length(output_linesearch)/2);
line_search =[];
ls_nodes=[];

fact=0;
for i=1: len_ls
    m= output_linesearch(i);
    if fact*stride <m
        line_search = [line_search, output_linesearch(len_ls +i)];
        ls_nodes=[ls_nodes,m];
        fact=fact+1;
    end
end

output_eqweight= eqweight_herding(maxtime,maxnodes,epsilon,n,dim,points,mu,c_var,meanarr,var,partion,functype,option);
len_eq=  int64(length(output_eqweight)/2);
eqweight =[];
eq_nodes=[];

fact=0;
for i=1: len_eq
    m= output_eqweight(i);
    if fact*stride < m
        eqweight = [eqweight, output_eqweight(len_eq +i)];
        eq_nodes=[eq_nodes,m];
        fact=fact+1;
    end
end


output_pmp= positive_matching_pursuit(maxtime,maxnodes,epsilon,n,dim,points,mu,c_var,meanarr,var,K,delta,partion,functype,option);
len_pmp=  int64(length(output_pmp)/2);
pmp =[];
pmp_nodes=[];

fact=0;
for i=1: len_pmp
    m= output_pmp(i);
    if fact*stride < m
        pmp = [pmp, output_pmp(len_pmp +i)];
        pmp_nodes=[pmp_nodes,m];
        fact=fact+1;
    end
end

output_fc_pmp= FC_PMP(maxtime,maxnodes,epsilon,n,dim,points,mu,c_var,meanarr,var,K,delta,partion,functype,option);
len_fc_pmp=  int64(length(output_fc_pmp)/2);
fc_pmp =[];
fc_pmp_nodes=[];

fact=0;
for i=1: len_fc_pmp
    m= output_fc_pmp(i);
    if fact*stride < m
        fc_pmp = [fc_pmp, output_fc_pmp(len_fc_pmp +i)];
        fc_pmp_nodes=[fc_pmp_nodes,m];
        fact=fact+1;
    end
end


output_gcos= greedy_cos(maxtime,maxnodes,epsilon,n,dim,points,mu,c_var,meanarr,var,K,delta,partion,functype,option);
len_gcos=  int64(length(output_gcos)/2);
gcos =[];
gcos_nodes=[];

fact=0;
for i=1: len_gcos
    m= output_gcos(i);
    if fact*stride < m
        gcos = [gcos, output_gcos(len_gcos +i)];
        gcos_nodes=[gcos_nodes,m];
        fact=fact+1;
    end
end


output_fc_gcos= FC_greedy_cos(maxtime,maxnodes,epsilon,n,dim,points,mu,c_var,meanarr,var,K,delta,partion,functype,option);
len_fc_gcos=  int64(length(output_fc_gcos)/2);
fc_gcos =[];
fc_gcos_nodes=[];

fact=0;
for i=1: len_fc_gcos
    m= output_fc_gcos(i);
    if fact*stride <m
        fc_gcos = [fc_gcos, output_fc_gcos(len_fc_gcos +i)];
        fc_gcos_nodes=[fc_gcos_nodes,m];
        fact=fact+1;
    end
end


output_FC= Fully_corrective(maxtime,maxnodes,epsilon,n,dim,points,c_var,mu,var,meanarr,partion,functype,option);
len_FC=  int64(length(output_FC)/2);
FC =[];
FC_nodes=[];

fact=0;
for i=1: len_FC
    m= output_FC(i);
    if fact*stride< m
        FC = [FC, output_FC(len_FC +i)];
        FC_nodes=[FC_nodes,m];
        fact=fact+1;
    end
end

output_SBQ=  SBQ(maxtime,maxnodes,epsilon,n,dim,points,mu,c_var,meanarr,var,partion,functype,option);
len_SBQ=  int64(length(output_SBQ)/2);
SBQ_err =[];
SBQ_nodes=[];

fact=0;
for i=1: len_SBQ
    m= output_SBQ(i);
    if fact*stride < m
        SBQ_err = [SBQ_err, output_SBQ(len_SBQ +i)];
        SBQ_nodes=[SBQ_nodes,m];
        fact=fact+1;
    end
end


%Gausslegendre
legendre_nodes=[];
legendre_error=[];

count_legendre = 3;
while(count_legendre^dim <maxnodes)
    error = GaussLegendre(count_legendre, epsilon, dim, functype, option)
    legendre_nodes = [legendre_nodes, count_legendre^dim]
    legendre_error = [legendre_error, error]
    count_legendre = count_legendre+1;
end


%Sphere
sphere_nodes =[];
sphere_error =[];

fileNameList = dir('SphereDesign/*.txt');
for ii = 1: length(fileNameList);
    barSplit = split(fileNameList(ii).name,"_");
    dotSplit = split(barSplit(2),".");
    numNodesStr = dotSplit(1);
    numNodes =  str2double(numNodesStr);
    if(numNodes > maxnodes)
        break
    else
        sphere_nodes = [sphere_nodes, numNodes];
        sphere_error = [sphere_error, SphericalDesign(fileNameList(ii).name, option)];
    end
end


order_23=[];
order_25=[];
order_exp=[];
order_sphere =[];

for ii=1:maxnodes
    m=double(ii);
    order_exp=[order_exp,exp(-double(m)^(1/dim))];
    order_23=[order_23,(1/double(m))^((dim+3)/(2*dim) )];
    order_25=[order_25,(1/double(m))^((dim+5)/(2*dim))];
    order_sphere=[order_sphere,(1/double(m))^(3/4)];
end
