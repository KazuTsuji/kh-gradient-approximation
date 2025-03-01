%compair the four methods on the convergence speed of the worst case errors for the number of nodes when dim=2

%dim=2 Sphere case

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

option=3;

maxtime=60;
max_time_copy=maxtime;
maxnodes=30000;
max_nodes_copy=maxnodes;
time_num=10;

%Gaussian dim =3 
%{
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
option=3;
maxtime=35;
max_time_copy=maxtime;
maxnodes=30000;
max_nodes_copy=maxnodes;
time_num=10;
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


stride=10;
for tt=1:time_num
    output_linesearch= linesearch(maxtime,maxnodes,epsilon,n,dim,points,c_var,mu,var,meanarr,partion,functype,option);
    len_ls=  int64(length(output_linesearch)/2);
    line_search =[];
    ls_nodes=[];

    fact=0;
    for i=1: len_ls
        if fact*stride < i
            line_search = [line_search, output_linesearch(len_ls +i)];
            ls_nodes=[ls_nodes,output_linesearch(i)];
            fact=fact+1;
        end
    end
    
    if tt==1
        line_search_ave=line_search;
        ls_nodes_ave=ls_nodes;
    else
        if(length(line_search) > length(line_search_ave))
            line_search = line_search(1:length(line_search_ave));
            ls_nodes = ls_nodes(1:length(line_search_ave));
        elseif(length(line_search) < length(line_search_ave))
            line_search_ave = line_search_ave(1:length(line_search));
            ls_nodes_ave = ls_nodes_ave(1:length(line_search));
        end
        line_search_ave=line_search_ave+line_search;
        ls_nodes_ave=ls_nodes_ave+ls_nodes;
    end
end
line_search=line_search_ave/time_num;
ls_nodes=ls_nodes_ave/time_num;
maxtime=max_time_copy;
maxnodes=max_nodes_copy;


stride=10;

for tt=1:time_num
    output_eqweight= eqweight_herding(maxtime,maxnodes,epsilon,n,dim,points,mu,c_var,meanarr,var,partion,functype,option);
    len_eq=  int64(length(output_eqweight)/2);
    eqweight =[];
    eq_nodes=[];

    fact=0;
    for i=1: len_eq
        if fact*stride < i
            eqweight = [eqweight, output_eqweight(len_eq +i)];
            eq_nodes=[eq_nodes,output_eqweight(i)];
            fact=fact+1;
        end
    end
    
    if tt==1
        eqweight_ave=eqweight;
        eq_nodes_ave=eq_nodes;
    else
        if(length(eqweight) > length(eqweight_ave))
            eqweight = eqweight(1:length(eqweight_ave));
            eq_nodes = eq_nodes(1:length(eqweight_ave));
        elseif(length(eqweight) < length(eqweight_ave))
            eqweight_ave = eqweight_ave(1:length(eqweight));
            eq_nodes_ave = eq_nodes_ave(1:length(eqweight));
        end
        eqweight_ave=eqweight_ave+eqweight;
        eq_nodes_ave=eq_nodes_ave+eq_nodes;
    end
end
eqweight=eqweight_ave/time_num;
eq_nodes=eq_nodes_ave/time_num;
maxtime=max_time_copy;
maxnodes=max_nodes_copy;



stride=10;
for tt=1:time_num
    output_pmp= positive_matching_pursuit(maxtime,maxnodes,epsilon,n,dim,points,mu,c_var,meanarr,var,K,delta,partion,functype,option);
    len_pmp=  int64(length(output_pmp)/2);
    pmp =[];
    pmp_nodes=[];

    fact=0;
    for i=1: len_pmp
        if fact*stride < i
            pmp = [pmp, output_pmp(len_pmp +i)];
            pmp_nodes=[pmp_nodes,output_pmp(i)];
            fact=fact+1;
        end
    end
    if tt==1
        pmp_ave=pmp;
        pmp_nodes_ave=pmp_nodes;
    else
        if(length(pmp) > length(pmp_ave))
            pmp = pmp(1:length(pmp_ave));
            pmp_nodes = pmp_nodes(1:length(pmp_ave));
        elseif(length(pmp) < length(pmp_ave))
            pmp_ave = pmp_ave(1:length(pmp));
            pmp_nodes_ave = pmp_nodes_ave(1:length(pmp));
        end
        pmp_ave = pmp_ave+pmp;
        pmp_nodes_ave=pmp_nodes_ave+pmp_nodes;
    end
end
pmp=pmp_ave/time_num;
pmp_nodes=pmp_nodes_ave/time_num;
maxtime=max_time_copy;
maxnodes=max_nodes_copy;


stride=10;

for tt=1:time_num
    output_fc_pmp= FC_PMP(maxtime,maxnodes,epsilon,n,dim,points,mu,c_var,meanarr,var,K,delta,partion,functype,option);
    len_fc_pmp=  int64(length(output_fc_pmp)/2);
    fc_pmp =[];
    fc_pmp_nodes=[];

    fact=0;
    for i=1: len_fc_pmp
        if fact*stride < i
            fc_pmp = [fc_pmp, output_fc_pmp(len_fc_pmp +i)];
            fc_pmp_nodes=[fc_pmp_nodes,output_fc_pmp(i)];
            fact=fact+1;
        end
    end
    
    if tt==1
        fc_pmp_ave=fc_pmp;
        fc_pmp_nodes_ave=fc_pmp_nodes;
    else
        if(length(fc_pmp) > length(fc_pmp_ave))
            fc_pmp = fc_pmp(1:length(fc_pmp_ave));
            fc_pmp_nodes = fc_pmp_nodes(1:length(fc_pmp_ave));
        elseif(length(fc_pmp) < length(fc_pmp_ave))
            fc_pmp_ave = fc_pmp_ave(1:length(fc_pmp) );
            fc_pmp_nodes_ave = fc_pmp_nodes_ave(1:length(fc_pmp));
        end
        fc_pmp_ave=fc_pmp_ave+fc_pmp;
        fc_pmp_nodes_ave=fc_pmp_nodes_ave+fc_pmp_nodes;
    end
end
fc_pmp=fc_pmp_ave/time_num;
fc_pmp_nodes=fc_pmp_nodes_ave/time_num;
maxtime=max_time_copy;
maxnodes=max_nodes_copy;


stride=10;
for tt=1:time_num
    output_gcos= greedy_cos(maxtime,maxnodes,epsilon,n,dim,points,mu,c_var,meanarr,var,K,delta,partion,functype,option);
    len_gcos=  int64(length(output_gcos)/2);
    gcos =[];
    gcos_nodes=[];

    fact=0;
    for i=1: len_gcos
        if fact*stride < i
            gcos = [gcos, output_gcos(len_gcos +i)];
            gcos_nodes=[gcos_nodes,output_gcos(i)];
            fact=fact+1;
        end
    end
    if tt==1
        gcos_ave=gcos;
        gcos_nodes_ave=gcos_nodes;
    else
        if(length(gcos) > length(gcos_ave))
            gcos = gcos(1:length(gcos_ave));
            gcos_nodes = gcos_nodes(1:length(gcos_ave));
        elseif(length(gcos) < length(gcos_ave))
            gcos_ave = gcos_ave(1:length(gcos));
            gcos_nodes_ave = gcos_nodes_ave(1:length(gcos));
        end
        gcos_ave=gcos_ave+gcos;
        gcos_nodes_ave=gcos_nodes_ave+gcos_nodes;
    end
end
gcos=gcos_ave/time_num;
gcos_nodes=gcos_nodes_ave/time_num;
maxtime=max_time_copy;
maxnodes=max_nodes_copy;


stride=10;

for tt=1:time_num
    output_fc_gcos= FC_greedy_cos(maxtime,maxnodes,epsilon,n,dim,points,mu,c_var,meanarr,var,K,delta,partion,functype,option);
    len_fc_gcos=  int64(length(output_fc_gcos)/2);
    fc_gcos =[];
    fc_gcos_nodes=[];

    fact=0;
    for i=1: len_fc_gcos
        if fact*stride < i
            fc_gcos = [fc_gcos, output_fc_gcos(len_fc_gcos +i)];
            fc_gcos_nodes=[fc_gcos_nodes,output_fc_gcos(i)];
            fact=fact+1;
        end
    end
    if tt==1
        fc_gcos_ave=fc_gcos;
        fc_gcos_nodes_ave=fc_gcos_nodes;
    else
        if(length(fc_gcos) > length(fc_gcos_ave))
            fc_gcos = fc_gcos(1:length(fc_gcos_ave));
            fc_gcos_nodes = fc_gcos_nodes(1:length(fc_gcos_ave));
        elseif(length(fc_gcos) < length(fc_gcos_ave))
            fc_gcos_ave = fc_gcos_ave(1:length(fc_gcos));
            fc_gcos_nodes_ave = fc_gcos_nodes_ave(1:length(fc_gcos));
        end
        fc_gcos_ave=fc_gcos_ave+fc_gcos;
        fc_gcos_nodes_ave=fc_gcos_nodes_ave+fc_gcos_nodes;
    end
end
fc_gcos=fc_gcos_ave/time_num;
fc_gcos_nodes=fc_gcos_nodes_ave/time_num;
maxtime=max_time_copy;
maxnodes=max_nodes_copy;


stride=10;
for tt=1:time_num
    output_FC= Fully_corrective(maxtime,maxnodes,epsilon,n,dim,points,c_var,mu,var,meanarr,partion,functype,option);
    len_FC=  int64(length(output_FC)/2);
    FC =[];
    FC_nodes=[];

    fact=0;
    for i=1: len_FC
        if fact*stride< i
            FC = [FC, output_FC(len_FC +i)];
            FC_nodes=[FC_nodes,output_FC(i)];
            fact=fact+1;
        end
    end
    if tt==1
        FC_ave=FC;
        FC_nodes_ave=FC_nodes;
    else
        if(length(FC) > length(FC_ave))
            FC = FC(1:length(FC_ave));
            FC_nodes = FC_nodes(1:length(FC_ave));
        elseif(length(FC) < length(FC_ave))
            FC_ave = FC_ave(1:length(FC));
            FC_nodes_ave = FC_nodes_ave(1:length(FC));
        end
        FC_ave=FC_ave+FC;
        FC_nodes_ave=FC_nodes_ave+FC_nodes;
    end
end
FC=FC_ave/time_num;
FC_nodes=FC_nodes_ave/time_num;
maxtime=max_time_copy;
maxnodes=max_nodes_copy;


stride=1;
for tt=1:time_num
    output_SBQ=  SBQ(maxtime,maxnodes,epsilon,n,dim,points,mu,c_var,meanarr,var,partion,functype,option);
    len_SBQ=  int64(length(output_SBQ)/2);
    SBQ_err =[];
    SBQ_nodes=[];

    fact=0;
    for i=1: len_SBQ
        if fact*stride < i
            SBQ_err = [SBQ_err, output_SBQ(len_SBQ +i)];
            SBQ_nodes=[SBQ_nodes,output_SBQ(i)];
            fact=fact+1;
        end
    end
    if tt==1
        SBQ_ave=SBQ_err;
        SBQ_nodes_ave=SBQ_nodes;
    else
         if(length(SBQ_err) > length(SBQ_ave))
            SBQ_err = SBQ_err(1:length(SBQ_ave));
            SBQ_nodes = SBQ_nodes(1:length(SBQ_ave));
        elseif(length(SBQ_err) < length(SBQ_ave))
            SBQ_ave = SBQ_ave(1:length(SBQ_err));
            SBQ_nodes_ave = SBQ_nodes_ave(1:length(SBQ_err));
        end
        SBQ_ave=SBQ_ave+SBQ_err;
        SBQ_nodes_ave=SBQ_nodes_ave+SBQ_nodes;
    end
end
SBQ_err=SBQ_ave/time_num;
SBQ_nodes=SBQ_nodes_ave/time_num;
maxtime=max_time_copy;
maxnodes=max_nodes_copy;
