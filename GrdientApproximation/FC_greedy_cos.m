%greedy maximization of cos 

function output = FC_greedy_cos(maxtime,T,epsilon,n,dim,points,mu,c_var,mean,var,K,delta,partion,functype,option)

tstart=tic;
t_calc=0.0;

time=[];


if functype==1
    func = @(x,y,epsilon) exp(-epsilon*abs(x-y));
    d_function = @(x,y,epsilon) exp(-epsilon*norm(x-y,1));
elseif functype==2
    func = @(x,y,epsilon) exp(-epsilon*(x-y)^2);
    d_function = @(x,y,epsilon) exp(-epsilon*(norm(x-y))^2);
elseif functype==3
    d_function = @(x,y,epsilon) (1+norm(x-y))*exp(-norm(x-y));
elseif functype==4
    d_function = @(x,y,epsilon) (1+norm(x-y)+(norm(x-y) ^2 )/3)*exp(-norm(x-y));
elseif functype==5
    d_function = @(x,y,epsilon) exp(-epsilon*norm(x-y));
elseif functype==6
    d_function = @(x,y,epsilon) 8/3 - norm([sin(0.5*pi* (1+x(1)))*cos(pi*(1+x(2))),sin(0.5*pi*(1+x(1)))*sin(pi*(1+x(2))),cos(0.5*pi*(1+x(1)))]-[sin(0.5*pi*(1+y(1)))*cos(pi*(1+y(2))),sin(0.5*pi*(1+y(1)))*sin(pi*(1+y(2))),cos(0.5*pi*(1+y(1)))] );   
else
    'error_arises'
end

%syms x;


%%%
p=zeros(T*K,dim);%array of nodes
ind=zeros(T*K,dim);%array of indexes of nodes
ind_num=[];
c=[];%array of weights
%%%

mu_norm=(derive_error(0,[],[],dim,[],epsilon,c_var,mean,var,mu,partion,functype))^2;%|mu|_k ^2

a=find(mu==max(mu));
ind_num=[ind_num,a(1)];
c=[1];
number=a(1)-1;

%add the index and node
for jj=1:dim;
    ind(1,jj)=idivide(int64(number),int64(n^(dim-jj)) )+1;
     p(1,jj)=points(ind(1,jj));
     number= rem(number,n^(dim-jj));
end
basevalue=d_function(p(1,:),p(1,:),epsilon);%K(x,x)のこと
x_t_norm=d_function(p(1,:),p(1,:),epsilon);%|x_t|^2
mu_x_t= mu(ind_num(1)); %<mu, x_t>
x_t_array=zeros(1,length(mu));%array of values of x_t at the candidate points

for j=0:n^(dim)-1;
    indexes=func_index(j,n,dim);
    points_ar=[];
    for l=1:dim;
        points_ar=[points_ar,points(indexes(l))];
    end
    value=d_function(points_ar,p(1,:),epsilon);
    x_t_array(j+1)=x_t_array(j+1)+ value;
end
%The initial point was added
%%%%%%%%

num_points=[1];
error_ar=[derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)];

nodes=[1];
error_nodes=[derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)];
iteration_t=[1];
error_t=[derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)];

time=[time,double((toc(tstart)-t_calc))];
tic
t_calc=double(t_calc)+double(toc);

same_count=0;%the number of times in which same points are selected from candidate set
t=1;

calctime=0;
while (t < (T+1)) &&(calctime < maxtime)
    %%%
    p_d= zeros(K,dim);
    coeff_d=[];
    ind_d=[];
    ind_d_num=[];
    
    flag=0;
    
     %memorize the inner products and norms at the previous iteration of k
     %d_k= added_points - sum(coeff_d)*x_t
     d_k_x_t=0;% inner product of added points and x_t <added points,x_t>
     d_k_mu=0;%inner product of added points and mu <added points,mu>
     d_k_norm=0;%norm of added points |added points|_k ^2
     d_k_array=zeros(1,length(mu));%values of d_k at each point of candidate set
    funcval_arr=zeros(K,length(mu));
    for k=1:K
       %pp=<mu-x_t, K(,y)-x_t> : p in our paper
       %qq=<mu-x_t,d_k> : q in our paper
       %xx=<d_k , K( ,y)-x_t> : beta in our paper
       %yy=|d_k|^2  : gamma in our paper
       %aa=|K(,y)-x_t|^2 : alpha in our paper
        
       %%%
       pp= mu- x_t_array - mu_x_t + x_t_norm;
       %pp= <mu-x_t,K( ,y)-x_t>
       
       if (length(coeff_d) ==0)
           xx= d_k_array;
       else
           xx= d_k_array - d_k_x_t -sum(coeff_d)*(x_t_array -x_t_norm ) ;
       end
       %xx= <d_k, K( ,y)-x_t> =<added_points - sum(coeff_d)*x_t, K( ,y)-x_t >
       
       
       aa= basevalue-2*x_t_array+x_t_norm;
       %aa=|K(,y)-x_t|^2
       %%%
       
       %yy=|d_k|^2=| added_points - sum(coeff_d)*x_t|^2
       if (length(coeff_d) ==0)
           yy= 0;
       else
           yy= d_k_norm -2*sum(coeff_d)*d_k_x_t+sum(coeff_d)^2*x_t_norm ;
       end
       
       %%qq=<mu-x_t,d_k> =<mu-x_t, added_points -sum(coeff_d)*x_t>
       if (length(coeff_d)==0)
           qq=0;
       else
           qq=d_k_mu - d_k_x_t -sum(coeff_d)*(mu_x_t - x_t_norm);
       end
       %%%$
      
       
       angval1= qq*xx - yy*pp;%qx-py
       angval2= pp.*xx - qq*aa;%px-aq
       
       
       denomi= aa.*(angval1).^2 + 2* xx.*(angval1).* (angval2) +yy*(angval2).^2;
       numerator= (angval1.*pp + qq*angval2).^2;
       
       if (length(coeff_d)==0)
           ang_arr= pp./(aa.^(1/2));
       else
           
           ang_arr =  (numerator)./denomi;
           
           %exclude <mu -x_t, d_k+c (K(, y)-x_t)> < 0
           signarr= qq+ ((angval1)./(angval2)).*pp;
      
           sign=find( (signarr <0)|((angval1)./(angval2) < 0) ); %exclude (cp+q)<0 or (c <0)
     
           if length(sign)>0
                for ii=1:length(sign)
                    ang_arr(sign(ii))=-1*ang_arr(sign(ii));
                end
           end
       end 
       %%%%
       
       ggg=find(ang_arr==max(ang_arr));%select the point that maximizes ang_arr(cos¥theta)
       ind_v=ggg(1); 
       number=ind_v-1;
       p_d_k=zeros(1,dim);
       for jj=1:dim;
            ind_dk_num=idivide(int64(number),int64(n^(dim-jj)) )+1;
            p_d_k(jj)=points(ind_dk_num);
            number= rem(number,n^(dim-jj));
       end
       
       %%%
       p_d_copy= p_d;
       coeff_d_copy= coeff_d;
       ind_d_copy= ind_d;
       ind_d_num_copy= ind_d_num;
       %%%
       
       if (length(coeff_d)==0)
           %lambda=mu(ind_v)-x_arr1(ind_v);
           lambda=1;
       else
           lambda=  angval1(ind_v)/angval2(ind_v) ;
       end
       
       bbb =find(ind_d_num == ind_v);
       if (numel(bbb)~=0)
            coeff_d(bbb(1))=coeff_d(bbb(1))+lambda;
            add_num_d=bbb(1);
       else
            coeff_d=[coeff_d,lambda];
            ind_d_num=[ind_d_num, ind_v];
            number= ind_v -1;
            for jj=1:dim;
                ind_d(length(coeff_d),jj)=idivide(int64(number),int64(n^(dim-jj)) )+1;
                p_d(length(coeff_d),jj)=points(ind_d(length(coeff_d),jj));
                number= rem(number,n^(dim-jj));
            end
            add_num_d=length(coeff_d);
            
            %関数値の配列を作る
            for j=0:n^(dim)-1;
               indexes=func_index(j,n,dim);
               at_point=[];
               for ll=1:dim;
                   at_point=[at_point,points(indexes(ll))];
               end
               value=d_function(at_point,p_d(length(coeff_d),:),epsilon);
               funcval_arr(add_num_d,j+1)=value;
           end
       end
       
       l_coeff=length(coeff_d);
       if  l_coeff >1
           K1=zeros(l_coeff,l_coeff);
           vec1=zeros(l_coeff,1);
           for ii=1:l_coeff
               vec1(ii)=mu(ind_d_num(ii))-x_t_array(ind_d_num(ii))-mu_x_t+x_t_norm;
               for jj=1:l_coeff
                  value=d_function(p_d(ii,:), p_d(jj,:),epsilon);
                  K1(ii,jj)=value+x_t_norm-x_t_array(ind_d_num(ii))-x_t_array(ind_d_num(jj));;
               end
           end
           options = optimoptions('quadprog','OptimalityTolerance',1e-10);
           [vector,funcvalue]=quadprog(K1,-vec1,[],[],[],[],zeros(l_coeff,1),[],[],options)
           %[vector,funcvalue]=FQ_optimization(K1, K2, 10^(-8),10^(-8),10^(-12),coeff_d_copy)
           coeff_d=transpose(vector);
       end
       
       d_k_array_copy=d_k_array;   
       d_k_x_t_copy= d_k_x_t;
       d_k_mu_copy= d_k_mu;
       d_k_norm_copy=d_k_norm;
       
       d_k_array=zeros(1,length(mu));
       for ii=1:length(coeff_d)
           for jj=1:n^(dim)
               d_k_array(jj)=d_k_array(jj)+ coeff_d(ii)*funcval_arr(ii,jj);
           end
       end
       d_k_x_t= inner_product(c,coeff_d,p,p_d,epsilon,dim,functype);
       d_k_mu=0;
       for ii=1:length(coeff_d)
           d_k_mu=d_k_mu+coeff_d(ii)*mu(ind_d_num(ii));
       end
       d_k_norm= inner_product(coeff_d,coeff_d,p_d,p_d,epsilon,dim,functype);

       
       
       %compute "align"
       norm_grad=sqrt(mu_norm -2*mu_x_t + x_t_norm);%norm of mu-x_t 
       
       %norm_d_1= norm of d_k' norm_d_2= norm of d_k　| added_points - sum(coeff_d)*x_t |
       norm_d_1= sqrt(d_k_norm - 2*sum(coeff_d)*d_k_x_t +sum(coeff_d)^2*x_t_norm);
       
       %inner_1=　<mu-x_t, d_k^'> = <mu-x_t, added_points　- sum(coeff_d)* x_t>
       inner_1 =d_k_mu -d_k_x_t -sum(coeff_d)*(mu_x_t-x_t_norm);
     
       if(k>1)
           norm_d_2= sqrt(d_k_norm_copy - 2*sum(coeff_d_copy)*d_k_x_t_copy +sum(coeff_d_copy)^2*x_t_norm);
            %inner_2=　<mu-x_t, d_k^'> = <mu-x_t, added_points　- sum(coeff_d_copy)*x_t>
           inner_2= d_k_mu_copy -d_k_x_t_copy -sum(coeff_d_copy)*(mu_x_t-x_t_norm);
       end
       
       align1= (inner_1)/(norm_d_1 *norm_grad );
       
       if (k==1)
           align2=-10;%negative large value
       else
           align2= (inner_2)/(norm_d_2 *norm_grad );
       end
      
       if (align1 - align2 > delta) && (flag==0)
           flag=0;
       else
           flag=1;
           
           p_d= p_d_copy;
           coeff_d= coeff_d_copy;
           ind_d= ind_d_copy;
           ind_d_num= ind_d_num_copy;
           
           
           d_k_array=d_k_array_copy;
           
           d_k_x_t= d_k_x_t_copy;
           d_k_mu= d_k_mu_copy;
           d_k_norm=d_k_norm_copy;
           
       end
       
       if (flag==1)
            break
       end
       
           
    end
    
    
    if(k==1)
        
        stepsize = min([1, ((inner_1 / (norm_d_1)^2))* sum(coeff_d)]); %stepsize= < mu- x_t, d_t -x_t>/|x_t -d_t|^2  **sum(coeff_d) (本来規格化する)
    else
        %stepsize = min([1, (inner_2 / (norm_d_2)^2)*sum(coeff_d) ]); %stepsizeは < mu- x_t, d_t -x_t>/|x_t -d_t|^2 *sum(coeff_d)
        stepsize= min([1,(inner_2 / (norm_d_2)^2)*sum(coeff_d)]);
    end
    
    
    c = (1- stepsize)*c  ;
    lenc= length(c);
    
    
    x_t_norm=(1- stepsize)^2*x_t_norm +(stepsize/sum(coeff_d))^2 *d_k_norm+ 2*((1- stepsize)*stepsize*d_k_x_t)/sum(coeff_d);%|x_t|^2 %|(1-stepsize)x_t+stepsize*(added points)/(coeff_d)|
    mu_x_t=(1- stepsize)*mu_x_t +(stepsize*d_k_mu)/sum(coeff_d); %<mu,(1-stepsize)*x_t+stepsize*added points/sum(coeff_d)>
    x_t_array= (1-stepsize)* x_t_array + (stepsize*d_k_array)/sum(coeff_d);%update x_t_array
    
    
    coeff_d= coeff_d/(sum(coeff_d));
    
     counts=0;   
     for i= lenc +1 : lenc+length(coeff_d)
        kkk=ind_d_num(i-lenc);
        bbb =find(ind_num == kkk );
        if (numel(bbb)~=0)
            c(bbb(1))= c(bbb(1))+ stepsize* coeff_d(i - lenc );
            same_count=same_count+1;
            counts=counts+1;
        else
            c =[c,stepsize* coeff_d(i - lenc )];
            ind_num = [ ind_num, ind_d_num(i-lenc) ];
            for k=1:dim
               p(i-counts,k)=p_d(i-lenc,k);
               ind(i-counts,k)=ind_d(i-lenc,k);
            end
        end
     end    
     

    %{
    num_points=[num_points,length(c)+same_count];
    error_ar = [error_ar,sqrt(mu_norm -2*mu_x_t + x_t_norm)];
    %}
     num_points=[num_points,length(c)+sum(same_count)];
     error_ar = [error_ar,sqrt(mu_norm -2*mu_x_t + x_t_norm)];
    
     
     t=t+K;
     nodes=[nodes,length(c)];
     error_nodes=[error_nodes,sqrt(mu_norm -2*mu_x_t + x_t_norm)];
    
     iteration_t=[iteration_t,t];
     error_t=[error_t,sqrt(mu_norm -2*mu_x_t + x_t_norm)];
     
      time=[time,double((toc(tstart)-t_calc))];
      tic
      t_calc=double(t_calc)+double(toc);
      
      calctime=double((toc(tstart)-t_calc));
      
end


if option==1
    output=[nodes,error_nodes];
elseif option==2
    output=[iteration_t,error_t];
elseif option==3
    output=[time,error_t];
elseif option==4
    output=[num_points, error_ar];
end

end