%kernel herding with linesearch

function output = Fully_corrective(maxtime,m,epsilon,n,dim,points,c_var,mu,var,mean,partion,functype,option)
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

syms x;

error_ar=[];
number_of_points=[];


p=zeros(m,dim);%array of nodes
ind=zeros(m,dim);%array of indexes of nodes
ind_num=[];
c=[];%array of weights
i=1;
same_count=0;

mu_norm=(derive_error(0,[],[],dim,[],epsilon,c_var,mean,var,mu,partion,functype))^2;%|mu|_k ^2

%memorize the inner products and norms at the previous iteration of t
x_t_norm=0;
mu_x_t =0;
x_t_array=zeros(1,length(mu));

nodes=[];
error_nodes=[];
iteration_t=[];
error_t=[];
funcval_arr=zeros(m,dim);
calctime=0;
while ( i < (m+1))&&(calctime < maxtime)
    if i==1;
        a=find(mu==max(mu));
        ind_num=[ind_num,a(1)];
        number=a(1)-1;
        
        for jj=1:dim;
            ind(i,jj)=idivide(int64(number),int64(n^(dim-jj)) )+1;
            p(i,jj)=points(ind(i,jj));
            number= rem(number,n^(dim-jj));
        end
       
    else
        b=x_t_array- mu;
        
        a=find(b==min(b));
        ind_num_i = a(1);
        p_i=zeros(1,dim);
        ind_i=zeros(1,dim);
        number=a(1)-1;
        
        for jj=1:dim;
            ind_i(jj)=idivide(int64(number),int64(n^(dim-jj)) )+1;
            p_i(jj)=points(ind_i(jj));
            number= rem(number,n^(dim-jj));
        end
    end
    %compute the step size
    if i==1;
        alpha=1;
        c=[c,alpha];
        x_t_norm= d_function(p(1,:),p(1,:),epsilon); %norm |x_t|
        mu_x_t= mu(a(1));%<mu,x_t>
        
        add_number=length(c);
        for j=0:n^(dim)-1;
            indexes=func_index(j,n,dim);
            points_ar=[];
            for l=1:dim;
                points_ar=[points_ar,points(indexes(l))];
            end
            value=d_function(points_ar,p(add_number,:),epsilon);
            x_t_array(j+1)=value;
            funcval_arr(1,j+1)=value;
        end
        
        nodes=[nodes,length(c)];
        error_nodes=[error_nodes,derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)];
    else
        bbb =find(ind_num == ind_num_i);
        if (numel(bbb)~=0)
            same_count=same_count+1;
            add_number=bbb(1);
        else
            ind_num= [ind_num, ind_num_i];
            for jj=1:dim
                ind(length(c)+1,jj)=ind_i(jj);
                p(length(c)+1,jj)=points(ind_i(jj));
                add_number=length(c)+1;
            end
            
           for j=0:n^(dim)-1;
               indexes=func_index(j,n,dim);
               points_ar=[];
               for ll=1:dim;
                   points_ar=[points_ar,points(indexes(ll))];
               end
               funcval_arr(length(c)+1,j+1)=d_function(points_ar,p(add_number,:),epsilon);
            end
        end
        %%%%%%%%%%%
        %最適化
        %%%%%%%%%%%%%
        H=zeros(length(c)+1,length(c)+1);
        ff=zeros(length(c)+1,1);
        for ii=1:length(c)+1
            ff(ii)=mu(ind_num(ii));
            for jj=1:length(c)+1
                H(ii,jj)=d_function(p(ii,:),p(jj,:),epsilon);
            end
        end
        A_eq=1+zeros(1,length(c)+1);
        b_eq=1;
        options = optimoptions('quadprog','OptimalityTolerance',1e-10);
        [vector,funcvalue]=quadprog(H,-ff,[],[],A_eq,b_eq,zeros(length(c)+1,1),[],[],options);
        
        c=vector;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x_t_norm=inner_product(c,c,p,p,epsilon,dim,functype);
        mu_x_t=0;
        for ii=1:length(c)
            mu_x_t=mu_x_t+c(ii)*mu(ind_num(ii));
        end
        
        x_t_array=zeros(1,length(mu));
        for ii=1:length(c)
           for jj=1:n^(dim)
               x_t_array(jj)=x_t_array(jj)+ c(ii)*funcval_arr(ii,jj);
           end
        end
    end
    
    nodes=[nodes,length(c)];
    error_nodes=[error_nodes,derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)];
    iteration_t=[iteration_t,i];
    error_t=[error_t, sqrt(mu_norm -2*mu_x_t + x_t_norm)];
    
    time=[time,double((toc(tstart)-t_calc))];
    tic
    t_calc=double(t_calc)+double(toc);
    
    error_ar=[error_ar,derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)];
    number_of_points=[number_of_points,length(c)+sum(same_count)];
    i=i+1;
    calctime=double((toc(tstart)-t_calc));
    
end
if option==1
    output=[nodes,error_nodes];
elseif option==2
    output=[iteration_t,error_t];
elseif option==3
    output=[time,error_t];
elseif option==4
    output=[number_of_points, error_ar];
end



end