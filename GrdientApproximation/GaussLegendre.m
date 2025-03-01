%kernel herding with equal weights

function output = GaussLegendre(m,epsilon,dim,functype,option)


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
else
    'error_arises'
    return
end

density_const= (1/ (erf(1)*sqrt(pi)) )^dim;

% Calculate Gauss Legendre nodes
[gl_node, gl_weight] = Gaulegwt(-1,1,m)

nodes=zeros(m^dim,dim);%array of nodes
quad_coeff=zeros(m^dim, 1);% array of weights

for ii=1:m^dim
    indexes = func_index(ii-1, m, dim)
    each_node = []
  
    factor =0.5^dim

    for jj=1:dim
        each_node =[each_node, gl_node(indexes(jj))]
        factor = factor* gl_weight(indexes(jj))
    end
    nodes(ii,:) = each_node
    quad_coeff(ii,:) = factor
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%




%mu_normの計算
if functype==1;
   val=1;
   for dd=1:dim
        val=val * integral2(@(x,y) exp(-epsilon*abs(x-y)).* exp(- x.^2).*exp(-y.^2),-1,1,-1,1);
    end
    mu_norm=val*density_const^2;
elseif functype==2;
    val=1;
    for dd=1:dim
        val=val *integral2(@(x,y) exp(-epsilon*(x-y).^2).*exp(- x.^2).*exp(- y.^2),-1,1,-1,1);
    end
    mu_norm=val*density_const^2
elseif functype==3
    if dim==2
        mu_norm=load('IntegralMatFiles/matern32_2d_val.mat').matern32_2dval;
    elseif dim==3
        mu_norm=load('IntegralMatFiles/matern32_3d_val.mat').matern32_3dval;
    end
elseif functype==4
    if dim==2
        mu_norm=load('IntegralMatFiles/matern52_2d_integvalue.mat').value2d;
    elseif dim==3
        mu_norm=load('IntegralMatFiles/matern52_3d.mat').value;
    end
elseif functype==5
    if dim==2
        mu_norm=load('IntegralMatFiles/exp_2d.mat').value;
    elseif dim==3
        mu_norm=load('IntegralMatFiles/exp_3d.mat').value;
    end
end




% mu_(x_i)の計算

mu_ar = zeros(m^dim, 1);

for ii=1:m^dim
    if functype == 1
        val = 1;
        for kk = 1:dim
            f = @(x) exp(-x.^2) .* exp(-epsilon * abs(x - nodes(ii,kk)));
            val = val * integral(f, -1, 1);
        end
        mu_ar(ii) = val * density_const;

    elseif functype == 2
        val = 1;
        for kk = 1:dim
            f = @(x) exp(-x.^2) .* exp(-epsilon * (x - nodes(ii,kk)).^2);
            val = val * integral(f, -1, 1);
        end
        mu_ar(ii) = val * density_const;
    
    elseif functype==3;
        p_ar=nodes(ii,:);
        if dim==2
            ff=@(x,y)(1+sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2 ) ) .*exp(-sqrt( (x-p_ar(1)).^2+(y-p_ar(2)).^2 ));
            mu_ar(ii)=(0.5)^dim *integral2(ff,-1,1,-1,1);
        elseif dim==3
            ff=@(x,y,z)(1+sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2 ) ) .*exp(-sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2));
            mu_ar(ii)=(0.5)^dim *integral3(ff,-1,1,-1,1,-1,1);
        end
    
    elseif functype==4;
        p_ar=nodes(ii,:);
        if dim==2
            ff=@(x,y)(1+sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2 )+((x-p_ar(1)).^2+(y-p_ar(2)).^2 )/3 ) .*exp(-sqrt( (x-p_ar(1)).^2+(y-p_ar(2)).^2 ));
            mu_ar(ii)=(0.5)^dim *integral2(ff,-1,1,-1,1);
        elseif dim==3
            ff=@(x,y,z)(1+sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2 )+((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2 )/3 ) .*exp(-sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2));
            mu_ar(ii)=(0.5)^dim *integral3(ff,-1,1,-1,1,-1,1);
        end
        
    elseif functype==5;
        p_ar=nodes(ii,:);
        if dim==2
            ff=@(x,y) exp(-epsilon*sqrt( (x-p_ar(1)).^2+(y-p_ar(2)).^2 ));
            mu_ar(ii)=(0.5)^dim *integral2(ff,-1,1,-1,1);
        elseif dim==3
            ff=@(x,y,z) exp(-epsilon*sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2));
            mu_ar(ii)=(0.5)^dim *integral3(ff,-1,1,-1,1,-1,1);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu_val_total=0;
quad_total=0;
    
%二次形式部分の計算

for ii=1: m^dim
    mu_val_total=mu_val_total+ quad_coeff(ii)* mu_ar(ii);
    for jj=1:m^dim
        quad_total=quad_total+ quad_coeff(ii) * quad_coeff(jj)* d_function(nodes(ii,:),nodes(jj,:),epsilon);
    end
end

error_value=mu_norm - 2*  mu_val_total + quad_total;
error_value=sqrt(error_value);



if option==4
    output= error_value;
end

end

