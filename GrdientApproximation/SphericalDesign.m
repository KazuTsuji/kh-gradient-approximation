%kernel herding with equal weights

function output = SphericalDesign(fileName,option)


d_function = @(x,y) 8/3 - norm(x-y) ;


% Load spherical nodes
nodes = load("SphereDesign/"+fileName)
num_node = size(nodes,1)
weight = 1/num_node
%%%%%%%%%%%%%%%%%%%%%%%%%

%mu_normの計算
mu_norm = 4/3
% mu_(x_i)の計算
mu_int =4/3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu_val_total=0;
quad_total=0;
    
%二次形式部分の計算

for ii=1: num_node
    mu_val_total=mu_val_total+ weight * mu_int;
    for jj=1:num_node
        quad_total=quad_total+ weight * weight* d_function(nodes(ii,:),nodes(jj,:));
    end
end

error_value=mu_norm - 2*  mu_val_total + quad_total;
error_value=sqrt(error_value);



if option==4
    output= error_value;
end

end

