function value= compute_integral_qmc(number,dim,m)
    %m=10^9
    %等間隔点で数値積分構成
    if number==1
        func=@(x) (1+norm(2*[x(1)-x(3),x(2)-x(4)]))*exp(-norm(2*[x(1)-x(3),x(2)-x(4)]));
    elseif number==2
        func=@(x) (1+norm(2*[x(1)-x(4),x(2)-x(5),x(3)-x(6)]))*exp(-norm(2*[x(1)-x(4),x(2)-x(5),x(3)-x(6)]));
    elseif number==3
        func=@(x) (1+norm(2*[x(1)-x(3),x(2)-x(4)])+(norm(2*[x(1)-x(3),x(2)-x(4)]) ^2 )/3)*exp(-norm(2*[x(1)-x(3),x(2)-x(4)]));
    elseif number==4
        func=@(x) (1+norm(2*[x(1)-x(4),x(2)-x(5),x(3)-x(6)])+(norm(2*[x(1)-x(4),x(2)-x(5),x(3)-x(6)]) ^2 )/3)*exp(-norm(2*[x(1)-x(4),x(2)-x(5),x(3)-x(6)]));
    end    
    array=sobolset(2*dim,'Leap',m);
    sizearray=size(array)
    %{
    index1=linspace(1,dim,dim);
    index2=linspace(dim+1,2*dim,dim);
    
    array1=array(:,index1);
    array2=array(:,index2);
    %}
    value=0;
    for ii=1:sizearray(1)
        value=value+func(array(ii,:));
    end
    value=value/sizearray(1);
end        