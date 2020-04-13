function ddxpdf=EMB_dfdx(f)
% EMB_DFDX computes symmetric discrete derivative of f across columns
    [~,Lx]=size(f);
    leftEdge=f(:,2)-f(:,1);
    rightEdge=f(:,Lx)-f(:,Lx-1);
    center = f(:,3:Lx)/2 - f(:,1:Lx-2)/2; %central derivative
    ddxpdf =[leftEdge,center,rightEdge];
end

