function [Rb,Rd,Ab,Vd,rl,rr,d] = EMB_data2geometry(data,fusionTime)
% EMB_GEOMETRY calculates geometrical measures from DIB length and droplet
% area measurements. in microns
% Rb: bilayer radius
% Rd: droplet radius
% Ab: bilayer area
% Vd: droplet volume
% rl: Rb/Rl bilayer radius normalized by radius of left droplet
% rr: Rb/Rr bilayer radius normalized by radius of right droplet
% d: distance between droplets

% exclude samples starting with nan between frame 1 and fusion time
a={data.r}; a=cat(3, a{:});
excluder=ones(size(a,3),1);
for n=1:length(data)
    if ~(sum(sum(isnan(a(1:fusionTime,:,n))))==0)
        excluder(n,1)=0;
    end
end
data=data(excluder==1);

%% compute geometry
scale=1.3; %micron/px

Rb=cat(3,data.DIB)/2*scale;
Rd=(cat(3,data.area)/pi).^0.5*scale;
        
Ab=Rb.^2*pi;
Vd=4/3*pi*Rd.^3;        

rl=Rb./Rd(:,1:5,:);
rr=Rb./Rd(:,2:6,:);

d=cat(3,data.d)*scale; d=diff(d,1,2);
end

