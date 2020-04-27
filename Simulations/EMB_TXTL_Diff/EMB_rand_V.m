function [t,y_R, r]=EMB_rand_V(a,b,c)

tic 

global mu sigmarel sigma ncompartment runs x1 y1 x2 y2 outcome
mu = a; % mean of distribution
sigmarel = b; % CV of distribution
sigma = a*b; % stdev of distribution
ncompartment = 6; % number of compartments
runs = c.*ncompartment; % number of runs with random values

r = mu + sigma.*randn(runs,1); % random values within normal distribution

for i=1:runs;
    if r(i) <= 0;
        while r(i) <=0;
            r(i) = mu + sigma.*randn(1,1);
        end
    end
end

for i=1:c;
    [t,outcome]=EMB_Diff_circuit_decay_V(r(ncompartment*i-ncompartment+1:ncompartment*i));
    [x1,y1]=size(outcome);
    for j=1:y1;
        y_R{j,1}(:,i)=outcome(:,j);
    end
end

[x2,y2]=size(y_R{1,1});

for j=1:y1;
    for i=1:x2;
        y_R{j,2}(i,1)=mean(y_R{j,1}(i,1:y2));
        y_R{j,3}(i,1)=std(y_R{j,1}(i,1:y2));
    end
end
    
t=t;
y_R=y_R;
r=r;

toc

end
