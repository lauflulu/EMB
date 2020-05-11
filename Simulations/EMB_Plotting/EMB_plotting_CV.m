function EMB_plotting_CV(t, y_R_1, y_R_2)

global x y nspecies subplot1 subplot2 subplot3 Color

[x,y]=size(y_R_1);
nspecies = 5;

Color{1}=[0.00392156885936856 0.294117659330368 0.266666680574417];
Color{2}=[0.0196078438311815 0.549019634723663 0.490196079015732];
Color{3}=[0.266666680574417 0.749019622802734 0.69803923368454];
Color{4}=[0.541176497936249 0.815686285495758 0.776470601558685];
Color{5}=[0.725490212440491 0.882352948188782 0.847058832645416];

% IPTG plot, variable 1
subplot1 = subplot(2,3,1, 'FontSize',14); hold on; axis square;
box(subplot1,'on');
xlabel('Time (min)');xlim(subplot1,[0 500]);ylabel('IPTG (M)');

for i=2:x./6+1;
    plot(t, y_R_1{nspecies.*i-4,4}(:,1), 'LineWidth', 1, 'Color', Color{i-1});
end

subplot4 = subplot(2,3,4, 'FontSize',14); hold on; axis square;
box(subplot4,'on');
xlabel('Time (min)');xlim(subplot4,[0 500]);ylabel('IPTG (M)');

for i=2:x./6+1;
    plot(t, y_R_2{nspecies.*i-4,4}(:,1), 'LineWidth', 1, 'Color', Color{i-1});
end

% active pLacO, ie inverse free LacI plot, variable 2
subplot2 = subplot(2,3,2, 'FontSize',14); hold on; axis square;

box(subplot2,'on');
xlabel('Time (min)');xlim(subplot2,[0 500]);ylabel('Ratio of active promoter');

for i=2:x./6+1;
    plot(t, y_R_1{nspecies.*i-2,4}(:,1), 'LineWidth', 1, 'Color', Color{i-1});
end

% active pLacO, ie inverse free LacI plot, variable 2
subplot5 = subplot(2,3,5, 'FontSize',14); hold on; axis square;

box(subplot5,'on');
xlabel('Time (min)');xlim(subplot5,[0 500]);ylabel('Ratio of active promoter');

for i=2:x./6+1;
    plot(t, y_R_2{nspecies.*i-2,4}(:,1), 'LineWidth', 1, 'Color', Color{i-1});
end

% YFP plot, variable 5
subplot3 = subplot(2,3,3, 'FontSize',14); hold on; axis square;

box(subplot3,'on');
xlabel('Time (min)');xlim(subplot3,[0 500]);ylabel('Reporter protein (a.u.)');

for i=2:x./6+1;
    plot(t, y_R_1{nspecies.*i,3}(:,1), 'LineWidth', 1, 'Color', Color{i-1});
end

subplot6 = subplot(2,3,6, 'FontSize',14); hold on; axis square;

box(subplot6,'on');
xlabel('Time (min)');xlim(subplot6,[0 500]);ylabel('Reporter protein (a.u.)');

for i=2:x./6+1;
    plot(t, y_R_2{nspecies.*i,3}(:,1), 'LineWidth', 1, 'Color', Color{i-1});
end
    