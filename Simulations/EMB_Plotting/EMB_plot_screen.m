function EMB_plot_screen(t, y_R);

global x y nscreen subplot1 subplot2 subplot3 Color

[x,y]=size(y_R);
nscreen = 5;

Color{1}=[0.00392156885936856 0.294117659330368 0.266666680574417];
Color{2}=[0.0196078438311815 0.549019634723663 0.490196079015732];
Color{3}=[0.266666680574417 0.749019622802734 0.69803923368454];
Color{4}=[0.541176497936249 0.815686285495758 0.776470601558685];
Color{5}=[0.725490212440491 0.882352948188782 0.847058832645416];

% IPTG plot, variable 1
subplot1 = subplot(1,3,1, 'FontSize',14); hold on; axis square;

box(subplot1,'on');
xlabel('Time (min)');xlim(subplot1,[0 500]);ylabel('IPTG (M)');

for i=1:nscreen;
    plot(t, [y_R{i,1}(:,6), y_R{i,1}(:,11), y_R{i,1}(:,16), y_R{i,1}(:,21), y_R{i,1}(:,26)], 'LineWidth', 1, 'Color', Color{nscreen-i+1});
end

% active pLacO, equivalent to variable 3
subplot2 = subplot(1,3,2, 'FontSize',14); hold on; axis square;

box(subplot2,'on');
xlabel('Time (min)');xlim(subplot2,[0 500]);ylabel('Ratio of active promoter');

for i=1:nscreen;
    plot(t, [y_R{i,1}(:,8), y_R{i,1}(:,13), y_R{i,1}(:,18), y_R{i,1}(:,23), y_R{i,1}(:,28)], 'LineWidth', 1, 'Color', Color{nscreen-i+1});
end

% reporter, variable 4
subplot3 = subplot(1,3,3, 'FontSize',14); hold on; axis square;

box(subplot3,'on');
xlabel('Time (min)');xlim(subplot3,[0 500]);ylabel('Reporter protein (a.u.)');

for i=1:nscreen;
    plot(t, [y_R{i,1}(:,10), y_R{i,1}(:,15), y_R{i,1}(:,20), y_R{i,1}(:,25), y_R{i,1}(:,30)], 'LineWidth', 1, 'Color', Color{nscreen-i+1});
end

end
