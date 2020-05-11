function EMB_plotting_circuit(t, y_R)

global x y nspecies ncomp subplot1 subplot2 subplot3 Color

[x,y]=size(y_R);
nspecies = 8;
ncomp = 6;

Color{1}=[0.00392156885936856 0.294117659330368 0.266666680574417];
Color{2}=[0.0196078438311815 0.549019634723663 0.490196079015732];
Color{3}=[0.266666680574417 0.749019622802734 0.69803923368454];
Color{4}=[0.541176497936249 0.815686285495758 0.776470601558685];
Color{5}=[0.725490212440491 0.882352948188782 0.847058832645416];

% IPTG plot, variable 1
subplot1 = subplot(1,3,1, 'FontSize',14); hold on; axis square;

box(subplot1,'on');
xlabel('Time (min)');xlim(subplot1,[0 500]);ylabel('IPTG (M)');

for i=2:ncomp;
    patch(t([1:end, end:-1:1]), [y_R{nspecies.*i-7,2}(1:end,1)+y_R{nspecies.*i-7,3}(1:end,1); y_R{nspecies.*i-7,2}(end:-1:1,1)-y_R{nspecies.*i-7,3}(end:-1:1,1)],Color{i-1},'FaceAlpha',0.3,'LineWidth',1,'LineStyle','none');
    plot(t, y_R{nspecies.*i-7,2}(:,1), 'LineWidth', 1, 'Color', Color{i-1});
end

% YFP, variable 8
subplot2 = subplot(1,3,2, 'FontSize',14); hold on; axis square;

box(subplot2,'on');
xlabel('Time (min)');xlim(subplot2,[0 500]);ylabel('YFP (M)');

for i=2:ncomp;
    patch(t([1:end, end:-1:1]), [y_R{nspecies.*i,2}(1:end,1)+y_R{nspecies.*i,3}(1:end,1); y_R{nspecies.*i,2}(end:-1:1,1)-y_R{nspecies.*i,3}(end:-1:1,1)],Color{i-1},'FaceAlpha',0.3,'LineWidth',1,'LineStyle','none');
    plot(t, y_R{nspecies.*i,2}(:,1), 'LineWidth', 1, 'Color', Color{i-1});
end

% RFP, variable 5
subplot3 = subplot(1,3,3, 'FontSize',14); hold on; axis square;

box(subplot3,'on');
xlabel('Time (min)');xlim(subplot3,[0 500]);ylabel('RFP (M)');

for i=2:ncomp;
    patch(t([1:end, end:-1:1]), [y_R{nspecies.*i-3,2}(1:end,1)+y_R{nspecies.*i-3,3}(1:end,1); y_R{nspecies.*i-3,2}(end:-1:1,1)-y_R{nspecies.*i-3,3}(end:-1:1,1)],Color{i-1}, 'FaceAlpha',0.3,'LineWidth',1,'LineStyle','none');
    plot(t, y_R{nspecies.*i-3,2}(:,1), 'LineWidth', 1, 'Color', Color{i-1});
end
    