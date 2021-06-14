
% load('sigClassTracesAll.mat')

tt = (1:25)*.5;
CC = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840];

figure, 

plot(tt, z1(sigAll(123,2:26)), 'Color', CC(1,:), 'LineWidth', 2), hold on
P = InterX([tt; z1(sigAll(123,2:26))], [tt; 0.2*ones(1,25)]);
plot(P(1,:),P(2,:), 'Color', CC(1,:), 'LineWidth', 1)
plot(P(1,:),P(2,:),'o',  'Color', CC(1,:)),hold on

plot(tt, z1(sigAll(398,311:335))-.03, 'Color', CC(2,:), 'LineWidth', 2), hold on
P = InterX([tt; z1(sigAll(398,311:335))], [tt; 0.2*ones(1,25)]);
plot(P(1,:)-.03,P(2,:)-.03, 'Color', CC(2,:), 'LineWidth', 1)
plot(P(1,:)-.03,P(2,:)-.03, 'o', 'Color', CC(2,:)),hold on

plot(tt, z1(sigAll(398,215:239))-.06, 'Color', CC(3,:), 'LineWidth', 2), hold on
P = InterX([tt; z1(sigAll(398,215:239))], [tt; 0.2*ones(1,25)]);
plot(P(1,:)-.06,P(2,:)-.06, 'Color', CC(3,:), 'LineWidth', 1)
plot(P(1,:)-.06,P(2,:)-.06, 'o',  'Color', CC(3,:)),

plot(tt, z1(sigAll(376,291:315))-.09, 'Color', CC(4,:), 'LineWidth', 2), hold on
P = InterX([tt; z1(sigAll(376,291:315))], [tt; 0.2*ones(1,25)]);
plot(P(1,:)-.09,P(2,:)-.09, 'Color', CC(4,:), 'LineWidth', 1)
plot(P(1,:)-.09,P(2,:)-.09,'o',  'Color', CC(4,:)), 

plot([0; 0], [-0.2; 0], '-k', 'LineWidth', 2);
plot([0; 2], [-.2; -.2], '-k', 'LineWidth', 2);
t1 = text(-1, -0.2, '20% \DeltaF/F','FontSize',20, 'Color', 'k', 'FontWeight','bold'); set(t1,'Rotation',90);
t2 = text(0, -.23, '2 sec','FontSize',20, 'Color', 'k', 'FontWeight','bold');

axis off
