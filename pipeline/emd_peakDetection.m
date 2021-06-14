
clear n;

F = F_tab(:,55)';

[imf,ort,nbits] = emd(F);

F_emd = F-imf(end,:)-imf(1,:)-imf(2,:);
n = findpeaks(F_emd);

% [raw, col] = find(n>std(n));
figure,plot([1:1:941], F_emd); hold on
for i=1:size(n,2)
    plot(n(i),F_emd(n(i)),'.r','MarkerSize',7);hold on
end

spks = EventDetection(F_tab(13,:));
spk = zeros(1,941);
spk(spks) = 1;
plot(tvec, F_tab(13,:)); hold on;
bar(tvec,spk*500,'EdgeColor','r');
clc