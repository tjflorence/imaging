len_function = [];


pvec = primes(100);


for ii = 1:10
    
    cprime = pvec(ii);
    
    grow_fact = cprime/2;
    
    len_function = [len_function grow_fact];
end

f1 = figure('color', 'w')

area_covered = 1:10;
plot(area_covered, 'r', 'LineWidth', 2)
hold on
plot(len_function, 'k', 'lineWidth', 2)

text(2, 18, 'area covered', 'color', 'r', 'Fontsize', 20)
text(2, 17, 'time per trial', 'FontSize', 20)

xlabel('number of bars', 'FontSize', 25)
ylabel('normalized increase', 'FontSize', 25)

set(gca, 'Xtick', [1 5 10], 'YTick', [0 10 20], 'FontSize', 20)
hold off

xlim([1 10])
ylim([0 20])

box off

export_fig('multibar_gain_limits.pdf', '-pdf')