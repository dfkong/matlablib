function subplot2(m,n)
figure;
for ii=1:n*(m-1)
    subplot(m,n,ii);
    set(gca, 'XMinorTick','on');
    set(gca, 'YMinorTick','on');
    myplot(m,n,ii);
end
for ii=n*(m-1)+1:m*n
    subplot(m,n,ii);
    set(gca, 'XMinorTick','on');
    set(gca, 'YMinorTick','on');
end
