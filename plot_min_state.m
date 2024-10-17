figure;
nb = 3;
ir = 0;
supercell = zeros(nsites*nb^2,2);
for ic = -(nb-1)/2 : (nb-1)/2
    for jc = -(nb-1)/2 : (nb-1)/2
        for is = 1 : nsites
            ir = ir + 1;
            supercell(ir,:) = mcell(is,:) + ic*ma1 + jc*ma2;
        end
    end
end

super_ics = zeros(nb^2*ncharges,1);
for in = 0 : nb^2-1
    super_ics(in*ncharges+1:ncharges*(in+1)) = ics_best(:) + in*nsites;
end


%if(add_disorder)
    disorder_colors = (disorder' - min(disorder))./max(disorder - min(disorder));
    sup_disorder_colors = repmat(disorder_colors,[nb^2,1]);
%     f = scatteredInterpolant(supercell(:,1),supercell(:,2),repmat(disorder_colors,[nb^2,1]),'natural');
%     z = f(xt',yt');
%     surf(xt',yt',z,'FaceAlpha',0.15);
%     shading interp
%     view(0,90)
%     colormap('hot')
%     grid off
%     colorbar

voronoi(supercell(:,1),supercell(:,2))
hold on
[v,c] = voronoin(supercell,{'Qbb'});
for i = 1:length(c)
fill(v(c{i},1),v(c{i},2),sup_disorder_colors(i),'FaceAlpha',0.25) ;
end
grid off
colormap('hot')
colorbar

%     hFills = h.FacePrims;
%     for i = 1:length(hFills)
%     % Have to set this. The default is 'truecolor' which ignores alpha.
%     h(i).ColorType = 'truecoloralpha';
%     % The 4th element is the 'alpha' value. First 3 are RGB. Note, the
%     % values expected are in range 0-255.
%     h(i).ColorData(4) = 0;
%     end 
%end
% [x,y] = meshgrid(-(nb-1)/2*L : (nb-1)/2*L,-(nb-1)/2*L : (nb-1)/2*L);
% xy = [x(:),y(:)];
% T = [a1;a2];
% xyt = xy*T;
% xt = reshape(xyt(:,1),size(x));
% yt = reshape(xyt(:,2),size(y));
% hold on
% plot(xt,yt,'k-')
% hold on
% plot(xt',yt','k-')
axis equal
hold on
scatter(supercell(:,1),supercell(:,2),100,zeros(nb^2*nsites,1),'filled','MarkerFaceColor','w','MarkerEdgeColor','k')
hold on
scatter(supercell(super_ics,1),supercell(super_ics,2),100,ones(nb^2*ncharges,1),'filled','MarkerFaceColor','k','MarkerEdgeColor','k')

axis([-3 3 -3 3]*alat)
%colormap(flipud(gray))
str_fill = num2str(filling);
caxis([0 1]);
set(gca,'FontSize',32)
xlabel(gca,'$x$ [\AA]','Interpreter','latex','FontSize',32)
ylabel(gca,'$y$ [\AA]','Interpreter','latex','FontSize',32)
title(gca,join(['$\nu=$',str_fill]),'Interpreter','latex')
saveas(gca,join([filename,'.fig']))