nb = 3;
ir = 0;
jr = 0;
ind_hex = zeros(nsites/3*2*nb^2,1);
ind_tri = zeros(nsites/3*nb^2,1);
supercell = zeros(nsites*nb^2,3);
for ic = -(nb-1)/2 : (nb-1)/2
    for jc = -(nb-1)/2 : (nb-1)/2
        for is = 1 : nsites
            ir = ir + 1;
            supercell(ir,:) = mcell(is,:) + ic*ma1 + jc*ma2;
            if(mod(ir,3)~=0)
                jr = jr+1;
                ind_hex(jr) = ir;
            end
        end
    end
end
ind_tri = setdiff(linspace(1,nsites*nb^2,nsites*nb^2),ind_hex);

super_ics = zeros(nb^2*ncharges,1);
for in = 0 : nb^2-1
    super_ics(in*ncharges+1:ncharges*(in+1)) = ics_best(:) + in*nsites;
end

% superlattice = supercell(2:2:end,:) - 1/3*(a1+a2);
% 
% voronoi(superlattice(:,1),superlattice(:,2))
% hold on
% if(add_disorder)
%     disorder_colors = (disorder' - min(disorder))./max(disorder - min(disorder));
%     sup_disorder_colors = repmat(disorder_colors,[nb^2,1]);
% %     f = scatteredInterpolant(supercell(:,1),supercell(:,2),repmat(disorder_colors,[nb^2,1]),'natural');
% %     z = f(xt',yt');
% %     surf(xt',yt',z,'FaceAlpha',0.15);
% %     shading interp
% %     view(0,90)
% %     colormap('hot')
% %     grid off
% %     colorbar
% 
% [v,c] = voronoin(superlattice,{'Qbb'});
% for i = 1:length(c)
% fill(v(c{i},1),v(c{i},2),sup_disorder_colors(i),'FaceAlpha',0.25) ;
% end
% grid off
% colormap('hot')
% colorbar

%     hFills = h.FacePrims;
%     for i = 1:length(hFills)
%     % Have to set this. The default is 'truecolor' which ignores alpha.
%     h(i).ColorType = 'truecoloralpha';
%     % The 4th element is the 'alpha' value. First 3 are RGB. Note, the
%     % values expected are in range 0-255.
%     h(i).ColorData(4) = 0;
%     end 
%end

[x,y] = meshgrid(-(nb-1)/2*L : (nb-1)/2*L,-(nb-1)/2*L : (nb-1)/2*L);
xy = [x(:),y(:)];
T = [a1;a2];
xyt = xy*T;
xt = reshape(xyt(:,1),size(x));
yt = reshape(xyt(:,2),size(y));
zt1 = 6.6*ones(size(xt));
zt2 = zeros(size(xt));

figure
plot3(xt,yt,zt1','r-','LineWidth',2)
hold on
plot3(xt',yt',zt1','r-','LineWidth',2)
axis equal
hold on
%,zeros(size(xt,1))
surf(xt,yt,6.6.*ones(size(xt,1)),'FaceAlpha',0.3)
shading interp
colormap([1 0 0])
view(25.855319653140974,29.427658393246325)
scatter3(supercell(ind_hex(1:2:end),1),supercell(ind_hex(1:2:end),2),supercell(ind_hex(1:2:end),3),100,zeros(size(ind_hex(1:2:end),1),1),'filled','MarkerFaceColor','w','MarkerEdgeColor','r','LineWidth',2)
scatter3(supercell(ind_hex(2:2:end),1),supercell(ind_hex(2:2:end),2),supercell(ind_hex(2:2:end),3),100,zeros(size(ind_hex(2:2:end),1),1),'filled','MarkerFaceColor','w','MarkerEdgeColor',[0.5 0 0.5],'LineWidth',2)
scatter3(supercell(super_ics,1),supercell(super_ics,2),supercell(super_ics,3),100,zeros(size(super_ics,1),1),'filled','MarkerFaceColor','r')%,'MarkerEdgeColor','k')
axis([-3*alat 3*alat -3*alat 3*alat 6 7])

%colormap(flipud(gray))
str_fill = num2str(filling);
caxis([0 1]);
set(gca,'FontSize',32)
xlabel(gca,'$x$ [\AA]','Interpreter','latex','FontSize',32)
ylabel(gca,'$y$ [\AA]','Interpreter','latex','FontSize',32)
title(gca,join(['$\nu=$',str_fill]),'Interpreter','latex')
saveas(gca,join([filename,'_Gamma','.fig']))


figure
plot3(xt,yt,zt2','b-','LineWidth',2)
hold on
plot3(xt',yt',zt2','b-','LineWidth',2)
axis equal
hold on
surf(xt,yt,zeros(size(xt,1)),'FaceAlpha',0.3)
shading interp
colormap([0 0 1])
view(25.855319653140974,29.427658393246325)
scatter3(supercell(ind_tri,1),supercell(ind_tri,2),supercell(ind_tri,3),100,zeros(size(ind_tri',1),1),'filled','MarkerFaceColor','w','MarkerEdgeColor','b','LineWidth',2)
scatter3(supercell(super_ics,1),supercell(super_ics,2),supercell(super_ics,3),100,zeros(size(super_ics,1),1),'filled','MarkerFaceColor','b')%,'MarkerEdgeColor','k')
axis([-3*alat 3*alat -3*alat 3*alat -0.5 0.5])


%colormap(flipud(gray))
str_fill = num2str(filling);
caxis([0 1]);
set(gca,'FontSize',32)
xlabel(gca,'$x$ [\AA]','Interpreter','latex','FontSize',32)
ylabel(gca,'$y$ [\AA]','Interpreter','latex','FontSize',32)
title(gca,join(['$\nu=$',str_fill]),'Interpreter','latex')
saveas(gca,join([filename,'_K','.fig']))