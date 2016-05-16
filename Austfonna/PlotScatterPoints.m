function PlotScatterPoints(X,Y,Z,Cmap,Cstep,c_units,clable,Msize)
zmin=min(Z);
zmax=max(Z);
map=colormap(Cmap);
color_steps=size(map,1);


for i=1:color_steps
    ind=find(Z<zmin+i*(zmax-zmin)/color_steps & Z>=zmin+(i-1)*(zmax-zmin)/color_steps);
    plot(X(ind),Y(ind),'.','MarkerSize',Msize,'Color',map(i,:));hold on;
end
%contourm(double(X),double(Y),double(Z),10);
if Cstep~=0
   Cstep
  h=colorbar;
  Yt = get(h,'YTick'); 
  set(h,'YTick',linspace(Yt(1),Yt(end),Cstep),...
    'YTickLabel',clable,...
    'FontSize',15);
  %h=xlabel(h,c_units);
  ylabel(h,c_units);
  
end

end