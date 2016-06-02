
function [arrOut, Xq, Yq]= interp2array(xs,xe,ys,ye,Data,stepV,stepD,lc)
% use stepV = 1 for full data, larger number to sub sample.
% lc lower cutoff

len = length(Data);
ss = Data(1:stepV:len,:);

Value = ss(:,3);
Value(Value>2000.) = nan;
Value(Value<lc)    = nan;

Xq  = linspace(xs,xe,floor((xe-xs)/stepD+1));
Yq  = linspace(ys,ye,floor((ye-ys)/stepD+1));

F = TriScatteredInterp(ss(:,1:2), Value,'natural');

[qy,qx] = meshgrid(Yq,Xq);
arrOut = F(qx,qy);
