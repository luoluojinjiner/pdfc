
function y = distance(a)
s = size(a,1);%Çó¾àÀë¾ØÕó
%for i=1:s
P_x = repmat(a(:,1),1,s);
P_y = repmat(a(:,2),1,s);
Pt_x = repmat(a(:,1)',s,1);
Pt_y = repmat(a(:,2)',s,1);
y = sqrt((P_x - Pt_x).^2 + (P_y - Pt_y).^2);
%end