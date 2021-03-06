gr1=[80,77,73,69;64,66,68,71;75,73,73,69;72,70,74,73;74,74,71,67;71,71,72,70;76,78,74,71;73,68,64,64;76,73,74,76;77,78,77,73];
gr2=[81,81,82,82;82,83,80,81;81,77,80,80;84,86,85,85;88,90,88,86;83,82,86,85;85,83,87,86;81,85,86,85;87,89,87,82;77,75,73,77];
gr3=[76,83,85,79;75,81,85,73;75,82,80,77;68,73,72,69;78,87,86,77;81,85,81,74;67,73,75,66;68,73,73,66;68,75,79,69;73,78,80,70];

sve=[gr1;gr2;gr3];
C(1,:)=mean(gr1,1);
C(2,:)=mean(gr2,1);
C(3,:)=mean(gr3,1);

plot(C')

lambda4(1,:)=ones(1,3);
lambda4(2:4,:)=-eye(3);

lambda3(1,:)=ones(1,2);
lambda3(2:3,:)=-eye(2);

lambda11(1,:)=ones(1,10);
lambda11(2:11,:)=-eye(10);

s1=cov(gr1);
s2=cov(gr2);
s3=cov(gr3);
e=9*s1+9*s2+9*s3;
s=e/27;
ss=lambda4'*s*lambda4;
ee=lambda4'*e*lambda4;

novaGr1=gr1*lambda4;
novaGr2=gr2*lambda4;
novaGr3=gr3*lambda4;
cv1=cov(novaGr1);
cv2=cov(novaGr2);
cv3=cov(novaGr3);
en=9*cv1+9*cv2+9*cv3;
sn=en/27;
X=diag([10,10,10]);

H=(lambda3'*C*lambda4)'*inv(lambda3'*inv(X)*lambda3)'*lambda3'*C*lambda4;
%st=det(H)
pog=det(ee+H)
det(ss)
st=det(ee)
%det(H+ee);
rez=st/pog
man=[novaGr1;novaGr2;novaGr3];
abc(1:10)=ones(1,10);
abc(11:20)=2*ones(1,10);
abc(21:30)=3*ones(1,10);
[D,p,stats]=manova1(man,abc,0.5)
%hottelingova stat - t^2 tr(g-1*p) f aproksimacija t2*e-1


k=lambda3'*C*lambda4;

