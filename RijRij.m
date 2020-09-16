
%====== Calculating RijRij-I3 and 3RijRij-I3 in interaction matrix A =====%
%=========================================================================%
function [rjkrjk1_I,rjkrjk2_I,rjkrjk3_I,rjkrjk4_I,rjkrjk5_I,rjkrjk6_I,rjkrjk31_I,...
    rjkrjk32_I,rjkrjk33_I,rjkrjk34_I,rjkrjk35_I,rjkrjk36_I,RJK]=RijRij(r_block)

rkj1=r_block(1,1)-r_block(:,1);
rkj2=r_block(1,2)-r_block(:,2);
rkj3=r_block(1,3)-r_block(:,3);

rk_to_rj=[rkj1  rkj2  rkj3];
rk_to_rj(1,:)=1;                % in order to scape from NAN Error
RJK=sqrt(rk_to_rj(:,1).^2+rk_to_rj(:,2).^2+rk_to_rj(:,3).^2);
rjkrjk=[rkj1./RJK  rkj2./RJK  rkj3./RJK];


rjkrjk1_I=rjkrjk(:,1).*rjkrjk(:,1)-1;
rjkrjk2_I=rjkrjk(:,1).*rjkrjk(:,2);
rjkrjk3_I=rjkrjk(:,1).*rjkrjk(:,3);
rjkrjk4_I=rjkrjk(:,2).*rjkrjk(:,2)-1;
rjkrjk5_I=rjkrjk(:,2).*rjkrjk(:,3);
rjkrjk6_I=rjkrjk(:,3).*rjkrjk(:,3)-1;

rjkrjk31_I=3*rjkrjk(:,1).*rjkrjk(:,1)-1;
rjkrjk32_I=3*rjkrjk(:,1).*rjkrjk(:,2);
rjkrjk33_I=3*rjkrjk(:,1).*rjkrjk(:,3);
rjkrjk34_I=3*rjkrjk(:,2).*rjkrjk(:,2)-1;
rjkrjk35_I=3*rjkrjk(:,2).*rjkrjk(:,3);
rjkrjk36_I=3*rjkrjk(:,3).*rjkrjk(:,3)-1;

rjkrjk1_I(1,1)=0;
rjkrjk2_I(1,1)=0;
rjkrjk3_I(1,1)=0;
rjkrjk4_I(1,1)=0;
rjkrjk5_I(1,1)=0;
rjkrjk6_I(1,1)=0;

rjkrjk31_I(1,1)=0;
rjkrjk32_I(1,1)=0;
rjkrjk33_I(1,1)=0;
rjkrjk34_I(1,1)=0;
rjkrjk35_I(1,1)=0;
rjkrjk36_I(1,1)=0;

end
%=========================================================================%