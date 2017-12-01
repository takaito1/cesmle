close all;
clear all;

suff={'001' '002' '009' '010' '011' '012' '013' '014' '015' '016' ...
      '017' '018' '019' '020' '021' '022' '023' '024' '025' '026' ...
      '027' '028' '029' '030' '031' '032' '033' '034' '035' '101' ...
      '102' '103' '104' '105'};
Nens=length(suff);

t1=[1970.5:2009.5]';
z  = ncread('/data/dataset/model/CESM-LE/mattlong/O2_001.nc','z_t');

for m=1:Nens

fn=['/data/dataset/model/CESM-LE/mattlong/O2_',suff{m},'.nc']

klev=find(z>=2e3*100&z<4e3*100);
core2D=ncread(char(fn),'O2',[1 1 klev(1) 51],[360 180 length(klev) 40]);
klev=find(z>=4e3*100);
core2A=ncread(char(fn),'O2',[1 1 klev(1) 51],[360 180 length(klev) 40]);

disp('loading masks');
masks = ones(360,180); masks(:,41:180)=0;

% cosine factor for global integration
y = -89.5:89.5;
x = .5:359.5;
[xx,yy]=meshgrid(x,y);
cosy=cos(yy'/180*pi);

% 1 = deep ocean; 2 = abyssal ocean
ND=size(core2D);
NA=size(core2A);

for i=1:2

   if i == 1
      mask = repmat(masks,[1 1 ND(3) ND(4)]);
      mask(mask==0)=NaN;
      mask(isnan(core2D))=NaN;
      cosy2 = repmat(cosy,[1 1 ND(3) ND(4)]);
      V(:,:,i)=squeeze(nansum(nansum(mask.*cosy2,1),2));
      o2le(:,:,i,m)=squeeze(nansum(nansum(core2D.*mask.*cosy2,1),2))./V(:,:,i);
   elseif i == 2
      mask = repmat(masks,[1 1 NA(3) NA(4)]);
      mask(mask==0)=NaN;
      mask(isnan(core2A))=NaN;
      cosy2 = repmat(cosy,[1 1 NA(3) NA(4)]);
      V(:,:,i)=squeeze(nansum(nansum(mask.*cosy2,1),2));
      o2le(:,:,i,m)=squeeze(nansum(nansum(core2A.*mask.*cosy2,1),2))./V(:,:,i);
   end

   Ndat = length(t1);
      tmp  = squeeze(nanmean(o2le(:,:,i,m),1))';
      tmp0=mean(tmp);
      tmp1=mean(t1);
      tmp2=mean(tmp.*t1);
      tmp3=tmp2-tmp1*tmp0;
      oxy_trend_le(i,m)=tmp3/var(t1);
      r2_le(i,m)=tmp3^2/var(t1)/var(tmp);
      tmp4=tmp0 - oxy_trend_le(i,m)*tmp1;
      tmp5=tmp - tmp4 - oxy_trend_le(i,m)*t1;
      tmp6=sum((t1 - tmp1).^2);
      tmp7=corrcoef(tmp(1:end-1),tmp(2:end));
      Neff(i,m)=(1-tmp7(1,2))/(1+tmp7(1,2))*Ndat;
      stderr_trend_le(i,m)=sqrt(sum(tmp5.^2)/(Neff(i,m)-2)/tmp6);
      tval95_le(i,m)=tinv(.975,Neff(i,m));
end

end
save SO_o2trend_le.mat oxy_trend_le o2le r2_le stderr_trend_le tval95_le;


