clear all;
close all;

OFn = zeros(5,12,4);
OFm = zeros(5,12,4);

for i = [1 2 3 4 5]
    for j = [4 6 8 10 12]
        for k = [0 1 2 3]
                MW(i,j,k+1) = 12.0*i+1.0*j+16.0*k;
                OFn(i,j,k+1) = 1.0*i+j/4.0-k/2.0;
                OFm(i,j,k+1) = OFn(i,j,k+1)*32.0/MW(i,j,k+1);
        end
    end
end

OFm_plot = zeros(4,5);

for o = [0 1 2 3]
    for c = [1 2 3 4 5]
        OFm_plot(o+1,c) = OFm(c,2*c+2,o+1);
        MW_plot(o+1,c) = MW(c,2*c+2,o+1);
    end
end

%enthalpy of formation of alkane oxygenates [kJ/mol]
h_fch4 = -74.87;
h_fch4o = -205.0;
% no data for ch4o2 (methandiol) or ch4o3
h_fc2h6 = -84.0;
h_fc2h6o = -234.0;
h_fc2h6o2 = -394.4;
% no data for c2h6o3
h_fc3h8 = -104.7;
h_fc3h8o = -256.0;
h_fc3h8o2 = -408.4;
h_fc3h8o3 = -577.9;
h_fc4h10 = -125.6;
h_fc4h10o = -277.0;
h_fc4h10o2 = -426.0;
% no data for c4h10o3
h_fc5h12 = -146.8;
h_fc5h12o = -298.0;
h_fc5h12o2 = -442.0;
% no data for c5h12o3

h_fco2 = -393.51;
h_fh2o = -241.83;
hf_fuels = [h_fch4 h_fc2h6 h_fc3h8 h_fc4h10 h_fc5h12; h_fch4o h_fc2h6o h_fc3h8o h_fc4h10o h_fc5h12o; NaN h_fc2h6o2 h_fc3h8o2 h_fc4h10o2 h_fc5h12o2];
nco2 = [1 2 3 4 5; 1 2 3 4 5; 1 2 3 4 5];
nh2o = (nco2.*2+2)./2;
deltaH = hf_fuels - nco2.*h_fco2 - nh2o.*h_fh2o;
LHV = deltaH./MW_plot(1:3,:);
%Assume kerosene is approx. C12H26
LHV_kerosene = 43.1;
OFm_kerosene = 18.5*32/(12.0*12+1.0*26);
%Assess single triol molecule with data (c3h8o3)
MW_c3h8o3 = 12.0*3 + 1.0*8 + 16.0*3;
deltaH_c3h8o3 = h_fc3h8o3 - 3.0*h_fco2 - 4*h_fh2o;
LHV_c3h8o3 = deltaH_c3h8o3/MW_c3h8o3;

%Relative increase/decrease in oxidizer mass for fixed energy content
%(relative to methane)
LHVratio = (1./LHV).*LHV(1,1);
LHVratio_kerosene = (1./LHV).*LHV_kerosene;
OFmratio = OFm_plot(1:3,:)./OFm_plot(1,1);
OFmratio_kerosene = OFm_plot(1:3,:)./OFm_kerosene;
RelOx = LHVratio.*OFmratio;
RelOx_kerosene = LHVratio_kerosene.*OFmratio_kerosene;
RelOx_c12h26 = (LHV(1,1)/LHV_kerosene)*(OFm_kerosene/OFm_plot(1,1));

fig = figure(1);
[C,h] = contourf([1 2 3 4 5], [0 1 2 3], OFm_plot(:,:));
h.LevelList = [0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.25 2.5 2.75 3 3.25 3.5 3.75];
colormap(jet);
ax = gca;
ax.XTick = [1 2 3 4 5];
ax.YTick = [0 1 2 3];
xlabel('carbon number');
ylabel('oxygen count');
ax.FontSize = 14;
ax.LabelFontSizeMultiplier = 1.1;
cb = colorbar;
cb.Label.String = 'Stoich O/F Mixture Ratio';
cb.Label.FontSize = 14;
hold on;
x_sp1 = [1 1 1 1];
x_sp2 = x_sp1.*2;
x_sp3 = x_sp1.*3;
x_sp4 = x_sp1.*4;
x_sp5 = x_sp1*5;
x_sp = cat(2,x_sp1,x_sp2,x_sp3,x_sp4,x_sp5);
y_sp = repmat([0 1 2 3],1,5);
s = scatter(x_sp,y_sp);
s.MarkerFaceColor = [0.5 0.5 0.5];
s.MarkerEdgeColor = 'black';
s.LineWidth = 1.0;
s.SizeData = 50;
fig.OuterPosition = [200 200 700 500];
pbaspect([1 1 1]);

fig = figure(2);
[C,h] = contourf([1 2 3 4 5], [0 1 2], LHV(:,:));
colormap(jet);
ax = gca;
ax.XTick = [1 2 3 4 5];
ax.YTick = [0 1 2];
xlabel('carbon number');
ylabel('oxygen count');
ax.FontSize = 14;
ax.LabelFontSizeMultiplier = 1.1;
cb = colorbar;
cb.Label.String = 'LHV [kJ/kg]';
cb.Label.FontSize = 14;
fig.OuterPosition = [200 200 700 500];
pbaspect([1 1 1]);

fig = figure(3);
[C,h] = contourf([1 2 3 4 5], [0 1 2], RelOx(:,:));
colormap(jet);
ax = gca;
ax.XTick = [1 2 3 4 5];
ax.YTick = [0 1 2];
xlabel('carbon number');
ylabel('oxygen count');
ax.FontSize = 14;
ax.LabelFontSizeMultiplier = 1.1;
cb = colorbar;
cb.Label.String = {'Oxidizer Mass Needed Relative to CH_{4}','for fixed energy content'};
cb.Label.FontSize = 14;
fig.OuterPosition = [200 200 700 500];
pbaspect([1 1 1]);
hold on;
s = scatter([1 2 3 4 5 1 2 3 4 5 1 2 3 4 5],[0 0 0 0 0 1 1 1 1 1 2 2 2 2 2]);
s.MarkerFaceColor = [0.5 0.5 0.5];
s.MarkerEdgeColor = 'black';
s.LineWidth = 1.0;
s.SizeData = 50;

% figure(4)
% [C,h] = contourf([1 2 3 4 5], [0 1 2], RelOx_kerosene(:,:));
% colormap(jet);
% ax = gca;
% ax.XTick = [1 2 3 4 5];
% ax.YTick = [0 1 2];
% xlabel('carbon number');
% ylabel('oxygen count');
% ax.FontSize = 14;
% ax.LabelFontSizeMultiplier = 1.1;
% cb = colorbar;
% cb.Label.String = 'Oxidizer Mass Needed Relative to RP-1';
% cb.Label.FontSize = 14;