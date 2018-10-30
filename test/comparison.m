clear all; close all; clc;

sections = build_sections_RCFT_3d();
data0 = build_data_from_sections_3d(sections);
data = data0(21);

benchmark = strength_interaction(data);
benchmark.axis = 'weak';
benchmark.echoOpenSeesOutput = true;

[surface1,surface2] = benchmark.surface_1_and_2;
surface3 = benchmark.surface_3;
surface4 = benchmark.surface_4;
surface5 = benchmark.surface_5(surface1);
% 
results_2d = benchmark.curve1_5_2d;



%% Figures
if strcmp(benchmark.axis,'weak')
    M1 = surface1.My;
    M2 = surface2.My;
    M3 = surface3.My;
    M4 = surface4.My;
    M5 = surface5.My;
    
elseif strcmp(benchmark.axis,'strong')
    M1 = surface1.Mz;
    M2 = surface2.Mz;
    M3 = surface3.Mz;
    M4 = surface4.Mz;
    M5 = surface5.Mz; 
end
figure
hold all
plot(M1,-surface1.P,'s-')
plot(results_2d.Curve1_M1,-results_2d.Curve1_P1,'o-')

title(sprintf('Curve 1-%s axis',benchmark.axis));    


figure
hold all
plot(M2,-surface2.P,'s-')
plot(results_2d.Curve2_M2,-results_2d.Curve2_P2,'o-')

title(sprintf('Curve 2-%s axis',benchmark.axis));    

figure
hold all
plot(M3,-surface3.P,'s-')
plot(results_2d.Curve3_M2,-results_2d.Curve3_P2,'o-')

title(sprintf('Curve 3-%s axis',benchmark.axis));    

figure
hold all
plot(M4,-surface4.P,'s-')
plot(results_2d.Curve4_M1,-results_2d.Curve4_P1,'o-')

title(sprintf('Curve 4-%s axis',benchmark.axis));    

figure
hold all
plot(M5,-surface5.P,'s-')
plot(results_2d.Curve5_M2,-results_2d.Curve5_P2,'o-')

title(sprintf('Curve 5-%s axis',benchmark.axis));    


