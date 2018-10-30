clc ; clear all,close all;

study = benchmark_study3d('RCFT');
% study.touch_data(1);
study.run_benchmark_study_3d; 
% check = study.check_data;
% study.unconservative_error_cat;
results = study.get_study_results;
% 
% hold all
% study.figure(1,'surface1')