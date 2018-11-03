clc ; clear all,close all;

study = benchmark_study3d('SRC');
% study.touch_data(1);
study.run_benchmark_study_3d; 
study.merge_results;
check = study.check_data;
results = study.get_study_results;
study.unconservative_error_cat;
error_accepted_as_real = 10;
study.run_benchmark_error_list(error_accepted_as_real);
% 
hold all
study.figure(14,'surface1')
study.figure(14,'surface4')
