classdef benchmark_study3d < handle
    
    properties
            
        section_type
        numErrors      = 9
        debug          = false 
    end
    methods
        % Constructor
        function obj = benchmark_study3d(section_type)
            obj.section_type = section_type;
        end
        
        function data = frames_data(obj)
            switch obj.section_type
                case 'SRC'
                    sections = build_sections_SRC_3d();
                case 'RCFT'
                    sections = build_sections_RCFT_3d();
            end
            data = build_data_from_sections_3d(sections);            

        end
        
        function num_frames = num_frames(obj)
            num_frames = size(frames_data(obj),2);
        end
        
        function frame_data  = frame_data(obj,frame_tag)
            frames_data   = obj.frames_data;
            frame_data = frames_data(frame_tag);
        end
        
        function frame_list = frame_list(obj)
            frame_list = 1:num_frames(obj);
        end

        function run_benchmark_study_3d(obj)

            parfor frame = frame_list(obj)
            
                benchmark   = strength_interaction(frame_data(obj,frame));
                
                results = get_results(obj,frame);
                if results.all_data_stored
                    continue
                end
                results = struct;
                results.frame_data = frame_data(obj,frame);                    
                if obj.debug
                    
                    [results.surface1,results.surface2] =...
                        benchmark.surface_1_and_2;
                    results.surface3 = benchmark.surface_3;
                    results.surface4 = benchmark.surface_4;
                    results.surface5 = benchmark.surface_5(results.surface1);
                    
                else
                    try
                        [results.surface1,results.surface2] =...
                            benchmark.surface_1_and_2;
                    catch
                        continue
                    end

                    results.surface3 = benchmark.surface_3;

                    try
                        results.surface4 = benchmark.surface_4;
                    catch
                        continue
                    end
                    try
                        results.surface5 = benchmark.surface_5(results.surface1);
                    catch
                        continue
                    end                     
                end
                
                results.first_error  = unconservative_error(obj,results.surface1,...
                    results.surface4);
     
                results.second_error = unconservative_error(obj,results.surface5,...
                    results.surface3);

                results.ind_rho = benchmark.ind_rho;
                results.ind_lambda = benchmark.ind_lambda...
                    (-max(results.surface1.P(1,:)));

                results.all_data_stored = true;
                
                obj.save_results(results,frame);

            end
        end
	
        function merge_results(obj)
            for frame = frame_list(obj)
                results_tmp = get_results(obj,frame);
                for fn = fieldnames(results_tmp)'
                   results.all(frame).(fn{1}) = results_tmp.(fn{1});
                end
            end
            save_study(obj,results);
        end
        
        function results = get_study_results(obj)
            if ~exist(study_file(obj),'file')
                results =struct;
                for frame = 1:num_frames(obj)
                    results.all(frame).all_data_stored = false;
                end
               save_study(obj,results)
            else
                
                load(study_file(obj),'results');
            end 
        end

         function results = get_results(obj,frame_tag)
            if ~exist(results_file(obj,frame_tag),'file')
                results = struct;
                results.all_data_stored = false;
                save_results(obj,results,frame_tag)
            else
                load(results_file(obj,frame_tag),'results');
            end 
         end        


        function val = check_data(obj)
            if ~exist(study_file(obj),'file')
                warning('No bunchmark study found for %s',obj.section_type)
                return 
            end
            results = get_study_results(obj);
            val = nan(size(obj.frame_list,2),1);
            for  frame = 1:num_frames(obj)
                if ~results.all(frame).all_data_stored
                    val(frame) = false;
                else
                    val(frame) = true;
                end
            end
        end

        function error=unconservative_error(obj,surface1,surface4)
            angle = linspace(0,pi/2,obj.numErrors);            
            numAngles = size(surface1.P,2);
            for i = 1:numAngles                                              

                M1(:,i) = sqrt(surface1.My(:,i).^2 + surface1.Mz(:,i).^2);
                M4(:,i) = sqrt(surface4.My(:,i).^2 + surface4.Mz(:,i).^2);

                id_1(i) = interactionDiagram2d(M1(:,i),abs(surface1.P(:,i)));
                id_4(i) = interactionDiagram2d(M4(:,i),abs(surface4.P(:,i)));

                for j=1:obj.numErrors
                    r_inelastic = id_1(i).radial_distance(angle(j));
                    r_elastic   = id_4(i).radial_distance(angle(j));   
                    error.all(j,i) = (r_inelastic-r_elastic)/r_inelastic*100;
                end                  

            end
            error.cat.minor = [min(min(min(error.all(          1                : round(obj.numErrors/3),  1),0)))...
                               min(min(min(error.all(  round(obj.numErrors/3)+1 : round(obj.numErrors*2/3),1),0)))...
                               min(min(min(error.all(round(obj.numErrors*2/3)+1 :       obj.numErrors,     1),0)))];
                           
            error.cat.major = [min(min(min(error.all(          1                : round(obj.numErrors/3),  numAngles),0)))...
                               min(min(min(error.all(  round(obj.numErrors/3)+1 : round(obj.numErrors*2/3),numAngles),0)))...
                               min(min(min(error.all(round(obj.numErrors*2/3)+1 :       obj.numErrors,     numAngles),0)))]; 
                           
            error.cat.biaxial= [min(min(min(error.all(          1                : round(obj.numErrors/3), 2:numAngles-1  ),0)))...
                               min(min(min(error.all(  round(obj.numErrors/3)+1 : round(obj.numErrors*2/3),2:numAngles-1  ),0)))...
                               min(min(min(error.all(round(obj.numErrors*2/3)+1 :       obj.numErrors,     2:numAngles-1  ),0)))]; 
                           
            error.cat.combined=[min(min(min(error.all(          1                : round(obj.numErrors/3),  :),0)))...
                               min(min(min(error.all(  round(obj.numErrors/3)+1 : round(obj.numErrors*2/3), :),0)))...
                               min(min(min(error.all(round(obj.numErrors*2/3)+1 :       obj.numErrors,      :),0)))];                
        end
        
        function unconservative_error_cat(obj)
            results = get_study_results(obj);
            results.cat.first_error   = cell(5,6);
            results.cat.second_error  = cell(5,6);
            for frame = frame_list(obj)
                
               if results.all(frame).all_data_stored
                   ind_lambda = results.all(frame).ind_lambda;
                   ind_rho = results.all(frame).ind_rho;
                   results.cat.first_error{ind_rho,ind_lambda} = vertcat(results.cat.first_error...
                       {ind_rho,ind_lambda},horzcat(results.all(frame).first_error.cat.biaxial,frame));
                   results.cat.second_error{ind_rho,ind_lambda} = vertcat(results.cat.second_error...
                       {ind_rho,ind_lambda},horzcat(results.all(frame).second_error.cat.biaxial,frame));
               else 
                   warning('unconservative_error_cat : Frame %i has not neccessary data stored',frame);
               end
            end
           save_study(obj,results);
        end

        function study_dir = study_dir(obj)
            study_dir = sprintf('%s/%s',pathOf.results_benchmark,obj.section_type);
        end

        function study_file = study_file(obj)
            study_file = [study_dir(obj) '/results.mat'];
        end

        function save_study(obj,results)
            save(study_file(obj),'results')
        end
        
        function results_dir = results_dir(obj)
            results_dir = sprintf('%s/%s',pathOf.results_study,obj.section_type);
        end        
        function results_file = results_file(obj,frame_tag)
            results_file = [results_dir(obj) sprintf('/results-%i.mat',frame_tag)];
        end
        
        function save_results(obj,results,frame_tag)
            save(results_file(obj,frame_tag),'results')
        end
        
        function touch_data(obj,frame_tag)
            results = get_results(obj,frame_tag);
            results.all_data_stored = false;
            save_results(obj,results)
        end
        
        function error_list = error_list(obj,accepted_error)
            results = get_study_results(obj);
            if ~isfield(results,'cat')
                unconservative_error_cat(obj)
                results = get_study_results(obj);
            end
            error_list = [0,0];
            for i = 1:size(results.cat.first_error,1)
                for j = 1:size(results.cat.first_error,2)
                    for k = 1:size(results.cat.first_error{i,j},1)
                        max_error = max(abs(results.cat.first_error{i,j}(k,1:3)));
                        if max_error >= accepted_error
                            frame_tag = results.cat.first_error{i,j}(k,4);
                            error_list = vertcat(error_list,[frame_tag,max_error]);
                        end
                    end
                end
            end
            error_list(1,:)=[];
            error_list = sortrows(error_list,1);
        end
        function run_benchmark_error_list(obj,accepted_error)
            
            frame_list = error_list(obj,accepted_error);
            parfor e = 1 : size(frame_list,1)
                frame = frame_list(e,1);
            
                benchmark = strength_interaction(frame_data(obj,frame));
                benchmark.extra_fine_DispStepSize = true;
                
                results = get_results(obj,frame);
                if obj.debug
                    
                    [results.surface1,results.surface2] =...
                        benchmark.surface_1_and_2;
                    results.surface5 = benchmark.surface_5(results.surface1);
                    
                else
                    try
                        [results.surface1,results.surface2] =...
                            benchmark.surface_1_and_2;
                    catch
                        continue
                    end

                    try
                        results.surface5 = benchmark.surface_5(results.surface1);
                    catch
                        continue
                    end                     
                end
                
                results.first_error  = unconservative_error(obj,results.surface1,...
                    results.surface4);
     
                results.second_error = unconservative_error(obj,results.surface5,...
                    results.surface3);

                delete(results_file(obj,frame));
                obj.save_results(results,frame);

            end
        end        
        
        
        
        
    
        function figure(obj,frame_tag,varargin)
            results = get_study_results(obj);
            num_figures = numel(varargin);
            hold all
            for i=1:num_figures
                if strcmp(varargin{i},'surface1')
                    surf(results.all(frame_tag).surface1.Mz,...
                        results.all(frame_tag).surface1.My,...
                        -results.all(frame_tag).surface1.P);

                elseif strcmp(varargin{i},'surface2')
                    surf(results.all(frame_tag).surface2.Mz,...
                        results.all(frame_tag).surface2.My,...
                        -results.all(frame_tag).surface2.P);                        
                elseif strcmp(varargin{i},'surface3')
                    surf(results.all(frame_tag).surface3.Mz,...
                        results.all(frame_tag).surface3.My,...
                        -results.all(frame_tag).surface3.P);                           
                elseif strcmp(varargin{i},'surface4')
                    surf(results.all(frame_tag).surface4.Mz,...
                        results.all(frame_tag).surface4.My,...
                        -results.all(frame_tag).surface4.P);                        
                elseif strcmp(varargin{i},'surface5')
                    surf(results.all(frame_tag).surface5.Mz,...
                        results.all(frame_tag).surface5.My,...
                        -results.all(frame_tag).surface5.P);                       
                end
            end
        end                   
    end
end
    
