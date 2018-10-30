classdef BenchmarkAnalysis3d_Driver
    
properties
    benchmark_frame_data
    % Fiber Section Definition Parameters
    FiberSectionDefinitonOptions
    FiberSectionDefinition_SectionID        = 1;
    FiberSectionDefinition_StartMatID       = 1;
    
    % Analysis Options 
    AnalysisOptions                         = struct;
    
    % Design Options 
    Delta0_over_L                           = 1/500;
    delta0_over_L                           = 1/1000;
    ElasticStiffnessType                    = 'ColumnStrength'
    
    % Elastic Analysis
    StiffnessReduction                      = 0.8;
    tau                                     = 0.8;    

end


methods
    %% Constructor
    function obj = BenchmarkAnalysis3d_Driver(data,analysisOptions)
        if nargin > 1
        obj.AnalysisOptions = analysisOptions;
        end        
        obj.benchmark_frame_data = data;
        % Set Fiber Section Definiton Options
        opts = struct;
        opts.includePackageDefinition = true;
        opts.nf1 = 30;
        opts.nf2 = 30;
        switch data.section_type
            case 'RCFT'
                opts.SteelMaterialType                  = 'ModifiedAbdelRahman';
                opts.AbdelRahmanResidualStressParameter = 0.75;
                opts.AbdelRahmanHardeningRatio          = 0.0;
                opts.ConcreteMaterialType               = 'ProposedForDesign';
    
            case 'SRC'
                opts.SteelMaterialType                  = 'ElasticPP';
                opts.ConcreteMaterialType               = 'ProposedForDesign';                    
                
            otherwise 
                error('Unknown section type: %s',data.section_type)
        end
        obj.FiberSectionDefinitonOptions = opts;
        
    end
    
    function BA = get3dAnalysisObject(obj,type,loading_angle,imperfection_angle)
        
        % Update data to match loading angle
        mydata = obj.benchmark_frame_data;
        switch mydata.frame_typeZ
            case 'Sidesway_Inhibited'
                mydata.betaTopZ = -cos(loading_angle); % The negative sign is used to make the moment deflection in the same direction of imperfection. Otherwise betaBotZ_over_betaTopZ should be negative in the build_data_from_section_3d 
                mydata.betaBotZ =  mydata.betaTopZ*mydata.betaBotZ_over_betaTopZ;
                
            case 'Sidesway_Uninhibited'
                if mydata.kqtopZ == 0 && mydata.kqbotZ == 0
                    error('Structure unstable')
                elseif mydata.kqtopZ == 0 || mydata.kqbotZ == 0
                    mydata.Hy = cos(loading_angle)/mydata.L;
                elseif mydata.kqtopZ == mydata.kqbotZ
                    mydata.Hy = cos(loading_angle)/(mydata.L/2);
                else
                    error('Not yet implemented for this condition');
                end                
                
            otherwise
                error('Unknown frame_typeZ: %s',mydata.frame_typeZ)
        end
        
        switch mydata.frame_typeY
            case 'Sidesway_Inhibited'
                mydata.betaTopY = sin(loading_angle);
                mydata.betaBotY = mydata.betaTopY*mydata.betaBotY_over_betaTopY;
                
            case 'Sidesway_Uninhibited'
                if mydata.kqtopY == 0 && mydata.kqbotY == 0
                    error('Structure unstable')
                elseif mydata.kqtopY == 0 || mydata.kqbotY == 0
                    mydata.Hz = sin(loading_angle)/mydata.L;
                elseif mydata.kqtopY == mydata.kqbotY
                    mydata.Hz = sin(loading_angle)/(mydata.L/2);
                else
                    error('Not yet implemented for this condition');
                end                
                
            otherwise
                error('Unknown frame_typeY: %s',mydata.frame_typeY)
        end
        
        % Update data for stiffness reduction
        switch type
            case 'elastic'
                switch mydata.frame_typeZ
                    case 'Sidesway_Inhibited'

                    case 'Sidesway_Uninhibited'
                        mydata.kqtopZ = obj.StiffnessReduction * mydata.kqtopZ; 
                        mydata.kqbotZ = obj.StiffnessReduction * mydata.kqbotZ;              
                    otherwise
                        error('Unknown frame_typeZ: %s',mydata.frame_typeZ)
                end

                switch mydata.frame_typeY
                    case 'Sidesway_Inhibited'

                    case 'Sidesway_Uninhibited'
                        mydata.kqtopY = obj.StiffnessReduction * mydata.kqtopY;
                        mydata.kqbotY = obj.StiffnessReduction * mydata.kqbotY;   
                    otherwise
                        error('Unknown frame_typeY: %s',mydata.frame_typeY)
                end
        end
        
        % Initial Imperfections
        switch type
            case 'elastic'
                mydata.delta0Y = 0;
                mydata.delta0Z = 0;
                mydata.Delta0Y = 0;
                mydata.Delta0Z = 0;
            case 'nonlinear'
                switch mydata.frame_typeZ
                    case 'Sidesway_Inhibited'
                        mydata.delta0Y = cos(imperfection_angle) * obj.delta0_over_L * mydata.L;
                        mydata.Delta0Y = 0;
                    case 'Sidesway_Uninhibited'
                        mydata.delta0Y = cos(imperfection_angle) * obj.delta0_over_L * mydata.L;
                        mydata.Delta0Y = cos(imperfection_angle) * obj.Delta0_over_L * mydata.L;
                    otherwise
                        error('Unknown frame_typeZ: %s',mydata.frame_typeZ)
                end
                switch mydata.frame_typeY
                    case 'Sidesway_Inhibited'
                        mydata.delta0Z = sin(imperfection_angle) * obj.delta0_over_L * mydata.L;
                    case 'Sidesway_Uninhibited'
                        mydata.delta0Z = sin(imperfection_angle) * obj.delta0_over_L * mydata.L;
                        mydata.Delta0Z = sin(imperfection_angle) * obj.Delta0_over_L * mydata.L;
                    otherwise
                        error('Unknown frame_typeY: %s',mydata.frame_typeY)
                end
        end        
        % Section Definition
        switch type
            case 'elastic'

                elastic_data = struct;
                [E,A,Iz,Iy,G,J] = mydata.section.sectionPropertiesForElasticAnalysis3d(obj.ElasticStiffnessType);
                elastic_data.sectionType = 'elastic';
                elastic_data.E  = E;
                elastic_data.A  = obj.StiffnessReduction * A;
                elastic_data.Iz = obj.StiffnessReduction * obj.tau * Iz;
                elastic_data.Iy = obj.StiffnessReduction * obj.tau * Iy;
                elastic_data.G  = G;    
                elastic_data.J  = J;                
                section_def = sprintf('section Elastic 1 %i %i %i %i %i %i \n',...
                       elastic_data.E, elastic_data.A, elastic_data.Iz,...
                       elastic_data.Iy, elastic_data.G,elastic_data.J);                

%                 section_def = FiberSectionDefinition(elastic_data,'3d',...
%                     obj.FiberSectionDefinition_SectionID,...
%                     obj.FiberSectionDefinition_StartMatID,struct);
                
            case 'nonlinear'
                section_def = FiberSectionDefinition(mydata.section,'3d',...
                    obj.FiberSectionDefinition_SectionID,...
                    obj.FiberSectionDefinition_StartMatID,...
                    obj.FiberSectionDefinitonOptions);
                
            otherwise
                error('Unknown type: %s',type)
        end
        
        BA = BenchmarkAnalysis3d_OpenSees(mydata,section_def,obj.AnalysisOptions);
    end
    
    function [BA,mydata] = get2dAnalysisObject(obj,type,axis)
        % @todo -check that this function is correct
        
        
        % Create data for a 2d object
        mydata = struct;
        mydata.L        = obj.benchmark_frame_data.L;
        mydata.section  = obj.benchmark_frame_data.section;
        mydata.axis     = axis;
        switch axis
            case {'z','strong'}
                mydata.frame_type = obj.benchmark_frame_data.frame_typeZ;
                switch mydata.frame_type
                    case 'Sidesway_Inhibited'
                        mydata.beta   = -1/obj.benchmark_frame_data.betaBotZ_over_betaTopZ;
                        mydata.delta0 = obj.delta0_over_L*mydata.L;
                    case 'Sidesway_Uninhibited'
                        mydata.kqtop  = obj.benchmark_frame_data.kqtopZ;
                        mydata.kqbot  = obj.benchmark_frame_data.kqbotZ;
                        mydata.gamma  = obj.benchmark_frame_data.gammaZ;
                        mydata.Delta0 = obj.Delta0_over_L*mydata.L;
                        mydata.delta0 = obj.delta0_over_L*mydata.L;
                    otherwise
                        error('Unknown frame_type: %s',mydata.frame_type)
                end                  
            case {'y','weak'}
                mydata.frame_type = obj.benchmark_frame_data.frame_typeY;
                switch mydata.frame_type
                    case 'Sidesway_Inhibited'
                        mydata.beta   = -1/obj.benchmark_frame_data.betaBotY_over_betaTopY;
                        mydata.delta0 = obj.delta0_over_L*mydata.L;
                    case 'Sidesway_Uninhibited'
                        mydata.kqtop  = obj.benchmark_frame_data.kqtopY;
                        mydata.kqbot  = obj.benchmark_frame_data.kqbotY;
                        mydata.gamma  = obj.benchmark_frame_data.gammaY;
                        mydata.Delta0 = obj.Delta0_over_L*mydata.L;
                        mydata.delta0 = obj.delta0_over_L*mydata.L;
                    otherwise
                        error('Unknown frame_type: %s',mydata.frame_type)
                end   
            otherwise
                error('Unknown axis: %s',axis)
        end
        
        
        
        % Section Definition
        switch type
            case 'elastic'
                elastic_data = struct;
                [E,A,I] = obj.benchmark_frame_data.section.sectionPropertiesForElasticAnalysis2d(axis,obj.ElasticStiffnessType);
                elastic_data.sectionType = 'elastic';
                elastic_data.E = E;
                elastic_data.A = A;
                elastic_data.I = I;
                
                section_def = FiberSectionDefinition(elastic_data,axis,...
                    obj.FiberSectionDefinition_SectionID,...
                    obj.FiberSectionDefinition_StartMatID,struct);
                
            case 'nonlinear'
                section_def = FiberSectionDefinition(mydata.section,axis,...
                    obj.FiberSectionDefinition_SectionID,...
                    obj.FiberSectionDefinition_StartMatID,...
                    obj.FiberSectionDefinitonOptions);
                
            otherwise
                error('Unknown type: %s',type)
        end        
        BA = BenchmarkAnalysis2d_OpenSees(mydata,section_def,obj.AnalysisOptions);
    end
end
end