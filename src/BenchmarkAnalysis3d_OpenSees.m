classdef BenchmarkAnalysis3d_OpenSees < OpenSeesAnalysis
    
    properties
        
        % Cross Section Properties
        section                                 = [];
        sectionDefinition                       = [];
        sectionTag                              = 1;
        
        % Frame Properties
        L                                       = [];
        frame_typeZ                             = [];
        frame_typeY                             = [];
        kqtopZ                                  = [];
        kqbotZ                                  = [];
        kqtopY                                  = [];
        kqbotY                                  = [];      
        gammaZ                                  = [];
        gammaY                                  = [];              
        betaTopY                                = [];
        betaTopZ                                = [];
        betaBotZ                                = [];
        betaBotY                                = [];                     
        Hy                                      = [];
        Hz                                      = [];        
        Delta0Z                                 = 0;
        delta0Z                                 = 0;
        Delta0Y                                 = 0;
        delta0Y                                 = 0;
        % theta0 - @todo eventually add a twist imperfection
        
        % Analysis Options
        eigenType                               = 'symmBandLapack'
        numElements                             = 6;
        numIntegrationPoints                    = 3;
        geomTransfType                          = 'Corotational';
        includeInitialGeometricImperfections    = true;
        eigenvalueLimitPointTolerance           = 1e-3;
        absoluteStrainLimit                     = 0.05;
        concreteCompressiveStrainLimit          = [];
        compressiveForceLimitType               = '';
        numStepsGravity                         = 100;
        numStepsLateral                         = 100;  
        limitPointWarning                       = false;
        runZeroAxialLoadAsSection               = true;
        reportMaximumLimitpoint                 = true; % if false report internal forces and deformations along the length
        extraResultsOutput                      = false;
        base_tolerance_force                    = 1e-3;
        controlled_dof_override                 = [];
        includeOpenSeesComposite                = true;
        convergence_test_output_flag            = 1
        
        % Leaning Column Load Properties
        leaning_column_E = 1.0;
        leaning_column_G = 1.0/(2*(1+0.3));
        leaning_column_A = 1.0e10;
        leaning_column_I = 1.0e10;
        leaning_column_J = 2.0e10;
    end
    
    methods
        function obj = BenchmarkAnalysis3d_OpenSees(data,sectionDefinition,analysisOptions)
            if nargin == 0
                % Called by superclass
                return
            end
            
            % Cross Section Properties
            obj.section = data.section;
            obj.sectionDefinition = sectionDefinition;

            % Frame Properties
            obj.L = data.L;
            obj.frame_typeZ = data.frame_typeZ;
            obj.frame_typeY = data.frame_typeY;

            % Initial geometric imperfections
            if isfield(data,'Delta0Z')
                obj.Delta0Z = data.Delta0Z;
            end
            if isfield(data,'Delta0Y')
                obj.Delta0Y = data.Delta0Y;
            end
            obj.delta0Z = data.delta0Z;
            obj.delta0Y = data.delta0Y;

            % Loading properties for bending about the Z axis 
            % (bending in the XY plane)
            switch obj.frame_typeZ
                case 'Sidesway_Inhibited'
                    obj.betaTopZ    = data.betaTopZ;
                    obj.betaBotZ    = data.betaBotZ;

                case 'Sidesway_Uninhibited'
                    obj.kqtopZ      = data.kqtopZ;
                    obj.kqbotZ      = data.kqbotZ;
                    obj.gammaZ      = data.gammaZ;
                    obj.Hy          = data.Hy;

                otherwise
                    error('Unknown frame_typeZ: %s',obj.frame_typeZ)
            end

            % Loading properties for bending about the Y axis 
            % (bending in the XZ plane)
            switch obj.frame_typeY
                case 'Sidesway_Inhibited'
                    obj.betaTopY    = data.betaTopY;
                    obj.betaBotY    = data.betaBotY;

                case 'Sidesway_Uninhibited'
                    obj.kqtopY      = data.kqtopY;
                    obj.kqbotY      = data.kqbotY;
                    obj.gammaY      = data.gammaY;
                    obj.Hz          = data.Hz;

                otherwise
                    error('Unknown frame_typeY: %s',obj.frame_typeY)
            end        
                
            % Analysis Options
            if nargin > 2
                names = fieldnames(analysisOptions);
                pnames = properties(obj);
                for i = 1:length(names)
                    if ismember(names{i},pnames)
                        obj.(names{i}) = analysisOptions.(names{i});
                    else
                        warning('Unknown option: %s',names{i})
                    end
                end
            end
        end
        
        function ang = loading_angle(obj,varargin)
            ang = atan(obj.M1y_over_X/obj.M1z_over_X);
            % @todo - as this is now, it will always return an angle
            % between 0 and 90 degrees. I think this is okay for now, but
            % we may want to correct it later.
        end
        function ang = imperfection_angle(obj)
            ang = atan(obj.delta0Z/obj.delta0Y);
        end
        function val = M1z_over_X(obj)
            switch obj.frame_typeZ
                case 'Sidesway_Inhibited'
                    val = max(abs(obj.betaBotZ),abs(obj.betaTopZ));
                case 'Sidesway_Uninhibited'
                    if obj.kqtopZ == 0 && obj.kqbotZ == 0
                        error('Structure unstable')
                    elseif obj.kqtopZ == 0 || obj.kqbotZ == 0
                        val = obj.Hy * obj.L ;
                    elseif obj.kqtopZ == obj.kqbotZ
                        val =  obj.Hy * obj.L/2 ;
                    else
                        error('M1z_over_X is not yet implemented for this condition');
                    end
                otherwise
                    error('Unknown frame_typeZ: %s',obj.frame_typeZ)
            end            
        end
        function val = M1y_over_X(obj)
            switch obj.frame_typeY
                case 'Sidesway_Inhibited'
                    val = max(abs(obj.betaTopY),abs(obj.betaBotY));
                case 'Sidesway_Uninhibited'
                    if obj.kqtopY == 0 && obj.kqbotY == 0
                        error('Structure unstable')
                    elseif obj.kqtopY == 0 || obj.kqbotY == 0
                        val = obj.Hz * obj.L ;
                    elseif obj.kqtopY == obj.kqbotY
                        val = obj.Hz * obj.L/2 ;
                    else
                        error('M1y_over_X is not yet implemented for this condition');
                    end
                otherwise
                    error('Unknown frame_typeY: %s',obj.frame_typeY)
            end            
        end
        
        function results = runAnalysis(obj,analysisType,X1,X2,tryNumber)
            % Compression load is negative

            if nargin < 5
                tryNumber = 1;
            end
            switch analysisType 
                case 'LimitPoint_Proportional'
                    X_over_P = X1; 
                    checkForLimitPoint = true;
                case 'LimitPoint_NonProportional'
                    P = X1;
                    checkForLimitPoint = true;
                    if P == 0 && obj.runZeroAxialLoadAsSection
                        results = obj.runSectionAnalysis(analysisType,X1,X2,tryNumber);
                        return
                    end
                case 'TargetForce_Proportional'
                    P = X1;
                    X = X2;
                    checkForLimitPoint = false;
                case  'TargetForce_NonProportional'
                    P = X1;
                    X = X2;
                    checkForLimitPoint = false;     
                % case 'TargetDisplacement_Proportional'            
                case 'TargetDrift_NonProportional'
                    P  = X1;
                    dX = X2; % @todo - how do we know what direction this is supposed to go in?
                    checkForLimitPoint = false; 
                otherwise
                    error('Unknown analysis type: %s',analysisType)
            end    
            
            % Filenames           
            inputFilename                     = obj.scratchFile('BenchmarkAnalysis3d_Driver.tcl');
            nodeDisplacementFilename          = obj.scratchFile('BenchmarkAnalysis3d_NodeDisp.out');
            elementForceFilename              = obj.scratchFile('BenchmarkAnalysis3d_EleForce.out');
            elementSectionDeformationFilename = obj.scratchFile('BenchmarkAnalysis3d_EleSectionDef.out');
            eigenFilename                     = obj.scratchFile('BenchmarkAnalysis3d_EigVals.out');

            %% Write Input File
            fid = fopen(inputFilename,'w');
            
            % Analysis Information
            fprintf(fid,'###  tryNumber %i ### \n',tryNumber);
            fprintf(fid,'###  Loading Angle - %g ###\n',loading_angle(obj));
            fprintf(fid,'###  Analysis Type - %s ###\n',analysisType);
            fprintf(fid,'###  frameTypeZ %s ### \n',obj.frame_typeZ);
            fprintf(fid,'###  frameTypeY %s ### \n',obj.frame_typeY);            
            
            % Initiate Model
            if obj.includeOpenSeesComposite
                fprintf(fid,'package require OpenSeesComposite \n');
                fprintf(fid,'namespace import OpenSeesComposite::* \n');
            end
            fprintf(fid,'model BasicBuilder -ndm 3 -ndf 6 \n');

            % Column Length
            fprintf(fid,'set columnSectionTag %i\n',obj.sectionTag);
            fprintf(fid,'set Lc %i \n',obj.L);
            
            % Information Specific to the Type of Frame
            if strcmp(obj.frame_typeZ,'Sidesway_Uninhibited')
                fprintf(fid,'set Hy     %g \n',obj.Hy);
                fprintf(fid,'set gammaZ %g \n',obj.gammaZ);                

            elseif strcmp(obj.frame_typeZ,'Sidesway_Inhibited')
                fprintf(fid,'set betaTopZ %g \n',obj.betaTopZ);
                fprintf(fid,'set betaBotZ %g \n',obj.betaBotZ);
                
            else
                error('Unknown frame type: %s',obj.frame_typeZ);
            end
            
            if strcmp(obj.frame_typeY,'Sidesway_Uninhibited')
                fprintf(fid,'set Hz     %g \n',obj.Hz);
                fprintf(fid,'set gammaY %g \n',obj.gammaY);

            elseif strcmp(obj.frame_typeY,'Sidesway_Inhibited')
                fprintf(fid,'set betaTopY %g \n',obj.betaTopY);
                fprintf(fid,'set betaBotY %g \n',obj.betaBotY);

            else
                error('Unknown frame type: %s',obj.frame_type);
            end            

            if (strcmp(obj.frame_typeZ,obj.frame_typeY)) 
                    loadStep = loading_steps(obj,obj.frame_typeZ,analysisType,tryNumber);
            elseif (strcmp(obj.frame_typeZ,'Sidesway_Uninhibited')) 

                if abs(imperfection_angle(obj)-pi/2)> 0.01
                    loadStep = loading_steps(obj,obj.frame_typeZ,analysisType,tryNumber);
                else
                    loadStep = loading_steps(obj,obj.frame_typeY,analysisType,tryNumber);
                end
            elseif (strcmp(obj.frame_typeY,'Sidesway_Uninhibited')) 

                if abs(imperfection_angle(obj)-0) > 0.01
                    loadStep = loading_steps(obj,obj.frame_typeY,analysisType,tryNumber);
                else
                    loadStep = loading_steps(obj,obj.frame_typeZ,analysisType,tryNumber);
                end
            else
                error('Not defined case')
            end
            
            % Information Specific to the Type of Loading
            switch analysisType
                case 'LimitPoint_Proportional'
                    fprintf(fid,'set coarseDispStepSize %g \n',loadStep.coarseDispStepSize);
                    fprintf(fid,'set fineDispStepSize %g \n',loadStep.fineDispStepSize);
                    fprintf(fid,'set maxDisp %g \n',loadStep.maxDisp);
                    fprintf(fid,'set maxAnalysisTime %i \n',loadStep.maxAnalysisTime);
                    fprintf(fid,'set baseForceTolerance %g \n',loadStep.baseForceTolerance);
                    fprintf(fid,'set baseDisplacementTolerance %g \n',loadStep.baseDisplacementTolerance);                    
                    fprintf(fid,'set X_over_P %g\n',X_over_P);

                case 'LimitPoint_NonProportional'
                    fprintf(fid,'set numStepsGravity %i\n',loadStep.numStepsGravity);
                    fprintf(fid,'set coarseDispStepSize %g \n',loadStep.coarseDispStepSize);
                    fprintf(fid,'set fineDispStepSize %g \n',loadStep.fineDispStepSize);
                    fprintf(fid,'set maxDisp %g \n',loadStep.maxDisp);
                    fprintf(fid,'set maxAnalysisTime %i \n',loadStep.maxAnalysisTime);
                    fprintf(fid,'set baseForceTolerance %g \n',loadStep.baseForceTolerance);
                    fprintf(fid,'set baseDisplacementTolerance %g \n',loadStep.baseDisplacementTolerance);                
                    fprintf(fid,'set P %g\n',P);     

                case 'TargetForce_Proportional'
                    fprintf(fid,'set numStepsGravity %i\n',loadStep.numStepsGravity);
                    fprintf(fid,'set numStepsLateral %i\n',loadStep.numStepsLateral);
                    fprintf(fid,'set baseForceTolerance %g \n',loadStep.baseForceTolerance);
                    fprintf(fid,'set baseDisplacementTolerance %g \n',loadStep.baseDisplacementTolerance);                                          
                    fprintf(fid,'set P %g\n',P);
                    fprintf(fid,'set X %g\n',X);
                
                case 'TargetForce_NonProportional'    
                    fprintf(fid,'set numStepsGravity %i\n',loadStep.numStepsGravity);
                    fprintf(fid,'set numStepsLateral %i\n',loadStep.numStepsLateral);
                    fprintf(fid,'set baseForceTolerance %g \n',loadStep.baseForceTolerance);
                    fprintf(fid,'set baseDisplacementTolerance %g \n',loadStep.baseDisplacementTolerance);                   
                    fprintf(fid,'set P %g\n',P);
                    fprintf(fid,'set X %g\n',X); 
                
                case 'TargetDrift_NonProportional'    
                    fprintf(fid,'set numStepsGravity %i\n',loadStep.numStepsGravity);
                    fprintf(fid,'set numStepsLateral %i\n',loadStep.numStepsLateral);
                    fprintf(fid,'set baseForceTolerance %g \n',loadStep.baseForceTolerance);
                    fprintf(fid,'set baseDisplacementTolerance %g \n',loadStep.baseDisplacementTolerance);
                    
                    fprintf(fid,'set P %g\n',P);
                    fprintf(fid,'set dX %g\n',dX);
                    
                otherwise
                    error('Unknown analysis type: %s',analysisType)
            end
            
            % Miscellaneous Analysis Parameters
            fprintf(fid,'set numEles %i\n',obj.numElements);
            fprintf(fid,'set midNode %i\n',0.5*obj.numElements+1); % @todo - we should add a check that numElements is an even number
            fprintf(fid,'set numNodes %i\n',obj.numElements+1);
            fprintf(fid,'set numIP %i\n',obj.numIntegrationPoints);
            fprintf(fid,'set startTime [clock seconds] \n');
            fprintf(fid,'set testOutputFlag %i \n',obj.convergence_test_output_flag);
            
            % Node and boundry conditions
            fprintf(fid,'\n# ############ BUILD MODEL ############ \n');
            fprintf(fid,'# Define nodes and constraints for column \n');
            for i=0:obj.numElements
                x = i/obj.numElements;
                if obj.includeInitialGeometricImperfections
                    imperfy = obj.delta0Y*sin(pi*x)+obj.Delta0Y*x;
                    imperfz = obj.delta0Z*sin(pi*x)+obj.Delta0Z*x;
                else
                    imperfy = 0;
                    imperfz = 0;
                end
                fprintf(fid,'node %i     %.04f      %.04f      %.04f \n',i+1,x*obj.L,imperfy,imperfz);
            end
            if obj.includeInitialGeometricImperfections
                imperfyt = obj.Delta0Y;
                imperfzt = obj.Delta0Z;
            else
                imperfyt = 0;
                imperfzt = 0;
            end
                
            for i=0:obj.numElements
                fprintf(fid,'mass %i   1.0   1.0   1.0   1.0   1.0   1.0 \n',i+1);
            end
            
            dof_bot   = [1 1 1 1 0 0]; % Initial fixity of the bottom node
            dof_top   = [0 0 0 0 0 0]; % Initial fixity of the top node
            
            fprintf(fid,'\n# ############ DEFINE BOUNDARY CONDITIONS ############ \n');
            if strcmp(obj.frame_typeZ,'Sidesway_Uninhibited')
                
                if obj.kqbotZ == 0 
                    % Free - nothing to do
                elseif obj.kqbotZ == Inf
                    % Fixed - change fixity
                    dof_bot(6) = 1;
                else
                    % Rotational Spring
                    fprintf(fid,' \n# Bottom Rotational Spring - Bending about the Z axis (in the XY plane) \n');
                    fprintf(fid,'node 200  0 0 0 \n');
                    fprintf(fid,'fix  200  1 1 1 1 1 1 \n');
                    fprintf(fid,'uniaxialMaterial Elastic 200 %g \n',obj.kqbotZ);
                    fprintf(fid,'element zeroLength 200 1 200 -mat 200 -dir 6 \n');                      
                end
                    
                if obj.kqtopZ == 0 
                    % Free - nothing to do
                elseif obj.kqtopZ == Inf
                    % Fixed - change fixity
                    dof_top(6) = 1;
                else
                    % Rotational Spring
                    fprintf(fid,' \n# Top Rotational Spring - Bending about the Z axis (in the XY plane) \n');
                    fprintf(fid,'node 201  $Lc %g %g \n',imperfyt,imperfzt);
                    fprintf(fid,'fix 201  1 1 1 1 1 1 \n');
                    fprintf(fid,'uniaxialMaterial Elastic 201 %g \n',obj.kqtopZ);
                    fprintf(fid,'element zeroLength 201 $numNodes 201 -mat 201 -dir 6 \n');                     
                end                
                
                if obj.gammaZ ~=0                 
                    fprintf(fid,'\n# Leaning Column - Bending about the Z axis (in the XY plane) \n');
                    fprintf(fid,'node 103  $Lc %g %g \n',imperfyt,imperfzt);
                    fprintf(fid,'node 102  0 0 0 \n');
                    fprintf(fid,'fix  103  0 0 1 0 0 0 \n');                    
                    fprintf(fid,'fix  102  1 1 1 1 0 0 \n');
                    fprintf(fid,'equalDOF $numNodes 103 2 \n');
                    fprintf(fid,'set geomTransfTagLeaningZ 3 \n');                                                    
                    fprintf(fid,'geomTransf %s $geomTransfTagLeaningZ 0 0 -1 \n',obj.geomTransfType);
                    fprintf(fid,'element elasticBeamColumn 102 102 103 %g %g %g %g %g %g $geomTransfTagLeaningZ \n',...
                        obj.leaning_column_A,obj.leaning_column_E,obj.leaning_column_G,...
                        obj.leaning_column_J,obj.leaning_column_I,obj.leaning_column_I);
                end
                
            elseif strcmp(obj.frame_typeZ,'Sidesway_Inhibited')
                dof_top(2) = 1;                    
            
            else
                error('frame_typeZ: %s not recognized',obj.frame_typeZ);
            end
            
            if strcmp(obj.frame_typeY,'Sidesway_Uninhibited')
                
                if obj.kqbotY == 0 
                    % Free - nothing to do
                elseif obj.kqbotY == Inf
                    % Fixed - change fixity
                    dof_bot(5) = 1;
                else
                    % Rotational Spring
                    fprintf(fid,' \n# Bottom Rotational Spring - Bending about the Y axis (in the XZ plane) \n');
                    fprintf(fid,'node 300  0 0 0 \n');
                    fprintf(fid,'fix  300  1 1 1 1 1 1 \n');
                    fprintf(fid,'uniaxialMaterial Elastic 300 %g \n',obj.kqbotY);
                    fprintf(fid,'element zeroLength 300 1 300 -mat 300 -dir 5 \n');
                end                
                
                if obj.kqtopY == 0 
                    % Free - nothing to do
                elseif obj.kqtopY == Inf
                    % Fixed - change fixity
                    dof_top(5) = 1;
                else
                    % Rotational Spring
                    fprintf(fid,' \n# Top Rotational Spring - Bending about the Y axis (in the XZ plane) \n');
                    fprintf(fid,'node 301  $Lc %g %g \n',imperfyt,imperfzt);
                    fprintf(fid,'fix 301  1 1 1 1 1 1 \n');
                    fprintf(fid,'uniaxialMaterial Elastic 301 %g \n',obj.kqtopY);                        
                    fprintf(fid,'element zeroLength 301 $numNodes 301 -mat 301 -dir 5 \n');                                         
                end                     
                
                if obj.gammaY ~=0
                    fprintf(fid,'\n# Leaning Column - Bending about the Y axis (in the XZ plane) \n');
                    fprintf(fid,'node 101  $Lc %g %g \n',imperfyt ,imperfzt);                        
                    fprintf(fid,'node 100  0 0 0 \n');
                    fprintf(fid,'fix  101  0 1 0 0 0 0 \n');                    
                    fprintf(fid,'fix  100  1 1 1 1 0 0 \n');
                    fprintf(fid,'equalDOF $numNodes 101 3 \n');
                    fprintf(fid,'set geomTransfTagLeaningY 2 \n');
                    fprintf(fid,'geomTransf %s $geomTransfTagLeaningY 0 0 -1 \n',obj.geomTransfType );                                                
                    fprintf(fid,'element elasticBeamColumn 101 100 101 %g %g %g %g %g %g $geomTransfTagLeaningY \n',...
                        obj.leaning_column_A,obj.leaning_column_E,obj.leaning_column_G,...
                        obj.leaning_column_J,obj.leaning_column_I,obj.leaning_column_I);
                end
                
            elseif strcmp(obj.frame_typeY,'Sidesway_Inhibited')
                dof_top(3) = 1;
            
            else
                error('frameType : %s not recognized',obj.frame_typeY);
            end
            
            fprintf(fid,'\n# Define constraints for bot and top nodes\n');
            fprintf(fid,'fix  1          %i  %i  %i  %i  %i  %i \n',...
                dof_bot(1),dof_bot(2),dof_bot(3),dof_bot(4),dof_bot(5),dof_bot(6));
            if any(dof_top)
                fprintf(fid,'fix  $numNodes  %i  %i  %i  %i  %i  %i \n',...
                    dof_top(1),dof_top(2),dof_top(3),dof_top(4),dof_top(5),dof_top(6));
            end
            
            fprintf(fid,'\n# ############ DEFINE PRIMARY BEAM-COLUMN ############ \n');
            
            fprintf(fid,'# Define Geometric Transformation\n');
            fprintf(fid,'set geomTransfTag 1 \n');                    
            fprintf(fid,'geomTransf %s $geomTransfTag 0.0 0.0 -1.0 \n',obj.geomTransfType);

            % Define the Section
            fprintf(fid,'\n# Define Section for the Beam-Column\n'); 
            if iscell(obj.sectionDefinition)
                for i = 1:length(obj.sectionDefinition)
                    fprintf(fid,'%s \n',obj.sectionDefinition{i});
                end
            else
                fprintf(fid,'%s \n',obj.sectionDefinition);
            end
            fprintf(fid,'\n# Define elements for column\n'); 
            for i =1:obj.numElements
                if tryNumber <=2
                    fprintf(fid,'element mixedBeamColumn3d %i %i %i $numIP $columnSectionTag $geomTransfTag \n',i,i,i+1);
                else
                    %fprintf(fid,'element mixedBeamColumn3d %i Xi %i $numIP $columnSectionTag $geomTransfTag -geomLinear \n',i,i,i+1);
                    fprintf(fid,'element mixedBeamColumn3d %i %i %i $numIP $columnSectionTag $geomTransfTag \n',i,i,i+1);
                end    
            end

            % Recorders
            fprintf(fid,'\n# ############ DEFINE RECORDERS ############ \n');
            fprintf(fid,'recorder Node -file %s -nodeRange 1 $numNodes -dof 1 2 3 4 5 6 disp \n',nodeDisplacementFilename);
            fprintf(fid,'recorder Element -file %s -eleRange 1 $numEles localForce \n',elementForceFilename);
            fprintf(fid,'recorder Element -file %s -eleRange 1 $numEles sectionDeformation_Force \n',elementSectionDeformationFilename);
            fprintf(fid,'set eigenFileId [open %s w] \n',eigenFilename);

            fprintf(fid,'\n# ############ LOADING AND ANALYSIS OPTIONS ############ \n');
            fprintf(fid,'# Analysis Options \n');         
            fprintf(fid,'system UmfPack \n');
            fprintf(fid,'constraints Transformation \n');
            fprintf(fid,'test NormUnbalance $baseForceTolerance 20 $testOutputFlag \n');                    
            fprintf(fid,'algorithm  Newton \n');
            fprintf(fid,'numberer Plain \n');

            % One Step with Nothing
            fprintf(fid,'\n# One step with nothing \n');
            fprintf(fid,'integrator LoadControl 0 \n');
            fprintf(fid,'analysis Static \n');
            fprintf(fid,'set ok [analyze 1] \n');
            fprintf(fid,'if {$ok == 0} {\n');
            fprintf(fid,'  %s \n',obj.lowestEigenValue_command);
            fprintf(fid,'  set firstEigenValue $lowestEigenValue \n');
            fprintf(fid,'  puts $eigenFileId "[getTime] $lowestEigenValue"  \n');
            fprintf(fid,'}\n');

            if strcmp(analysisType,'LimitPoint_Proportional')

                fprintf(fid,'\n# ############ DISPLACEMENT CONTROL / AXIAL ONLY LOADING ############ \n');                
                fprintf(fid,'# Apply gravity loads \n');                            
                
                
                if (strcmp(obj.frame_typeZ,obj.frame_typeY)) 
                    frame_type = obj.frame_typeZ;

                    if imperfection_angle(obj) <= pi/4
                        if isempty(obj.controlled_dof_override); controlled_dof=2; else; controlled_dof=1;end
                    else
                        if isempty(obj.controlled_dof_override); controlled_dof=3; else; controlled_dof=1;end
                    end
                    if (strcmp(frame_type,'Sidesway_Inhibited'))
                        if (abs(obj.betaBotY) > abs(obj.betaTopY))
                            fprintf(fid,'set controlledNode [expr $midNode-1] \n');
                        elseif (abs(obj.betaBotY) < abs(obj.betaTopY))
                            fprintf(fid,'set controlledNode [expr $midNode+1] \n'); 
                        else
                            fprintf(fid,'set controlledNode [expr $midNode] \n');                        
                        end
                    else
                        fprintf(fid,'set controlledNode $numNodes \n');                
                    end

                elseif (strcmp(obj.frame_typeZ,'Sidesway_Uninhibited')) 

                    if abs(imperfection_angle(obj)-pi/2)> 0.01
                        if isempty(obj.controlled_dof_override); controlled_dof=2; else; controlled_dof=1;end
                        fprintf(fid,'set controlledNode $numNodes \n');   
                    else
                        if isempty(obj.controlled_dof_override); controlled_dof=3; else; controlled_dof=1;end
                        if (abs(obj.betaBotY) > abs(obj.betaTopY))
                            fprintf(fid,'set controlledNode [expr $midNode-1] \n');
                        elseif (abs(obj.betaBotY) < abs(obj.betaTopY))
                            fprintf(fid,'set controlledNode [expr $midNode+1] \n'); 
                        else
                            fprintf(fid,'set controlledNode [expr $midNode] \n');                        
                        end
                    end
                elseif (strcmp(obj.frame_typeY,'Sidesway_Uninhibited')) 

                    if abs(imperfection_angle(obj)-0) > 0.01
                        if isempty(obj.controlled_dof_override); controlled_dof=3; else; controlled_dof=1;end
                        fprintf(fid,'set controlledNode $numNodes \n');   
                    else
                        if isempty(obj.controlled_dof_override); controlled_dof=2; else; controlled_dof=1;end
                        if (abs(obj.betaBotZ) > abs(obj.betaTopZ))
                            fprintf(fid,'set controlledNode [expr $midNode-1] \n');                            
                        elseif (abs(obj.betaBotZ) < abs(obj.betaTopZ))
                            fprintf(fid,'set controlledNode [expr $midNode+1] \n'); 
                        else
                            fprintf(fid,'set controlledNode [expr $midNode] \n');                        
                        end
                    end
                else
                    error('Not defined case')
                end
                fprintf(fid,'pattern Plain 1 Linear { \n');
                if strcmp(obj.frame_typeZ,'Sidesway_Inhibited') &&...
                   strcmp(obj.frame_typeY,'Sidesway_Inhibited')

                    fprintf(fid,'  load $numNodes -1.0 0.0 0.0 0.0 [expr $X_over_P*$betaTopY]  [expr $X_over_P*$betaTopZ] \n');
                    fprintf(fid,'  load 1          0.0 0.0 0.0 0.0 [expr $X_over_P*$betaBotY]  [expr $X_over_P*$betaBotZ] \n');                

                elseif strcmp(obj.frame_typeY,'Sidesway_Uninhibited')&&...
                   strcmp(obj.frame_typeZ,'Sidesway_Inhibited')
                    fprintf(fid,'  load $numNodes -1.0 0.0 [expr $X_over_P*$Hz] 0.0 0.0 [expr $X_over_P*$betaTopZ] \n'); 
                    fprintf(fid,'  load 1          0.0 0.0 0.0                  0.0 0.0 [expr $X_over_P*$betaBotZ]  \n');                

                    if obj.gammaY ~= 0 
                        fprintf(fid,'  load  101 [expr -1*$gammaY] 0.0 0.0 0.0 0.0 0.0 \n');  
                    end                    
                elseif strcmp(obj.frame_typeZ,'Sidesway_Uninhibited')&&...
                   strcmp(obj.frame_typeY,'Sidesway_Inhibited')
                    fprintf(fid,'  load $numNodes -1.0 [expr $X_over_P*$Hy] 0.0 0.0 [expr $X_over_P*$betaTopY] 0.0 \n');
                    fprintf(fid,'  load 1          0.0 0.0                  0.0 0.0 [expr $X_over_P*$betaBotY] 0.0 \n');      
                    if obj.gammaZ ~= 0
                        fprintf(fid,'  load  103 [expr -1*$gammaZ] 0.0 0.0 0.0 0.0 0.0 \n');  
                    end
                    
                elseif strcmp(obj.frame_typeY,'Sidesway_Uninhibited')&&...
                   strcmp(obj.frame_typeZ,'Sidesway_Uninhibited')
                    fprintf(fid,'  load $numNodes  -1.0 [expr $X_over_P*$Hy] [expr $X_over_P*$Hz] 0.0 0.0 0.0 \n');
                    
                    if obj.gammaY ~= 0 
                        fprintf(fid,'  load  101 [expr -1*$gammaY] 0.0 0.0 0.0 0.0 0.0 \n');  
                    end                
                    if obj.gammaZ ~= 0
                        fprintf(fid,'  load  103 [expr -1*$gammaZ] 0.0 0.0 0.0 0.0 0.0 \n');  
                    end
                    
                else
                    error('frameType not recgonized');
                end
                
                fprintf(fid,'} \n');
                fprintf(fid,'set ok 0 \n');
                fprintf(fid,'set iStep 0 \n');             
              
                fprintf(fid,'set dispStepSize $coarseDispStepSize \n');                                       
                fprintf(fid,'while {$ok == 0} { \n');
                fprintf(fid,'    set previousTime [getTime] \n');                
                fprintf(fid,'    test NormUnbalance $baseForceTolerance 10 $testOutputFlag \n');                                      
                fprintf(fid,'    integrator DisplacementControl $controlledNode %i $dispStepSize \n',controlled_dof);                                    
                fprintf(fid,'    algorithm  Newton \n');
                fprintf(fid,'    set ok [analyze 1] \n');
                fprintf(fid,'    if {$ok != 0} {  \n');                    
                fprintf(fid,'       integrator DisplacementControl $controlledNode %i  [expr 0.1*$dispStepSize] \n',controlled_dof);
                fprintf(fid,'       test NormUnbalance $baseForceTolerance 10 $testOutputFlag \n');                                      
                fprintf(fid,'       algorithm Newton \n');
                fprintf(fid,'       set ok [analyze 1] \n');
                fprintf(fid,'       } \n');
                fprintf(fid,'    if {$ok != 0} { \n');
                fprintf(fid,'       integrator DisplacementControl $controlledNode %i [expr 0.01*$dispStepSize] \n',controlled_dof);
                fprintf(fid,'       test NormUnbalance $baseForceTolerance 10 $testOutputFlag \n');                                      
                fprintf(fid,'       algorithm Newton \n');
                fprintf(fid,'       set ok [analyze 1] \n');
                fprintf(fid,'       } \n');
                fprintf(fid,'    if {$ok != 0} { \n');
                fprintf(fid,'       integrator DisplacementControl $controlledNode %i [expr 0.1*$dispStepSize] \n',controlled_dof);
                fprintf(fid,'       test NormDispIncr $baseDisplacementTolerance 10 $testOutputFlag \n');                    
                fprintf(fid,'       algorithm  KrylovNewton \n');
                fprintf(fid,'       set ok [analyze 1] \n');
                fprintf(fid,'       } \n');
                fprintf(fid,'    if {$ok != 0} { \n');
                fprintf(fid,'       integrator DisplacementControl $controlledNode %i [expr 0.01*$dispStepSize] \n',controlled_dof);
                fprintf(fid,'       test NormDispIncr $baseDisplacementTolerance 10 $testOutputFlag \n');                    
                fprintf(fid,'       algorithm Newton \n');
                fprintf(fid,'       set ok [analyze 1] \n');
                fprintf(fid,'       } \n');
                fprintf(fid,'    if {$ok == 0} { \n');
                fprintf(fid,'       %s \n',obj.lowestEigenValue_command);
                fprintf(fid,'       puts $eigenFileId "[getTime] $lowestEigenValue" \n');
                fprintf(fid,'       set currentTime [getTime] \n');
                fprintf(fid,'       if {[expr $currentTime - $previousTime] < 0.0} {exit 1} \n');  
                fprintf(fid,'       if {$lowestEigenValue < 0.0} { exit 1 } \n');
                fprintf(fid,'       # Analysis completed sucessfully \n'); 
                fprintf(fid,'       if { [nodeDisp $controlledNode 1] > $maxDisp } { \n');
                fprintf(fid,'       # Maximum deformation reached \n');
                fprintf(fid,'       exit 7 \n');
                fprintf(fid,'       } \n');                    
                fprintf(fid,'       incr iStep \n');     
                fprintf(fid,'       if { [expr $lowestEigenValue/$firstEigenValue] > 0.25 } { \n');
                fprintf(fid,'          set dispStepSize $coarseDispStepSize \n');
                fprintf(fid,'          } else { \n');
                fprintf(fid,'          set dispStepSize $fineDispStepSize \n');
                fprintf(fid,'          } \n');
                fprintf(fid,'       } else { \n');
                fprintf(fid,'       exit 4 \n');
                fprintf(fid,'       } \n');
                fprintf(fid,'   if { [expr [clock seconds]-$startTime] > $maxAnalysisTime } { \n');
                fprintf(fid,'   # Reached Maximum Time \n');
                fprintf(fid,'       exit 8 \n');
                fprintf(fid,'   } \n');              
                fprintf(fid,'} \n');
                fprintf(fid,'# Step limit reached \n');
                fprintf(fid,'exit 2 \n');        

            elseif strcmp(analysisType,'LimitPoint_NonProportional')

                fprintf(fid,'\n# ############ DISPLACEMENT CONTROL / NON-PROPORTIONAL LOADING ############ \n');            
                % Axial Loading
                if (strcmp(obj.frame_typeZ,obj.frame_typeY)) 
                    frame_type = obj.frame_typeZ;
                    if loading_angle(obj) <= pi/4
                        controlled_dof=2;
                    else
                        controlled_dof=3; 
                    end
                    if (strcmp(frame_type,'Sidesway_Inhibited'))
                        if (abs(obj.betaBotY) > abs(obj.betaTopY))
                            fprintf(fid,'set controlledNode [expr $midNode-1] \n');
                        elseif (abs(obj.betaBotY) < abs(obj.betaTopY))
                            fprintf(fid,'set controlledNode [expr $midNode+1] \n'); 
                        else
                            fprintf(fid,'set controlledNode [expr $midNode] \n');                        
                        end
                    else
                        fprintf(fid,'set controlledNode $numNodes \n');                
                    end

                elseif (strcmp(obj.frame_typeZ,'Sidesway_Uninhibited')) 

                    if abs(loading_angle(obj)-pi/2) > 0.01
                        controlled_dof=2;
                        fprintf(fid,'set controlledNode $numNodes \n');   
                    else
                        controlled_dof=3;
                        if (abs(obj.betaBotZ) > abs(obj.betaTopZ))
                            fprintf(fid,'set controlledNode [expr $midNode-1] \n');
                        elseif (abs(obj.betaBotZ) < abs(obj.betaTopZ))
                            fprintf(fid,'set controlledNode [expr $midNode+1] \n'); 
                        else
                            fprintf(fid,'set controlledNode [expr $midNode] \n');                        
                        end
                    end
                elseif (strcmp(obj.frame_typeY,'Sidesway_Uninhibited')) 

                    if (loading_angle(obj)-0) > 0.01
                        controlled_dof=3;
                        fprintf(fid,'set controlledNode $numNodes \n');   
                    else
                        controlled_dof=2;
                        if (abs(obj.betaBotY) > abs(obj.betaTopY))
                            fprintf(fid,'set controlledNode [expr $midNode-1] \n');
                        elseif (abs(obj.betaBotY) < abs(obj.betaTopY))
                            fprintf(fid,'set controlledNode [expr $midNode+1] \n'); 
                        else
                            fprintf(fid,'set controlledNode [expr $midNode] \n');                        
                        end
                    end
                else
                    error('Not defined case')
                end 
                if P ~= 0     
                    fprintf(fid,'\n# Apply gravity loads \n');
                    fprintf(fid,'pattern Plain 1 Linear { \n');
                    fprintf(fid,'  load $numNodes $P 0.0 0.0 0.0 0.0 0.0 \n');
                    if strcmp(obj.frame_typeZ,'Sidesway_Uninhibited')
                        if obj.gammaZ ~= 0 
                            fprintf(fid,'  load  103 [expr $P*$gammaZ] 0.0 0.0 0.0 0.0 0.0 \n');  
                        end
                    end
                    if strcmp(obj.frame_typeY,'Sidesway_Uninhibited')
                        if obj.gammaY ~= 0 
                            fprintf(fid,'  load  101 [expr $P*$gammaY] 0.0 0.0 0.0 0.0 0.0 \n');  
                        end
                    end
                    fprintf(fid,'} \n');               

                    fprintf(fid,'set ok 0 \n');
                    fprintf(fid,'set iStep 0 \n');                     
                    fprintf(fid,'integrator LoadControl [expr 1.0/$numStepsGravity] \n');
                    fprintf(fid,'test NormDispIncr $baseDisplacementTolerance 10 $testOutputFlag \n');                
                    fprintf(fid,'analysis Static \n');
                    fprintf(fid,'while {$numStepsGravity > $iStep && $ok == 0} { \n');
                    fprintf(fid,'  set ok [analyze 1] \n');                    
                    fprintf(fid,'  if {$ok != 0} { \n');
                    fprintf(fid,'     test NormDispIncr [expr 100*$baseDisplacementTolerance] 10 $testOutputFlag \n');                                      
                    fprintf(fid,'     set ok [analyze 1] \n');
                    fprintf(fid,'     } \n');
                    fprintf(fid,'  if {$ok != 0} { \n');
                    fprintf(fid,'     test NormDispIncr [expr 1000*$baseDisplacementTolerance] 10 $testOutputFlag \n');                                      
                    fprintf(fid,'     set ok [analyze 1] \n');
                    fprintf(fid,'     } \n');                    
                    fprintf(fid,'  if {$ok == 0} {  \n');
                    fprintf(fid,'     %s \n',obj.lowestEigenValue_command);
                    fprintf(fid,'     puts $eigenFileId "[getTime] $lowestEigenValue"  \n');
                    fprintf(fid,'     incr iStep  \n');
                    fprintf(fid,'     } \n');
                    fprintf(fid,'} \n');

                    fprintf(fid,'if {$ok != 0} {  \n');
                    fprintf(fid,'   # Analysis failed in gravity loading \n');
                    fprintf(fid,'   exit 3 \n');                
                    fprintf(fid,'   } \n');
                    fprintf(fid,'loadConst -time 0.0 \n');
                end

                % Lateral Loading
                fprintf(fid,'\n# Apply lateral loads \n');

                fprintf(fid,'pattern Plain 2 Linear { \n');                                             
                if strcmp(obj.frame_typeZ,'Sidesway_Inhibited')&&...
                    strcmp(obj.frame_typeY,'Sidesway_Inhibited')
                        
                    fprintf(fid,'  load $numNodes  0.0 0.0 0.0 0.0 $betaTopY $betaTopZ  \n');
                    fprintf(fid,'  load 1          0.0 0.0 0.0 0.0 $betaBotY $betaBotZ  \n');        

                elseif strcmp(obj.frame_typeY,'Sidesway_Uninhibited')&&...
                    strcmp(obj.frame_typeZ,'Sidesway_Inhibited')
                    fprintf(fid,'  load $numNodes  0.0 0.0 $Hz 0.0 0.0 $betaTopZ \n');
                    fprintf(fid,'  load 1          0.0 0.0 0.0 0.0 0.0 $betaBotZ  \n');
                    
                elseif strcmp(obj.frame_typeZ,'Sidesway_Uninhibited')&&...
                    strcmp(obj.frame_typeY,'Sidesway_Inhibited')
                    fprintf(fid,'  load $numNodes  0.0 $Hy 0.0 0.0 $betaTopY 0.0  \n');
                    fprintf(fid,'  load 1          0.0 0.0 0.0 0.0 $betaBotY 0.0   \n');
                elseif strcmp(obj.frame_typeY,'Sidesway_Uninhibited')&&...
                    strcmp(obj.frame_typeZ,'Sidesway_Uninhibited')
                    fprintf(fid,'  load $numNodes  0.0 $Hy $Hz 0.0 0.0 0.0  \n');
                else
                    error('frameType not recgonized');
                end
               
                
                fprintf(fid,'} \n');

                fprintf(fid,'set ok 0 \n');
                fprintf(fid,'set iStep 0 \n');            
                fprintf(fid,'set controlledDOF %i \n',controlled_dof);  
                fprintf(fid,'set dispStepSize $coarseDispStepSize \n');                                                
                fprintf(fid,'while {$ok == 0} { \n');
                fprintf(fid,'    set previousTime [getTime] \n'); 
                fprintf(fid,'    test NormDispIncr $baseDisplacementTolerance 10 $testOutputFlag \n');                    
                fprintf(fid,'    algorithm  Newton \n');
                fprintf(fid,'    integrator DisplacementControl $controlledNode %i [expr $dispStepSize] \n',controlled_dof);                                    
                fprintf(fid,'    set ok [analyze 1] \n');
                fprintf(fid,'    if {$ok != 0} {  \n');                    
                fprintf(fid,'       integrator DisplacementControl $controlledNode %i  [expr 0.1*$dispStepSize] \n',controlled_dof);
                fprintf(fid,'       test NormDispIncr $baseDisplacementTolerance 10 $testOutputFlag \n');                    
                fprintf(fid,'       algorithm  Newton \n');
                fprintf(fid,'       set ok [analyze 1] \n');
                fprintf(fid,'       } \n');
                fprintf(fid,'    if {$ok != 0} { \n');
                fprintf(fid,'       integrator DisplacementControl $controlledNode %i [expr 0.01*$dispStepSize] \n',controlled_dof);
                fprintf(fid,'       test NormDispIncr $baseDisplacementTolerance 10 $testOutputFlag \n');                    
                fprintf(fid,'       algorithm  Newton \n');
                fprintf(fid,'       set ok [analyze 1] \n');
                fprintf(fid,'       } \n');
                fprintf(fid,'    if {$ok != 0} { \n');
                fprintf(fid,'       integrator DisplacementControl $controlledNode %i [expr 0.001*$dispStepSize] \n',controlled_dof);
                fprintf(fid,'       test NormDispIncr $baseDisplacementTolerance 10 $testOutputFlag \n');                    
                fprintf(fid,'       algorithm  Newton \n');
                fprintf(fid,'       set ok [analyze 1] \n');
                fprintf(fid,'       } \n');
                fprintf(fid,'    if {$ok != 0} { \n');
                fprintf(fid,'       integrator DisplacementControl $controlledNode %i [expr 0.001*$dispStepSize] \n',controlled_dof);
                fprintf(fid,'       test NormDispIncr [expr 1e1*$baseDisplacementTolerance] 10 $testOutputFlag \n');                    
                fprintf(fid,'       algorithm  Newton \n');
                fprintf(fid,'       set ok [analyze 1] \n');
                fprintf(fid,'       } \n');
                fprintf(fid,'    if {$ok != 0} {  \n');                    
                fprintf(fid,'       integrator DisplacementControl $controlledNode %i  [expr 0.001*$dispStepSize] \n',controlled_dof);
                fprintf(fid,'       test NormDispIncr [expr 1e2*$baseDisplacementTolerance] 10 $testOutputFlag \n');                    
                fprintf(fid,'       algorithm  Newton \n');
                fprintf(fid,'       set ok [analyze 1] \n');
                fprintf(fid,'       } \n');
                fprintf(fid,'    if {$ok != 0} { \n');
                fprintf(fid,'       integrator DisplacementControl $controlledNode %i [expr 0.0001*$dispStepSize] \n',controlled_dof);
                fprintf(fid,'       test NormDispIncr [expr 1e1*$baseDisplacementTolerance] 10 $testOutputFlag \n');                    
                fprintf(fid,'       algorithm  Newton \n');
                fprintf(fid,'       set ok [analyze 1] \n');
                fprintf(fid,'       } \n');
                fprintf(fid,'    if {$ok != 0} { \n');
                fprintf(fid,'       integrator DisplacementControl $controlledNode %i [expr 0.0001*$dispStepSize] \n',controlled_dof);
                fprintf(fid,'       test NormDispIncr [expr 1e2*$baseDisplacementTolerance] 10 $testOutputFlag \n');                    
                fprintf(fid,'       algorithm  Newton \n');
                fprintf(fid,'       set ok [analyze 1] \n');
                fprintf(fid,'       } \n');                                              
                fprintf(fid,'    if {$ok == 0} { \n');
                fprintf(fid,'       %s \n',obj.lowestEigenValue_command);
                fprintf(fid,'       puts $eigenFileId "[getTime] $lowestEigenValue" \n');
                fprintf(fid,'       set currentTime [getTime] \n');
                fprintf(fid,'       if {[expr $currentTime - $previousTime] < 0.0} {exit 1} \n');                                                            
                fprintf(fid,'       if {$lowestEigenValue < 0.0} { exit 1 } \n');
                fprintf(fid,'       # Analysis completed sucessfully \n'); 
                fprintf(fid,'       if { [expr abs([nodeDisp $controlledNode $controlledDOF])] > $maxDisp } { \n');
                fprintf(fid,'          # Maximum deformation reached \n');
                fprintf(fid,'          exit 7 \n');
                fprintf(fid,'          } \n');                       
                fprintf(fid,'       incr iStep \n');                         
                fprintf(fid,'       if { [expr $lowestEigenValue/$firstEigenValue] > 0.25 } { \n');
                fprintf(fid,'          set dispStepSize $coarseDispStepSize \n');
                fprintf(fid,'          } else { \n');
                fprintf(fid,'          set dispStepSize $fineDispStepSize \n');                   
                fprintf(fid,'          } \n');
                fprintf(fid,'       } else { \n');
                fprintf(fid,'       exit 4 \n');
                fprintf(fid,'       } \n');
                fprintf(fid,'       if { [expr [clock seconds]-$startTime] > $maxAnalysisTime } { \n');
                fprintf(fid,'           # Reached Maximum Time \n');
                fprintf(fid,'           exit 8 \n');
                fprintf(fid,'          } \n');             

                fprintf(fid,'} \n');
                fprintf(fid,'# Step limit reached \n');            
                fprintf(fid,'exit 2 \n');              

            elseif strcmp(analysisType,'TargetForce_Proportional')
                fprintf(fid,'\n# ############ LOAD CONTROL / PROPORTIONAL LOADING ############ \n');
                fprintf(fid,'# Apply gravity loads \n');                            
                fprintf(fid,'pattern Plain 1 Linear { \n');
                if strcmp(obj.frame_typeZ,'Sidesway_Inhibited') &&...
                   strcmp(obj.frame_typeY,'Sidesway_Inhibited')
               
                    fprintf(fid,'  load $numNodes  $P  0.0 0.0 0.0 [expr $X*$betaTopY] [expr $X*$betaTopZ] \n');
                    fprintf(fid,'  load 1          0.0 0.0 0.0 0.0 [expr $X*$betaBotY] [expr $X*$betaBotZ] \n');  

                elseif strcmp(obj.frame_typeY,'Sidesway_Uninhibited') &&...
                   strcmp(obj.frame_typeZ,'Sidesway_Inhibited') 
                    fprintf(fid,'  load $numNodes  $P  0.0 $Hz 0.0 0.0 [expr $X*$betaTopZ] \n'); 
                    fprintf(fid,'  load 1          0.0 0.0 0.0 0.0 0.0 [expr $X*$betaBotZ] \n');  
                    if obj.gammaY ~= 0 
                        fprintf(fid,'  load  101 [expr $P*$gammaY] 0.0 0.0 0.0 0.0 0.0 \n');  
                    end 
                    
                elseif strcmp(obj.frame_typeZ,'Sidesway_Uninhibited') &&...
                   strcmp(obj.frame_typeY,'Sidesway_Inhibited') 
                    fprintf(fid,'  load $numNodes  $P  $Hy 0.0 0.0 [expr $X*$betaTopY] 0.0  \n'); 
                    fprintf(fid,'  load 1          0.0 0.0 0.0 0.0 [expr $X*$betaBotY] 0.0  \n');  
                    if obj.gammaZ ~= 0
                        fprintf(fid,'  load  103 [expr $P*$gammaZ] 0.0 0.0 0.0 0.0 0.0 \n');  
                    end
                elseif strcmp(obj.frame_typeY,'Sidesway_Uninhibited') &&...
                   strcmp(obj.frame_typeZ,'Sidesway_Uninhibited') 
                    fprintf(fid,'  load $numNodes  $P  $Hy $Hz 0.0 0.0 0.0  \n'); 
                    if obj.gammaZ ~= 0
                        fprintf(fid,'  load  103 [expr $P*$gammaZ] 0.0 0.0 0.0 0.0 0.0 \n');  
                    end
                    if obj.gammaY ~= 0 
                        fprintf(fid,'  load  101 [expr $P*$gammaY] 0.0 0.0 0.0 0.0 0.0 \n');  
                    end                    
                else
                    error('frameType not recgonized');
                end
                               
                fprintf(fid,'} \n');            
                fprintf(fid,'set ok 0 \n');
                fprintf(fid,'set iStep 0 \n');             
                fprintf(fid,'set loadStep [expr 1.0/$numStepsLateral] \n');             
                fprintf(fid,'integrator LoadControl $loadStep \n');             
                fprintf(fid,'test NormUnbalance $baseForceTolerance 30 $testOutputFlag \n');
                fprintf(fid,'analysis Static \n');
                fprintf(fid,'while {$numStepsLateral > $iStep && $ok == 0} { \n');
                fprintf(fid,'    set ok [analyze 1] \n');
                fprintf(fid,'    if {$ok != 0} { \n');
                fprintf(fid,'       test NormUnbalance $baseForceTolerance 30 $testOutputFlag \n');
                fprintf(fid,'       algorithm NewtonLineSearch .8 \n');
                fprintf(fid,'       set ok [analyze 1] \n');
                fprintf(fid,'       test NormUnbalance $baseForceTolerance 30 $testOutputFlag \n');
                fprintf(fid,'       algorithm Newton \n');
                fprintf(fid,'       } \n');
                fprintf(fid,'   if {$ok != 0} { \n');
                fprintf(fid,'      test NormUnbalance $baseForceTolerance 30 $testOutputFlag \n');
                fprintf(fid,'      algorithm ModifiedNewton -initial \n');
                fprintf(fid,'      set ok [analyze 1] \n');
                fprintf(fid,'      test NormUnbalance $baseForceTolerance 30 $testOutputFlag \n');
                fprintf(fid,'      algorithm Newton \n');
                fprintf(fid,'      } \n');
                fprintf(fid,'    if {$ok == 0} { \n');
                fprintf(fid,'       set lowestEigenValue [eigen -standard -symmBandLapack 1] \n');
                fprintf(fid,'       puts $eigenFileId "[getTime] $lowestEigenValue" \n');
                fprintf(fid,'       incr iStep \n');
                fprintf(fid,'       } \n');
                fprintf(fid,'} \n');
                fprintf(fid,'if {$ok == 0} { \n');
                fprintf(fid,'   # Analysis completed sucessfully \n');
                fprintf(fid,'   exit 5 \n');
                fprintf(fid,'   } else { \n');
                fprintf(fid,'   # Analysis failed in lateral loading \n');
                fprintf(fid,'   exit 6 \n');
                fprintf(fid,'   } \n');
            
            elseif strcmp(analysisType,'TargetForce_NonProportional')
                fprintf(fid,'\n# ############ LOAD CONTROL / NON-PROPORTIONAL LOADING ############ \n');
                % Axial Loading
                if P ~= 0     
                    fprintf(fid,'# Apply gravity loads \n');
                    fprintf(fid,'pattern Plain 1 Linear { \n');
                    fprintf(fid,'  load $numNodes $P 0.0 0.0 0.0 0.0 0.0 \n');

                    if strcmp(obj.frame_typeZ,'Sidesway_Uninhibited')
                        if obj.gammaZ ~= 0 
                            fprintf(fid,'  load  103 [expr $P*$gammaZ] 0.0 0.0 0.0 0.0 0.0 \n');  
                        end
                    end
                    if strcmp(obj.frame_typeY,'Sidesway_Uninhibited')
                        if obj.gammaY ~= 0 
                            fprintf(fid,'  load  101 [expr $P*$gammaY] 0.0 0.0 0.0 0.0 0.0 \n');  
                        end
                    end
                    fprintf(fid,'} \n');               

                    fprintf(fid,'set ok 0 \n');
                    fprintf(fid,'set iStep 0 \n');                     
                    fprintf(fid,'integrator LoadControl [expr 1.0/$numStepsGravity] \n');
                    fprintf(fid,'test NormUnbalance $baseForceTolerance 30 $testOutputFlag \n');                
                    fprintf(fid,'analysis Static \n');
                    fprintf(fid,'while {$numStepsGravity > $iStep && $ok == 0} { \n');
                    fprintf(fid,'      set ok [analyze 1] \n');
                    fprintf(fid,'      if {$ok == 0} {  \n');
                    fprintf(fid,'          %s \n',obj.lowestEigenValue_command);
                    fprintf(fid,'          puts $eigenFileId "[getTime] $lowestEigenValue"  \n');
                    fprintf(fid,'          incr iStep  \n');
                    fprintf(fid,'          } \n');
                    fprintf(fid,'      } \n');

                    fprintf(fid,'if {$ok != 0} {  \n');
                    fprintf(fid,'   # Analysis failed in gravity loading \n');
                    fprintf(fid,'   exit 3 \n');                
                    fprintf(fid,'   } \n');
                    fprintf(fid,'loadConst -time 0.0 \n');
                end   
                % Lateral Loading
                fprintf(fid,'# Apply lateral loads \n');

                fprintf(fid,'pattern Plain 2 Linear { \n');                                             

                if strcmp(obj.frame_typeY,'Sidesway_Inhibited') && ...
                   strcmp(obj.frame_typeZ,'Sidesway_Inhibited')
                                        
                   fprintf(fid,'  load $numNodes  0.0 0.0 0.0 0.0 [expr $X*$betaTopY] [expr $X*$betaTopZ] \n');
                   fprintf(fid,'  load 1          0.0 0.0 0.0 0.0 [expr $X*$betaBotY] [expr $X*$betaBotZ] \n');                        
                  
                elseif strcmp(obj.frame_typeZ,'Sidesway_Uninhibited')&& ...
                       strcmp(obj.frame_typeY,'Sidesway_Inhibited')

                    fprintf(fid,'  load $numNodes  0.0 $Hy 0.0 0.0 [expr $X*$betaTopY] 0.0 \n');     
                    fprintf(fid,'  load 1          0.0 0.0 0.0 0.0 [expr $X*$betaBotY] 0.0 \n');
                    
                elseif strcmp(obj.frame_typeZ,'Sidesway_Inhibited')&& ...
                       strcmp(obj.frame_typeY,'Sidesway_Uninhibited')

                    fprintf(fid,'  load $numNodes  0.0 0.0 $Hz 0.0 0.0 [expr $X*$betaTopZ] \n');     
                    fprintf(fid,'  load 1          0.0 0.0 0.0 0.0 0.0 [expr $X*$betaBotZ] \n');
                elseif strcmp(obj.frame_typeY,'Sidesway_Uninhibited')&& ...
                       strcmp(obj.frame_typeZ,'Sidesway_Uninhibited')

                    fprintf(fid,'  load $numNodes  0.0 $Hy $Hz 0.0 0.0 0.0 \n');     
                else
                    error('frameType not recgonized');
                end
     
                load_step = 1/loadStep.numStepsLateral;
                fprintf(fid,'  set loadStep %g \n',load_step);
                
                fprintf(fid,'} \n');

                fprintf(fid,'set ok 0 \n');
                fprintf(fid,'set iStep 0 \n');            
                fprintf(fid,'integrator LoadControl $loadStep \n');
                fprintf(fid,'test NormUnbalance $baseForceTolerance 30 $testOutputFlag \n');                
                fprintf(fid,'analysis Static \n');
                fprintf(fid,'while {$numStepsLateral > $iStep && $ok == 0} { \n');
                fprintf(fid,'  set ok [analyze 1] \n'); 
                fprintf(fid,'  if {$ok != 0} { \n');
                fprintf(fid,'      test NormUnbalance $baseForceTolerance 30 $testOutputFlag \n');
                fprintf(fid,'      algorithm NewtonLineSearch .8 \n');
                fprintf(fid,'      set ok [analyze 1] \n');
                fprintf(fid,'      test NormUnbalance $baseForceTolerance 30 $testOutputFlag \n');
                fprintf(fid,'      algorithm Newton \n');
                fprintf(fid,'      } \n');
                fprintf(fid,'  if {$ok != 0} { \n');
                fprintf(fid,'      test NormUnbalance $baseForceTolerance 30 $testOutputFlag \n');
                fprintf(fid,'      algorithm ModifiedNewton -initial \n');
                fprintf(fid,'      set ok [analyze 1] \n');
                fprintf(fid,'      test NormUnbalance $baseForceTolerance 30 $testOutputFlag \n');
                fprintf(fid,'      algorithm Newton \n');
                fprintf(fid,'      } \n'); 
                fprintf(fid,'  if {$ok != 0} { \n');
                fprintf(fid,'      test NormUnbalance [expr 100*$baseForceTolerance] 30 $testOutputFlag \n');
                fprintf(fid,'      algorithm ModifiedNewton -initial \n');
                fprintf(fid,'      set ok [analyze 1] \n');
                fprintf(fid,'      test NormUnbalance $baseForceTolerance 30 $testOutputFlag \n');
                fprintf(fid,'      algorithm Newton \n');
                fprintf(fid,'      } \n');                 
                fprintf(fid,'  if {$ok == 0} { \n');
                fprintf(fid,'     %s \n',obj.lowestEigenValue_command);
                fprintf(fid,'     puts $eigenFileId "[getTime] $lowestEigenValue" \n');
                fprintf(fid,'     incr iStep \n');                         
                fprintf(fid,'     } \n');
                fprintf(fid,'  } \n');
                fprintf(fid,'  if {$ok == 0} { \n');
                fprintf(fid,'      # Analysis completed sucessfully \n');
                fprintf(fid,'      exit 5 \n');
                fprintf(fid,'  } else { \n');
                fprintf(fid,'      # Analysis failed in lateral loading \n');
                fprintf(fid,'      exit 6 \n');
                fprintf(fid,'  } \n');
                
            elseif strcmp(analysisType,'TargetDrift_NonProportional')
                fprintf(fid,'\n# ############ LOAD CONTROL / NON-PROPORTIONAL LOADING ############ \n');
                error('CHECK BEFORE USE')
                % Axial Loading
                if P ~= 0     
                    fprintf(fid,'# Apply gravity loads \n');
                    fprintf(fid,'pattern Plain 1 Linear { \n');
                    fprintf(fid,'  load $numNodes $P 0.0 0.0 0.0 0.0 0.0 \n');
                    if strcmp(obj.frame_typeZ,'Sidesway_Uninhibited')
                        if obj.gammaZ ~= 0
                            fprintf(fid,'  load  103 [expr $P*$gammaZ] 0.0 0.0 0.0 0.0 0.0 \n');  
                        end
                    end
                    if strcmp(obj.frame_typeY,'Sidesway_Uninhibited')
                        if obj.gammaY ~= 0 
                            fprintf(fid,'  load  101 [expr $P*$gammaY] 0.0 0.0 0.0 0.0 0.0 \n');  
                        end
                    end                    
                    fprintf(fid,'} \n');               
                    fprintf(fid,'set ok 0 \n');
                    fprintf(fid,'set iStep 0 \n');                     
                    fprintf(fid,'integrator LoadControl [expr 1.0/$numStepsGravity] \n');
                    fprintf(fid,'test NormUnbalance $baseForceTolerance 30 $testOutputFlag \n');                
                    fprintf(fid,'analysis Static \n');
                    fprintf(fid,'while {$numStepsGravity > $iStep && $ok == 0} { \n');
                    fprintf(fid,'      set ok [analyze 1] \n');
                    fprintf(fid,'   if {$ok == 0} {  \n');
                    fprintf(fid,'      %s \n',obj.lowestEigenValue_command);
                    fprintf(fid,'      puts $eigenFileId "[getTime] $lowestEigenValue"  \n');
                    fprintf(fid,'      incr iStep  \n');
                    fprintf(fid,'      } \n');
                    fprintf(fid,'} \n');

                    fprintf(fid,'if {$ok != 0} {  \n');
                    fprintf(fid,'   # Analysis failed in gravity loading \n');
                    fprintf(fid,'   exit 3 \n');                
                    fprintf(fid,'   } \n');
                    fprintf(fid,'loadConst -time 0.0 \n');
                end                        
                % Lateral Loading
                fprintf(fid,'# Apply lateral loads \n');
                fprintf(fid,'pattern Plain 2 Linear { \n');                                             
                if strcmp(obj.frame_typeY,'Sidesway_Inhibited')
                    if obj.betaTopZ ~= 0 
                        fprintf(fid,'  load $numNodes  0.0 0.0 0.0 0.0 0.0 $betaTopZ \n');
                    end
                    fprintf(fid,'  load 1          0.0 0.0 0.0 0.0 0.0 $betaBotZ \n');                                            
                elseif strcmp(obj.frame_typeY,'Sidesway_Uninhibited')

                    fprintf(fid,'  load $numNodes  0.0 $Hy 0.0 0.0 0.0 0.0 \n');                         
                else
                    error('frameType not recgonized');
                end
    
                if strcmp(obj.frame_typeZ,'Sidesway_Inhibited')
                    if obj.betaTopY ~= 0 
                        fprintf(fid,'  load $numNodes  0.0 0.0 0.0 0.0 $betaTopY 0.0 \n');
                    end
                    fprintf(fid,'  load 1          0.0 0.0 0.0 0.0 $betaBotY 0.0 \n');        

                elseif strcmp(obj.frame_typeZ,'Sidesway_Uninhibited')

                    fprintf(fid,'  load $numNodes  0.0 0.0 $Hz 0.0 0.0 0.0 \n');               
                else
                    error('frameType not recgonized');
                end                                 
                fprintf(fid,'} \n');

                fprintf(fid,'set ok 0 \n');
                fprintf(fid,'set iStep 0 \n');            

                if strcmp(obj.frame_type,'Sidesway_Inhibited')

                   if (obj.betaBot > obj.betaTop)
                    fprintf(fid,'set controlledNode [expr $midNode-1] \n');
                   elseif (obj.betaBot < obj.betaTop)
                    fprintf(fid,'set controlledNode [expr $midNode+1] \n'); 
                   else
                    fprintf(fid,'set controlledNode [expr $midNode] \n');                        
                   end

                elseif strcmp(obj.frame_type,'Sidesway_Uninhibited')
                    fprintf(fid,'set controlledNode $numNodes \n');        
                else
                    error('FrameType not recgonized');
                end
                fprintf(fid,'set dispStepSize [expr $dX/double($numStepsLateral)] \n');                                                
                fprintf(fid,'integrator DisplacementControl $controlledNode $controlledDOF $dispStepSize \n');  
                fprintf(fid,'test NormUnbalance $baseForceTolerance 30 $testOutputFlag \n');                
                fprintf(fid,'analysis Static \n');   
                fprintf(fid,'while {$numStepsLateral > $iStep && $ok == 0} { \n');            
                fprintf(fid,'  set ok [analyze 1] \n'); 
                fprintf(fid,'  if {$ok != 0} { \n');
                fprintf(fid,'      test NormUnbalance $baseForceTolerance 30 $testOutputFlag \n');
                fprintf(fid,'      algorithm NewtonLineSearch .8 \n');
                fprintf(fid,'      set ok [analyze 1] \n');
                fprintf(fid,'      test NormUnbalance $baseForceTolerance 30 $testOutputFlag \n');
                fprintf(fid,'      algorithm Newton \n');
                fprintf(fid,'  } \n');
                fprintf(fid,'  if {$ok != 0} { \n');
                fprintf(fid,'      test NormUnbalance $baseForceTolerance 30 $testOutputFlag \n');
                fprintf(fid,'      algorithm ModifiedNewton -initial \n');
                fprintf(fid,'      set ok [analyze 1] \n');
                fprintf(fid,'      test NormUnbalance $baseForceTolerance 30 $testOutputFlag \n');
                fprintf(fid,'      algorithm Newton \n');
                fprintf(fid,'  } \n');                      
                fprintf(fid,'  if {$ok == 0} { \n');
                fprintf(fid,'      %s \n',obj.lowestEigenValue_command);
                fprintf(fid,'      puts $eigenFileId "[getTime] $lowestEigenValue" \n');
                fprintf(fid,'      incr iStep \n');                         
                fprintf(fid,'  } \n');
                fprintf(fid,'} \n');
                fprintf(fid,'if {$ok == 0} { \n');
                fprintf(fid,'    # Analysis completed sucessfully \n');
                fprintf(fid,'    exit 5 \n');
                fprintf(fid,'} else { \n');
                fprintf(fid,'    # Analysis failed in lateral loading \n');
                fprintf(fid,'    exit 6 \n');
                fprintf(fid,'} \n');
            else
                error('analysisType: %s not recgonized',analysisType);
            end
            
            % Close File
            fprintf(fid,'close $eigenFileId \n');
            fclose(fid);
            
            
            %% Run Analysis
            [status, result] = obj.runOpenSees(inputFilename);
            
            
            %% Store Analysis Data
            results = struct;
            results.member_type = 'Beam-Column';
            results.textOutput = result;
            switch status
                case 1
                    results.exitStatus = 'Peak Obtained';
                case 2
                    results.exitStatus = 'Step Limits Reached';
                case 3
                    fprintf('%s\n',result);
                    warning('Analysis failed in gravity loading');
                    if obj.numStepsGravity < 10000
                        warning('Repeating the analysis with larger numStepsGravity');
                        obj.numStepsGravity = obj.numStepsGravity * 2;
                        results = runAnalysis(obj,analysisType,X1,X2,X3,tryNumber); % @todo - update this
                        return
                    else
                        error('Analysis failed in gravity loading');
                    end
                case 4
                    results.exitStatus = 'Analysis Failed In Main Loading';
                case 5
                    results.exitStatus = 'Loads Reached';
                case 6
                    fprintf('%s\n',result);
                    warning('Analysis Failed In Lateral Loading');
                    results.exitStatus = 'Analysis Failed In Lateral Loading';
                    % if obj.numStepsLateral < 10000                    
                    %     warning('Repeating the analysis with larger numStepsLateral');
                    %     obj.numStepsLateral = obj.numStepsLateral * 2;
                    %     results = runAnalysis(obj,analysisType,X1,X2,X3,tryNumber);
                    %     return
                    % end
                case 7
                    results.exitStatus = 'Deformation Limits Reached';
                case 8
                    results.exitStatus = 'Analysis Timed Out';
                otherwise
                    if strcmp(obj.eigenType,'genBandArpack') 
                        warning('Repeating the analysis with symmBandLapack');
                        obj.eigenType = 'symmBandLapack';
                        results = runAnalysis(obj,analysisType,X1,X2,X3,tryNumber);
                        obj.eigenType = 'genBandArpack';
                        return
                    else
                        fprintf('%s\n',result);
                        error('Analysis ended in an unknown manner, exit code = %i',status);
                    end                    

            end


            %% Read Results
            temp = dlmread(nodeDisplacementFilename);
            nodeDisp.x         = temp(:,1:6:end);
            nodeDisp.y         = temp(:,2:6:end);
            nodeDisp.z         = temp(:,3:6:end);
            nodeDisp.lateral   = sqrt(nodeDisp.z.^2 + nodeDisp.y.^2);
            
            nodeDisp.qx        = temp(:,4:6:end);
            nodeDisp.qy        = temp(:,5:6:end);
            nodeDisp.qz        = temp(:,6:6:end);
            nodeDisp.q_lateral = sqrt(nodeDisp.qz.^2 + nodeDisp.qy.^2);

            temp = dlmread(elementForceFilename);
            eleForce.Pi  = temp(:,1:12:end);
            eleForce.Vyi = temp(:,2:12:end);
            eleForce.Vzi = temp(:,3:12:end);            
            eleForce.Ti  = temp(:,4:12:end);
            eleForce.Myi = temp(:,5:12:end);
            eleForce.Mzi = temp(:,6:12:end);
            
            eleForce.Pj  = temp(:,7:12:end);
            eleForce.Vyj = temp(:,8:12:end);
            eleForce.Vzj = temp(:,9:12:end);            
            eleForce.Tj  = temp(:,10:12:end);
            eleForce.Myj = temp(:,11:12:end);
            eleForce.Mzj = temp(:,12:12:end);            
            
            eleForce.P_node  = horzcat(-eleForce.Pi(:,1),eleForce.Pj);
            eleForce.Mz_node = abs(horzcat(eleForce.Mzi(:,1),eleForce.Mzj));
            eleForce.My_node = abs(horzcat(eleForce.Myi(:,1),eleForce.Myj));
            eleForce.M_node  = (eleForce.Mz_node.^2 + eleForce.My_node.^2).^0.5;
            
            temp = dlmread(elementSectionDeformationFilename);
            sectionDeformation.axial      = temp(:,1:3:end);
            sectionDeformation.curvatureY = temp(:,2:3:end);
            sectionDeformation.curvatureZ = temp(:,3:3:end);
            
            
            temp = dlmread(eigenFilename);
            time = temp(:,1);
            lowestEigenvalue = temp(:,2);                


            %% Post-Process Results

            % Applied Loads
            switch analysisType
                case 'LimitPoint_Proportional'
                    ind_grav = [];
                    ind_main = 1:numel(time);
                    
                    P1 = -time;
                    x  = time;
                    X  = X_over_P;
                case 'TargetForce_Proportional'
                    ind_grav = [];
                    ind_main = 1:numel(time);

                    P1 = time * P;
                    x  = time;
                case {'LimitPoint_NonProportional','TargetDrift_NonProportional'}
                    if P == 0
                        ind_grav = [];
                        ind_main = 1:numel(time);

                        P1 = zeros(size(time));
                        x  = time;
                    else
                        ind_grav = 1:(obj.numStepsGravity+1);
                        ind_main = (obj.numStepsGravity+1):numel(time);

                        time_1 = time(1:(obj.numStepsGravity+1));
                        P1_1 = time_1*P;
                        x_1  = zeros(size(P1_1));

                        time_2 = time((obj.numStepsGravity+2):end);
                        P1_2 = ones(size(time_2))*P;
                        x_2  = time_2;

                        P1 = vertcat(P1_1,P1_2);
                        x  = vertcat(x_1,x_2);
                    end
                    X = 1;
                case 'TargetForce_NonProportional'
                    if P == 0
                        ind_grav = [];
                        ind_main = 1:numel(time);

                        P1 = zeros(size(time));
                        x  = time;
                    else
                        ind_grav = 1:(obj.numStepsGravity+1);
                        ind_main = (obj.numStepsGravity+1):numel(time);

                        time_1 = time(1:(obj.numStepsGravity+1));
                        P1_1 = time_1*P;
                        x_1  = zeros(size(P1_1));

                        time_2 = time((obj.numStepsGravity+2):end);
                        P1_2 = ones(size(time_2))*P;
                        x_2  = time_2;

                        P1 = vertcat(P1_1,P1_2);
                        x  = vertcat(x_1,x_2);
                    end
                otherwise
                    error('Unknown analysis type: %s',analysisType)
            end

            % First Order Moment
            data = struct;
            data_fields = fields(obj);
            for i = 1:length(data_fields)
                data.(data_fields{i}) = obj.(data_fields{i}); % @todo - what is this for? is "data" used anywhere?
            end 

            M1y = obj.M1y_over_X()*x * X;
            M1z = obj.M1z_over_X()*x * X;
            M1  = sqrt(M1y.^2+M1z.^2);
            
            
            % Internal Loads
            if obj.reportMaximumLimitpoint
                [~,ind_M] = max(eleForce.M_node,[],2);
                M2z = eleForce.Mz_node(:,ind_M(end));                        
                M2y = eleForce.My_node(:,ind_M(end));            
                P2  = eleForce.P_node(:,ind_M(end));
            else
                M2z = eleForce.Mz_node;                        
                M2y = eleForce.My_node;            
                P2  = eleForce.P_node; 
            end 
            % Deformations            
            if obj.reportMaximumLimitpoint
                [~,ind_lateral] = max(nodeDisp.lateral,[],2);
                def.z = nodeDisp.z(:,ind_lateral(end));                        
                def.y = nodeDisp.y(:,ind_lateral(end));            
                def.x = max(abs(nodeDisp.x),[],2);
            else
                def.z = nodeDisp.z;                        
                def.y = nodeDisp.y;            
                def.x = nodeDisp.x; 
            end
            
            % Rotations            
            if obj.reportMaximumLimitpoint
                [~,ind_qlateral] = max(nodeDisp.q_lateral,[],2);
                def.qz = nodeDisp.qz(:,ind_qlateral(end));                        
                def.qy = nodeDisp.qy(:,ind_qlateral(end));            
                def.qx = max(abs(nodeDisp.qx),[],2);
            else
                def.qz = nodeDisp.qz;                        
                def.qy = nodeDisp.qy;            
                def.qx = nodeDisp.qx; 
            end
            
            % Strains
            maxAbsoluteStrain = ...
                obj.section.longitudinalStrain3d(...
                sectionDeformation.axial,sectionDeformation.curvatureY,...
                sectionDeformation.curvatureZ,'MaxAbsolute');
            maxAbsoluteStrain = max(maxAbsoluteStrain,[],2);

            if obj.section.hasConcrete
                maxConcreteCompressiveStrain = ...
                    obj.section.longitudinalStrain3d(...
                    sectionDeformation.axial,sectionDeformation.curvatureY,...
                    sectionDeformation.curvatureZ,'MaxConcreteCompressive');
                maxConcreteCompressiveStrain = max(maxConcreteCompressiveStrain,[],2);
            end

            % Store Path Results
            results.path.eigen = lowestEigenvalue;
            results.path.P1    = P1;
            results.path.M1y   = M1y;
            results.path.M1z   = M1z;
            results.path.M1    = M1;
            results.path.P2    = P2;
            results.path.M2y   = M2y;
            results.path.M2z   = M2z;
            results.path.def   = def;
            results.path.maxAbsoluteStrain = maxAbsoluteStrain;
            if obj.section.hasConcrete
                results.path.maxConcreteCompressiveStrain = maxConcreteCompressiveStrain;
            else
                results.path.maxConcreteCompressiveStrain = [];
            end
            results.path.topDisplacement = nodeDisp.x(:,end);
            results.path.botRotation.z = -nodeDisp.qz(:,1);
            results.path.botRotation.y = -nodeDisp.qy(:,1);

            if ~isempty(ind_main)
                results.mainPath.eigen  = lowestEigenvalue(ind_main);
                results.mainPath.P1     = P1(ind_main);
                results.mainPath.M1y    = M1y(ind_main);
                results.mainPath.M1z    = M1z(ind_main);
                results.mainPath.M1     = M1(ind_main);                
                results.mainPath.P2     = P2(ind_main,:);
                results.mainPath.M2y    = M2y(ind_main,:);
                results.mainPath.M2z    = M2z(ind_main,:);                
                results.mainPath.def.z  = def.z(ind_main);
                results.mainPath.def.y  = def.y(ind_main);
                results.mainPath.def.x  = def.x(ind_main);
                results.mainPath.def.qz  = def.z(ind_main);
                results.mainPath.def.qy  = def.y(ind_main);
                results.mainPath.def.qx  = def.x(ind_main);
                
                results.mainPath.maxAbsoluteStrain = maxAbsoluteStrain(ind_main);
                if obj.section.hasConcrete
                    results.mainPath.maxConcreteCompressiveStrain = maxConcreteCompressiveStrain(ind_main);
                else
                    results.mainPath.maxConcreteCompressiveStrain = [];
                end
                results.mainPath.topDisplacement = nodeDisp.x(ind_main,end);
                results.mainPath.botRotation.z = -nodeDisp.qz(ind_main,1);
                results.mainPath.botRotation.y = -nodeDisp.qy(ind_main,1);
            end

            if ~isempty(ind_grav)
                results.gravPath.eigen  = lowestEigenvalue(ind_grav);
                results.gravPath.P1     = P1(ind_grav);
                results.gravPath.M1y    = M1y(ind_grav);
                results.gravPath.M1z    = M1z(ind_grav);
                results.gravPath.M1     = M1(ind_grav);                
                results.gravPath.P2     = P2(ind_grav,:);
                results.gravPath.M2y    = M2y(ind_grav,:);
                results.gravPath.M2z    = M2z(ind_grav,:); 
                results.gravPath.def.z  = def.z(ind_grav);
                results.gravPath.def.y  = def.y(ind_grav);
                results.gravPath.def.x  = def.x(ind_grav);
                results.gravPath.def.qz  = def.qz(ind_grav);
                results.gravPath.def.qy  = def.qy(ind_grav);
                results.gravPath.def.qx  = def.qx(ind_grav);
                
                results.gravPath.maxAbsoluteStrain = maxAbsoluteStrain(ind_grav);
                if obj.section.hasConcrete
                    results.gravPath.maxConcreteCompressiveStrain = maxConcreteCompressiveStrain(ind_grav);
                else
                    results.gravPath.maxConcreteCompressiveStrain = [];
                end
                results.gravPath.topDisplacement = nodeDisp.x(ind_grav,end);
                results.gravPath.botRotation.z = -nodeDisp.qz(ind_grav,1);
                results.gravPath.botRotation.y = -nodeDisp.qy(ind_grav,1);
            end


            % Limit Point
            if checkForLimitPoint
                results = obj.findLimitPoint(results);
            end

            % Extra Results
            if obj.extraResultsOutput
                nEL = obj.numElements;
                nIP = obj.numIntegrationPoints;
                Lele = obj.L/nEL;

                % Distribution along length by elements
                xi = nan(1,2*nEL);
                temp = linspace(0,obj.L,nEL+1);
                xi(1:2:2*nEL) = temp(1:(end-1));
                xi(2:2:2*nEL) = temp(2:end);
                
                P = nan(length(time),2*nEL);
                Mz = nan(length(time),2*nEL);
                My = nan(length(time),2*nEL);                
                Vz = nan(length(time),2*nEL);
                Vy = nan(length(time),2*nEL);
                
                P(:,1:2:2*nEL) =  eleForce.Pi;
                P(:,2:2:2*nEL) = -eleForce.Pj;
                Mz(:,1:2:2*nEL) =  eleForce.Mzi;
                Mz(:,2:2:2*nEL) = -eleForce.Mzj;
                
                My(:,1:2:2*nEL) =  eleForce.Myi;
                My(:,2:2:2*nEL) = -eleForce.Myj;
                
                Vz(:,1:2:2*nEL) =  eleForce.Vzi;
                Vz(:,2:2:2*nEL) = -eleForce.Vzj;
                
                Vy(:,1:2:2*nEL) =  eleForce.Vyi;
                Vy(:,2:2:2*nEL) = -eleForce.Vyj;

                % Store results
                results.distE.xi = xi;
                results.distE.P = P;
                results.distE.Mz = Mz;
                results.distE.My = My;
                results.distE.Vz = Vz;
                results.distE.Vy = Vy;

                % Distribution along length by sections
                [ixi,iwt] = lobatto_quadrature(nIP,Lele);
                xi = nan(1,nEL*nIP);
                wt = nan(1,nEL*nIP);
                for i = 1:nEL
                    xi(1+(i-1)*nIP:i*nIP) = ixi + (i-1)*Lele;
                    wt(1+(i-1)*nIP:i*nIP) = iwt;
                end

                % Store Results
                results.distS.xi = xi;
                results.distS.wt = wt;
                results.distS.axialStrain = sectionDeformation.axial;
                results.distS.curvatureZ = sectionDeformation.curvatureZ;
                results.distS.curvatureY = sectionDeformation.curvatureY; 
            end

            % Clean Folder
            if obj.deleteFilesAfterAnalysis
                delete(inputFilename,nodeDisplacementFilename,...
                    elementForceFilename,elementSectionDeformationFilename,...
                    eigenFilename)
            end         
        end        

        function results = runSectionAnalysis(obj,analysisType,X1,X2,tryNumber)
            if nargin < 5
                tryNumber = 1;
            end
            switch analysisType
                case 'LimitPoint_Axial'
                    checkForLimitPoint = true;
                case 'LimitPoint_Proportional'
                    e = X1;
                    checkForLimitPoint = true;
                case 'LimitPoint_Proportional2'
                    e = X1;
                    checkForLimitPoint = true;                    
                case 'LimitPoint_NonProportional'
                    P = X1;
                    checkForLimitPoint = true;
                case 'TargetForce_Proportional'
                    P  = X1;
                    Mz = X2;
                    
                    checkForLimitPoint = false;
                case 'TargetForce_NonProportional'
                    P  = X1;
                    Mz = X2;
                    checkForLimitPoint = false;                  
                % case 'TargetDisplacement_Proportional'
                % case 'TargetDisplacement_NonProportional'
                case 'TargetCurvature_NonProportional'
                    P = X1;
                    kappaZ = X2;
                    checkForLimitPoint = false;
                otherwise
                    error('Unknown analysis type: %s',analysisType)
            end
            
            % Run Analysis
            si = SectionAnalysis(obj.sectionDefinition);
            si.scratchPath              = obj.scratchPath;
            si.deleteFilesAfterAnalysis = obj.deleteFilesAfterAnalysis;
            
            % Step Sizes
            switch tryNumber
                case 1
                    deformationStep_Axial = -0.0005;
                    maxNumSteps_Axial = 1000;
                    deformationStep_Bending = 0.00001/obj.section.depth('strong');
                    maxNumSteps_Bending = round(1/deformationStep_Bending);
                case 2
                    deformationStep_Axial = -0.0005;
                    maxNumSteps_Axial = 1000;
                    deformationStep_Bending = 0.000001/obj.section.depth('strong');
                    maxNumSteps_Bending = round(1/deformationStep_Bending);
                otherwise
                    error('Undefined for tryNumber: %i',tryNumber)
            end
            
            switch analysisType
                case 'LimitPoint_Axial'
                    siResults = si.runSectionToPeak3d(loading_angle(obj,'SectionAnalysis'),'AxialOnly',...
                        deformationStep_Axial,maxNumSteps_Axial);
                case 'LimitPoint_Proportional'
                    siResults = si.runSectionToPeak3d('Proportional',e,...
                        deformationStep_Axial,maxNumSteps_Axial);
                case 'LimitPoint_Proportional2'
                    siResults = si.runSectionToPeak3d('Proportional2',e,...
                        deformationStep_Axial,maxNumSteps_Axial);
                case 'LimitPoint_NonProportional'
                    siResults = si.runSectionToPeak3d(loading_angle(obj,'SectionAnalysis'),'NonProportional',...
                        deformationStep_Bending,maxNumSteps_Bending,P,10);
                case 'TargetForce_Proportional'
                    error('Not implemented');
                    siResults = si.runSectionAnalysis([2 0 2],[0 0; P M],[P/10 M/200]);
                case 'TargetForce_NonProportional'
                    error('Not implemented');
                    siResults = si.runSectionAnalysis([2 0 2],[0 0; P 0; P M],[P/10 M/200]);
                case 'TargetCurvature_NonProportional'
                    error('Not implemented');
                    siResults = si.runSectionAnalysis([2 0 1],[0 0; P 0; P kappa],[P/10 kappa/200]);
                otherwise
                    error('Unknown loadingType');
            end
                    
            
            % Store Analysis Data
            results = struct;
            results.member_type = 'Section';
            results.textOutput  = siResults.textOutput;
            results.exitStatus  = siResults.status;
            
            % Store Path Results
            results.path.P1             = siResults.axialForce;
            results.path.M1z            = siResults.momentz;
            results.path.M1y            = siResults.momenty;
            results.path.M1             = sqrt(results.path.M1z.^2+results.path.M1y.^2);            
            results.path.P2             = results.path.P1;
            results.path.M2z            = results.path.M1z;
            results.path.M2y            = results.path.M1y;            
            results.path.def.z          = siResults.dispz;
            results.path.def.y          = siResults.dispy;
            results.path.axialStrain    = siResults.axialStrain;
            results.path.curvature.z    = siResults.curvaturez;
            results.path.curvature.y    = siResults.curvaturey;
            
            
            results.path.maxAbsoluteStrain = obj.section.longitudinalStrain3d(...
                results.path.axialStrain,results.path.curvature.z,results.path.curvature.y,...
                'MaxAbsolute');
            
            if obj.section.hasConcrete
                results.path.maxConcreteCompressiveStrain = obj.section.longitudinalStrain3d(...
                    results.path.axialStrain,results.path.curvature.z,results.path.curvature.y,...
                    'MaxConcreteCompressive');
            end
            
            
            % Limit Point
            if checkForLimitPoint
                results.path.eigen = siResults.eigen;
                results = obj.findLimitPoint(results);
            end 
            
        end
        
        function results = findLimitPoint(obj,results)
            
            % Find Stability Limit Point
            [ind,x] = find_limit_point_in_vector(results.path.eigen,0.0);
            
            if ~isempty(ind)
                stabilityLimitPoint.limit_type = 'Stability Limit (Reached)';
                stabilityLimitPoint.good = true;
                stabilityLimitPoint.time = ind+x;
            else
                [ind,x] = find_limit_point_in_vector(...
                    results.path.eigen/results.path.eigen(1),...
                    obj.eigenvalueLimitPointTolerance);
                
                if ~isempty(ind)
                    stabilityLimitPoint.limit_type = 'Stability Limit (Tolerance Reached)';
                    stabilityLimitPoint.good = true;
                    stabilityLimitPoint.time = ind+x;
                else
                    ind = results.path.M1(1:end)-vertcat(0.0,results.path.M1(1:end-1)) < 0;
                    
                    if ~isempty(find(ind,1))
                        stabilityLimitPoint.limit_type = 'Stability Limit (Max Time Reached-M)';
                        stabilityLimitPoint.good = true;
                        stabilityLimitPoint.time = find(ind,1)-1;
                    else
                        ind = abs(results.path.P1(1:end))-abs(vertcat(0.0,results.path.P1(1:end-1))) < 0;
                        if ~isempty(find(ind,1))
                            stabilityLimitPoint.limit_type = 'Stability Limit (Max Time Reached-P)';
                            stabilityLimitPoint.good = true;
                            stabilityLimitPoint.time = find(ind,1)-1;
                        else                        
                            stabilityLimitPoint.limit_type = 'Stability Limit (Not Reached)';
                            stabilityLimitPoint.good = false;
                            stabilityLimitPoint.time = Inf;
                        end
                    end                    
                end
            end
            
            
            % Find Strain Limit
            [ind,x] = find_limit_point_in_vector(...
                results.path.maxAbsoluteStrain,obj.absoluteStrainLimit);
            
            if isempty(ind)
                strainLimitPoint.limit_type = 'Strain Limit (Not Reached)';
                strainLimitPoint.good = false;
                strainLimitPoint.time = Inf;
            else
                strainLimitPoint.limit_type = 'Strain Limit (Reached)';
                strainLimitPoint.good = true;
                strainLimitPoint.time = ind+x;
            end            
 
            % Find Concrete Compressive Strain Limit Point
            if ~isempty(obj.concreteCompressiveStrainLimit) && obj.section.hasConcrete
                [ind,x] = find_limit_point_in_vector(...
                    results.path.maxConcreteCompressiveStrain,...
                    obj.concreteCompressiveStrainLimit);
                
                if isempty(ind)
                    concreteCompressiveStrainLimitPoint.limit_type = 'Concrete Compressive Strain Limit (Not Reached)';
                    concreteCompressiveStrainLimitPoint.good = false;
                    concreteCompressiveStrainLimitPoint.time = Inf;
                else
                    concreteCompressiveStrainLimitPoint.limit_type = 'Concrete Compressive Strain Limit (Reached)';
                    concreteCompressiveStrainLimitPoint.good = true;
                    concreteCompressiveStrainLimitPoint.time = ind+x;
                end
            else
                concreteCompressiveStrainLimitPoint.limit_type = 'Concrete Compressive Strain Limit (Not Checked)';
                concreteCompressiveStrainLimitPoint.good = false;
                concreteCompressiveStrainLimitPoint.time = Inf;
            end
            
            % Find Compressive Force Limit Point
            if ~isempty(obj.compressiveForceLimitType) && ~strcmpi(obj.compressiveForceLimitType,'none')
                
                switch lower(obj.compressiveForceLimitType)
                    case 'pnco'
                        compressiveForceLimit = obj.section.Pnco;
                    otherwise
                        error('Unknown compressiveForceLimitType')
                end
                
                [ind,x] = find_limit_point_in_vector(results.path.P2,compressiveForceLimit);
                
                if isempty(ind)
                    compressiveForceLimitPoint.limit_type = 'Compressive Force Limit (Not Reached)';
                    compressiveForceLimitPoint.good = false;
                    compressiveForceLimitPoint.time = Inf;
                else
                    compressiveForceLimitPoint.limit_type = 'Compressive Force Limit (Reached)';
                    compressiveForceLimitPoint.good = true;
                    compressiveForceLimitPoint.time = ind+x;
                end
            else
                compressiveForceLimitPoint.limit_type = 'Compressive Force Limit (Not Checked)';
                compressiveForceLimitPoint.good = false;
                compressiveForceLimitPoint.time = Inf;
            end
            
            % Determine Controling Limit
            limitPointTimes = [
                stabilityLimitPoint.time
                strainLimitPoint.time
                concreteCompressiveStrainLimitPoint.time
                compressiveForceLimitPoint.time];
            
            controllingLimitPointTime = min(limitPointTimes);
            
            if stabilityLimitPoint.time == controllingLimitPointTime
                results.limitPoint = stabilityLimitPoint;                
            elseif strainLimitPoint.time == controllingLimitPointTime
                results.limitPoint = strainLimitPoint;
            elseif concreteCompressiveStrainLimitPoint.time == controllingLimitPointTime
                results.limitPoint = concreteCompressiveStrainLimitPoint;
            elseif compressiveForceLimitPoint.time == controllingLimitPointTime
                results.limitPoint = compressiveForceLimitPoint;
            else
                error('error')
            end
            
            % Warn if limit point could not be determined
            if obj.limitPointWarning && ~results.limitPoint.good
                warning('Limit point could not be determined, \n   Exit Status: %s \n   Eigenvalue Ratio: %g (limit: %g)\n   Max Strain: %g (limit: %g)',...
                    results.exitStatus,...
                    min(results.path.eigen/results.path.eigen(1)),obj.eigenvalueLimitPointTolerance,...
                    max(results.path.maxAbsoluteStrain),obj.absoluteStrainLimit);
            end
            
            % Compute and Store Limit Point Data
            if results.limitPoint.time ~= Inf
                ind = floor(results.limitPoint.time);
                x = results.limitPoint.time - ind;
                results.limitPoint.P1  = interpolate_vector(results.path.P1,ind,x);
                results.limitPoint.M1y = interpolate_vector(results.path.M1y,ind,x);
                results.limitPoint.M1z = interpolate_vector(results.path.M1z,ind,x); 

                if obj.reportMaximumLimitpoint
                    results.limitPoint.P2  = interpolate_vector(results.path.P2,ind,x);
                    results.limitPoint.M2y = interpolate_vector(results.path.M2y,ind,x);
                    results.limitPoint.M2z = interpolate_vector(results.path.M2z,ind,x);
                    results.limitPoint.def.y = interpolate_vector(results.path.def.y,ind,x);
                    results.limitPoint.def.z = interpolate_vector(results.path.def.z,ind,x);
                    
                else
                    for i=1:obj.numElements+1
                        results.limitPoint.P2(i)  = interpolate_vector(results.path.P2(:,i),ind,x);
                        results.limitPoint.M2y(i) = interpolate_vector(results.path.M2y(:,i),ind,x);
                        results.limitPoint.M2z(i) = interpolate_vector(results.path.M2z(:,i),ind,x);
                        results.limitPoint.def.y(i) = interpolate_vector(results.path.def.y(:,i),ind,x);
                        results.limitPoint.def.z(i) = interpolate_vector(results.path.def.z(:,i),ind,x);
                    end
                end
            else
                results.limitPoint.P1  = 0.0;
                results.limitPoint.M1y = 0.0;
                results.limitPoint.M1z = 0.0; 
                if obj.reportMaximumLimitpoint
                    results.limitPoint.P2  = 0.0;
                    results.limitPoint.M2y = 0.0;
                    results.limitPoint.M2z = 0.0;
                    results.limitPoint.def = 0.0;
                    
                else
                    results.limitPoint.P2  = zeros(1,obj.numElements+1);
                    results.limitPoint.M2y = zeros(1,obj.numElements+1);
                    results.limitPoint.M2z = zeros(1,obj.numElements+1);
                    results.limitPoint.def.y = zeros(1,obj.numElements+1);
                    results.limitPoint.def.z = zeros(1,obj.numElements+1);
                    
                end
            end
            
        end
        
        % Command for 
        function str = lowestEigenValue_command(obj)
            switch obj.eigenType
                case 'genBandArpack'
                    str = 'set lowestEigenValue [eigen 1]';
                case 'symmBandLapack'
                    str = 'set lowestEigenValue [eigen  -standard -symmBandLapack 1]';
                otherwise
                    error('Unknown eigenType: %s',obj.eigenType);
            end
        end
            
        % Loading Steps based on frame type and analysis type
        function loadingSteps = loading_steps(obj,frame_type,analysisType,tryNumber)
            
            if obj.L < 100
                e = 0.8;
            elseif obj.L > 100 && obj.L < 300
                e = 1;
            elseif obj.L > 300 && obj.L < 700
                e = 1.2;
            elseif obj.L > 700 
                e = 1.4;
            end

            if strcmp(frame_type,'Sidesway_Uninhibited')
                % Loading Steps
                switch analysisType
                    case 'LimitPoint_Proportional' %'Sidesway_Uninhibited'                      
                        if isempty(obj.controlled_dof_override)
                            direction = 1;
                        else
                            direction = -1;
                        end

                        switch tryNumber
                            case 1
                                loadingSteps.coarseDispStepSize        = obj.L/5000*direction;
                                loadingSteps.fineDispStepSize          = obj.L/50000*direction;
                                loadingSteps.maxDisp                   = 0.10*obj.L;
                                loadingSteps.maxAnalysisTime           = 240;
                                loadingSteps.baseForceTolerance        = obj.base_tolerance_force;
                                loadingSteps.baseDisplacementTolerance = 1e-7*obj.L;
                            case 2
                                loadingSteps.coarseDispStepSize        = obj.L/100000*direction;
                                loadingSteps.fineDispStepSize          = obj.L/1000000*direction;
                                loadingSteps.maxDisp                   = 0.10*obj.L;
                                loadingSteps.maxAnalysisTime           = 240;
                                loadingSteps.baseForceTolerance        = obj.base_tolerance_force;
                                loadingSteps.baseDisplacementTolerance = 1e-7*obj.L;

                            otherwise
                                error('Undefined for tryNumber: %i',tryNumber)
                        end

                    case 'LimitPoint_NonProportional'%'Sidesway_Uninhibited'
                        switch tryNumber
                            case 1
                                loadingSteps.numStepsGravity           = obj.numStepsGravity;
                                loadingSteps.coarseDispStepSize        = obj.L^e/2000;
                                loadingSteps.fineDispStepSize          = obj.L^e/20000;
                                loadingSteps.maxDisp                   = 0.15*obj.L;
                                loadingSteps.maxAnalysisTime           = 240;
                                loadingSteps.baseForceTolerance        = obj.base_tolerance_force;
                                loadingSteps.baseDisplacementTolerance = 1e-7*obj.L;                                
                            case 2
                                loadingSteps.numStepsGravity           = obj.numStepsGravity;
                                loadingSteps.coarseDispStepSize        = obj.L^e/100000;
                                loadingSteps.fineDispStepSize          = obj.L^e/1000000;
                                loadingSteps.maxDisp                   = 0.15*obj.L;
                                loadingSteps.maxAnalysisTime           = 600;
                                loadingSteps.baseForceTolerance        = obj.base_tolerance_force;
                                loadingSteps.baseDisplacementTolerance = 1e-7*obj.L;                                 
                            case 3
                                loadingSteps.numStepsGravity           = obj.numStepsGravity;
                                loadingSteps.coarseDispStepSize        = obj.L/1000000;
                                loadingSteps.fineDispStepSize          = obj.L/10000000;
                                loadingSteps.maxDisp                   = 0.15*obj.L;
                                loadingSteps.maxAnalysisTime           = 600;
                                loadingSteps.baseForceTolerance        = obj.base_tolerance_force;
                                loadingSteps.baseDisplacementTolerance = 1e-7*obj.L;   
                            otherwise
                                error('Undefined for tryNumber: %i',tryNumber)
                        end       

                    case 'TargetForce_Proportional'
                        
                        loadingSteps.numStepsGravity           = obj.numStepsGravity;
                        loadingSteps.numStepsLateral           = obj.numStepsLateral;
                        loadingSteps.baseForceTolerance        = obj.base_tolerance_force;
                        loadingSteps.baseDisplacementTolerance = 1e-8*obj.L;                                              

                    case  'TargetForce_NonProportional'
                        
                        loadingSteps.numStepsGravity           = obj.numStepsGravity;
                        loadingSteps.numStepsLateral           = obj.numStepsLateral;
                        loadingSteps.baseForceTolerance        = obj.base_tolerance_force;
                        loadingSteps.baseDisplacementTolerance = 1e-8*obj.L;
                        
                    case 'TargetDrift_NonProportional'
                        loadingSteps.numStepsGravity           = obj.numStepsGravity;
                        loadingSteps.numStepsLateral           = obj.numStepsLateral;
                        loadingSteps.baseForceTolerance        = obj.base_tolerance_force;
                        loadingSteps.baseDisplacementTolerance = 1e-8*obj.L;

                    otherwise
                        error('Unknown analysis type: %s',analysisType)
                end                

            elseif strcmp(frame_type,'Sidesway_Inhibited')

                % Loading Steps
                switch analysisType
                    case 'LimitPoint_Proportional'

                        if isempty(obj.controlled_dof_override)
                            direction = 1;
                        else
                            direction = -1;
                        end

                        switch tryNumber
                            case 1
                                loadingSteps.coarseDispStepSize        = obj.L/100000*direction;
                                loadingSteps.fineDispStepSize          = obj.L/1000000*direction;
                                loadingSteps.maxDisp                   = 0.05*obj.L;
                                loadingSteps.maxAnalysisTime           = 240;
                                loadingSteps.baseForceTolerance        = obj.base_tolerance_force;
                                loadingSteps.baseDisplacementTolerance = 1e-7*obj.L;                                

                            case 2
                                loadingSteps.coarseDispStepSize        = obj.L/600000*direction;
                                loadingSteps.fineDispStepSize          = obj.L/6000000*direction;
                                loadingSteps.maxDisp                   = 0.05*obj.L;
                                loadingSteps.maxAnalysisTime           = 600;
                                loadingSteps.baseForceTolerance        = obj.base_tolerance_force;
                                loadingSteps.baseDisplacementTolerance = 1e-7*obj.L;                                 

                            case 3
                                loadingSteps.coarseDispStepSize        = obj.L/200000*direction;
                                loadingSteps.fineDispStepSize          = obj.L/2000000*direction;
                                loadingSteps.maxDisp                   = 0.05*obj.L;
                                loadingSteps.maxAnalysisTime           = 600;
                                loadingSteps.baseForceTolerance        = obj.base_tolerance_force;
                                loadingSteps.baseDisplacementTolerance = 1e-7*obj.L; 
                                
                            otherwise
                                error('Undefined for tryNumber: %i',tryNumber)
                        end

                    case 'LimitPoint_NonProportional'
                        switch tryNumber
                            case 1
                                
                                loadingSteps.numStepsGravity           = obj.numStepsGravity;
                                loadingSteps.coarseDispStepSize        = obj.L/2000;
                                loadingSteps.fineDispStepSize          = obj.L/20000;
                                loadingSteps.maxDisp                   = 0.10*obj.L;
                                loadingSteps.maxAnalysisTime           = 300;
                                loadingSteps.baseForceTolerance        = obj.base_tolerance_force;
                                loadingSteps.baseDisplacementTolerance = 1e-7*obj.L;                                                                 

                            case 2
                                
                                loadingSteps.numStepsGravity           = obj.numStepsGravity;
                                loadingSteps.coarseDispStepSize        = obj.L/100000;
                                loadingSteps.fineDispStepSize          = obj.L/1000000;
                                loadingSteps.maxDisp                   = 0.10*obj.L;
                                loadingSteps.maxAnalysisTime           = 600;
                                loadingSteps.baseForceTolerance        = obj.base_tolerance_force;
                                loadingSteps.baseDisplacementTolerance = 1e-7*obj.L;                                     
                            case 3
                                
                                loadingSteps.numStepsGravity           = obj.numStepsGravity;
                                loadingSteps.coarseDispStepSize        = obj.L/1000000;
                                loadingSteps.fineDispStepSize          = obj.L/10000000;
                                loadingSteps.maxDisp                   = 0.10*obj.L;
                                loadingSteps.maxAnalysisTime           = 600;
                                loadingSteps.baseForceTolerance        = obj.base_tolerance_force;
                                loadingSteps.baseDisplacementTolerance = 1e-7*obj.L;  
                            otherwise
                                error('Undefined for tryNumber: %i',tryNumber)
                        end
                        
                    case 'TargetForce_Proportional'
                        
                        loadingSteps.numStepsGravity           = obj.numStepsGravity;
                        loadingSteps.numStepsLateral           = obj.numStepsLateral;
                        loadingSteps.baseForceTolerance        = obj.base_tolerance_force;
                        loadingSteps.baseDisplacementTolerance = 1e-8*obj.L;                                              

                    case  'TargetForce_NonProportional'
                        
                        loadingSteps.numStepsGravity           = obj.numStepsGravity;
                        loadingSteps.numStepsLateral           = obj.numStepsLateral;
                        loadingSteps.baseForceTolerance        = obj.base_tolerance_force;
                        loadingSteps.baseDisplacementTolerance = 1e-8*obj.L;                  

                    case 'TargetDrift_NonProportional'
                        loadingSteps.numStepsGravity           = obj.numStepsGravity;
                        loadingSteps.numStepsLateral           = obj.numStepsLateral;
                        loadingSteps.baseForceTolerance        = obj.base_tolerance_force;
                        loadingSteps.baseDisplacementTolerance = 1e-8*obj.L;                        
                        
                    otherwise
                        error('Unknown analysis type: %s',analysisType)
                end                     

            else
                error('Unknown frame type: %s',obj.frame_type);
            end
        end        
 
    end
    
end

