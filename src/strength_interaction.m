classdef strength_interaction < BenchmarkAnalysis3d_OpenSees

% This object gets the frame and section data of a steel-concrete
% composite beam-column and outputs surface 1 through 5 for 
% maximum or each element internal forces along the length of a beam-column
%{
surface 1 : Maximum Applied Loads Determined by Second-Order Inelastic
Analysis
surface 2 : Maximum Internal forces Determined by Second-Order Inelastic
Analysis
surface 3 : Proposed/psd_ACB Nominal Beam-column Strength 
surface 4 : Maximum Applied Loads Permitted by the Proposed/psd_ACB
Methodology
surface 5 : Maximum Internal Forces Determined by Second-Order Elastic
Analy
%}

properties
    % Cross Section and frame Properties
    frame_data
    % Analysis Options
    axis                                    = 'all'; % 'all' or 'weak' or 'strong' or 'biaxial'
    loading_angle                           = [];   % Just for 'biaxial' axis
    check_for_errors_in_surface1            = false;
    extra_fine_DispStepSize                 = false;
    % Elastic Analysis
    applyNotionalLoad                       = [];
    % Interaction Options        
    numPoints = 9;
    numAngles = 5;
end

    methods
        %% Constructor
        function obj = strength_interaction(data)
            obj.frame_data = data;         
        end
        function alpha = alpha(obj)
            switch obj.frame_data.section_type
                case 'RCFT'
                    alpha = 1.5;
                case 'SRC'
                    alpha = 1.25;
            end        
        end

        function val = num_angles(obj)
            if strcmp(obj.axis,'strong') || strcmp(obj.axis,'weak') || strcmp(obj.axis,'biaxial')
                val = 1;
            elseif strcmp(obj.axis,'all')
                val = obj.numAngles;
            else
                error('undefined axis');
            end
        end
        function angle = angle(obj,i)
            if strcmp(obj.axis,'strong')
                angle = 0;
            elseif strcmp(obj.axis,'weak') 
                angle = pi/2;
            elseif strcmp(obj.axis,'all')
                angles = linspace(0,pi/2,num_angles(obj));
                angle = angles(i);
            elseif strcmp(obj.axis,'biaxial')
                angle = obj.loading_angle;        
            else
                error('undefined axis');
            end        
        end

        function  val = psd_ACB(obj) 
            val.PNO = abs(obj.frame_data.section.Pnco);
            val.MNOz=abs(obj.frame_data.section.Mno('strong'));
            val.MNOy=abs(obj.frame_data.section.Mno('weak'));

            [Pa_weak,~] = obj.frame_data.section.pointA('weak');
            [Pc_weak,val.Mcy] = obj.frame_data.section.pointC('weak');

            [Pa_strong,~] = obj.frame_data.section.pointA('strong');
            [Pc_strong,val.Mcz] = obj.frame_data.section.pointC('strong');
            if strcmp(obj.axis,'strong')
                xi=obj.frame_data.section.stabilityReduction('strong',val.PNO);
                val.xPa = abs(xi*Pa_strong);
                val.xPc = abs(xi*Pc_strong);  
            elseif strcmp(obj.axis,'weak')
                xi=obj.frame_data.section.stabilityReduction('weak',val.PNO);
                val.xPa = abs(xi*Pa_weak);
                val.xPc = abs(xi*Pc_weak);
            elseif strcmp(obj.axis,'biaxial') || strcmp(obj.axis,'all')
                xi=obj.frame_data.section.stabilityReduction('min',val.PNO);
                val.xPa = abs(min(xi*Pa_strong,xi*Pa_weak));
                val.xPc = abs(min(xi*Pc_strong,xi*Pc_weak));
            end
        end

        function B2 = driftRatio(obj,P)
            stiffness_reduction = 0.8;
            tau = 0.8;
            
            EIeffY = stiffness_reduction * tau * obj.frame_data.section.EIeff('weak');
            EIeffZ = stiffness_reduction * tau * obj.frame_data.section.EIeff('strong');
            kqtopY = stiffness_reduction *  obj.frame_data.kqtopY;
            kqbotY = stiffness_reduction *  obj.frame_data.kqbotY;
            kqtopZ = stiffness_reduction *  obj.frame_data.kqtopZ;
            kqbotZ = stiffness_reduction *  obj.frame_data.kqbotZ;          

            if strcmp(obj.frame_data.frame_typeY,'Sidesway_Uninhibited')

                elastic_weak_object = BenchmarkAnalysis2d_Elastic_Sidesway_Uninhibited...
                    (EIeffY,obj.frame_data.L,kqtopY,kqbotY,...
                    obj.frame_data.gammaY,0,0);
                elastic_weak_object.includeInitialGeometricImperfections = false;
                B2_weak = elastic_weak_object.driftRatio(P);
            elseif strcmp(obj.frame_data.frame_typeY,'Sidesway_Inhibited')
                B2_weak = 1;
            end
            if strcmp(obj.frame_data.frame_typeZ,'Sidesway_Uninhibited')

                elastic_strong_object = BenchmarkAnalysis2d_Elastic_Sidesway_Uninhibited...
                    (EIeffZ,obj.frame_data.L,kqtopZ,kqbotZ,...
                    obj.frame_data.gammaZ,0,0);
                elastic_strong_object.includeInitialGeometricImperfections = false;
                B2_strong = elastic_strong_object.driftRatio(P);

            elseif strcmp(obj.frame_data.frame_typeZ,'Sidesway_Inhibited')
                B2_strong = 1;
            end
            B2 = max(B2_weak,B2_strong);    
        end

        function [total_load,notional_load] = HwithNotionalLoad(obj,axis,H,P,B2,angle)
            if strcmp(axis,'z')
                frame_type = obj.frame_data.frame_typeZ;
            elseif strcmp(axis,'y')
                frame_type = obj.frame_data.frame_typeY;
            end         
            if strcmp(frame_type, 'Sidesway_Uninhibited') && obj.applyNotionalLoad
                if strcmp(axis,'z')
                    gamma = obj.frame_data.gammaZ;
                    dir = cos(angle);
                elseif strcmp(axis,'y')
                    gamma = obj.frame_data.gammaY;
                    dir = sin(angle);
                end            
                if B2 > 1.7
                    notional_load = 0.002 * dir * (1+gamma) * abs(P);
                    total_load   =  H + notional_load;
                else

                    if H == 0 && P ~=0
                        nl = 0.002 * dir * (1+gamma) * abs(P);
                        total_load   =  nl ;
                        notional_load = nl;
                    else
                        total_load = H;
                        notional_load = 0;
                    end
                end
            elseif strcmp(frame_type, 'Sidesway_Uninhibited') && ~obj.applyNotionalLoad 
                notional_load = 0;
                total_load   =  H;
            elseif strcmp(frame_type, 'Sidesway_Inhibited')
                notional_load = 0;
                total_load   =  H;
            else
                error('Not defined type')
            end
        end

        function [initGeomImperfAng,num] = initGeomImperfAng(obj,angle,i)

            if strcmp(obj.frame_data.frame_typeY,'Sidesway_Uninhibited') ...
                && strcmp(obj.frame_data.frame_typeZ,'Sidesway_Uninhibited')

                ind_imp = obj.frame_data.kqbotY ~= 0 | obj.frame_data.kqbotZ == 0;

            elseif strcmp(obj.frame_data.frame_typeY,'Sidesway_Inhibited') ...
                && strcmp(obj.frame_data.frame_typeZ,'Sidesway_Inhibited') 

                ind_imp = true;
            elseif strcmp(obj.frame_data.frame_typeY,'Sidesway_Inhibited') ...
                && strcmp(obj.frame_data.frame_typeZ,'Sidesway_Uninhibited') 

                ind_imp = true;            
            elseif strcmp(obj.frame_data.frame_typeY,'Sidesway_Uninhibited') ...
                && strcmp(obj.frame_data.frame_typeZ,'Sidesway_Inhibited') 

                ind_imp = false;
            else 
                error('Not implemented yet');
            end
            if angle == pi/2 || ind_imp
                initGeomImperfAngs = unique([0,angle]);
            else
                initGeomImperfAngs = unique([0,angle,pi/2]);
            end
            num = length(initGeomImperfAngs);
            if  exist('i','var')
                initGeomImperfAng = initGeomImperfAngs(i);
            else
                initGeomImperfAng = initGeomImperfAngs;   
            end
        end

        function [notionalLoadAng,num] = notionalLoadAng(obj,angle,i)

            if strcmp(obj.frame_data.frame_typeY,'Sidesway_Inhibited') ...
                && strcmp(obj.frame_data.frame_typeZ,'Sidesway_Inhibited') 

                notionalLoadAng = nan;
                num = 1;
                return
            elseif strcmp(obj.frame_data.frame_typeY,'Sidesway_Uninhibited') ...
                && strcmp(obj.frame_data.frame_typeZ,'Sidesway_Inhibited') 

                notionalLoadAng = pi/2;
                num = 1;
                return
            elseif strcmp(obj.frame_data.frame_typeZ,'Sidesway_Uninhibited') ...
                && strcmp(obj.frame_data.frame_typeY,'Sidesway_Inhibited') 

                notionalLoadAng = 0;
                num = 1;
                return

            elseif strcmp(obj.frame_data.frame_typeY,'Sidesway_Uninhibited') ...
                && strcmp(obj.frame_data.frame_typeZ,'Sidesway_Uninhibited')

                ind_nl = obj.frame_data.kqbotY ~= 0 | obj.frame_data.kqbotZ == 0;
            else 
                error('Not implemented yet');
            end

            if angle == pi/2 || ind_nl
                notionalLoadAngs = unique([0,angle]);
            else
                notionalLoadAngs = unique([0,angle,pi/2]);
            end
            num = length(notionalLoadAngs);
            if  exist('i','var')
                notionalLoadAng = notionalLoadAngs(i);
            else
                notionalLoadAng = notionalLoadAngs;   
            end
        end        
        %% surface 1 and 2
        function [surface1,surface2] = surface_1_and_2(obj)

            analysisOptions = struct;
            analysisOptions.controlled_dof_override = [];
            analysisOptions.includeInitialGeometricImperfections = true;
            analysisOptions.reportMaximumLimitpoint = true;    
            if  strcmp(obj.frame_data.section_type,'SRC')
              analysisOptions.eigenType = 'symmBandLapack';
            end

            for i=1:num_angles(obj)
                [~,M2y,M2z,P0(i,1)] = second_order_analysis(obj,'LimitPoint_Proportional'...
                    ,angle(obj,i),analysisOptions);    
                if obj.reportMaximumLimitpoint
                    M0y(i,1) = M2y;
                    M0z(i,1) = M2z;                                
                else
                    M0y{i,1} = M2y(1,:);
                    M0z{i,1} = M2z(1,:);                   
                end
            end
            if obj.reportMaximumLimitpoint
                P1 = nan(obj.numPoints,num_angles(obj));
                M1y = nan(obj.numPoints,num_angles(obj));
                M1z = nan(obj.numPoints,num_angles(obj));

                P2 = nan(obj.numPoints,num_angles(obj));
                M2y = nan(obj.numPoints,num_angles(obj));
                M2z = nan(obj.numPoints,num_angles(obj));
                status = cell(obj.numPoints,num_angles(obj));                    
            else
                P2 = cell(obj.numPoints,num_angles(obj));
                M2y = cell(obj.numPoints,num_angles(obj));
                M2z = cell(obj.numPoints,num_angles(obj));              
                status = cell(obj.numPoints,num_angles(obj));                    
            end            
            for i = 1:num_angles(obj)           
                    [P2(:,i),M2y(:,i),M2z(:,i),P1(:,i),M1y(:,i),M1z(:,i),...
                        status(:,i)] = second_order_analysis(obj,'LimitPoint_NonProportional',...
                        angle(obj,i),analysisOptions,P0(i,:),M0y(i,:),M0z(i,:));  
            end       
            M1 = (M1z.^2+M1y.^2).^(0.5);
            if obj.check_for_errors_in_surface1
                for i = 2:obj.numPoints-1
                    for j = 1:num_angles(obj)
                        ratio = M1(i,j)/mean([M1(i-1,j),M1(i+1,j)]);
                        if ratio < 0.90
                            [P2(i,j),M2y(i,j),M2z(i,j),P1(i,j),M1y(i,j),M1z(i,j),...
                                ~] = second_order_analysis(obj,...
                                'LimitPoint_NonProportional_2',...
                                angle(obj,j),analysisOptions,P1(i,j));                    
                        end
                    end
                end
                for i = 1:obj.numPoints
                    for j = 2:num_angles(obj)-1
                        ratio = M1(i,j)/mean([M1(i,j-1),M1(i,j+1)]);
                        if ratio < 0.85
                            [P2(i,j),M2y(i,j),M2z(i,j),P1(i,j),M1y(i,j),M1z(i,j),...
                                ~] = second_order_analysis(obj,...
                                'LimitPoint_NonProportional_2',...
                                angle(obj,j),analysisOptions,P1(i,j));                    
                        end
                    end
                end
            end
            if obj.reportMaximumLimitpoint
                surface2 = struct;
                surface2.P = P2;
                surface2.My = M2y;
                surface2.Mz = M2z;         
            else
                surface2 = struct;
                for i =1:obj.numPoints
                    for j=1:num_angles(obj)
                        for k=1:obj.numElements+1
                            surface2(k).Mz(i,j)=M2z{i,j}(1,k);
                            surface2(k).My(i,j)=M2y{i,j}(1,k);
                            surface2(k).P(i,j)=P2{i,j}(1,k);                      
                        end
                    end
                end
            end         

            surface1 = struct;
            surface1.P = P1;
            surface1.My = M1y;
            surface1.Mz = M1z;
            surface1.status = status;
        end

        %% surface 3
        function [surface3] = surface_3(obj)
            P3  = nan(obj.numPoints,num_angles(obj));
            M3y = nan(obj.numPoints,num_angles(obj));
            M3z = nan(obj.numPoints,num_angles(obj));
            M3  = nan(obj.numPoints,num_angles(obj));  
            Ps1  = linspace(obj.psd_ACB.xPa,0,obj.numPoints-1);
            ind = Ps1 < obj.psd_ACB.xPc;
            Ps = horzcat(Ps1(1:find(ind,1)-1),obj.psd_ACB.xPc,Ps1(find(ind,1):end));
            for ii = 1:num_angles(obj)
                for i = 1:length(Ps)
                    if Ps(i) >= obj.psd_ACB.xPc
                        M3(i) = (1-(Ps(i)-obj.psd_ACB.xPc)/(obj.psd_ACB.xPa-obj.psd_ACB.xPc))/...
                            ((sin(angle(obj,ii))/obj.psd_ACB.Mcy)^alpha(obj) +...
                            (cos(angle(obj,ii))/obj.psd_ACB.Mcz)^alpha(obj))^(1/alpha(obj));
                        M3z(i,ii) = M3(i)*cos(angle(obj,ii));
                        M3y(i,ii) = M3(i)*sin(angle(obj,ii));
                        P3(i,ii) = -Ps(i);
                    else               
                        M3(i) = 1/((sin(angle(obj,ii))/obj.psd_ACB.Mcy)^alpha(obj)...
                            + (cos(angle(obj,ii))/obj.psd_ACB.Mcz)^alpha(obj))^(1/alpha(obj));
                        M3z(i,ii) = M3(i)*cos(angle(obj,ii));
                        M3y(i,ii) = M3(i)*sin(angle(obj,ii));
                        P3(i,ii) = -Ps(i); 
                    end  
                end  
            end
            surface3 = struct;
            surface3.P = P3;
            surface3.My = M3y;
            surface3.Mz = M3z;        
        end

        %% surface 4
        function [surface4] = surface_4(obj)

            analysisOptionsElastic = struct;
            analysisOptionsElastic.controlled_dof_override = 1;
            analysisOptionsElastic.includeInitialGeometricImperfections = false;
            analysisOptionsElastic.reportMaximumLimitpoint = false;  
            obj.applyNotionalLoad = true;
            analysisOptionsElastic.geomTransfType = 'PDelta';
            P4 = nan(obj.numPoints,num_angles(obj));
            M4y = nan(obj.numPoints,num_angles(obj));
            M4z = nan(obj.numPoints,num_angles(obj));
            status = cell(obj.numPoints,num_angles(obj));                    


            % Elsatic buckling load    
            Peb = second_order_analysis(obj,'LimitPoint_Proportional'...
                ,angle(obj,1),analysisOptionsElastic); 
            P0 = -min(obj.psd_ACB.xPa,min(abs(0.99*Peb)));
            P1 = transpose(linspace(P0,0,obj.numPoints));     

            M1y = ones(obj.numPoints,1) * obj.psd_ACB.Mcy/20;
            M1z = ones(obj.numPoints,1) * 0;

            for i = 1:num_angles(obj)
                [~,~,~,P4(:,i),M4y(:,i),M4z(:,i),status(:,i)] = ...
                    second_order_analysis(obj,'TargetForce_NonProportional',angle(obj,i),...
                    analysisOptionsElastic,P1(:,1),M1y(:,1),M1z(:,1));                         
            end        

            P4(1,:) = -min(-P4(1,:));%%%%%%%%

            surface4        = struct;
            surface4.P      = P4;
            surface4.My     = M4y;
            surface4.Mz     = M4z;
            surface4.status = status;

        end   
        %% surface 5
        function surface5 = surface_5(obj,surface1)

            analysisOptionsElastic = struct;
            analysisOptionsElastic.controlled_dof_override = 1;
            analysisOptionsElastic.includeInitialGeometricImperfections = false;
            analysisOptionsElastic.reportMaximumLimitpoint = true;
            analysisOptionsElastic.geomTransfType = 'PDelta';
            obj.applyNotionalLoad = true;
            if obj.reportMaximumLimitpoint
                P5 = nan(obj.numPoints,num_angles(obj));
                M5y = nan(obj.numPoints,num_angles(obj));
                M5z = nan(obj.numPoints,num_angles(obj));
                status = cell(obj.numPoints,num_angles(obj));                    
            else
                P5 = cell(obj.numPoints,num_angles(obj));
                M5y = cell(obj.numPoints,num_angles(obj));
                M5z = cell(obj.numPoints,num_angles(obj));                
                status = cell(obj.numPoints,num_angles(obj));                    
            end        

            % Elsatic buckling load    
            Peb = second_order_analysis(obj,'LimitPoint_Proportional'...
                ,angle(obj,1),analysisOptionsElastic); 

            for i = 1:num_angles(obj)
                [P5(:,i),M5y(:,i),M5z(:,i),~,~,~,status(:,i)] = ...
                    second_order_analysis(obj,'TargetForce_Proportional',angle(obj,i),...
                    analysisOptionsElastic,surface1.P(:,i),surface1.My(:,i),surface1.Mz(:,i),Peb);                         
            end        

            if obj.reportMaximumLimitpoint

                surface5 = struct;
                surface5.P = P5;
                surface5.My = M5y;
                surface5.Mz = M5z;            
                surface5.status = status;            
            else
                surface5 = struct;            
                for i =1:obj.numPoints
                    for j=1:num_angles(obj)
                        for k=1:obj.numElements+1                        
                            surface5(k).Mz(i,j)=M5z{i,j}(1,k);
                            surface5(k).My(i,j)=M5y{i,j}(1,k);
                            surface5(k).P(i,j)=P5{i,j}(1,k);                        
                        end
                    end
                end
            end         
        end


        %% Second order analysis for each angle
        function [P2,M2y,M2z,P1,M1y,M1z,status] = second_order_analysis(...
                obj,analysisType,angle,analysisOptions,P0,M0y,M0z,Peb)



            BA_Driver = BenchmarkAnalysis3d_Driver(obj.frame_data,analysisOptions);


            if isempty(analysisOptions.controlled_dof_override)
                [~,num_imp_angles] = obj.initGeomImperfAng(angle);
                for i=1:num_imp_angles
                    be(i) = BA_Driver.get3dAnalysisObject('nonlinear',angle,obj.initGeomImperfAng(angle,i));
                    scratchFolder = sprintf('%s-%s',obj.frame_data.frame_name,obj.frame_data.section_name);
                    be(i).scratchPath = [pathOf.scratch '/' scratchFolder];
                end
            else       
                [~,num_nl_angles] = obj.notionalLoadAng(angle);
                num_imp_angles = num_nl_angles;
                for i=1:num_nl_angles
                    be(i) = BA_Driver.get3dAnalysisObject('elastic',angle,nan);
                    scratchFolder = sprintf('%s-%s-Elastic',obj.frame_data.frame_name,obj.frame_data.section_name);
                    be(i).scratchPath = [pathOf.scratch '/' scratchFolder];
                end
            end

            mkdir(pathOf.scratch,scratchFolder);
            if strcmp(analysisType,'LimitPoint_Proportional')

                % Initilize results
                status = cell(1,1);
                % Run axial only analysis to get Pn        
                P1_comp = inf(1,num_imp_angles);
                for i = 1:num_imp_angles
                    if obj.extra_fine_DispStepSize
                        try_number = 2;
                    else
                        try_number = 1;
                    end
                    results(i) = be(i).runAnalysis(analysisType,0,[],try_number);
                    if ~results(i).limitPoint.good
                        results(i) = be(i).runAnalysis(analysisType,0,[],try_number+1);
                    end
                    if ~results(i).limitPoint.good
                        results(i) = be(i).runAnalysis(analysisType,0,[],try_number+2);
                    end                    
                    if ~results(i).limitPoint.good
                        error('Limit Point Not Obtained: Axial Only');
                    end
                    P1_comp(i) = results(i).limitPoint.P1;
                end
                [~,ind] = min(abs(P1_comp));
                status{1} = {results(ind).limitPoint.limit_type};
                M2y = results(ind).limitPoint.M2y;
                M2z = results(ind).limitPoint.M2z;
                P1  = results(ind).limitPoint.P1;
                P2  = results(ind).limitPoint.P2;

                % Run nonproportional analyses at different axial loads
            elseif strcmp(analysisType,'LimitPoint_NonProportional')

                % Initilize results
                if obj.reportMaximumLimitpoint
                    M2y = zeros(obj.numPoints,1);
                    M2z = zeros(obj.numPoints,1);
                    P2 = zeros(obj.numPoints,1);
                    P2(1)  = P0;
                    M2y(1) = M0y;
                    M2z(1) = M0z;
                else
                    M2y = cell(obj.numPoints,1);
                    M2z = cell(obj.numPoints,1);
                    P2  = cell(obj.numPoints,1);
                    P2{1}    = ones(1,7)*P0;
                    M2y{1,1} = M0y{1,1};
                    M2z{1,1} = M0z{1,1};
                end
                P1 = zeros(obj.numPoints,1);            
                P1(1) = P0;
                M1y = zeros(obj.numPoints,1);
                M1z = zeros(obj.numPoints,1);
                status = cell(obj.numPoints,1);
                Ps = linspace(P0,0,obj.numPoints);
                for i = 2:length(Ps)
                    % Run LimitPoint NonProportional analysis                       
                    M1_comp = inf(1,num_imp_angles);
                    for j = 1:num_imp_angles           
                        be(j).numStepsGravity = 100;
                        if obj.extra_fine_DispStepSize
                            try_number = 2;
                        else
                            try_number = 1;
                        end
                            
                        results(j) = be(j).runAnalysis(analysisType,Ps(i),0,try_number);
                        if ~results(j).limitPoint.good
                            results(j) = be(j).runAnalysis(analysisType,Ps(i),0,try_number+1);
                        end
                        if ~results(j).limitPoint.good
                            error('Limit Point Not Obtained: Lateral Loading');
                        end
                        M1_comp(j) = sqrt(results(j).limitPoint.M1z^2 + results(j).limitPoint.M1y^2);
                    end                
                    [~,ind] = min(abs(M1_comp));
                    if obj.reportMaximumLimitpoint
                        M2y(i,1) = results(ind).limitPoint.M2y;
                        M2z(i,1) = results(ind).limitPoint.M2z;
                        P2(i,1)  = results(ind).limitPoint.P2;
                    else
                        M2y{i,1} = results(ind).limitPoint.M2y;
                        M2z{i,1} = results(ind).limitPoint.M2z;
                        P2{i,1}  = results(ind).limitPoint.P2;
                    end                    

                    status{i} = results(ind).limitPoint.limit_type;
                    P1(i) = Ps(i);
                    M1y(i,1) = results(ind).limitPoint.M1y;
                    M1z(i,1) = results(ind).limitPoint.M1z;
                    clear results;
                end
            elseif strcmp(analysisType,'LimitPoint_NonProportional_2')

                    % Run LimitPoint NonProportional analysis 
                    M1_comp = inf(1,num_imp_angles);
                    for j = 1:num_imp_angles

                        be(j).numStepsGravity = 100;
                        results(j) = be(j).runAnalysis('LimitPoint_NonProportional',P0,0,2);
                        if ~results(j).limitPoint.good
                            error('Limit Point Not Obtained: Lateral Loading');
                        end
                        M1_comp(j) = sqrt(results(j).limitPoint.M1z^2 + results(j).limitPoint.M1y^2);
                    end                
                    [~,ind] = min(abs(M1_comp));                

                    if obj.reportMaximumLimitpoint
                        M2y = results(ind).limitPoint.M2y;
                        M2z = results(ind).limitPoint.M2z;
                        P2  = results(ind).limitPoint.P2;
                    else
                        M2y = results(ind).limitPoint.M2y;
                        M2z = results(ind).limitPoint.M2z;
                        P2  = results(ind).limitPoint.P2;
                    end                    

                    status{1,1} = results(ind).limitPoint.limit_type;
                    P1 = P0;
                    M1y = results(ind).limitPoint.M1y;
                    M1z = results(ind).limitPoint.M1z; 

            elseif strcmp(analysisType,'TargetForce_Proportional')            

                % Initilize results
                    M1y = zeros(obj.numPoints,1);
                    M1z = zeros(obj.numPoints,1);
                    P1 = zeros(obj.numPoints,1);

                if obj.reportMaximumLimitpoint

                    M2y = zeros(obj.numPoints,1);
                    M2z = zeros(obj.numPoints,1);
                    P2 = zeros(obj.numPoints,1);

                else
                    M2y = cell(obj.numPoints,1);
                    M2z = cell(obj.numPoints,1);
                    P2  = cell(obj.numPoints,1);

                end
                status = cell(obj.numPoints,1);
                M1_comp = inf(1,num_nl_angles);
                for i = 1:obj.numPoints
                    if  abs(P0(i,1)) > abs(Peb)                    
                        if obj.reportMaximumLimitpoint
                            P2(i,1)  = nan;
                            M2y(i,1) = nan;
                            M2z(i,1) =  nan;
                        else
                            P2{i,1}  = ones(1,obj.numElements+1) * nan;
                            M2y{i,1} = ones(1,obj.numElements+1) * nan;
                            M2z{i,1} = ones(1,obj.numElements+1) * nan;
                        end                    
                        continue
                    else

                        % Run Target Force NonProportional analysis       
                        for j = 1:num_nl_angles

                            B2 = obj.driftRatio(P0(i,1));
                            
                            if strcmp(obj.frame_data.frame_typeZ,'Sidesway_Inhibited')  
                                be(j).betaTopZ   =  M0z(i);
                                be(j).betaBotZ   =  be(j).betaTopZ * obj.frame_data.betaBotZ_over_betaTopZ;
                            else
                                H0y(i) = M0z(i)/obj.lever_arm_z;
                                [Y(i,j),~] = obj.HwithNotionalLoad('z',...
                                    H0y(i),P0(i,1),B2,obj.notionalLoadAng(angle,j));                             
                                be(j).Hy   = Y(i,j);
                            end 
                            if strcmp(obj.frame_data.frame_typeY,'Sidesway_Inhibited')
                                be(j).betaTopY   =  M0y(i); 
                                be(j).betaBotY   =  be(j).betaTopY * obj.frame_data.betaBotY_over_betaTopY;                
                            else
                                H0z(i) = M0y(i)/obj.lever_arm_y;
                                [Z(i,j),~] = obj.HwithNotionalLoad('y',...
                                    H0z(i),P0(i,1),B2,obj.notionalLoadAng(angle,j));                            
                                be(j).Hz   = Z(i,j);
                            end

                            results_elastic(j) = be(j).runAnalysis(analysisType...
                              ,P0(i,1),1,1);
                            M1_comp(1,j) = sqrt(results_elastic(j).path.M1y(end)^2 +...
                                results_elastic(j).path.M1z(end)^2);
                        end
                        [~,ind] = min(abs(M1_comp(1,:))); 

                        if obj.reportMaximumLimitpoint

                            M2y(i,1) = results_elastic(ind).path.M2y(end);
                            M2z(i,1) = results_elastic(ind).path.M2z(end);
                            P2(i,1)  = results_elastic(ind).path.P2(end);                    
                        else
                            M2y{i,1} = results_elastic(ind).path.M2y(end,:);
                            M2z{i,1} = results_elastic(ind).path.M2z(end,:);
                            P2{i,1}  = results_elastic(ind).path.P2(end,:);
                        end                    

                        status{i} = results_elastic(ind).exitStatus;
                    end
                end   

            elseif strcmp(analysisType,'TargetForce_NonProportional')

                % Initilize results
                M1y = zeros(obj.numPoints,1);
                M1z = zeros(obj.numPoints,1);
                P1 = zeros(obj.numPoints,1);

                M2y = zeros(obj.numPoints,1);
                M2z = zeros(obj.numPoints,1);
                P2 = zeros(obj.numPoints,1);
                M0 = M0y;
                status = cell(obj.numPoints,1);            
                % Run Target Force NonProportional analysis

                for i = 1:obj.numPoints
                    if i == 1
                        P1_comp = inf(1,num_nl_angles);
                        B2 = obj.driftRatio(P0(i,1));
                        for j = 1:num_nl_angles              

                            be(j).numStepsLateral = 100;
                            be(j).numStepsGravity = 100;
                            if strcmp(obj.frame_data.frame_typeZ,'Sidesway_Inhibited')        
                                be(j).betaTopZ   =  0;
                                be(j).betaBotZ   =  be(j).betaTopZ * obj.frame_data.betaBotZ_over_betaTopZ;               
                            else
                                [Y(j),notionalY(j)] = obj.HwithNotionalLoad('z',...
                                    0,P0(i,1),B2,obj.notionalLoadAng(angle,j));
                                be(j).Hy   = Y(i,j);
                            end 
                            if strcmp(obj.frame_data.frame_typeY,'Sidesway_Inhibited')        
                                be(j).betaTopY   =  0;
                                be(j).betaBotY   =  be(j).betaTopY * obj.frame_data.betaBotY_over_betaTopY;
                            else
                                [Z(j),notionalZ(j)] = obj.HwithNotionalLoad('y',...
                                    0,P0(i,1),B2,obj.notionalLoadAng(angle,j));                              
                                be(j).Hz   = Z(i,j);
                            end

                            results_elastic(j) = be(j).runAnalysis('TargetForce_Proportional'...
                              ,P0(i,1),1,1);                

                            [ratioLimitPoint,results_elastic(j)] = ratio_limitPoint(obj,results_elastic(j),P0(1,1));
                            [min_time,ind] = min(ratioLimitPoint.time);                
                            ratioLimitPointFlag = false;
                            if  strcmp(ratioLimitPoint.limit_type(ind),'Main')
                                ratioLimitPointFlag = true;
                                status{i} = ratioLimitPoint.limit_type{ind};
                                ind_floor(j) = floor(min_time);
                                x(j) = min_time - ind_floor(j);
                                P1_comp(1,j)  = interpolate_vector(results_elastic(j).mainPath.P1,ind_floor(j),x(j));

                            end
                        end
                        if ratioLimitPointFlag
                            [~,ind] = min(abs(P1_comp)); 
                            P1(i,1) = P1_comp(i,ind);

                            if strcmp(obj.frame_data.frame_typeZ,'Sidesway_Inhibited')        
                                M1z(i,1) = interpolate_vector(results_elastic(ind).mainPath.M1z,ind_floor(ind),x(ind));              
                            else
                                M1z(i,1) = interpolate_vector(results_elastic(ind).mainPath.M1z,ind_floor(ind),x(ind))...
                                    - notionalY(ind) * obj.lever_arm_z;                             
                            end 
                            if strcmp(obj.frame_data.frame_typeY,'Sidesway_Inhibited')        
                                M1y(i,1) = interpolate_vector(results_elastic(ind).mainPath.M1y,ind_floor(ind),x(ind));
                            else
                                M1y(i,1) = interpolate_vector(results_elastic(ind).mainPath.M1y,ind_floor(ind),x(ind))...
                                    - notionalZ(ind) * obj.lever_arm_y;                            
                            end                        

                            P0 = transpose(linspace(P1(i,1),0,obj.numPoints));
    %                         ind = abs(Ps1) < obj.psd_ACB.xPc;
    %                         P0 = vertcat(Ps1(1:find(ind,1)-1),-obj.psd_ACB.xPc,Ps1(find(ind,1):end));
                            clear results_elastic
                            continue  
                        end
                    end                                 
                    try_number = 1;
                    while try_number < 400   
                        if strcmp(obj.frame_data.frame_typeZ,'Sidesway_Inhibited')        
                            be(1).betaTopZ   =  try_number * M0(i) * cos(angle);
                            be(1).betaBotZ   =  be(1).betaTopZ * obj.frame_data.betaBotZ_over_betaTopZ;                
                        else 
                            be(1).Hy   = try_number * M0(i) * cos(angle)/obj.lever_arm_z;
                        end 
                        if strcmp(obj.frame_data.frame_typeY,'Sidesway_Inhibited')        
                            be(1).betaTopY   =  try_number * M0(i) * sin(angle);
                            be(1).betaBotY   =  be(1).betaTopY * obj.frame_data.betaBotY_over_betaTopY;
                        else
                            be(1).Hz   = try_number * M0(i) * sin(angle)/obj.lever_arm_y;
                        end                    
                        results_elastic = be(1).runAnalysis(analysisType...
                          ,P0(i,1),1,1);

                        [ratioLimitPoint,results_elastic] = ratio_limitPoint(obj,results_elastic,P0(i,1));

                        [min_time,ind] = min(ratioLimitPoint.time);

                        if strcmp(ratioLimitPoint.limit_type(ind),'Main')
                            status{i} = ratioLimitPoint.limit_type{ind};
                            B2 = obj.driftRatio(P0(i,1));
                            ind = floor(min_time);
                            x = min_time - ind;
                            M1_comp = inf(1,num_nl_angles);
                            for j = 1:num_nl_angles 
                                P1_comp(1,j)  = interpolate_vector(results_elastic.mainPath.P1,ind,x);
                                if strcmp(obj.frame_data.frame_typeZ,'Sidesway_Inhibited')        

                                    M1z_comp(1,j) = interpolate_vector(results_elastic.mainPath.M1z,ind,x);               

                                else
                                    [~,notionalY(j)] = obj.HwithNotionalLoad('z',...
                                        be(1).Hy,P0(i,1),B2,obj.notionalLoadAng(angle,j));
                                    M1z_comp(1,j) = interpolate_vector(results_elastic.mainPath.M1z,ind,x)...
                                        - notionalY(j) * obj.lever_arm_z;                        
                                end 
                                if strcmp(obj.frame_data.frame_typeY,'Sidesway_Inhibited')        
                                    M1y_comp(1,j) = interpolate_vector(results_elastic.mainPath.M1y,ind,x);

                                else
                                    [~,notionalZ(j)] = obj.HwithNotionalLoad('y',...
                                        be(1).Hz,P0(i,1),B2,obj.notionalLoadAng(angle,j)); 
                                    M1y_comp(1,j) = interpolate_vector(results_elastic.mainPath.M1y,ind,x)...
                                        - notionalZ(j) * obj.lever_arm_y;                              
                                end                            

                                M1_comp(1,j) = sqrt(M1y_comp(1,j)^2 + M1z_comp(1,j)^2);

                            end
                            [~,ind] = min(M1_comp(1,:));                
                            M1y(i,1) = M1y_comp(1,ind);
                            M1z(i,1) = M1z_comp(1,ind);                        
                            P1(i,1)  = P1_comp(1,ind);
                            break
                        elseif strcmp(ratioLimitPoint.limit_type(ind),'Gravity')
                            status{i} = ratioLimitPoint.limit_type{ind};
                            B2 = obj.driftRatio(P0(i,1));
                            ind = floor(min_time);
                            x = min_time - ind;
                            for j = 1:num_nl_angles                      
                                P1_comp(1,j)  = interpolate_vector(results_elastic.gravPath.P1,ind,x);
                                if strcmp(obj.frame_data.frame_typeZ,'Sidesway_Inhibited')        

                                    M1z_comp(1,j) = interpolate_vector(results_elastic.gravPath.M1z,ind,x);               

                                else
                                    [~,notionalY(j)] = obj.HwithNotionalLoad('z',...
                                        be(1).Hy,P0(i,1),B2,obj.notionalLoadAng(angle,j));
                                    M1z_comp(1,j) = interpolate_vector(results_elastic.gravPath.M1z,ind,x)...
                                        - notionalY(j) * obj.lever_arm_z;                        
                                end 
                                if strcmp(obj.frame_data.frame_typeY,'Sidesway_Inhibited')        
                                    M1y_comp(1,j) = interpolate_vector(results_elastic.gravPath.M1y,ind,x);

                                else
                                    [~,notionalZ(j)] = obj.HwithNotionalLoad('y',...
                                        be(1).Hz,P0(i,1),B2,obj.notionalLoadAng(angle,j)); 
                                    M1y_comp(1,j) = interpolate_vector(results_elastic.gravPath.M1y,ind,x)...
                                        - notionalZ(j) * obj.lever_arm_y;                              
                                end

                                M1_comp(1,j) = sqrt(M1y_comp(1,j)^2 + M1z_comp(1,j)^2);

                            end
                            [~,ind] = min(M1_comp(1,:));                
                            M1y(i,1) = M1y_comp(1,ind);
                            M1z(i,1) = M1z_comp(1,ind);                        
                            P1(i,1)  = P1_comp(1,ind);
                            break
                        elseif strcmp(ratioLimitPoint.limit_type(ind),'Boundary')
                            status{i} = ratioLimitPoint.limit_type{ind};
                            B2 = obj.driftRatio(P0(i,1));
                            ind = floor(min_time);
                            x = min_time - ind;
                            for j = 1:num_nl_angles  
                                P1_comp(1,j)  = interpolate_vector(results_elastic.path.P1,ind,x);
                                if strcmp(obj.frame_data.frame_typeZ,'Sidesway_Inhibited')        

                                    M1z_comp(1,j) = interpolate_vector(results_elastic.path.M1z,ind,x);               

                                else
                                    [~,notionalY(j)] = obj.HwithNotionalLoad('z',...
                                        be(1).Hy,P0(i,1),B2,obj.notionalLoadAng(angle,j));
                                    M1z_comp(1,j) = interpolate_vector(results_elastic.path.M1z,ind,x)...
                                        - notionalY(j) * obj.lever_arm_z;                        
                                end 
                                if strcmp(obj.frame_data.frame_typeY,'Sidesway_Inhibited')        
                                    M1y_comp(1,j) = interpolate_vector(results_elastic.path.M1y,ind,x);

                                else
                                    [~,notionalZ(j)] = obj.HwithNotionalLoad('y',...
                                        be(1).Hz,P0(i,1),B2,obj.notionalLoadAng(angle,j)); 
                                    M1y_comp(1,j) = interpolate_vector(results_elastic.path.M1y,ind,x)...
                                        - notionalZ(j) * obj.lever_arm_y;                              
                                end
                                M1_comp(1,j) = sqrt(M1y_comp(1,j)^2 + M1z_comp(1,j)^2);
                            end
                            [~,ind] = min(M1_comp(i,:));                
                            M1y(i,1) = M1y_comp(1,ind);
                            M1z(i,1) = M1z_comp(1,ind);                        
                            P1(i,1)  = P1_comp(1,ind);
                            break
                        else
                            try_number = try_number * 2;
                            M1y(i,1) = nan;
                            M1z(i,1) = nan;
                            P1(i,1)  = nan;                      
                        end
                    end
                end
            end
            % Delete scratch folder
            try
                rmdir(be(1).scratchPath);
            end
        end    

        function [ratioLimitPoint,results_elastic] = ratio_limitPoint(obj,results_elastic,P0)

            for ii = 1:obj.numElements+1                

                if isfield(results_elastic,'gravPath') && P0 ~= 0
                    Pr_grav =  results_elastic.gravPath.P2(:,ii);

                    % Calculate ratio in gravity loading
                    ind_Pc_grav = abs(Pr_grav) < obj.psd_ACB.xPc;
                    Mzr = results_elastic.gravPath.M2z(ind_Pc_grav,ii);
                    Myr = results_elastic.gravPath.M2y(ind_Pc_grav,ii);

                    results_elastic.gravPath.ratio(ind_Pc_grav,ii) = ((Myr/obj.psd_ACB.Mcy).^alpha(obj)...
                        + (Mzr/obj.psd_ACB.Mcz).^alpha(obj)).^(1/alpha(obj));
                    ind_Pc_grav = abs(Pr_grav) >= obj.psd_ACB.xPc;

                    Mzr = results_elastic.gravPath.M2z(ind_Pc_grav,ii);
                    Myr = results_elastic.gravPath.M2y(ind_Pc_grav,ii);                
                    Pr =  abs(results_elastic.gravPath.P2(ind_Pc_grav,ii));

                    results_elastic.gravPath.ratio(ind_Pc_grav,ii) = (Pr-obj.psd_ACB.xPc)/...
                        (obj.psd_ACB.xPa-obj.psd_ACB.xPc)+((Myr/obj.psd_ACB.Mcy).^alpha(obj) +...
                        (Mzr/obj.psd_ACB.Mcz).^alpha(obj)).^(1/alpha(obj));
                end

                % Calculate ratio in main loading

                Pr_main =  results_elastic.mainPath.P2(:,ii);


                ind_Pc_main = abs(Pr_main) < obj.psd_ACB.xPc;
                Mzr = results_elastic.mainPath.M2z(ind_Pc_main,ii);
                Myr = results_elastic.mainPath.M2y(ind_Pc_main,ii);

                results_elastic.mainPath.ratio(ind_Pc_main,ii) = ((Myr/obj.psd_ACB.Mcy).^alpha(obj)...
                    + (Mzr/obj.psd_ACB.Mcz).^alpha(obj)).^(1/alpha(obj));

                ind_Pc_main = abs(Pr_main) >= obj.psd_ACB.xPc;

                Mzr = results_elastic.mainPath.M2z(ind_Pc_main,ii);
                Myr = results_elastic.mainPath.M2y(ind_Pc_main,ii);                
                Pr =  abs(results_elastic.mainPath.P2(ind_Pc_main,ii));

                results_elastic.mainPath.ratio(ind_Pc_main,ii) = (Pr-obj.psd_ACB.xPc)/...
                    (obj.psd_ACB.xPa-obj.psd_ACB.xPc)+((Myr/obj.psd_ACB.Mcy).^alpha(obj) +...
                    (Mzr/obj.psd_ACB.Mcz).^alpha(obj)).^(1/alpha(obj));                


                % Calculate ratio in combined loading

                Pr_all =  results_elastic.path.P2(:,ii);
                ind_Pc = abs(Pr_all) < obj.psd_ACB.xPc;
                Mzr = results_elastic.path.M2z(ind_Pc,ii);
                Myr = results_elastic.path.M2y(ind_Pc,ii);

                results_elastic.path.ratio(ind_Pc,ii) = ((Myr/obj.psd_ACB.Mcy).^alpha(obj)...
                    + (Mzr/obj.psd_ACB.Mcz).^alpha(obj)).^(1/alpha(obj));

                ind_Pc = abs(Pr_all) >= obj.psd_ACB.xPc;

                Mzr = results_elastic.path.M2z(ind_Pc,ii);
                Myr = results_elastic.path.M2y(ind_Pc,ii);                
                Pr =  abs(results_elastic.path.P2(ind_Pc,ii));

                results_elastic.path.ratio(ind_Pc,ii) = (Pr-obj.psd_ACB.xPc)/...
                    (obj.psd_ACB.xPa-obj.psd_ACB.xPc)+((Myr/obj.psd_ACB.Mcy).^alpha(obj) +...
                    (Mzr/obj.psd_ACB.Mcz).^alpha(obj)).^(1/alpha(obj));                

                if isfield(results_elastic,'gravPath') && P0 ~= 0
                    [ind_grav,x_grav] = find_limit_point_in_vector(results_elastic.gravPath.ratio(:,ii),1.0);
                else
                    ind_grav = [];
                end
                [ind_main,x_main] = find_limit_point_in_vector(results_elastic.mainPath.ratio(:,ii),1.0);
                [ind_all,x_all] = find_limit_point_in_vector(results_elastic.path.ratio(:,ii),1.0);

                if ~isempty(ind_grav)
                    ratioLimitPoint.limit_type{1,ii} = 'Gravity';
                    ratioLimitPoint.good(1,ii) = true;
                    ratioLimitPoint.time(1,ii) = ind_grav + x_grav;
                elseif ~isempty(ind_main)
                    ratioLimitPoint.limit_type{1,ii} = 'Main';
                    ratioLimitPoint.good(1,ii) = true;
                    ratioLimitPoint.time(1,ii) = ind_main + x_main;
                elseif ~isempty(ind_all)
                    ratioLimitPoint.limit_type{1,ii} = 'Boundary';
                    ratioLimitPoint.good(1,ii) = true;
                    ratioLimitPoint.time(1,ii) = ind_all + x_all;                
                else
                    ratioLimitPoint.limit_type{1,ii} = 'Not Reached';
                    ratioLimitPoint.good(1,ii) = false;
                    ratioLimitPoint.time(1,ii) = inf;                

                end

            end  
        end
        function lever_arm_y = lever_arm_y(obj)
            if strcmp(obj.frame_data.frame_typeY,'Sidesway_Inhibited')                      
                lever_arm_y = 1;
            elseif strcmp(obj.frame_data.frame_typeY,'Sidesway_Uninhibited') 
                if (obj.frame_data.kqbotY ~= 0 && obj.frame_data.kqtopY ~= 0 && ...
                    obj.frame_data.kqtopY == obj.frame_data.kqbotY)
                  lever_arm_y = obj.frame_data.L/2 ;
                elseif (obj.frame_data.kqtopY == 0 || obj.frame_data.kqbotY == 0)
                  lever_arm_y = obj.frame_data.L ;
                else
                  error('lever_arm_y is not defined');
                end  
            end
        end
        function lever_arm_z = lever_arm_z(obj)

            if strcmp(obj.frame_data.frame_typeZ,'Sidesway_Inhibited')        
                lever_arm_z = 1;
            elseif strcmp(obj.frame_data.frame_typeZ,'Sidesway_Uninhibited')                         
                if (obj.frame_data.kqbotZ ~= 0 && obj.frame_data.kqtopZ ~= 0 &&...
                    obj.frame_data.kqtopZ == obj.frame_data.kqbotZ)
                  lever_arm_z = obj.frame_data.L/2 ;
                elseif (obj.frame_data.kqtopZ == 0 || obj.frame_data.kqbotZ == 0)
                  lever_arm_z = obj.frame_data.L ;
                else
                  error('lever_arm_z is not defined');
                end
            end
        end    


        function ind_lambda = ind_lambda(obj,P1max)
            PnOverPo=min(abs(P1max/obj.psd_ACB.PNO),1);        
            lambda_oe = sqrt(AISC_inverse_column_curve(PnOverPo));

            if (lambda_oe <= 0.25)
                        ind_lambda=1;
                  elseif ( lambda_oe > 0.25 ) && (lambda_oe <= 0.5)
                        ind_lambda=2;
                  elseif (lambda_oe > 0.5) && (lambda_oe <= 0.75)
                        ind_lambda=3;
                  elseif (lambda_oe > 0.75) && (lambda_oe <= 1.0)
                        ind_lambda=4;
                  elseif (lambda_oe > 1.0) && (lambda_oe <= 1.5)
                        ind_lambda=5;
                  elseif lambda_oe > 1.5
                        ind_lambda=6;
            else 
                        ind_lambda= nan;
            end        
        end
        function ind_rho = ind_rho(obj)
            rho_s =  round(100*obj.frame_data.section.getSectionData('steelratio'))/100;
            switch obj.frame_data.section_type  
                case 'SRC'
                    rho_s_range = [0.01,0.04,0.09,0.12];
                switch rho_s
                      case rho_s_range(1)
                            ind_rho=1;
                      case rho_s_range(2)
                            ind_rho=2;
                      case rho_s_range(3)
                            ind_rho=3;
                      case rho_s_range(4)
                            ind_rho=4;
                    otherwise 
                            ind_rho = nan;
                end   
                case 'RCFT'
                    rho_s_range = [0.03,0.05,0.10,0.11,0.19,0.195,0.28];%%%%%%%%%%%%%%%%%%%%%%%%%%
                switch rho_s
                      case rho_s_range(1)
                            ind_rho=1;
                      case rho_s_range(2)
                            ind_rho=2;
                      case rho_s_range(3)
                            ind_rho=3;
                      case rho_s_range(4)
                            ind_rho=3;                            
                      case rho_s_range(5)
                            ind_rho=4;
                      case rho_s_range(6)
                            ind_rho=4;                            
                      case rho_s_range(7)
                            ind_rho=5;
                      otherwise 
                           error('ro error')              
                end  
            end            
        end
        function results = curve1_5_2d(obj)

            % OpenSees Analysis Options
            numPointsCurve12 = obj.numPoints;

            analysisOptions = struct;
            analysisOptions.absoluteStrainLimit = 0.05;


            % psd_ACB Method Options
            designStrengthType          = 'AISC';
            numPointsCurve4             = obj.numPoints;
            notionalLoadObject          = notional_load(0.000,0.002,1.7);
            columnStiffnessReduction    = 0.8;
            beamStiffnessReduction      = 0.8;
            tauType                     = 'Composite';
            Py_tau                      = nan;
            effectiveLengthFactorType   = 'one';
            EI                          = obj.frame_data.section.EIeff(obj.axis);

            %% Curves 1 and 2 (OpenSees Analysis)
            BA_Driver = BenchmarkAnalysis3d_Driver(obj.frame_data,analysisOptions);
            [ba,data] = BA_Driver.get2dAnalysisObject('nonlinear',obj.axis);
            ba.deleteFilesAfterAnalysis = true;

            % Initilize results
            Curve1_limit_type = cell(numPointsCurve12,1);
            results.Curve1_M1   = nan(numPointsCurve12,1);
            results.Curve1_P1   = nan(numPointsCurve12,1);
            results.Curve2_M2   = nan(numPointsCurve12,1);
            results.Curve2_P2   = nan(numPointsCurve12,1);
            results.Curve1_def  = nan(numPointsCurve12,1);

            % Run axial only analysis to get Pn
            iResults = ba.runAnalysis('LimitPoint_Proportional',0,[],1);

            if ~iResults.limitPoint.good
                iResults = ba.runAnalysis('LimitPoint_Proportional',0,[],2);
            end

            if ~iResults.limitPoint.good
                error('Limit Point Not Obtained: Axial Only');
            end

            Curve1_limit_type{1}   = iResults.limitPoint.limit_type;
            results.Curve1_M1(1)   = iResults.limitPoint.M1;
            results.Curve1_P1(1)   = iResults.limitPoint.P1;
            results.Curve2_M2(1)   = iResults.limitPoint.M2;
            results.Curve2_P2(1)   = iResults.limitPoint.P2;
            results.Curve1_def(1)  = iResults.limitPoint.def;

            % Run nonproportional analyses at different axial loads
            Ps = linspace(results.Curve1_P1(1),0,numPointsCurve12);
            for i = 2:length(Ps)
                P = Ps(i);
                if P ~= 0
                    iResults = ba.runAnalysis('LimitPoint_NonProportional',P,[],1);

                    if ~iResults.limitPoint.good
                        iResults = ba.runAnalysis('LimitPoint_NonProportional',P,[],2);
                    end

                    if ~iResults.limitPoint.good
                        iResults = ba.runAnalysis('LimitPoint_NonProportional',P,[],3);
                    end

                    if ~iResults.limitPoint.good
                        error('Limit Point Not Obtained: Axial Load Level %i (%g)',i,P);
                    end

                    Curve1_limit_type{i} = iResults.limitPoint.limit_type;
                    results.Curve1_M1(i)   = iResults.limitPoint.M1;
                    results.Curve1_P1(i)   = iResults.limitPoint.P1;
                    results.Curve2_M2(i)   = iResults.limitPoint.M2;
                    results.Curve2_P2(i)   = iResults.limitPoint.P2;
                    results.Curve1_def(i)  = iResults.limitPoint.def;
                else
                    iResults = ba.runSectionAnalysis('LimitPoint_NonProportional',0,[],1);

                    if ~iResults.limitPoint.good
                        error('Limit Point Not Obtained: Axial Load Level %i (%g)',i,P);
                    end

                    Curve1_limit_type{i}   = iResults.limitPoint.limit_type;
                    results.Curve1_M1(i)   = iResults.limitPoint.M1;
                    results.Curve1_P1(i)   = 0;
                    results.Curve2_M2(i)   = iResults.limitPoint.M2;
                    results.Curve2_P2(i)   = 0;
                    results.Curve1_def(i)  = NaN;
                end
            end

            %% Curve 3
            [results.Curve3_P2,results.Curve3_M2] = obj.frame_data.section.beamColumnInteraction2d(obj.axis,designStrengthType,'CompPos');

            %% Curve 4
            try
                data.beta = 1/data.beta;
            end
            elastic_frame = BenchmarkAnalysis2d_Elastic(data,EI);  
            elastic_frame.includeInitialGeometricImperfections = false;

            [results.Curve4_P1,results.Curve4_M1,results.Curve3_P2b,results.Curve3_M2b] = elastic_frame.designInteraction(...
                obj.frame_data.section,obj.axis,...
                designStrengthType,numPointsCurve4,...
                notionalLoadObject,...
                effectiveLengthFactorType,...
                columnStiffnessReduction,beamStiffnessReduction,tauType);

            %% Curve 5
            [results.Curve5_P2,results.Curve5_M2] = elastic_frame.impliedInteraction(results.Curve1_P1,results.Curve1_M1,...
                notionalLoadObject,...
                columnStiffnessReduction,beamStiffnessReduction,tauType,Py_tau);
        end
    end    
end









