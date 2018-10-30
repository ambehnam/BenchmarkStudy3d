function data = build_data_from_sections_3d(sections)


frames = build_frames;
numSections = length(sections);
numFrames   = length(frames);
numData     = numSections*numFrames;


data(numData) = struct;
for iSection = 1:numSections
    EIeffZ  = sections(iSection).section.EIeff('Strong');
    EIeffY  = sections(iSection).section.EIeff('Weak');
    
    for iFrame = 1:numFrames
        
        H = sections(iSection).section.depth('Strong');
        B = sections(iSection).section.depth('Weak');
        L = frames(iFrame).lambda*sqrt(H*B);
        
        iData = iFrame + (iSection-1)*numFrames;
        
        section_fields = fields(sections(iSection));
        for i = 1:length(section_fields)
            data(iData).(section_fields{i}) = sections(iSection).(section_fields{i});
        end
        data(iData).section_id  = iSection;
        data(iData).frame_id    = iFrame;
        data(iData).frame_name  = frames(iFrame).frame_name;
        data(iData).L           = L;
        
        data(iData).frame_typeZ = frames(iFrame).frame_typeZ;
        switch data(iData).frame_typeZ
            case 'Sidesway_Uninhibited'
                data(iData).kqtopZ      = (6*EIeffZ)/(frames(iFrame).GtopZ*L);
                data(iData).kqbotZ      = (6*EIeffZ)/(frames(iFrame).GbotZ*L);
                data(iData).gammaZ      = frames(iFrame).gammaZ;
                data(iData).Hy          = 1;
                
            case 'Sidesway_Inhibited'
                data(iData).betaBotZ_over_betaTopZ = frames(iFrame).betaBotZ_over_betaTopZ;
                
            otherwise
                error('Unknown frame type: %s',data(iData).frame_typeZ)
        end        
        
        data(iData).frame_typeY = frames(iFrame).frame_typeY;
        switch data(iData).frame_typeY
            case 'Sidesway_Uninhibited'
                data(iData).kqtopY      = (6*EIeffY)/(frames(iFrame).GtopY*L);
                data(iData).kqbotY      = (6*EIeffY)/(frames(iFrame).GbotY*L);
                data(iData).gammaY      = frames(iFrame).gammaY;
                data(iData).Hz          = 1;
                
            case 'Sidesway_Inhibited'
                data(iData).betaBotY_over_betaTopY = frames(iFrame).betaBotY_over_betaTopY;
                
            otherwise
                error('Unknown frame type: %s',data(iData).frame_typeY)
        end              

        % Set length within section object
        data(iData).section.Lx = L;
        data(iData).section.Ly = L;
    end    
end

end


function frames = build_frames
iFrame = 1;

% Group 1
lambda = [12 24 36];
betaBotZ_over_betaTopZ = [0.5 -1.0];
betaBotY_over_betaTopY = [0.5 -1.0];
% Define Frames
for i = 1:length(lambda)
    for j = 1:length(betaBotZ_over_betaTopZ)
        for k = 1:length(betaBotY_over_betaTopY)
            frames(iFrame).lambda      = lambda(i);            
            
            frames(iFrame).frame_typeZ = 'Sidesway_Inhibited';
            frames(iFrame).betaBotZ_over_betaTopZ = betaBotZ_over_betaTopZ(j);

            frames(iFrame).frame_typeY = 'Sidesway_Inhibited';            
            frames(iFrame).betaBotY_over_betaTopY = betaBotY_over_betaTopY(k);
            
            frames(iFrame).frame_name = sprintf('G1(A)-L_%i-%i',...
                frames(iFrame).lambda,iFrame);
            
            iFrame = iFrame+1;
        end
    end
end



% Group 2
lambda = [6 12 24];
EndRestraintPairZ = 'ABCD';
EndRestraintPairY = 'ABCD';
gamma = 0;
% Define Frames
for i = 1:length(lambda)
    for j = 1:length(EndRestraintPairZ)
        for k = 1:length(EndRestraintPairY)
            frames(iFrame).lambda       = lambda(i);            
            
            frames(iFrame).frame_typeZ  = 'Sidesway_Uninhibited';
            [Gtop,Gbot] = EndRestraintPair(EndRestraintPairZ(j),gamma);
            frames(iFrame).GtopZ        = Gtop;
            frames(iFrame).GbotZ        = Gbot;
            frames(iFrame).gammaZ       = gamma;

            frames(iFrame).frame_typeY  = 'Sidesway_Uninhibited';
            [Gtop,Gbot] = EndRestraintPair(EndRestraintPairY(k),gamma);
            frames(iFrame).GtopY        = Gtop;
            frames(iFrame).GbotY        = Gbot;
            frames(iFrame).gammaY       = gamma;
            
            frames(iFrame).frame_name = sprintf('G2(B)-L_%i-%i',...
                frames(iFrame).lambda,iFrame);
            
            iFrame = iFrame+1;
        end
    end
end



% Group 3
lambda = [6 12];
EndRestraintPairZ = 'ABCD';
EndRestraintPairY = 'ABCD';
gamma = [1 3];
% Define Frames
for i = 1:length(lambda)
    for j = 1:length(EndRestraintPairZ)
        for k = 1:length(EndRestraintPairY)
            for l = 1:length(gamma)
                frames(iFrame).lambda       = lambda(i);            

                frames(iFrame).frame_typeZ  = 'Sidesway_Uninhibited';
                [Gtop,Gbot] = EndRestraintPair(EndRestraintPairZ(j),gamma(l));
                frames(iFrame).GtopZ        = Gtop;
                frames(iFrame).GbotZ        = Gbot;
                frames(iFrame).gammaZ       = gamma(l);

                frames(iFrame).frame_typeY  = 'Sidesway_Uninhibited';
                [Gtop,Gbot] = EndRestraintPair(EndRestraintPairY(k),gamma(l));
                frames(iFrame).GtopY        = Gtop;
                frames(iFrame).GbotY        = Gbot;
                frames(iFrame).gammaY       = gamma(l);

                frames(iFrame).frame_name = sprintf('G3(B)-L_%i-%i',...
                    frames(iFrame).lambda,iFrame);

                iFrame = iFrame+1;
            end
        end
    end
end


% Group 4
lambda = 24;
EndRestraintPairZ = 'AB';
EndRestraintPairY = 'AB';
gamma = [1 3];
% Define Frames
for i = 1:length(lambda)
    for j = 1:length(EndRestraintPairZ)
        for k = 1:length(EndRestraintPairY)
            for l = 1:length(gamma)
                frames(iFrame).lambda       = lambda(i);            

                frames(iFrame).frame_typeZ  = 'Sidesway_Uninhibited';
                [Gtop,Gbot] = EndRestraintPair(EndRestraintPairZ(j),gamma(l));
                frames(iFrame).GtopZ        = Gtop;
                frames(iFrame).GbotZ        = Gbot;
                frames(iFrame).gammaZ       = gamma(l);

                frames(iFrame).frame_typeY  = 'Sidesway_Uninhibited';
                [Gtop,Gbot] = EndRestraintPair(EndRestraintPairY(k),gamma(l));
                frames(iFrame).GtopY        = Gtop;
                frames(iFrame).GbotY        = Gbot;
                frames(iFrame).gammaY       = gamma(l);

                frames(iFrame).frame_name = sprintf('G4(B)-L_%i-%i',...
                    frames(iFrame).lambda,iFrame);

                iFrame = iFrame+1;
            end
        end
    end
end



% Group 5
lambda = [6 12 24];
EndRestraintPairZ = 'ABCD';
gammaZ = [0];
betaBotY_over_betaTopY = [0.5 -1.0];
% Define Frames
for i = 1:length(lambda)
    for j = 1:length(EndRestraintPairZ)
        for k = 1:length(gammaZ)
            for l = 1:length(betaBotY_over_betaTopY)
                frames(iFrame).lambda       = lambda(i);            

                frames(iFrame).frame_typeZ  = 'Sidesway_Uninhibited';
                [Gtop,Gbot] = EndRestraintPair(EndRestraintPairZ(j),gammaZ(k));
                frames(iFrame).GtopZ        = Gtop;
                frames(iFrame).GbotZ        = Gbot;
                frames(iFrame).gammaZ       = gammaZ(k);

                frames(iFrame).frame_typeY = 'Sidesway_Inhibited';            
                frames(iFrame).betaBotY_over_betaTopY = betaBotY_over_betaTopY(l);

                frames(iFrame).frame_name = sprintf('G5(C)-L_%i-%i',...
                    frames(iFrame).lambda,iFrame);

                iFrame = iFrame+1;
            end
        end
    end
end


% Group 6
lambda = [6 12 24];
EndRestraintPairZ = 'A';
gammaZ = [1 3];
betaBotY_over_betaTopY = [0.5 -1.0];
% Define Frames
for i = 1:length(lambda)
    for j = 1:length(EndRestraintPairZ)
        for k = 1:length(gammaZ)
            for l = 1:length(betaBotY_over_betaTopY)
                frames(iFrame).lambda       = lambda(i);            

                frames(iFrame).frame_typeZ  = 'Sidesway_Uninhibited';
                [Gtop,Gbot] = EndRestraintPair(EndRestraintPairZ(j),gammaZ(k));
                frames(iFrame).GtopZ        = Gtop;
                frames(iFrame).GbotZ        = Gbot;
                frames(iFrame).gammaZ       = gammaZ(k);

                frames(iFrame).frame_typeY = 'Sidesway_Inhibited';            
                frames(iFrame).betaBotY_over_betaTopY = betaBotY_over_betaTopY(l);

                frames(iFrame).frame_name = sprintf('G6(C)-L_%i-%i',...
                    frames(iFrame).lambda,iFrame);

                iFrame = iFrame+1;
            end
        end
    end
end

end

function [Gtop,Gbot] = EndRestraintPair(pair,gamma)
switch pair
    case 'A'
        Gtop = 0;
        Gbot = 0;
    case 'B'
        if gamma == 0 
            Gtop = 3;
            Gbot = 3;        
        else
            Gtop = 1;
            Gbot = 1;        
        end
    case 'C'
        Gtop = 0;
        Gbot = Inf;
    case 'D'
        if gamma == 0 
            Gtop = 3;
        else
            Gtop = 1;
        end
        Gbot = Inf;
    otherwise
        error('Unknown pair: %s',pair)
end

end