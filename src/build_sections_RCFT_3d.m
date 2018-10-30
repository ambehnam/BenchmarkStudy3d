function sections = build_sections_RCFT_3d()

iSection = 1;
fc = [4 10];

Fy = 46;
Fu = 58;
chosenSections = [
    6.0 6.0 1/2
    9.0 9.0 1/2
    8.0 8.0 1/4
    9.0 9.0 1/8
    14.0 14.0 1/8
    10.0 5.0 3/8
    10.0 5.0 3/16];
D = chosenSections(:,1);
B = chosenSections(:,2);
t = 0.93*chosenSections(:,3);

sections(length(D)*length(fc)) = struct;
for i = 1:length(D)
    for j = 1:length(fc)
        sections(iSection).section = RCFT(D(i),B(i),t(i),Fy,fc(j),'US');
        sections(iSection).section.Fu = Fu;
        sections(iSection).section.neglectLocalBuckling = true;
        sections(iSection).section.Ec = 57*sqrt(1000*fc(j));
        sections(iSection).section_type = 'RCFT';
        sections(iSection).section_name = ...
            sprintf('RCFT-%s-%i',upper(listLetter(i)),fc(j));
        sections(iSection).units        = 'US';
        iSection = iSection+1;
    end
end

end