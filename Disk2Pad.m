% CONVERT DISK TO PAD
function padCf = Disk2Pad(Cf,elecName,elecArray,Nelec,numTargs)
padCf = zeros(Nelec,numTargs);
for e = 1:Nelec
    curElec = elecName{e};
    pos = curElec(regexp(curElec,'(?!z)\D'));
    if regexp(curElec,'z'); num = 0; address = {'1' 'z' '2'}; else; ...
            num = str2double(curElec(regexp(curElec,'\d')));
        if num == 1; address = {'z' '1' '3'}; elseif num == 2; address = {'z' '2' '4'}; ...
        else; address = (num-2):2:num+2; end
    end
    if strcmp(pos,'Fp')
        archy = [elecArray{find(cell2mat(cellfun(@sum,cellfun(@(x) strcmp(x,pos),elecArray,'uni',0),'uni',0))):find(cell2mat(cellfun(@sum,cellfun(@(x) strcmp(x,pos),elecArray,'uni',0),'uni',0)))+1}];
    elseif strcmp(pos,'O')
        archy = [elecArray{find(cell2mat(cellfun(@sum,cellfun(@(x) strcmp(x,pos),elecArray,'uni',0),'uni',0)))-1:find(cell2mat(cellfun(@sum,cellfun(@(x) strcmp(x,pos),elecArray,'uni',0),'uni',0)))}];
    else
        archy = [elecArray{find(cell2mat(cellfun(@sum,cellfun(@(x) strcmp(x,pos),elecArray,'uni',0),'uni',0)))-1:find(cell2mat(cellfun(@sum,cellfun(@(x) strcmp(x,pos),elecArray,'uni',0),'uni',0)))+1}];
    end
    streets = elecName(ismember(cellfun(@(x) x(regexp(x,'(?!z)\D')),elecName,'uni',0),archy));
    if ismember('z',address)
        neighbors = streets(ismember(cellfun(@(x) x(regexp(x,'[z]|\d')),streets,'uni',0),address));
    else
        neighbors = streets(ismember(str2double(cellfun(@(x) x(regexp(x,'\d')),streets,'uni',0)),address));
    end
    padIdx = ismember(elecName,neighbors);
    padCf(e,:) = sum(Cf(padIdx,:));
    disp(['Created Pad Electrode for: ' curElec ' with ' num2str(sum(padIdx)) ' electrodes ...'])
end