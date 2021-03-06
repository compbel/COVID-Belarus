% The inputs are:
% treeFile: the file with ML tree
% metaFile: the file with metadata from GISAID
% tripsFile: the file with the human mobility statistics provided by 
% European Commission Knowledge Center on Migration and Demography (KCMD)
% Available at KCMD Dynamic Data Hub
% cntCodeFile: the file with the decoding of 3-letter country codes from
% the previous file. Available at KCMD Dynamic Data Hub

clear;
treeFile = 'tree.nwk';
metaFile = 'nextstrain__metadata.xlsx';
tripsFile = 'KCMD_DDH_data_KCMD-EUI GMP_ Estimated trips.xlsx';
cntCodeFile = 'KCMD_DDH_meta_KCMD-EUI GMP_ Estimated trips.xlsx';

metadata = readtable(metaFile);
codes = readtable(cntCodeFile);
trips = readtable(tripsFile);

names = cellstr(metadata.strain);
dates = datetime(metadata.date);

%% 
% extract country names for sequences and unify the names of some geographical locations
% for GISAID and KCMD Dynamic Data Hub
n = length(names);
countries = cell(1,n);
for i = 1:n
    C = split(names{i},'/');
    if lower(names{i}(1)) == names{i}(1) || startsWith(names{i},'ENV')
        cnt = C{2};
    else
        cnt = C{1};
    end
    if strcmp(cnt,'USA') || strcmp(cnt,'DRC')
        countries{i} = cnt;
        continue;
    end
    countries{i} = cnt(1);
    for j = 2:length(cnt)
        if lower(cnt(j)) == cnt(j)
            countries{i} = [countries{i} cnt(j)];
        else
            countries{i} = [countries{i} ' ' cnt(j)];
        end
    end
    if strcmp(countries{i},'England') || strcmp(countries{i},'Scotland') || strcmp(countries{i},'Wales') || strcmp(countries{i},'Northern Ireland')
        countries{i} = 'United Kingdom';
    end
    if strcmp(countries{i},'Zhejiang') || strcmp(countries{i},'Jiangsu') || strcmp(countries{i},'Hefei') || strcmp(countries{i},'Hangzhou') || strcmp(countries{i},'Anhui') || strcmp(countries{i},'Fujian') || strcmp(countries{i},'Shulan') || strcmp(countries{i},'Urumqi') || strcmp(countries{i},'Wuhan')
        countries{i} = 'China';
    end
    if strcmp(countries{i},'Guangdong') || strcmp(countries{i},'Foshan') || strcmp(countries{i},'Jiangxi') || strcmp(countries{i},'Beijing') || strcmp(countries{i},'Sichuan') || strcmp(countries{i},'Yunnan') || strcmp(countries{i},'Qingdao')  || strcmp(countries{i},'Shanghai') 
        countries{i} = 'China';
    end
    if strcmp(countries{i},'Chongqing') || strcmp(countries{i},'Shenzhen') || strcmp(countries{i},'Shandong') || strcmp(countries{i},'Jingzhou') || strcmp(countries{i},'Tianmen') || strcmp(countries{i},'Guangzhou') || strcmp(countries{i},'Harbin') || strcmp(countries{i},'Liaoning')
        countries{i} = 'China';
    end
    if strcmp(countries{i},'Bahrain')
        countries{i} = 'Bahrein';
    end
    if strcmp(countries{i},'Bucuresti')
        countries{i} = 'Romania';
    end
end

cntList = unique(countries);
nCnt = length(cntList);

% find index of each country code in cntList
codeInd = zeros(1,size(codes,1));
for i = 1:size(codes,1)
    ind = find(strcmp(cntList,codes.country{i}));
    if ~isempty(ind)
        codeInd(i) = ind;
    end
end
%% 
% calculate transition matrix for traits
Q = zeros(nCnt,nCnt);
for i = 1:size(trips,1)
    i
    code1 = trips.reportingCountry(i);
    code2 = trips.secondaryCountry(i);
    ind1 = codeInd(find(strcmp(codes.code,code1)));
    ind2 = codeInd(find(strcmp(codes.code,code2)));
    if (ind1 ~= 0) && (ind2 ~= 0)
        Q(ind1,ind2) = Q(ind1,ind2) + trips.value(i);
    end
end
yearsObs = 6;
Q = Q/yearsObs;
Q = Q./(sum(Q,2));
Q(isnan(Q)) = 0;
ind = find(Q == 0);
eps = 0.1*min(Q(Q > 0));
Q(ind) = eps;
[centr,val] = eigs(Q',1);
centr = centr'/sum(centr);
Q = Q - diag(diag(Q));
Q(1:nCnt+1:end) = -sum(Q,2);
Q = Q';
%% 
% Estimate maximum joint likelihood traits for internal nodes of the ML
% phylogeny
tree = phytreeread(treeFile);
[AM,WM] = phytree2graph(tree,get(tree,'NumNodes'),dates);
[L,C,labels] = pupko_log(AM,WM,cntList,countries,Q,[],centr);

% output: labels - countries assigned to nodes of the tree
% L and C: dynamic programming likelihoods and labels calculated by the
% algorithm from Pupko et.al., A  fast  algorithm  for  joint  reconstruction  
%of  ancestral  amino  acid  sequences. Molecular biology and evolution,49617(6):890?896, 2000.
