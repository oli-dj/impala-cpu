function [ list ] = multigrid_list( TI, tau, gridlevel)
%MULTIGRID_LIST
%
%   TI: Training image (2D or 3D)
%   tau: data template
%   gridlevel: (integer) (2^n) n = 0,1,2...
%

cat = unique(TI(:))';
num_cat = length(cat);
list = [];

%2D or 3D
switch length(size(TI))
    case 2 %2D
        tau = tau(:,1:2);
        template_length = size(tau,1);
        
        if gridlevel == 1
            list = populate_impala_list( TI, tau);
        else
            for x = 1:gridlevel
                for y = 1:gridlevel
                    list = [list ; populate_impala_list(...
                        TI(x:gridlevel:end,y:gridlevel:end),tau)];
                end
            end
        end
    case 3 %3D
        tau = tau(:,1:3);
        template_length = size(tau,1);
        if gridlevel == 1
            list = populate_impala_list( TI, tau);
        else
            for x = 1:gridlevel
                for y = 1:gridlevel
                    for z = 1:gridlevel
                        list = [list ; populate_impala_list(...
                            TI(x:gridlevel:end,y:gridlevel:end,...
                            z:gridlevel:end),tau)];
                    end
                end
            end
        end
end

%Find unique patterns
[d,~,Id] = unique(cell2mat(list(:,1)),'rows');
list_length = length(d);

% Create count matrix
C = cell2mat(list(:,2));

%Preallocate final count matrix
c = zeros(list_length, num_cat);

%For each unique pattern
for i = 1:list_length
    %Sum counts for each facies
    c(i,:) = sum(C(Id == i,:),1);
    %Had forgotten ",1) " in sum! :D
end

list = mat2cell([d c],ones(1,list_length),[template_length num_cat]);