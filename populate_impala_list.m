function [ list ] = populate_impala_list( TI, tau )
%POPULATE_IMPALA_LIST Create IMPALA style pattern library for MPS
%
%   Inputs:
%       TI: Training image 2D or 3D
%       tau: data template
%
%   Ouputs:
%       list: cell array with d and c vectors
%           d: data event
%           c: counts
%
%   Óli D. Jóhannsson 2017


template_length = size(tau,1);
cat = unique(TI(:))';
num_cat = length(cat);
list_counter = 0;

% 2D or 3D
switch length(size(TI))
    
    case 2 % 2D
        % only take x and y components of tau
        tau = tau(:,1:2);
        
        [num_x,num_y] = size(TI);
        min_x = 1 + abs(min(tau(:,1)));
        max_x = num_x - max(tau(:,1));
        min_y = 1 + abs(min(tau(:,2)));
        max_y = num_y - max(tau(:,2));
        
        list_length = (1+max_x-min_x)*(1+max_y-min_y);
        list = cell(list_length,2);
        
        for i = min_x : max_x
            for j = min_y : max_y
                list_counter = list_counter + 1;
                
                d = zeros(1,template_length);
                
                for h = 1:template_length
                    d(h) = TI(i+tau(h,1),j+tau(h,2));
                    
                end
                
                c = 1.*(cat == TI(i,j));
                
                list{list_counter,1} = d;
                list{list_counter,2} = c;
            end
        end
        
        %Find unique patterns
        [d,~,Id] = unique(cell2mat(list(:,1)),'rows');
        final_list_length = size(d,1);
        
        % Create count matrix
        C = cell2mat(list(:,2));
        
        %Preallocate final count matrix
        c = zeros(final_list_length, num_cat);
        
        %For each unique pattern
        for i = 1:final_list_length
            %Sum counts for each facies
            c(i,:) = sum(C(Id == i,:),1);
            %Had forgotten ",1) " in sum! :D
        end
        
    case 3 % 3D
        tau = tau(:,1:3);
        [num_x,num_y,num_z] = size(TI);
        min_x = 1 + abs(min(tau(:,1)));
        max_x = num_x - abs(max(tau(:,1)));
        min_y = 1 + abs(min(tau(:,2)));
        max_y = num_y - abs(max(tau(:,2)));
        min_z = 1 + abs(min(tau(:,3)));
        max_z = num_z - abs(max(tau(:,3)));
        
        for i = min_x : max_x
            for j = min_y : max_y
                for k = min_z : max_z
                    %list_counter = list_counter + 1;
                    %list_counter = 1 +...
                    %    (i-min_x)*(1+max_y-min_y)*(1+max_z-min_z) +...
                    %    (j-min_y)*(1+max_z-min_z) + (k-min_z);
                    
                    d = zeros(1,template_length);
                    
                    for h = 1:template_length
                        d(h) = TI(i+tau(h,1),j+tau(h,2),k+tau(h,3));
                        
                    end
                    
                    c = 1.*(cat == TI(i,j,k));
                    
                    %list{list_counter,1} = d;
                    %list{list_counter,2} = c;
                    D(i,j,k,:) = d;
                    C(i,j,k,:) = c;
                    
                end
            end
        end
        
        D = reshape(D,[max_x*max_y*max_z,size(tau,1)]);
        C = reshape(C,[max_x*max_y*max_z,length(cat)]);
        
        tic
        %Find unique patterns
        %[d,~,Id] = unique(D,'rows','sorted');
        %[d,ID,Id] = unique(D,'rows','sorted');
        [d,~,Id] = unique(D,'rows');
        final_list_length = size(d,1);
        toc
        
        %Preallocate final count matrix
        c = zeros(final_list_length, num_cat);
        
        %For each unique pattern
        % Note: Parfor is >2x faster on a 4 core for a large list.
        %       There maybe a faster way to do this though.
        tic
        for i = 1:final_list_length
            %Sum counts for each facies
            c(i,:) = sum(C(Id == i,:),1);
        end
        toc
        %c
%         tic
%         c3 = zeros(final_list_length, num_cat);
%         for i = 1:final_list_length
%             for j = 1:length(C)
%                 if Id(j) == i
%                     c3(i,:) = c3(i,:) + C(j,:);
%                 end
%             end
%         end
%         toc
%         %c3
%         isequal(c,c3) 

end
%whos d
%whos c
%num_cat
%final_list_length
list = mat2cell([d c],ones(1,final_list_length),[template_length num_cat]);


