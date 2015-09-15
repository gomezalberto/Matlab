function patches = findPatches(number_of_patches, region, border, varargin)
% findPatches Function that extracts a given number of patches from an image region
%
% This version should match the c++ version
% region is an ImageType object
%
% Author: Alberto Gomez, Biomedical Engineering, KCL, 2011

patches = [];
global npatches;
npatches=1;

blank_border =0;

for i=1:size(varargin,2)
    if (strcmp(varargin{i},'blanks'))
        blank_border = varargin{i+1};
        i=i+1;
    end
end


InputDimension = numel(region.spacing);

%%%%%%%%%%%%%%%%%%%%%%%%% Find first patches %%%%%%%%%%%%%%%%%%%
maxPatchSize = ceil( region.size(1:InputDimension)./ number_of_patches);

if isa(region,'VectorImageType')
    region_patch = VectorPatchType(region);
    region_patch.datax = region.data;
    region_patch.datax = region.datax;
    region_patch.datax = region.datay;
    region_patch.datax = region.dataz;
elseif isa(region,'ImageType') % this will also yield true for VectortImageType
    region_patch = PatchType(region);
    region_patch.data = region.data;
else
    region_patch = region;
end
    

patches = patchRegion(maxPatchSize, region_patch, patches, border);

if blank_border>0
    % remove 1 node at the end only if it is not an end patch
    for i=1:numel(patches)
        end_borders = region.borderIndex2-patches(i).borderIndex2;
        removable_dimensions = find(end_borders);
        for j=removable_dimensions
            new_size = patches(i).size;
            new_size(j) = new_size(j)-blank_border;
        end
        
        patch_ = region.extractPatch(patches(i).index, new_size, border );
        patches(i)=patch_;
    end
end

end

function patches = patchRegion(maxPatchSize,region, patches, border)

    InputDimension = numel(region.spacing);


    npatches = floor(region.size(1:InputDimension) ./ maxPatchSize); % Patches of size (maxPatchSize)
    remaining_nodes=region.size(1:InputDimension)-npatches.*maxPatchSize;
    total_npatches = npatches + double(remaining_nodes>0); % Patches of size (maxPatchSize) and smaller

    total_new_patches = prod(total_npatches);
            
    for i=1:total_new_patches

        str = [];
        str2 =[];
        for ii=1:numel(region.spacing)
            str = [str 'i' num2str(ii) ','];
            str2 = [str2 'i' num2str(ii) ' '];
        end
        eval(['['  str(1:end-1) '] = ind2sub(total_npatches'',i); index = [' str2(1:end-1) ']'';' ]);
        
        patch_index= (index-1).*maxPatchSize+region.index;

        current_size = zeros(numel(region.spacing),1);
        for ii=1:numel(region.spacing)
            in_patches = find( index- npatches<=0);
            out_patches = find( index- npatches>0);
            current_size(in_patches) = maxPatchSize(in_patches);
            current_size(out_patches) = remaining_nodes(out_patches);
        end

        current_patch = region.extractPatch(patch_index, current_size, border);
        patches = [patches ; current_patch];
                                    
    end % End further subdivision
end
                                    
                                    
                                    
                                    
                                    
             