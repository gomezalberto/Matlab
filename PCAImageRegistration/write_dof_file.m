function [ ] = write_dof_file( dof_filename, dof )
%write_dof_file Writes to file filename the dof in matrix dof
%   Input:  -filename   name of the file
%           -dof        matrix Nx3 with dof
% 
% Author: Devis Peressutti, Biomedical Engineering, KCL, 2011

% radiants to degrees
r2d = 180/pi;
% Convert radiants to degrees
dof(4:end) = r2d*dof(4:end);
mdof = zeros(numel(dof),3);
mdof(:,end) = dof;
% Open file
fid = fopen(dof_filename,'w');
% Check if can open file
if (fid~=-1)
    display(sprintf('Writing dof to %s',dof_filename));
    fprintf(fid,'DOF: %d\n',numel(dof));
    fprintf(fid,'%.1f\t%.1f\t%.2f\n',transpose(mdof));
    fclose(fid);
else
    error('Can''t write to %s',filename);
end

end

