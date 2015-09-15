function [ dof ] = read_dof_file( dof_filename, ndof )
%read_dof_file Read in rigid DOF file
%   This function reads in the ndof parameters stored in the
%   dof_filename file

% degrees to radiants
d2r = pi/180;
% Open file
fid = fopen(dof_filename, 'r');
% Format of file is compliant with IRTK (first line is DOF: ndofs)
C = textscan(fid, '%f', ndof*3, 'delimiter', ',', 'HeaderLines', 1)';
c = C{1};
% Close file
fclose(fid);
% Reshape matrix - first two columns are filled with zeros
tmp = reshape(c, 3,ndof)';
% Retain only last column
dof = tmp(:,3);
% Convert degrees to radiants
dof(4:end) = d2r*dof(4:end);

end

