function M = rigidMatrixFromParameters(x)
% RIGIDMATRIXFROMPARAMETERS Computes the rigid matrix from the input rigid paramrters
% Data from the maxima script
% 
% Author: Alberto Gomez, Biomedical Engineering, KCL, 2013

    tx=x(1);
    ty=x(2);
    tz=x(3);
    a=x(4);
    b=x(5);
    c=x(6);
    
    
    M = [cos(b)*cos(c)-sin(a)*sin(b)*sin(c) -cos(a)*sin(c) sin(a)*cos(b)*sin(c)+sin(b)*cos(c) tx
        cos(b)*sin(c)+sin(a)*sin(b)*cos(c) cos(a)*cos(c) sin(b)*sin(c)-sin(a)*cos(b)*cos(c) ty
        -cos(a)*sin(b) sin(a) cos(a)*cos(b) tz
        0 0 0 1];
        
end