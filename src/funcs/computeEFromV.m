function EField = computeEFromV(V, x, y, z)
    % computeEFieldFromPotential Computes the electric field from a symbolic potential
    %
    % Inputs:
    %   V - Symbolic potential as a function of x, y, and z
    %   x, y, z - Symbolic variables representing spatial coordinates
    %
    % Outputs:
    %   Ex, Ey, Ez - Symbolic expressions for the electric field components along x, y, and z

    % Compute partial derivatives of the potential (gradient components)
    Ex = -diff(V, x); % Electric field component along x
    Ey = -diff(V, y); % Electric field component along y
    Ez = -diff(V, z); % Electric field component along z
    EField = {Ex, Ey, Ez};
end