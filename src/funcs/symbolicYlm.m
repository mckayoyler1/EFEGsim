function Ylm = symbolicYlm(l, m)
    % Function to calculate the symbolic spherical harmonic Y_lm(theta, phi)
    % Inputs:
    %   l     - Degree (non-negative integer)
    %   m     - Order (integer satisfying -l <= m <= l)
    %   theta - Symbolic or numeric polar angle (0 <= theta <= pi)
    %   phi   - Symbolic or numeric azimuthal angle (0 <= phi < 2*pi)
    % Output:
    %   Ylm   - Symbolic spherical harmonic Y_lm(theta, phi)
    
    % Ensure symbolic variables for theta and phi
    syms theta phi real
    syms x real

    % Define the Legendre polynomial P_l(x) using Rodrigues' formula
    P_l = (1/(2^l * factorial(l))) * diff((x^2 - 1)^l, x, l);
    
    % Define the Associated Legendre polynomial P_l^m(x)
    P_lm = (-1)^m * (1 - x^2)^(m/2) * diff(P_l, x, m);
    
    % Substitute x = cos(theta) to express P_l^m in terms of theta
    P_lm_cos_theta = subs(P_lm, x, cos(theta));
    
    % Define the normalization constant for the spherical harmonic
    N_lm = sqrt((2*l + 1)/(4*sym(pi)) * factorial(l - m)/factorial(l + m));
    % Define the spherical harmonic Y_lm(theta, phi)
    Ylm = N_lm * P_lm_cos_theta * exp(1i * m * phi);
end