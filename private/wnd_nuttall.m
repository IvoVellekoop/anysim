function f = wnd_nuttall(L)
%WND_NUTTALL Nutall window function for anti-reflection boundary

% This function computes the integral of the ordinary Nuttall window
%   Nuttall, Albert. "Some windows with very good sidelobe behavior."
%   IEEE Transactions on Acoustics, Speech, and Signal Processing 29.1
%   (1981): 84-91.
%

    % original formulation:
    %f = cumsum(nuttallwin(L+1));
    %f = f(1:end-1)/f(end);

    % exact formulation (normalized coefficients from page 89 of ref) 
    x = (0:L-1)'/(L-1);
    a2 = [-0.4891775 0.1365995/2 -0.0106411/3].' / (0.3635819 * 2 * pi);
    f = sin(x* [1 2 3] *2*pi) * a2 + x;
end

