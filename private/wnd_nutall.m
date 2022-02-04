function f = wnd_nutall(L)
%WND_NUTALL Nutall window function for anti-reflection boundary
    f = cumsum(nuttallwin(L+1));
    f = f(1:end-1)/f(end);
end

