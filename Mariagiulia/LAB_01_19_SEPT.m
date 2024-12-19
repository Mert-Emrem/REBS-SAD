% 19 sept 2024

% INTRODUCTION 

% strings and char
h = 'ciao' ; % s i mp l e s t r i n g
i = char('pippo') ; % a n o t h e r s t r i n g o f c h a r a c t e r s
l = strcat([h,'',i]) ; % h o r i z o n t a l c o n c a t e n a t i o n o f s t r i n g s

% dynamic allocation for matrix is very inefficient

% es 1
f = @(x) 5.*cos(1/2 .* x) + 10.*log10(3.*x);
v1 = [0:0.001:1];
v2 = [0:0.001:100];
figure
plot(v1, f(v1));
figure
plot(v2, f(v2));

% es 3
eps = -pi/4;
b = [cos(eps), sin(eps)];
ang_rad = atan2(b(2),b(1)); % è intelligente e così ti da l'angolo nel quadrante gisuto mentre atan non sa darti gli angoli in II e III quadrante
ang_deg = rad2deg(ang_rad);
ang_arcsec = ang_deg*3600;


