function [phi]=MC(th), phi=max(0,min([(1+th)/2 2 2*th]));   end   % Flux limiters