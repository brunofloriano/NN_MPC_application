if strcmp(save_mainname,'dataBalloon ')

    Ecc = -beta_cc * ( L{s}(1,1) + gamma_cc*( sum(diag(L{s})) -  L{s}(1,1)) );
    Eac = beta_ac * (1 - virtual_area_covered);%sum(sum(virtual_encompassed_area))/sum(sum(total_circle_area)));%virtual_area_covered/total_area);
    Et = 0;%-beta_t * (exp(alpha_t*t(virtual_time_counter)) - 1);
    
    E = [Ecc; Eac; Et];

    Jinst(virtual_time_counter,s) = norm(E);
end