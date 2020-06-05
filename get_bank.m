function phi_target = get_bank(phi_max, phi_o, V_1, phi_1, V)

phi_target = phi_max - (phi_max-phi_o)*exp(log((phi_max-phi_1)/(phi_max-phi_o))/V_1*V);