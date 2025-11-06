function [F_drag, M_drag] = atmosferic_drag(r_cp, v_rel, S, C_d, p)
    %Atmosferic Drag
    %F_D​=−ρ * C_d * A * v_rel^2​ * v^rel
    %    -------------------------------
    %                  2
    %𝜌: densidad atmosférica [kg/m³],
    %C_d: coeficiente de arrastre (≈ 2.2 para satélites pequeños),
    %A: área efectiva expuesta al flujo [m²],
    %v_rel: velocidad del satélite respecto a la atmósfera,
    %v^rel: dirección unitaria de la velocidad relativa.
    
    v_norm = v_rel/norm(v_rel)
    ​F_drag = p*C_d*S*(v_rel^2)*v_norm
    M_drag = cross(r_cp, F_drag)
end