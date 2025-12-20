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
    v_hat = v_rel/norm(v_rel);
    v_norm = norm(v_rel);
    F_drag = -v_hat*(v_norm^2)*p*C_d*S/2;
    M_drag = cross(r_cp, F_drag);
end