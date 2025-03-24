% Helper function to find which shell contains a given radius
function idx = find_shell_index(r, shells)
    for i = 1:length(shells)-1
        if r >= shells(i).inner_radius && r <= shells(i).outer_radius
            idx = i;
            return;
        end
    end
    % If we get here, must be in outermost shell
    idx = length(shells);
end