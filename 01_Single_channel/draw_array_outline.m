function draw_array_outline(valid)

[nr,nc] = size(valid);

for i = 1:nr
    for j = 1:nc

        if ~valid(i,j)
            continue
        end

        % sopra
        if i == 1 || ~valid(i-1,j)
            plot([j-0.5 j+0.5],[i-0.5 i-0.5],'k','LineWidth',1)
        end

        % sotto
        if i == nr || ~valid(i+1,j)
            plot([j-0.5 j+0.5],[i+0.5 i+0.5],'k','LineWidth',1)
        end

        % sinistra
        if j == 1 || ~valid(i,j-1)
            plot([j-0.5 j-0.5],[i-0.5 i+0.5],'k','LineWidth',1)
        end

        % destra
        if j == nc || ~valid(i,j+1)
            plot([j+0.5 j+0.5],[i-0.5 i+0.5],'k','LineWidth',1)
        end
    end
end

end