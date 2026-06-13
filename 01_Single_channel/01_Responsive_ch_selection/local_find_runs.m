function runs = local_find_runs(mask)
    mask = mask(:)';
    dmask = diff([false, mask, false]);
    starts = find(dmask == 1);
    ends   = find(dmask == -1) - 1;
    runs = [starts(:), ends(:)];
end