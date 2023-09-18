function savefolderlocal = save_figure(gcf,name,save_folder,save_mainname,time_now)

    savefolderlocal = [save_folder '/' save_mainname time_now];
    saveas(gcf,[savefolderlocal '/' save_mainname  time_now '-' name '.eps']);
    saveas(gcf,[savefolderlocal '/' save_mainname  time_now '-' name '.fig']);