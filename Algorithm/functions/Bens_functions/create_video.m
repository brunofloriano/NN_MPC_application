function create_video(F,title)
    
    
    writerObj = VideoWriter(title, 'MPEG-4');
    writerObj.FrameRate = 10;
    open(writerObj);
    for i=1:length(F)
        frame = F(i);
        writeVideo(writerObj, frame);
    end

    close(writerObj);
end