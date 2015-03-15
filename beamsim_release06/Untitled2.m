fname=@(x)sprintf('probe_pitchsweep_%2.3f_',x*1e3);
fname2=@(x)sprintf('%s\\%sbeam_xy.png',fname(x),fname(x));
avi = avifile( 'animation.avi');
avi.Fps=5;
for idx_sweep=1:length(sweep_pitch)
    f1=fname2(sweep_pitch(idx_sweep));
    frame=imread(f1);
    avi=addframe(avi,frame);
end
avi=close(avi);
    