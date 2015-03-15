%myaa My anti-aliased printout procedure.
function myaa(simulation,fname)
if simulation.verbose
    txt=sprintf('storing %s',fname);
    fprintf('%s',txt);
end
tempfile='tempfile.png';
set(gcf, 'PaperUnits','centimeters');
%width=15; % centimeters
width=simulation.papersize;
aspect=4/3; % aspect ratio
set(gcf, 'PaperSize', [width width/aspect]);
myfiguresize = [0, 0, width, width/aspect];
set(gcf, 'PaperPosition', myfiguresize)
print(gcf,'-r1200','-dpng', tempfile);
raw_hires=imread(tempfile);
%size(raw_hires)
raw_lowres = single(imresize(raw_hires,1/4,'bicubic'))/256;
imwrite(raw_lowres,fname,'png');
try
    deletefile(tempfile);
catch E
    % silent fallover
end

if simulation.verbose
    for ii=1:length(txt);
        fprintf('\b');
    end
end

