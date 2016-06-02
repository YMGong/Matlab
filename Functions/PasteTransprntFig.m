function PasteTransprntFig(filespec)
    %% Draw Plot
    figure
    x = linspace(0,2*pi);
    y = rand(1,100);
    polar(x,y)

    %% Adjust Figure Settings
    %All that's needed is gcf Color settings and editmenufcn() operation.
    set(gcf, 'PaperPosition', [0.25 0.25 8.0 6.0])
    set(gcf, 'Position' , [1 1 459 362])
    set(gcf, 'Color' , 'None')

    set(gca, 'Position' , [0 0 1 1])

    editmenufcn(gcf, 'EditCopyFigure')

    %% PowerPoint Commands (pulled from saveppt.m)
    ppt = actxserver('PowerPoint.Application'); %Start an ActiveX session with PowerPoint.

    if ~exist(filespec, 'file')
        op = invoke(ppt.Presentations, 'Add'); %Create new presentation.
    else
        op = invoke(ppt.Presentations, 'Open', filespec, [], [], 0); %Open existing presentation.
    end

    slide_count = get(op.Slides, 'Count'); %Get current number of slides.
    slide_count = int32(double(slide_count) + 1);
    new_slide = invoke(op.Slides, 'Add', slide_count, 11); %Add a new slide (with title object).

    pic1 = invoke(new_slide.Shapes, 'Paste'); %Paste the contents of the Clipboard.

    %% Exiting Commands
    if ~exist(filespec,'file')
        invoke(op, 'SaveAs', filespec,1); %Save file as new.
    else
        invoke(op, 'Save'); %Save existing file.
    end
    invoke(op, 'Close'); %Close the presentation window.
    invoke(ppt, 'Quit'); %Quit PowerPoint.
end