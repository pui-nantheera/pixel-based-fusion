%pointclick Graphical input from mouse.
%ginput by Mathworks, but modified to pointclick by Patrick Baker 05/05/2000
%   This is a hacked version of ginput which allows you to zoom in using
%   the mouse.  The following actions are valid:
%   SPACE  selects the current point
%   'f'    gives NaN as the value for the current point
%   RETURN ends the point selection
%   '-'    zooms out completely
%
%   [X,Y] = pointclick(N) gets N points from the current axes and returns 
%   the X- and Y-coordinates in length N vectors X and Y.  The cursor
%   can be positioned using a mouse (or by using the Arrow Keys on some 
%   systems).  Data points are entered by pressing a mouse button
%   or any key on the keyboard except carriage return, which terminates
%   the input before N points are entered.
%
%   [X,Y] = pointclick gathers an unlimited number of points until the
%   return key is pressed.
% 
%   [X,Y,BUTTON] = pointclick(N) returns a third result, BUTTON, that 
%   contains a vector of integers specifying which mouse button was
%   used (1,2,3 from left) or ASCII numbers if a key on the keyboard
%   was used.
function [out1,out2,out3] = pointclick(arg1)

%   Copyright (c) 1984-96 by The MathWorks, Inc.
%   $Revision: 5.18 $  $Date: 1996/11/10 17:48:08 $

zoom on
out1 = []; out2 = []; out3 = []; y = [];
c = computer;
if ~strcmp(c(1:2),'PC') & ~strcmp(c(1:2),'MA')
    tp = get(0,'TerminalProtocol');
else
    tp = 'micro';
end

if ~strcmp(tp,'none') & ~strcmp(tp,'x') & ~strcmp(tp,'micro'),
    if nargout == 1,
        if nargin == 1,
            eval('out1 = trmginput(arg1);');
        else
            eval('out1 = trmginput;');
        end
    elseif nargout == 2 | nargout == 0,
        if nargin == 1,
            eval('[out1,out2] = trmginput(arg1);');
        else
            eval('[out1,out2] = trmginput;');
        end
        if  nargout == 0
            out1 = [ out1 out2 ];
        end
    elseif nargout == 3,
        if nargin == 1,
            eval('[out1,out2,out3] = trmginput(arg1);');
        else
            eval('[out1,out2,out3] = trmginput;');
        end
    end
else

    fig = gcf;
    figure(gcf);

    if nargin == 0
        how_many = -1;
        b = [];
    else
        how_many = arg1;
        b = [];
        if  isstr(how_many) ...
            | size(how_many,1) ~= 1 | size(how_many,2) ~= 1 ...
            | ~(fix(how_many) == how_many) ...
            | how_many < 0
            error('Requires a positive integer.')
        end
        if how_many == 0
            ptr_fig = 0;
            while(ptr_fig ~= fig)
                ptr_fig = get(0,'PointerWindow');
            end
            scrn_pt = get(0,'PointerLocation');
            loc = get(fig,'Position');
            pt = [scrn_pt(1) - loc(1), scrn_pt(2) - loc(2)];
            out1 = pt(1); y = pt(2);
        elseif how_many < 0
            error('Argument must be a positive integer.')
        end
    end

    pointer = get(gcf,'pointer');
    set(gcf,'pointer','crosshair');
    fig_units = get(fig,'units');
    char = 0;
    begin_how_many = how_many;
    while how_many ~= 0
        % Use no-side effect WAITFORBUTTONPRESS
        waserr = 0;
%        eval('keydown = wfbp;', 'waserr = 1;');
        keydown=waitforbuttonpress;
        if(waserr == 1)
      if(ishandle(fig))
        set(fig,'pointer',pointer,'units',fig_units);
        error('Interrupted');
      else
        error('Interrupted by figure deletion');
      end
        end
      
        ptr_fig = get(0,'CurrentFigure');
        if(ptr_fig == fig)
            if keydown
                char = get(fig, 'CurrentCharacter');
                button = abs(get(fig, 'CurrentCharacter'));
                scrn_pt = get(0, 'PointerLocation');
                set(fig,'units','pixels')
                loc = get(fig, 'Position');
                pt = [scrn_pt(1) - loc(1), scrn_pt(2) - loc(2)];
                set(fig,'CurrentPoint',pt);
                if (char==13)
                    break;
                end
                switch char
                    case '+',
                        zoom(2)
                    case {'-', 'a'},
                        zoom out
                    case ' ',	
                        pt = get(gca, 'CurrentPoint');
%                        disp([how_many pt(1,1:2)]);
                        how_many = 0;
                        r=gca;
                        out1 = [out1;gca,pt(1,1)];
                        y = [y;pt(1,2)];
                        b = [b;button];
                        zoom out
                    case 'f',
                        pt = [NaN NaN];
                        disp([how_many pt(1,1:2)]);
                        how_many = how_many - 1;
                        out1 = [out1;pt(1,1)];
                        y = [y;pt(1,2)];
                        b = [b;button];
                    case 'p',
                        how_many=begin_how_many;
                        out1=[];
                        y=[];
                        b=[];
                        disp('RESET!');
                end
            end
        end
    end

    set(fig,'pointer',pointer,'units',fig_units);

    if nargout > 1
        out2 = y;
        if nargout > 2
            out3 = b;
        end
    else
        out1 = [out1 y];
    end

end
idxes=find(sum(out1<0,2));
if length(idxes) > 0,
  out1(idxes,:)=NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = wfbp
%WFBP   Replacement for WAITFORBUTTONPRESS that has no side effects.

% Remove figure button functions
fprops = {'windowbuttonupfcn','buttondownfcn', ...
          'windowbuttondownfcn','windowbuttonmotionfcn'};
fig = gcf;
fvals = get(fig,fprops);
set(fig,fprops,{'','','',''})

% Remove all other buttondown functions
ax = findobj(fig,'type','axes');
if isempty(ax)
    ch = {};
else
    ch = get(ax,{'Children'});
end
for i=1:length(ch),
    ch{i} = ch{i}(:)';
end
h = [ax(:)',ch{:}];
vals = get(h,{'buttondownfcn'});
mt = repmat({''},size(vals));
set(h,{'buttondownfcn'},mt);

% Now wait for that buttonpress, and check for error conditions
waserr = 0;
eval(['if nargout==0,', ...
      '   waitforbuttonpress,', ...
      'else,', ...
      '   keydown = waitforbuttonpress;',...
      'end' ], 'waserr = 1;');

% Put everything back
if(ishandle(fig))
  set(fig,fprops,fvals)
  set(h,{'buttondownfcn'},vals)
end

if(waserr == 1)
  error('Interrupted');
end

if nargout>0, key = keydown; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
