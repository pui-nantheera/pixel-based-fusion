%function weight = Find_Weighting(Im, overlay, RegMap)
%% Program to allow for the choosing of the weighting of regions



temp = cell(1,1);
temp{1} = '1.0';
c = 0;
%% For joint
%if size(overlay,3) ==1
    for i = 1:size(Im,3)
%        overlay(:,:,i) = max((overlay(:,:,1) == 255)*255, Im(:,:,i));
        RegMap(:,:,i) = RegMap(:,:,1);
    end;
%end;

uiwait(msgbox('In order to weight a region, click on a region and then press Return.  Pressing Backspace or Delete removes the previously selected point.','Adding Weightings','help','modal'));

for i = 1:size(overlay,3)
    hfig = figure;
    imshow(overlay(:,:,i),[0 255]);
    button = 'Yes';
    while button == 'Yes';
        [x,y] = getpts;
        val = inputdlg('Enter the weight (suggested range 0.0 - 2.0)', 'Enter Weight');
        valnum = str2double(val);
        while isnan(valnum)
            uiwait(errordlg('Value entered is not a number', 'Input Error!'));
            val = inputdlg('Enter the weight (suggested range 0.0 - 2.0)', 'Enter Weight');
            valnum = str2double(val);
        end;
        c = c+1;
        weight(c,:) = [i RegMap(round(y), round(x), i) valnum];
        button = questdlg('Do you want to add more weightings for this image?','Continue Operation','Yes','No','Yes');
        if length(button) ~= 3
            button(3) = ' ';
        end;
    end;
    close(hfig);
end;