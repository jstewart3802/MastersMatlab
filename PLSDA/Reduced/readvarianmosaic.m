function [wavenumbers, data, width, height, filename, acqdate]=readvarianmosaic(filename, keepme)

%edited by James Stewart Nov 2020 to remove some things not needed by our
%program and to make it easier to read and understand



%   Function: readvarianmosaic
%   Usage: [wavenumbers, data, width, height, filename, acqdate] = readvarianmosaic_v4_1(filename, keepme);
%   Usage: [wavenumbers, data, width, height, filename, acqdate] = readvarianmosaic_v4_1();
%           (second version prompts for a filename)
%
%   Extracts the spectra from a Varian .dmt/.dms/.dmd file combination.
%   Plots an image of the total signal.
%
%   input:
%   'filename' string containing the full path to the .dms file (optional)
%   'keepme' vector of wavenumber values in pairs indicating the limits of regions to retain (optional)
% 
%   output:
%   'wavenumbers' is a list of the wavenumbers related to the data
%   'data' is a 3D cube of the data in the file (128xX x 128xY x wavenumbers)
%   'width' is width in pixels of the entire mosaic
%   'height' is height in pixels of the entire mosaic
%   'filename' is a string containing the full path to the .dms file
%   'acqdate' is a string containing the date and time of acquisition
%
%                     *******Caution******* 
%   This code is a hack of the Varian format and the location of the data
%   within the file may vary. Always check the output to make sure it is
%   sensible. If you have a file that doesn't work, please contact Alex. 
%
%   Copyright (c) 2011 - 2015, Alex Henderson 
%   Contact email: alex.henderson@manchester.ac.uk
%   Licenced under the GNU General Public License (GPL) version 3
%   http://www.gnu.org/copyleft/gpl.html
%   Other licensing options are available, please contact Alex for details
%   If you use this file in your work, please acknowledge the author(s) in
%   your publications. 
%
%       version 4.1 March 2015


%       version 4.1 March 2015 Alex Henderson, Small change to fix
%       filename case sensitivity issues on Linux. No other functional
%       changes. 
%       version 4.0 January 2014 Alex Henderson, Moved the data allocation
%       outside the file read loop. Also moved the wavenumber truncation
%       for keepme scenarios outside the file read loop.
%       version 3.0 November 2013 Alex Henderson, Added 'keepme' to allow
%       the data to have spectral regions removed during import
%       version 2.0 August 2012 Alex Henderson, replaced loop through
%       matrix with permute and flipdim
%       version 1.3 December 2011 Alex Henderson, added GPL licence and
%       incorporated the getfilename function
%       version 1.2 November 2011 Alex Henderson, the dmt filename is all
%       lowercase while the other filenames are of mixed case. Now we ask
%       for the .dms filename instead and build a lowercase .dmt filename
%       from that. This was only a problem in Linux. 
%       version 1.1 October 2011 Alex Henderson, the dms file only matches
%       the image if the number of tiles is small. Now we read each tile
%       separately. 
%       version 1.0 October 2011 Alex Henderson, initial release, based on
%       readvarian v2

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% get the filename if not supplied

if (exist('filename', 'var') == 0)
    filename=getfilename('*.dms', 'Varian Mosaic Files (*.dms)');

    if (isfloat(filename) && (filename==0))
        return;
    end
    
    filename=filename{1};
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% determine the mosaic dimensions
[pathstr, name, ext] = fileparts(filename); 
tiles_in_x_dir = xtiles(fullfile(pathstr,name));
tiles_in_y_dir = ytiles(fullfile(pathstr,name));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% extract the wavenumbers and date from the dmt file

[pathstr, name, ext] = fileparts(filename); 

% The Varian software stores all files with mixed case filenames except the
% dmt file which is all lowercase. Therefore we build this from the dms
% filename. 
name = lower(name);

dmtfilename = fullfile(pathstr,[name '.dmt']);

[fid, message] = fopen(dmtfilename, 'r', 'a');
if(fid == -1) 
    disp(['reading dmt file: ', dmtfilename]);
    error(message); 
end;

% wavenumbers
    status = fseek(fid, 2228, 'bof');
    if(status == -1), message = ferror(fid, 'clear'); error(message); end;
    startwavenumber = double(fread(fid, 1, 'int32'));
    
    status = fseek(fid, 2236, 'bof');
    if(status == -1), message = ferror(fid, 'clear'); error(message); end;
    numberofpoints = double(fread(fid, 1, 'int32'));
    
    status = fseek(fid, 2216, 'bof');
    if(status == -1), message = ferror(fid, 'clear'); error(message); end;
    wavenumberstep = fread(fid, 1, 'double');
    
    % some validation
        if(startwavenumber < 0)
            error('Start wavenumber is negative. Cannot read this file.');
        end
        
%         if(numberofpoints ~= vals)
%             error('Number of data points does not match between .dms and .dmt files. Cannot read this file.');
%         end            
    
    wavenumbers = 1:(numberofpoints+startwavenumber-1);
    wavenumbers = wavenumbers * wavenumberstep;
    wavenumbers = wavenumbers(startwavenumber:end);

% date
    % Longest date is: Wednesday, September 30, 2011 00:00:00
    status = fseek(fid, 0, 'bof');
    if(status == -1), message = ferror(fid, 'clear'); error(message); end;

    str = fread(fid, inf, '*char')';
    expr = 'Time Stamp.{44}\w+, (\w+) (\d\d), (\d\d\d\d) (\d\d):(\d\d):(\d\d)';

    [start_idx, end_idx, extents, matches, tokens, names, splits] = regexp(str, expr);

    monthword=tokens{1,1}{1};
    day=tokens{1,1}{2};
    year=tokens{1,1}{3};
    hours=tokens{1,1}{4};
    minutes=tokens{1,1}{5};
    seconds=tokens{1,1}{6};
    
    switch monthword
       case 'January'
          month='01';
       case 'February'
          month='02';
       case 'March'
          month='03';
       case 'April'
          month='04';
       case 'May'
          month='05';
       case 'June'
          month='06';
       case 'July'
          month='07';
       case 'August'
          month='08';
       case 'September'
          month='09';
       case 'October'
          month='10';
       case 'November'
          month='11';
       case 'December'
          month='12';
       otherwise
          month='99';
    end;

    acqdate = [day, ' ', monthword, ' ', year, ', ', hours, ':', minutes,':',seconds];
    
fclose(fid);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Handle sections we want to keep, if specified

if (exist('keepme', 'var'))

   if(isempty(keepme))
       error('keepme limits list is empty');
   end
   
   keepme=columnvector(keepme);

   if(rem(length(keepme),2) ~= 0)
       error('keepme limits list is not an even number');
   end

   keepme=sort(keepme);           
   keepme=reshape(keepme,2,[]);
   keepme=keepme';
   numberofsections=size(keepme,1);
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Allocate memory for output
% 
limits=[];
if (exist('keepme', 'var'))
    keptwavenumbers=[];
    for section=1:numberofsections
        limits=indicesfromvalues(keepme(section,:), wavenumbers);
        keptwavenumbers=cat(2,keptwavenumbers,wavenumbers(limits(1):limits(2)));
    end
    wavenumbers=keptwavenumbers;
end

data=zeros(128*tiles_in_y_dir, 128*tiles_in_x_dir, length(wavenumbers));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% read each .dmd file
[pathstr, name, ext] = fileparts(filename); 

for y = 1:tiles_in_y_dir
    for x = 1:tiles_in_x_dir

        current_extn = sprintf('_%04d_%04d.dmd', x-1, y-1);
        tempfilename = fullfile(pathstr,[name, current_extn]);
        [fid, message] = fopen(tempfilename, 'r', 'a');
        if(fid == -1) 
            disp(['reading dmd file: ', tempfilename]);
            error(message); 
        end;

        tempdata = fread(fid, inf, '*float32');
        fclose(fid);
        
        tempdata=tempdata(256:end);
        tempdata=reshape(tempdata,128,128,[]);

        % remove spectral ranges during file read stage
        if (exist('keepme', 'var'))
            keptdata=[];
            for section=1:numberofsections
                keptdata=cat(3,keptdata,tempdata(:,:,limits(1):limits(2)));
            end
            tempdata=keptdata;
        end
        
        % rotate the image to match the spectrometer's output
        tempdata=permute(tempdata,[2,1,3]);
        tempdata=flipdim(tempdata,1);
        
        % insert this tile into the image
        data((1+((y-1)*128)) : (y*128), (1+((x-1)*128)) : (x*128), :) = tempdata;
    end
end
clear tempdata;

data=double(data);
[height, width, vals]=size(data);

end % of main function
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function tiles_in_x_dir = xtiles(basefilename)

% count tiles in x dimension

tiles_in_x_dir = 1;
[pathstr, name, ext] = fileparts(basefilename); 
finished = 0;
counter = 0;
while (~finished)
    current_extn = sprintf('_%04d_0000.dmd', counter);
    tempfilename = [basefilename, current_extn];
    [fid] = fopen(tempfilename, 'r', 'a');
    if(fid == -1) 
        tiles_in_x_dir = counter;
        finished = 1;
    else
        fclose(fid);
    end
    counter = counter + 1;
end

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function tiles_in_y_dir = ytiles(basefilename)

% count tiles in y dimension

tiles_in_y_dir = 1;
[pathstr, name, ext] = fileparts(basefilename); 
finished = 0;
counter = 0;
while (~finished)
    current_extn = sprintf('_0000_%04d.dmd', counter);
    tempfilename = [basefilename, current_extn];
    [fid] = fopen(tempfilename, 'r', 'a');
    if(fid == -1) 
        tiles_in_y_dir = counter;
        finished = 1;
    else
        fclose(fid);
    end
    counter = counter + 1;
end

end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function filename=getfilename(filter, filtername)

%   Function: getfilename
%   Usage: [filename] = getfilename(filter, filtername);
%
%   Collects a single filename from the user.
%   For multiple filenames use getfilenames.m
%
%   'filter' and 'filtername' are strings of the form...
%       filter = '*.mat'
%       filtername = 'MAT Files (*.mat)'
%   'filename' is a char array containing the name of the file including
%   the path
%
%   (c) May 2011, Alex Henderson
%

% Mostly based on getfilenames and tweaked to only accept a single filename

filetypes = {   filter, filtername; ...
                '*.*',    'All Files (*.*)'};

% example...            
%filetypes = {   '*.mat',  'MAT Files (*.mat)'; ...
%                '*.*',    'All Files (*.*)'};

setappdata(0,'UseNativeSystemDialogs',false);

[filenames, pathname] = uigetfile(filetypes, 'Select file...', 'MultiSelect', 'off');

if (isfloat(filenames) && (filenames==0))
    disp('Error: No filename selected');
    filename=0;
    return;
end

if(iscell(filenames))
    % change from a row of filenames to a column of filenames
    % if only one file is selected we have a single string (not a cell
    % array)
    filenames = filenames';
else
    % convert the filename to a cell array (with one entry)
    filenames=cellstr(filenames);
end

for i=1:size(filenames,1)
    filenames{i,1}=[pathname,filenames{i,1}];
end
filename=filenames;

end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
