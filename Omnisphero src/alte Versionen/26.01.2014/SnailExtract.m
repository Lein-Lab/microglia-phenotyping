%fprintf (1,'Hello World');
%Write CSV file into an array.
datafile = 'Example.csv'
delimiter = ';'
lineCounter = 0;
if exist(datafile) > 0
    filePointer = fopen(datafile,'r');
    while (~feof(filePointer))
        lineCounter = lineCounter+1;
        line = fgetline(filePointer);
        cells = regexp(line,delimiter,'split');
        for i=1:length(cells)
            csv(i,lineCounter) = cells(i);
        end            
    end
else
    warning(['Could not find file ',datafile]);
    return;
end