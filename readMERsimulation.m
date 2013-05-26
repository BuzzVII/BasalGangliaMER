function [header,DataMER,ppheader,DataPP]=readMERsimulation(filename)
%function [header,DataMER,ppheader,DataPP]=readMERsimulation(filename)
%This function takes a simulation output file that has used
%writeMERsimulation.m and extracts the MER data, the PP data and then
%creates a structure file, header, that contains the simulation parameters
    
    fid=fopen(filename);
    if (fid < 0)
        error('Could not open %s', file_name);
    end
    
    %If error occures while reading file catch and close the file
    try
        %read in the full file as characters to find string markers
        fullfile=fread(fid,inf,'char').';
    
        %Set string markers that indicate the start of each data set
        PPstartstring='%START_PP_DATA%';
        MERstartstring='%START_DATA%';
    
        %initialize finding inications variables to not found
        fileindex=0;
        PPfound=0;
        MERfound=0;
    
        %repeat until the PP and MER data start position is found, or the eof
        %is reached.
        while (~PPfound || ~MERfound) && fileindex<length(fullfile)
            fileindex=fileindex+1;
            if ~MERfound && (fileindex+12)<length(fullfile)%check if MER start string is present and location if found
                if strcmp(char(fullfile(fileindex:fileindex+11)),MERstartstring)
                    MERfound=1;
                    MERstart=(fileindex+11);
                end
            end
            if ~PPfound && (fileindex+15)<length(fullfile)%check if PP start string is present and location if found
                if strcmp(char(fullfile(fileindex:fileindex+14)),PPstartstring)
                    PPfound=1;
                    PPstart=(fileindex+14);
                end
            end
        end
    
        %Check that both MER and PP position are found, if not warn user and
        %avoid reading data from file.
        if ~PPfound && ~MERfound
            warning('One or more datasets missing in file');
            header=[];
            ppheader=[];
            DataMER=[];
            DataPP=[];
        else
        
            %Extract the header from the data then process header
            headertext=char(fullfile(1:PPstart));
            [header,ppheader]=ParseHeader(headertext);
        
            %calculate the length of the pp data
            ppsize=(MERstart-(PPstart+13))*8/64;
        
            %find the location of the PP data and read it out
            ppread=fseek(fid,PPstart,'bof');
            DataPP=fread(fid,ppsize,'int64');
    
            %find the location of the MER data and read it out
            merread=fseek(fid,MERstart,'bof');
           DataMER=fread(fid,inf,'float64');
        
            %Make sure that both MER and PP data are sucesfully read, if not
            %return warning.
           if ppread || merread 
                warning('Data not succesfully read');
           end
        end
        fclose(fid);            % close the file
        
    catch exception
        fclose(fid);            % close the file
        rethrow(exception);     % rethrow the error
    end
end
function [header,ppheader]=ParseHeader(headertext)
%function header=GetHeader(headertext)
%This function takes in the header from a simulation file in the format of
%a string. This string is then broken down into the individual paramater 
%values, which are stored in the structure file header and returned. 
      header = struct('date', [], ...
                'patient' , 'Sim', ...         % simulation patient
                'fs', [], ...                  % in Hz
                'ts', [], ...                  % recording time is s
                'samples', [], ...             % the number of recorded samples
                'num_neurons', [], ...
                'neuron_model', [], ...
                'noise_type', [], ...
                'noise_params', [], ...
                'low_filter_type', [], ...
                'low_filter_params', [], ...   % e.g. order then cutoff
                'high_filter_type', [], ...
                'high_filter_params', [], ...  % e.g. order then cutoff
                'firing_distribution', [], ...
                'firing_params', [], ...
                'data_format', [] ...
                );
      ppheader = struct('fs', [], ...          % in Hz
                'ts', [], ...                  % recording time is s
                'samples', [], ...             % the number of recorded samples
                'num_neurons', [],...
                'firing_distribution', [], ...
                'firing_params', [], ...
                'data_format', [] ...
                );            
            
      %set initial conditions for continuing to read
      pphead=0;
      datastart=0;
      
      %seperate strings found in headertext sepearated by : or ;
      parseheader=textscan(headertext,'%s %s','delimiter',[';',':']);
      tokens=parseheader{1}; %assign the strings seperated by ; into field tokens
      values=parseheader{2}; %and strings seperated : into field values

      %Sort all the values into the appropriate header token
      for tokenindex=1:length(tokens)     
           token=tokens{tokenindex}; %assign current token
           
           %check is token indicates start of PP header or Data
           if strcmp(token,'%start_PP_Header%')
                    pphead=1;
                    fprintf('finished reading MER header\n')
           end
           if strcmp(token,'%START_PP_DATA%')
                    datastart=1;
                    fprintf('finished reading PP header\n Header parsing complete!\n')
           end
           
           %assign header fild values while in MER header text
           if ~pphead && ~datastart
               value=values{tokenindex};
                if strcmpi(token,'date')
                    header.date=value;
                    fprintf('date found\n')
                elseif strcmpi(token,'fs')
                    header.fs=str2num(value);
                    ppheader.fs=str2num(value);
                    fprintf('fs found\n')
                elseif strcmpi(token,'neuron_model')
                    header.neuron_model=value;
                    fprintf('neuron model found\n')
                elseif strcmpi(token,'noise_model')
                    header.noise_type=value;
                    fprintf('noise model found\n')
                elseif strcmpi(token,'noise_para')
                    header.noise_params=value;
                    fprintf('noise model found\n')
                elseif strcmpi(token,'post_low_filter_params')
                    header.low_filter_params=value;
                    fprintf('low pass filter parameters found\n')
                elseif strcmpi(token,'post_low_filter_type')
                    header.high_filter_type=value;
                    fprintf('low pass filter type found\n')
                elseif strcmpi(token,'post_high_filter_type')
                    header.low_filter_type=value;
                    fprintf('high pass filter type found\n')
                elseif strcmpi(token,'post_high_filter_params')
                    header.high_filter_params=value;
                    fprintf('high pass filter parameters found\n')
                elseif strcmpi(token,'distribution')
                    header.firing_distribution=value;
                    ppheader.firing_distribution=value;
                    fprintf('isi distribution found\n')
                elseif strcmpi(token,'distribution_parameters')
                    header.firing_params=str2num(value);
                    ppheader.firing_params=str2num(value);
                    fprintf('distribution parameters found\n')
                elseif strcmpi(token,'num_neurons')
                    header.num_neurons=str2num(value);
                    fprintf('Number of neurons in MER found\n')
                elseif strcmpi(token,'data_endianity')
                    header.data_format.endianity=value;
                    fprintf('MER data endianity found\n')
                elseif strcmpi(token,'data_type')
                    header.data_format.type=value;
                    fprintf('MER data type found\n')
                elseif strcmpi(token,'data_width')
                    header.data_format.width=value;
                    fprintf('MER data width found\n')
                else
                    warning('****Invalid header token: %s****\n',token)
                end            
                
           %assign pp header field values while in PP header text   
           elseif pphead && ~datastart
               value=values{tokenindex};
                if strcmp(token,'%start_PP_Header%')
                    pphead=1;
                    fprintf('Preparing PP header\n')
                elseif strcmp(token,'pp_length')
                    ppheader.ts=str2num(value)*ppheader.fs;
                    header.ts=str2num(value)*ppheader.fs;
                    ppheader.samples=str2num(value);
                    header.samples=str2num(value);
                    fprintf('PP data length found\n');
                elseif strcmp(token,'pp_num_neurons')
                    ppheader.num_neurons=str2num(value);
                    fprintf('PP number of neurons found\n')
                elseif strcmp(token,'pp_data_endianity')
                    ppheader.data_format.endianity=value;
                    fprintf('PP data endianity found\n')
                elseif strcmp(token,'pp_data_type')
                    ppheader.data_format.type=value;
                    fprintf('PP data type found\n')
                elseif strcmp(token,'pp_data_width')        
                    ppheader.data_format.width=value;
                    fprintf('PP data width found\n')
                else
                    warning('****Invalid PP header token: %s****\n',token)
                end 
            end
      end  
end