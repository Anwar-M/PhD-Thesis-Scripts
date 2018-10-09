function [numerical_inp,array_input,array_input2,psd_file,geo_file,sources]= read_input(BB)

%============================================================================================
%Function reads input data from a text file

%Input  - BB , string : name of the input file

%Output - numerical_inp, 1D array : all data that is not an array
%		- array_input, 2D array   :  arrays of type (x,3)
%		- array_input2, 2D array  : arrays of type (x,2)
%		- psd_file, string 		  : name of the PSD vs freq input file
%		- geo_file, string		  : name of the geometry file

%NOTE: As input parameters are added, they can be included by addding positions in the arrays

%============================================================================================

%Open file
infile = fopen(BB,'r');
if(infile == -1)
	error('File can''t be read/found...');
end

%Read data from input file
while ~feof(infile)
	
	line=fgets(infile);
	aux_str = split(line);
	
	if(strcmp('ref_pressure',aux_str(1))==1)
		numerical_inp(1) = str2double(aux_str(2));
	end
	
	if(strcmp('speed_sound',aux_str(1))==1)
		numerical_inp(2) = str2double(aux_str(2));
	end
	
	if(strcmp('AC_position',aux_str(1))==1)
		array_input(1,:) = str2double(aux_str(2:4));
	end
	
	if(strcmp('humidity',aux_str(1))==1)
		numerical_inp(3) = str2double(aux_str(2));
	end
	
	if(strcmp('temperature',aux_str(1))==1)
		numerical_inp(4) = str2double(aux_str(2));
	end
	
	if(strcmp('AC_velocity',aux_str(1))==1)
		array_input(2,:) = str2double(aux_str(2:4));
	end
	
	if(strcmp('time_flight',aux_str(1))==1)
		numerical_inp(5) = str2double(aux_str(2));
	end
	
	if(strcmp('delta_t',aux_str(1))==1)
		numerical_inp(6) = str2double(aux_str(2));
	end
	
	if(strcmp('source_pos',aux_str(1))==1)
		array_input(3,:) = str2double(aux_str(2:4));
	end
	
	if(strcmp('source_file',aux_str(1))==1)
		aux = str2double(aux_str(2));
		stringg = aux_str(3);
		psd_file{aux} = strjoin(stringg);
	end
	
	if(strcmp('geo_file',aux_str(1))==1)
		aux = str2double(aux_str(2));
		stringg = aux_str(3);
		geo_file{aux} = strjoin(stringg);
	end
	
	if(strcmp('Nx',aux_str(1))==1)
		numerical_inp(7) = str2double(aux_str(2));
	end
	
	if(strcmp('Ny',aux_str(1))==1)
		numerical_inp(8) = str2double(aux_str(2));
	end
	
	if(strcmp('P1',aux_str(1))==1)
		array_input2(1,:) = str2double(aux_str(2:3));
	end
	
	if(strcmp('P2',aux_str(1))==1)
		array_input2(2,:) = str2double(aux_str(2:3));
	end
	
	if(strcmp('P3',aux_str(1))==1)
		array_input2(3,:) = str2double(aux_str(2:3));
	end
	
	if(strcmp('P4',aux_str(1))==1)
		array_input2(4,:) = str2double(aux_str(2:3));
	end
	
	if(strcmp('N_sources',aux_str(1))==1)
		numerical_inp(9) = str2double(aux_str(2));
	end
	
	if(strcmp('source_pos',aux_str(1))==1)
		aux = str2double(aux_str(2));
		sources(aux,:) = str2double(aux_str(3:5));
	end
	
end

%Check some parameters
if(size(sources,1)~=numerical_inp(9))
	display('Number of sources and source positions not equal');
	return;
end
if(size(psd_file,2)~=numerical_inp(9))
	display('Number of sources and source strength not equal');
	return;
end
if(size(geo_file,2)~=numerical_inp(9))
	display('Number of sources and geometry not equal');
	return;
end


%Close file
fclose(infile);