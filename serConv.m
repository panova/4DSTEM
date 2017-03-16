%serConv() function
%	Description:
%		Convert all .ser image files in one folder to .tif (8-bit and 16-bit).
%	Requires:
%		serReader.m
%
%Written by Peter Ercius 2009/08/06
%
%---------------------- NO WARRANTY ------------------
%THIS PROGRAM IS PROVIDED AS-IS WITH ABSOLUTELY NO WARRANTY
%OR GUARANTEE OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
%INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
%MERCHANABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
%IN NO EVENT SHALL THE AUTHOR BE LIABLE
%FOR DAMAGES RESULTING FROM THE USE OR INABILITY TO USE THIS
%PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA
%BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR
%THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH
%ANY OTHER PROGRAM). 
%------------------------------------------------------------------------

function serConv(dirname)

    if nargin < 1
       dirname = uigetdir();
    end
    
    if dirname == 0
        disp('No directory selected. Exiting...')
        return;
    end
    binfilenames = dir(fullfile(dirname,'*.ser'));

    for ii=1:length(binfilenames)
        %FID = fopen(fullfile(dirname,binfilenames(ii).name),'rb');

        disp(fullfile(dirname,binfilenames(ii).name))

        Ser = serReader(fullfile(dirname,binfilenames(ii).name));

        if isfield(Ser,'image') && ~isempty(Ser.image)
            [status,message,messageid] = mkdir(dirname,'8bit'); %make directory for 8-bit mages
            imageData = imrotate(Ser.image,90);

			imageData_8 = localBytscl(imageData);
			
            %Convert to 16-bit integers
            imageData = uint16(imageData);

            %Write out the image in 8-bit tif format
            nn = strfind(binfilenames(ii).name,'.ser'); %get the filename only
            imwrite(imageData_8,fullfile(dirname,'8bit',[binfilenames(ii).name(1:nn-1) '.tif']),'tif');

            %Write out the 16-bit tif image
            imwrite(imageData,fullfile(dirname,[binfilenames(ii).name(1:nn-1) '_16bit.tif']),'tif');
        end

        %Multiple spectra or image files
        if isfield(Ser,'data')
            nn = strfind(binfilenames(ii).name,'.ser'); %get the filename only
            imageSize = size(Ser.data{1});
            dataSize = length(Ser.data);
            %Check for normal images
            if imageSize(1) == imageSize(2)
                [status,message,messageid] = mkdir(dirname,'8bit'); %make directory for 8-bit mages
                imwrite(localBytscl(Ser.data{1}),fullfile(dirname,'8bit',[binfilenames(ii).name(1:nn-1) '.tif']),'tif')
                if dataSize > 1
                    curFname = fullfile(dirname,'8bit',[binfilenames(ii).name(1:nn-1) '.tif']);
                    for jj = 2:dataSize
                        t = Tiff(curFname,'a');
                        t.setTag('Photometric',Tiff.Photometric.MinIsBlack);
                        t.setTag('Compression',Tiff.Compression.None);
                        t.setTag('BitsPerSample',8);
                        t.setTag('SamplesPerPixel',1);
                        t.setTag('SampleFormat',Tiff.SampleFormat.UInt);
                        t.setTag('ExtraSamples',Tiff.ExtraSamples.Unspecified);
                        t.setTag('ImageLength',imageSize(1));
                        t.setTag('ImageWidth',imageSize(2));
                        t.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
                        t.write(localBytscl(Ser.data{jj}));
                        t.close();
                    end
                    disp('Only 8-bit series written')
                end

            %Check for normal spectra
            elseif( (imageSize(1) == 2048) || (imageSize(2) == 2048) )
                specS = size(Ser.data{1});
                out = [Ser.calibration{1}' cell2mat(Ser.data)]; 
                save(fullfile(dirname,[binfilenames(ii).name(1:nn-1) '.txt']), 'out','-ascii')
            else
                disp('Unknown 3D data set type. Unable to extract the data.')
            end
        end
    end
end

function imageData_8 = localBytscl(x)
	
    arInfo = size(x);
    
    xx = double(x(:));
    
    minval = min(xx(:));
    xNorm = xx(:) - minval; %set min to 0
    maxval = max(xNorm);
    xNorm = xNorm ./ maxval; %normlize values to 0 - 1

    xNorm(xNorm < 0) = 0;
    xNorm(xNorm > 1) = 1;

    imageData_8 = uint8(reshape(1+fix(xNorm .* 255),arInfo(1),arInfo(2))); %use 1+fix(...) because image(y) expects values 1 - 256 (Matlab values start at 1 not 0)
end

% function writeTiffStackSlice(image,fname)
%     dataSize = size(image);
%     
%     t = Tiff(fname,'a');
%     t.setTag('Photometric',Tiff.Photometric.MinIsBlack);
%     t.setTag('Compression',Tiff.Compression.None);
%     t.setTag('BitsPerSample',8);
%     t.setTag('SamplesPerPixel',1);
%     t.setTag('SampleFormat',Tiff.SampleFormat.UInt);
%     t.setTag('ExtraSamples',Tiff.ExtraSamples.Unspecified);
%     t.setTag('ImageLength',dataSize(1));
%     t.setTag('ImageWidth',dataSize(2));
%     t.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
%     t.write(image);
%     t.close();
% end