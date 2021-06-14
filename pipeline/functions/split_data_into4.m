function split_data_into4(data3D, outputPath, bitDepth)

    data3D_1 = data3D(1:size(data3D,1)/2, 1:size(data3D,2)/2, :);
    data3D_2 = data3D(1:size(data3D,1)/2, size(data3D,2)/2:size(data3D,2), :);   
    data3D_3 = data3D(size(data3D,1)/2:size(data3D,1), 1:size(data3D,2)/2, :);    
    data3D_4 = data3D(size(data3D,1)/2:size(data3D,1), size(data3D,2)/2:size(data3D,2), :);
    
    if bitDepth == 16
         for i = 1 : size(data3D,3)
            imwrite(uint16(data3D_1(:,:,i)), [outputPath '_crop1.tif'], 'WriteMode', 'append');
            imwrite(uint16(data3D_2(:,:,i)), [outputPath '_crop2.tif'], 'WriteMode', 'append');
            imwrite(uint16(data3D_3(:,:,i)), [outputPath '_crop3.tif'], 'WriteMode', 'append');
            imwrite(uint16(data3D_4(:,:,i)), [outputPath '_crop4.tif'], 'WriteMode', 'append');
         end
    else
        for i = 1 : size(data3D,3)
            imwrite(uint8(data3D_1(:,:,i)), [outputPath '_crop1.tif'], 'WriteMode', 'append');
            imwrite(uint8(data3D_2(:,:,i)), [outputPath '_crop2.tif'], 'WriteMode', 'append');
            imwrite(uint8(data3D_3(:,:,i)), [outputPath '_crop3.tif'], 'WriteMode', 'append');
            imwrite(uint8(data3D_4(:,:,i)), [outputPath '_crop4.tif'], 'WriteMode', 'append');
         end
    end
    
    
end