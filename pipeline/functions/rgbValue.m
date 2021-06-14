function rgbV=rgbValue(col, exp)
    if col<0
        rgbV=[ 0 0 1];
    elseif col<1/4
        rgbV=[ 0 (col*4)^exp 1];
    elseif col<2/4
        rgbV=[ 0 1 (2-col*4)^exp];
    elseif col<3/4
        rgbV=[(col*4-2)^exp 1 0];
    elseif col<1
        rgbV=[1 (4-col*4)^exp 0];
    else
        rgbV=[1 0 0];
    end
end