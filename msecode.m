function err = msecode(x, y)


validateattributes(x,{'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', ...
    'single','double'},{'nonsparse'},mfilename,'A',1);
validateattributes(y,{'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', ...
    'single','double'},{'nonsparse'},mfilename,'B',1);

if ~isa(x,class(y))
    error(message('images:validate:differentClassMatrices','A','B'));
end
    
if ~isequal(size(x),size(y))
    error(message('images:validate:unequalSizeMatrices','A','B'));
end

if isempty(x) % If x is empty, y must also be empty
    err = [];
    return;
end

if isinteger(x)     
    x = double(x);
    y = double(y);
end

err = ((norm(x(:)-y(:),2).^2)/numel(x))/100;
