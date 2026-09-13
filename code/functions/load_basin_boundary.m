function boundary = load_basin_boundary(filePath)
%LOAD_BASIN_BOUNDARY Load a 2 x N [lon;lat] basin-boundary array.
D = load(filePath);
fields = fieldnames(D);
boundary = [];
for i = 1:numel(fields)
    v = D.(fields{i});
    if isnumeric(v) && ndims(v)==2 && size(v,1)==2 && size(v,2)>2
        boundary = double(v);
        return
    end
end
error('No 2 x N basin-boundary array found in %s',filePath);
end
