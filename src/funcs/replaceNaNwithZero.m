function result = replaceNaNwithZero(data)
    % This function replaces all NaN values in the input with zeros
    %
    % Inputs:
    %   data - Can be a scalar, vector, matrix, or multi-dimensional array
    %
    % Outputs:
    %   result - Same size and type as input, but with all NaNs replaced by 0
    
    % Create a copy of the input data
    result = data;
    
    % Find all NaN values
    nanIndices = isnan(result);
    
    % Replace NaN values with zeros
    result(nanIndices) = 0;
    
    % Display information about replacements
    numReplaced = sum(nanIndices(:));
    if numReplaced > 0
        fprintf('Replaced %d NaN values with zeros.\n', numReplaced);
    else
        fprintf('No NaN values found in the data.\n');
    end
end