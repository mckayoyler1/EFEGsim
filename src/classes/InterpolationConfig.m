classdef InterpolationConfig
    properties
        Method          % Interpolation method (e.g., 'linear', 'nearest')
        Extrapolation   % Extrapolation method (e.g., 'none', 'linear')
    end

    methods
        function obj = InterpolationConfig(method, extrapolation)
            % Constructor with default values
            if nargin < 1 || isempty(method)
                method = 'linear'; % Default interpolation method
            end

            if nargin < 2 || isempty(extrapolation)
                extrapolation = 'none'; % Default extrapolation method
            end

            obj.Method = method;
            obj.Extrapolation = extrapolation;
        end

        function validateConfig(obj)
            % Validate the interpolation configuration
            validMethods = {'linear', 'nearest', 'natural', 'cubic', 'v4'};
            validExtrapolations = {'none', 'linear', 'nearest'};

            if ~ismember(obj.Method, validMethods)
                error("Invalid interpolation method '%s'. Must be one of: %s.", obj.Method, strjoin(validMethods, ', '));
            end

            if ~ismember(obj.Extrapolation, validExtrapolations)
                error("Invalid extrapolation method '%s'. Must be one of: %s.", obj.Extrapolation, strjoin(validExtrapolations, ', '));
            end
        end
    end
end