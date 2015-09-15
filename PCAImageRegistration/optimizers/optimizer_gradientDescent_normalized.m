function [ xout, output ] = optimizer_gradientDescent_normalized( cost_function, x0, options)
%OPTIMIZER_GRADIENTDESCENT_NORMALIZED  Minimization using gradient descent with normalized gradient.
% 
% Gradient must be provided by the cost function
%
% Input arguments:
%   options: structure with the following elements:
%       options.DisplayIter - true/false;
%       options.th_fvalue - value
%       options.initialStepNorm - value
%       options.lambda_decrease_rate - value (0,1)
%
% Author: Alberto Gomez, Biomedical Engineering, KCL, 2013

[current_value, current_gradient_ ] = cost_function(x0);
current_lambda =options.initialStepNorm;
current_x = x0;

iterations=0;
output.f_values = zeros(options.max_iterations,1);
output.Af = zeros(options.max_iterations,1);
output.iteration_parameters = zeros(options.max_iterations,numel(x0));
output.Aparams = zeros(options.max_iterations,1);
output.direction_change = zeros(options.max_iterations,1);

if ~isfield(options, 'Scaling')
    options.Scaling = ones(size(x0));
end

if ~isfield(options, 'use_gradient')
    options.use_gradient = true;
end
if ~isfield(options, 'force_iterations')
    options.force_iterations = false;
end

if ~numel(current_gradient_)
    current_gradient = zeros(size(x0'));
    for i=1:numel(current_gradient)
        deltas = zeros(size(x0));
        deltas(i)=1;
        % negative value
        current_x_m = current_x-options.Scaling(i)*deltas*current_lambda/2;
        fm = cost_function(current_x_m);    
        current_x_M = current_x+options.Scaling(i)*deltas*current_lambda/2;
        fM = cost_function(current_x_M);
        
        current_gradient(i) = (fM-fm)/(current_x_M(i)-current_x_m(i));
        %current_gradient(i) = (fM-current_value);
    end
else
    current_gradient = current_gradient_;
end

while true
    iterations = iterations+1;
 
    previous_value = current_value;
    previous_x = current_x;
    current_gradient = current_gradient.*options.Scaling';
    normalized_gradient = current_gradient/norm(current_gradient);
    
    
    %current_gradient_ = current_gradient_.*options.Scaling';
    %normalized_gradient_ = current_gradient_/norm(current_gradient_);
    
    
    current_x = previous_x + current_lambda*(-normalized_gradient').*options.Scaling;
    [current_value, current_gradient_ ] = cost_function(current_x);
    
    if ~numel(current_gradient_) || ~options.use_gradient
        % estimate gradient
        current_gradient = zeros(size(x0'));
        for i=1:numel(current_gradient)
            deltas = zeros(size(x0));
            deltas(i)=1;
            % negative value
            current_x_m = current_x-options.Scaling(i)*deltas*current_lambda/2;
            fm = cost_function(current_x_m);
            current_x_M = current_x+options.Scaling(i)*deltas*current_lambda/2;
            fM = cost_function(current_x_M);
            
            current_gradient(i) = (fM-fm)/(current_x_M(i)-current_x_m(i));
            %current_gradient(i) = (fM-current_value);
        end
    else
        current_gradient = current_gradient_;
    end
    
    % check if gradient changed directions. If so, divide lambda by two.
    if sign(normalized_gradient*current_gradient') <0
        current_lambda = current_lambda*options.lambda_decrease_rate;
        output.direction_change(iterations)=1;
    end
    
    Af = (current_value - previous_value)^2/previous_value^2;
    Aparams = norm(current_x - previous_x)/norm(previous_x);
    
    % save current iteration
    output.f_values(iterations)=current_value;
    output.iteration_parameters(iterations,:)=current_x;
    output.Aparams(iterations)=Aparams;
    output.Af(iterations)=Af;
    if options.DisplayIter
        disp(['# ' num2str(iterations) '- Value: ' num2str(current_value)  ', Residual: ' num2str(Af) ', Aparams: ' num2str(Aparams) ])
    end
    
    
    % See if we have finished
    if Af < options.th_Afvalue && ~options.force_iterations
        if options.DisplayIter
            disp(['Function decrease no longer meaningful: ' num2str(Af) '<' num2str(options.th_Afvalue)])
        end
        break;
    end
    
    if Aparams < options.th_Aparams && ~options.force_iterations
        if options.DisplayIter
            disp(['Convergence attained Aparams: ' num2str(Aparams) '<'  num2str(options.th_Aparams)])
        end
        break;
    end
    
    if iterations > options.max_iterations
        if options.DisplayIter
            disp(['Maximum number of iterations attained: ' num2str(iterations) ])
        end
        break;
    end
    
end

xout = current_x;
output.f_value = current_value;
output.f_gradient = current_gradient;
output.number_of_iterations =iterations ;


end

