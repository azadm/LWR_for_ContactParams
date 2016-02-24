function lwrModel = LWRModel(X_train, Y_train, sigma, cutoff)
% learned mapping from jointAngles(motor commands) to center of pressure 
% using LWR
% computes the forward mapping
% and the Jacobian (1st order derivative) for control 
%
% Notation: rows are samples, columns are dimensions
%
%E. Rueckert, Sept 2015, rueckert.elmar@gmail.com

lwrModel.flgNormalizeData = 1;

if lwrModel.flgNormalizeData
    [lwrModel.Inputs, lwrModel.PS] = normalize(X_train);
else
    lwrModel.Inputs = X_train;
    lwrModel.PS = [];
end
lwrModel.Outputs = Y_train;

lwrModel.numInputDimensions = size(X_train,2);
lwrModel.numPredictors = size(Y_train,2);
lwrModel.numSamples = size(X_train,1);

if nargin < 3
    sigma = 0.1;
end

sigmaDimOne = 1/sigma*ones(1,lwrModel.numInputDimensions); 
if nargin < 4
    cutoff = 0;
end
cutoffInit = cutoff;
MinNumPts = 2;
lwrModel.flgUseMex = 1;

%Forward model prediction and Jacobian (optional)
lwrModel.predict = @(queryPoint_)predictWithFeatureLocal(queryPoint_, lwrModel, sigmaDimOne, cutoffInit, MinNumPts);
            
lwrModel.normalizeQueryPts = @(queryPoints_)normalizeQueryPoint(queryPoints_, lwrModel.PS, lwrModel.flgNormalizeData);

end

function [y,PS] = normalize(x)
    [yTrans, PS] = mapstd(x');
    y = yTrans';
end

function [y] = denormalize(x,PS)
    y = mapstd('reverse',x',PS)';
end

function [y] = normalizeQueryPoint(x, PS, flgNormalizeData)
    if flgNormalizeData
        y = mapstd('apply',x',PS)';
    else
        y = x;
    end
end

function [predictions, Jacobians, betas] = predictWithFeatureLocal(queryPoint, lwrModel, sigmaDimOne, cutoffInit, MinNumPts)

queryPoint = normalizeQueryPoint(queryPoint, lwrModel.PS, lwrModel.flgNormalizeData);

if lwrModel.flgUseMex
    fnEval = @MexLocalLWR;
else
    fnEval = @MatlabLocalLWR;
end
flgComputeJacobian =  nargout > 1;

predictions = zeros(1,lwrModel.numPredictors);
if flgComputeJacobian
    Jacobians = nan(lwrModel.numPredictors, lwrModel.numInputDimensions);
    betas = nan(lwrModel.numPredictors, 1+lwrModel.numInputDimensions);
end
X = lwrModel.Inputs;
Y = lwrModel.Outputs;
% queryPoint

 for predictorId=1:lwrModel.numPredictors

    y = Y(:,predictorId);
    if flgComputeJacobian
        [prediction, beta] = fnEval(X, y, queryPoint, sigmaDimOne, cutoffInit, MinNumPts);
%         prediction - (beta(1) + sum(beta(2:end).*queryPoint))
%         beta(1) 
%         keyboard
        predictions(1,predictorId) = prediction;
        Jacobians(predictorId,:) = beta(2:end);
        betas(predictorId,:) = beta;
    else
        prediction = fnEval(X, y, queryPoint, sigmaDimOne, cutoffInit, MinNumPts);
        predictions(1,predictorId) = prediction;
    end
    
 end
% betas
end

function [prediction, beta] = MatlabLocalLWR(X, y, queryPoint, sigmaDimOne, cutoff, MinNumPts)

    prediction = nan;

    Dprecision = diag(sigmaDimOne);
    
    xDiff = bsxfun(@minus, X, queryPoint);
    w = exp(-0.5*diag(xDiff*Dprecision*xDiff'));
    
%     figure;
%     hist(w)
    
    valid_ids = find(w > cutoff);
    if 1==1
        whileIts = 0;
        while isempty(valid_ids) || length(valid_ids) < MinNumPts 
            cutoff = cutoff*0.5;
            valid_ids = find(w > cutoff);
            whileIts = whileIts + 1;
            if whileIts > 30
                break
            end
        end
    end
    
    numValid = length(valid_ids);
%     [numValid cutoff]
    if numValid < 2
        warning('not neough points found');
        return;
    end
    
%     keyboard
    XSubSet = [ones(length(valid_ids),1) X(valid_ids,:)];
    ySubSet = y(valid_ids);

    W = diag(w(valid_ids));
    beta = ((1e-6 + sqrt(W)*XSubSet)\(sqrt(W)*ySubSet));
    %beta = (XSubSet'*sqrt(W)*XSubSet) \ (XSubSet'*sqrt(W)*ySubSet);

    prediction = beta(1) + queryPoint*beta(2:end); 
   
end
