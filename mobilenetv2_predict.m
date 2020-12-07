
% mobilenetv2 function
function out = mobilenetv2_predict(in) 
persistent mynet;

if isempty(mynet)
    mynet = coder.loadDeepLearningNetwork('mobilenetv2','mobilenetv2');
end
 
out = mynet.predict(in);
