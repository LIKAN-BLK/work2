function chans = read_param_struct(param, parameters)

chans = zeros(1, size(param, 2));

for i = 1:length(param)
    chans(i) = find(cellfun(@(x) ~isempty(x) && x==1 , strfind(parameters.ChannelNames.Value,char(param(i)))));
end