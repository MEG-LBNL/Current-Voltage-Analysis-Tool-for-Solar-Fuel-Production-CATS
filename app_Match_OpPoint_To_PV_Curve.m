function [intensity, shading] = app_Match_OpPoint_To_PV_Curve(OpVolt, ...
    OpCurrent, Rs, Rsh, n, Isc_initial, Voc_initial, TempPV)
%Calculates the light intensity (shading of PV) by water vapor as it
%matches to operating point to a correspondent PV I-V curve

%% Find curve for operating point -> determine effective light intensity
intensity = zeros(size(OpVolt,1),1);
    
options = optimset('TolFun', 1e-15,'Display','off');
D = parallel.pool.DataQueue;
h = waitbar(0, 'Please wait ...'); %Start waitbar
afterEach(D, @nUpdateWaitbar); %When data is sent, update waitbar
p = 1; %Progress value in waitbar

parfor iPoint = 1:size(OpVolt,1)
    Operating_point = [OpVolt(iPoint) -OpCurrent(iPoint)];%[V],[A]
    x0 = Operating_point(2)/Isc_initial; %Start the search at x0 sun
    %Calc. the difference between a curve at x light int. and the point
    intensity(iPoint) = fsolve(@(x) ...
        abs(app_calculate_IV_PV(Operating_point(1), Isc_initial, ...
        Voc_initial, x, Rs, Rsh, n, TempPV) -...
        Operating_point(2)),x0, options);
    if ~mod(iPoint, 100), send(D, iPoint); end %Update waitbar
end
close(h) %Close waitbar

shading = 1 - intensity;

function nUpdateWaitbar(~)
    waitbar(p/size(OpVolt,1), h); %Update progress
    p = p + 100; %Increase by 100 [unit?]
end
end