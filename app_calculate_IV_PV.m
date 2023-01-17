function [CurrentMod, I0] = app_calculate_IV_PV(Voltage, Isc_initial,...
    Voc_initial, light_concentration, Rs, Rsh, n, TempPV)
%Outputs the calculated current for a given voltage which can be a whole IV
%curve or one operating point, uses the characteristic PV parameters such
%as Rs, Rsh and n and the light concentration (can be lower than 1 sun)

kB = 1.380649 * 10^-23; %Boltzmann constant [J/K]
q = 1.602176634 * 10^-19; %Elementary charge [C]

%Short circuit current [A] is dir. proportional to light intensity
Isc = Isc_initial.*light_concentration;
%Open circuit potential [V] has a log. prop. to light intensity
%Should only use this equation if the temperature is constant
Voc = Voc_initial + n*kB.*TempPV/q.*log(light_concentration);

%Reverse sat. current [A], can be calculated from Voc (V = Voc and I = 0)
%high quality cell: IL = Isc
I0 = (Isc-Voc/Rsh)./(exp(q.*Voc/n/kB./TempPV)-1);

CurrentMod = zeros(size(Voltage,1),1);
CurrentGuess = Isc; %[A]

%Evaluate single-diode solar cell equation
for i = 1:size(Voltage,1)
    funDiode = @(iC) ...
        (Isc - I0*(exp((Voltage(i)+iC*Rs)/(n*kB*TempPV/q))-1))-...
        ((Voltage(i) + iC*Rs)/Rsh)-iC;
    if(false && license('test', 'optimization_toolbox'))
        options = optimoptions('lsqnonlin','Display','off');
        CurrentMod(i) = lsqnonlin(funDiode,CurrentGuess,0,12,options);
    else
        CurrentMod(i) = fzero(funDiode, CurrentGuess);
    end
    CurrentGuess = CurrentMod(i); %Make last calculated value next guess
end

end