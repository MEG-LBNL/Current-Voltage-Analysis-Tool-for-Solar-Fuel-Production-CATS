classdef CATS_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        CATSUIFigure                  matlab.ui.Figure
        GridLayout                    matlab.ui.container.GridLayout
        LeftPanel                     matlab.ui.container.Panel
        v01Label                      matlab.ui.control.Label
        TabGroup2                     matlab.ui.container.TabGroup
        GeneralTab                    matlab.ui.container.Tab
        ParallelCheck                 matlab.ui.control.CheckBox
        CalculationoptionsPanel       matlab.ui.container.Panel
        FactortoincreasecalculationspeedEditField  matlab.ui.control.NumericEditField
        FactortoincreasecalculationspeedEditFieldLabel  matlab.ui.control.Label
        E0EditField                   matlab.ui.control.NumericEditField
        MinimumpotentialdifferenceforfullcellreactionofinterestLabel  matlab.ui.control.Label
        EndsEditField                 matlab.ui.control.NumericEditField
        EndsEditFieldLabel            matlab.ui.control.Label
        StartsEditField               matlab.ui.control.NumericEditField
        StartsEditFieldLabel          matlab.ui.control.Label
        StartandendtimeforcurrentvoltageanalysisLabel  matlab.ui.control.Label
        MethodoptionsPanel            matlab.ui.container.Panel
        StablePVparametersCheckBox    matlab.ui.control.CheckBox
        ViasimpleareaestimationCheckBox  matlab.ui.control.CheckBox
        CanECIVcurvesbeestimatedCheckBox  matlab.ui.control.CheckBox
        TimedependentECIVcurvesavailableCheckBox  matlab.ui.control.CheckBox
        ImportsettingsPanel           matlab.ui.container.Panel
        UnitsLabel                    matlab.ui.control.Label
        HeaderrowsEditField           matlab.ui.control.NumericEditField
        HeaderrowsEditFieldLabel      matlab.ui.control.Label
        TimeDropDown                  matlab.ui.control.DropDown
        TimeDropDownLabel             matlab.ui.control.Label
        VoltageDropDown               matlab.ui.control.DropDown
        VoltageDropDownLabel          matlab.ui.control.Label
        CurrentDropDown               matlab.ui.control.DropDown
        CurrentDropDownLabel          matlab.ui.control.Label
        DelimiterDropDown             matlab.ui.control.DropDown
        DelimiterDropDownLabel        matlab.ui.control.Label
        PVfittingTab_2                matlab.ui.container.Tab
        ValuesforcalculationPanel     matlab.ui.container.Panel
        MinimumvoltagePVParameter     matlab.ui.control.NumericEditField
        MinimumvoltageforPVparameterestimationVLabel  matlab.ui.control.Label
        SpeedFactor                   matlab.ui.control.NumericEditField
        FactortoincreasespeedforPVparameterestimationEditFieldLabel  matlab.ui.control.Label
        PVtemperatureCEditField       matlab.ui.control.NumericEditField
        PVtemperatureCLabel           matlab.ui.control.Label
        PVsizecm2EditField            matlab.ui.control.NumericEditField
        PVsizecm2EditFieldLabel       matlab.ui.control.Label
        CalculatedparametersPanel     matlab.ui.container.Panel
        I0AEditField                  matlab.ui.control.NumericEditField
        I0AEditFieldLabel             matlab.ui.control.Label
        nEditField                    matlab.ui.control.NumericEditField
        nEditFieldLabel               matlab.ui.control.Label
        RshOhmEditField               matlab.ui.control.NumericEditField
        RshOhmLabel                   matlab.ui.control.Label
        RsOhmEditField                matlab.ui.control.NumericEditField
        RsOhmEditFieldLabel           matlab.ui.control.Label
        VocnoshadingVEditField        matlab.ui.control.NumericEditField
        VocnoshadingVEditFieldLabel   matlab.ui.control.Label
        IscnoshadingAEditField        matlab.ui.control.NumericEditField
        IscnoshadingAEditFieldLabel   matlab.ui.control.Label
        IVsyncTab                     matlab.ui.container.Tab
        ColumnnumberLabel             matlab.ui.control.Label
        DelimiterLabel                matlab.ui.control.Label
        VoltageDropDown_2             matlab.ui.control.DropDown
        VoltageDropDown_2Label        matlab.ui.control.Label
        CurrentDropDown_2             matlab.ui.control.DropDown
        CurrentDropDown_2Label        matlab.ui.control.Label
        TimeEditField                 matlab.ui.control.NumericEditField
        TimeEditFieldLabel            matlab.ui.control.Label
        CurrentEditField_2            matlab.ui.control.NumericEditField
        CurrentEditField_2Label       matlab.ui.control.Label
        VoltageEditField_2            matlab.ui.control.NumericEditField
        VoltageEditField_2Label       matlab.ui.control.Label
        Starttimesformatyyyymmdd_HHMMSSLabel  matlab.ui.control.Label
        CurrentEditField              matlab.ui.control.EditField
        CurrentEditFieldLabel         matlab.ui.control.Label
        VoltageEditField              matlab.ui.control.EditField
        VoltageEditFieldLabel         matlab.ui.control.Label
        CenterPanel                   matlab.ui.container.Panel
        TabGroup                      matlab.ui.container.TabGroup
        ImporteddataTab               matlab.ui.container.Tab
        ChangetohButton               matlab.ui.control.Button
        ChangetomAButton              matlab.ui.control.Button
        UIAxes                        matlab.ui.control.UIAxes
        PVfittingTab                  matlab.ui.container.Tab
        UIAxes2_2                     matlab.ui.control.UIAxes
        UIAxes2                       matlab.ui.control.UIAxes
        CalculatedresultsTab          matlab.ui.container.Tab
        ChangetohButton_2             matlab.ui.control.Button
        ChangetomAButton_2            matlab.ui.control.Button
        UIAxes_3                      matlab.ui.control.UIAxes
        UIAxes_2                      matlab.ui.control.UIAxes
        RightPanel                    matlab.ui.container.Panel
        SavecalculatedresultsButton   matlab.ui.control.Button
        SynccurrentvoltagedataButton  matlab.ui.control.Button
        ImporttimedependentcurrentvoltagedataButton  matlab.ui.control.Button
        StartcalculationButton        matlab.ui.control.Button
        EstimatePVparametersButton    matlab.ui.control.Button
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
        twoPanelWidth = 768;
    end

    
    properties (Access = private)
        RawDataIV % Saves imported, time-dependent current-voltage data
        PlotUnits % Units for imported data
        PlotUnits_D % Units for deconvoluted data
        RawDataPV % IV curves of PV for parameter estimation
        PVparameters % Characteristic PV parameters
        DeconvLoss % Deconvoluted losses (and gains) calculated
    end
    
    methods (Access = private)
        
        function RawDataIV = importIV(app, strMsg)
            %function to import the data
            %Choose specific message for import if desired
            switch nargin 
                case 1
                    strImportMsg = 'Choose IV data';
                case 2
                    strImportMsg = strMsg;
            end
            Delimit = app.DelimiterDropDown.Value;
            switch Delimit
                case 'Comma'
                    strDelimit = ',';
                case 'Tab'
                    strDelimit = '\t';
                case 'Space'
                    strDelimit = ' ';
            end
            CurrentUnits = app.CurrentDropDown.Value;
            VoltageUnits = app.VoltageDropDown.Value;
            TimeUnits = app.TimeDropDown.Value;
            iHeaderLines = app.HeaderrowsEditField.Value;

            %Select file(s)
            [filename,pathName] = uigetfile({'*.txt; *.csv; *.mpt',...
                'All files (*.txt; *.csv; *.mpt)'},strImportMsg,'',...
                'MultiSelect','on');
            if(size(filename,2) == 1)
                RawDataIV = {}; %If no file selected, stop program
                return;
            end
            
            if(iscell(filename))
                filename = sort(filename); %Sort alphabetically if list of files
            end
            
            if iscell(filename) %For a list of files
                FullfileNameTxt = fullfile(pathName,filename);
            else %Only one file
                filename = cellstr(filename);
                FullfileNameTxt = cellstr(fullfile(pathName,filename));
            end
            IVData = cell(size(filename,2),2);
            
            %Load file(s)
            for iSize = 1:size(filename,2)
                    Time = 1;
                    Voltage = 2;
                    Current = 3;
                    Columns = [Time Voltage Current];
                    IVData{iSize,1} = ...
                        importdata(FullfileNameTxt{1,iSize},strDelimit,iHeaderLines);
                    IVData{iSize,1}.data = IVData{iSize,1}.data(:,Columns);
                    IVData{iSize,1}.colheaders = IVData{iSize,1}.textdata;
                    IVData{iSize,2} = filename{1,iSize};
            end
            
            for iSize = 1:size(filename,2)
                if(strncmp(TimeUnits,'h',1))
                    IVData{iSize,1}.data(:,1) = IVData{iSize,1}.data(:,1)*3600; %Convert h to s
                elseif(strncmp(TimeUnits,'min',3))
                    IVData{iSize,1}.data(:,1) = IVData{iSize,1}.data(:,1)*60; %Convert min to s
                end
                if(strncmp(VoltageUnits,'mV',2))
                    IVData{iSize,1}.data(:,2) = IVData{iSize,1}.data(:,2)/1000; %Convert mV to V
                end
                if(strncmp(CurrentUnits,'mA',2))
                    IVData{iSize,1}.data(:,3) = IVData{iSize,1}.data(:,3)/1000; %Convert mA to A
                end
            end
            
            RawDataIV = IVData;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plotIVt(app, RawDataIV, PlotUnits)
            %Plots measured current and voltage over time
            TimeData = RawDataIV{1,1}.data(:,1);
            VoltageData = RawDataIV{1,1}.data(:,2);
            CurrentData = RawDataIV{1,1}.data(:,3);
            
            if(strncmp(PlotUnits{1,1},'h',1))
                TimeData = TimeData/3600; %s to h
                xAxisLab = 'Time [h]';
            else
                xAxisLab = 'Time [s]';
            end
            
            yyaxis(app.UIAxes,'left')
            app.UIAxes.XLabel.String = xAxisLab;
            app.UIAxes.YLabel.String ='Voltage [V]';
            plot(app.UIAxes, TimeData, VoltageData)
            
            yyaxis(app.UIAxes,'right')
            if(strncmp(PlotUnits{3,1},'mA',2))
                CurrentData = CurrentData*1000; %A to mA
                app.UIAxes.YLabel.String='Current [mA]';
            else
                app.UIAxes.YLabel.String='Current [A]';
            end
            plot(app.UIAxes,TimeData, CurrentData)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function PVparameters = estimate_PV_Param(app)
            %Estimates series resistance, shunt resistance and ideality factor from PV
            %I-V curves
            PVparameters = table();
            
            %CellSize = app.PVsizecm2EditField.Value; %PV size [cm2]
            TempPV = 273 + app.PVtemperatureCEditField.Value; %[°C] -> [K]  
            
            %% Import current and voltage data
            app.RawDataPV = importIV(app);
            NumberOfFiles = size(app.RawDataPV,1);
            if(NumberOfFiles == 0), return; end
            
            IVData = cell(1,NumberOfFiles);
            
            for FileNumber = 1:NumberOfFiles
                mask = true(length(app.RawDataPV{FileNumber,1}.data(:,1)),1);
                iFirstHalf = sum(mask)/2; %Divide by 2 if only positive or neg. part
                iOnes = 0;
                for iPoint = 1:size(mask,1)
                    if(mask(iPoint))
                        iOnes = iOnes + 1;
                    end
                    if(iOnes > iFirstHalf) %> for increasing curve, < for decr. curve
                       mask(iPoint) = 0; 
                    end
                end
                IVData{1,FileNumber}(:,1) = app.RawDataPV{FileNumber,1}.data(mask,2); %Voltage
                IVData{1,FileNumber}(:,2) = app.RawDataPV{FileNumber,1}.data(mask,3); %Current
            end
            
            %% Starting values for iteration
            Rs = 15; %Series resistance [Ohm]
            Rsh = 100000; %Shunt resistance [Ohm]
            n = 3; %Ideality factor
            
            %% Calculate short circuit curr.(Isc)and open circuit volt.(Voc)from curves
            Voc_V = zeros(size(IVData,2),1);
            Isc_V = zeros(size(IVData,2),1);
            for iCurve = 1:size(IVData,2)
                [~, Voc_V(iCurve), Isc_V(iCurve), ~, ~, ~] = calc_FF(app, ...
                    [IVData{1,iCurve}(:,1) abs(IVData{1,iCurve}(:,2))], 1, 2, false);
                vVoltage = IVData{1,iCurve}(:,1);
                vCurrent = IVData{1,iCurve}(:,2);
                vCurrent(vVoltage > Voc_V(iCurve)) = []; %Remove values beyond Voc
                vVoltage(vVoltage > Voc_V(iCurve)) = []; %Remove values beyond Voc
                IVData{1,iCurve} = [];
                IVData{1,iCurve}(:,1) = vVoltage;
                IVData{1,iCurve}(:,2) = abs(vCurrent); %Make all values positive
            end
            
            %This assumes that the curve at 1 sun (no shading) is provided
            %and no curves above 1 sun are recorded
            %light_concentration_V = Isc_V./max(Isc_V);
            
            light_concentration_V = ones(length(Isc_V),1);
            
            %% Calculate I-V curve
            %Loop through voltage [V] to obtain I-V curve
            xValues = zeros(size(IVData,2),3); %Rs, Rsh, n
            x0 = [Rs, Rsh, n]; %Starting values for iteration
            nthVal = app.SpeedFactor.Value; %Increase speed, use only every nth value of the curve
            minVo = app.MinimumvoltagePVParameter.Value; %Consider only values in the operating voltage range [V]
            MaxError = zeros(length(IVData),1);
            TotError = zeros(length(IVData),1);
            %Res_norm = zeros(length(IVData),1);
            if(size(IVData,2) > 1 && app.ParallelCheck.Value)
              parforArg = Inf; %Run regular parfor loop
            else
              parforArg = 0; %Runs essentially like a for loop
            end

            parfor (iIter = 1:size(IVData,2), parforArg) %Extract parameters
                %Get voltage data for curve fit [V]
                Voltage = IVData{1,iIter}(:,1); 
                %Get experimental current data for curve fit [A]
                CurrentExp = IVData{1,iIter}(:,2);
                
                %Remove NaN values
                Voltage(isnan(Voltage))=[]; 
                CurrentExp(isnan(CurrentExp))=[];
                
                %Consider only values in the operating voltage range [V]
                maxVo = Voc_V(iIter);
                CurrentExp(Voltage > maxVo) = []; %Consider only values until maxVo
                Voltage(Voltage > maxVo) = []; %Values > maxVo impossible dur.oper.
                CurrentExp(Voltage < minVo) = []; %Consider only values above minVo
                Voltage(Voltage < minVo) = []; %Values < minVo uncommon dur.oper.
                Isc_initial = Isc_V(iIter);
                Voc_initial = Voc_V(iIter);
                
                %Reduce the amount of values taken into account
                Voltage = Voltage(1:nthVal:floor(length(Voltage)/nthVal)*nthVal);
                CurrentExp = ...
                    CurrentExp(1:nthVal:floor(length(CurrentExp)/nthVal)*nthVal);
        
                %Function that calculates difference between exp. and model data
                fun = @(x,Voltage)(app_calculate_IV_PV(Voltage, ...
                    Isc_initial, Voc_initial, light_concentration_V(iIter),...
                    x(1), x(2), x(3), TempPV));
                    
                %Minimize difference between exp. and model
                options = optimoptions('lsqcurvefit','TolFun',1e-12,...
                    'Display','off','TolX',1e-10);
                [x,~,residual,~,~] = lsqcurvefit(fun,x0,...
                    Voltage,CurrentExp,[0 0 1],[100 10^7 100],options);
                %Save values
                %Res_norm(iIter) = resnorm;
                TotError(iIter) = sum(abs(residual));
                MaxError(iIter) = max(abs(residual));
                xValues(iIter,:) = x;
                
                fprintf('%f\t%f\t%f\n',x(1),x(2),x(3)) %Results for Rs, Rsh and n
            end
            
            %% plot results
            I0 = zeros(size(light_concentration_V,1),1);
            for iPlot = 1:size(light_concentration_V,1)
                Voltage = IVData{1,iPlot}(:,1);
                Voltage(isnan(Voltage))=[];
                CurrentExp = IVData{1,iPlot}(:,2);
                CurrentExp(isnan(CurrentExp))=[];
                Isc_initial = Isc_V(iPlot);
                Voc_initial = Voc_V(iPlot);
                
                %Calculate modeled curves
                [CurrentMod, I0(iPlot)] = app_calculate_IV_PV(Voltage, Isc_initial, ...
                    Voc_initial,light_concentration_V(iPlot),...
                    xValues(iPlot,1),xValues(iPlot,2), xValues(iPlot,3), TempPV);
                
                plot(app.UIAxes2, Voltage, CurrentExp*1000, 'LineWidth', 1)
                hold(app.UIAxes2, 'on')
                plot(app.UIAxes2, Voltage, CurrentMod*1000, 'Color', [0.8 0.8 0.8], 'LineWidth', 1)
            end
            xlim(app.UIAxes2, [minVo 2.5])
            ylim(app.UIAxes2, [0 ceil(max(Isc_V)*1000)])
            legend(app.UIAxes2, 'Experimental', 'Modeled')
            
            %% plot results with averaged values
            I0Ave = zeros(size(light_concentration_V,1),1);
            MaxErrorAve = zeros(length(light_concentration_V),1);
            TotErrorAve = zeros(length(light_concentration_V),1);
            ColorOrd = get(0,'defaultAxesColorOrder');
            
            %Plot original and modeled curves
            for iPlot = 1:size(light_concentration_V,1)
                Voltage = IVData{1,iPlot}(:,1);
                CurrentExp = IVData{1,iPlot}(:,2);
                
                %Remove NaN values
                CurrentExp(isnan(CurrentExp))=[];
                Voltage(isnan(Voltage))=[];
                
                plot(app.UIAxes2_2, Voltage, CurrentExp*1000, 'Color', ColorOrd(iPlot,:),...
                    'LineWidth', 2), hold(app.UIAxes2_2, 'on')
                
                Isc_initial = Isc_V(iPlot);
                Voc_initial = Voc_V(iPlot);
                
                %Calculate values 
                [CurrentMod, I0Ave(iPlot)] = app_calculate_IV_PV(Voltage, Isc_initial, ...
                    Voc_initial,light_concentration_V(iPlot),...
                    mean(xValues(:,1)),mean(xValues(:,2)), mean(xValues(:,3)), TempPV);
                
                plot(app.UIAxes2_2, Voltage(1:10:floor(length(Voltage)/10)*10), ...
                    CurrentMod(1:10:floor(length(Voltage)/10)*10)*1000, ...
                    'o', 'Color', ColorOrd(iPlot,:), 'LineWidth', 1)
                
                %Calc. error, consider only values in the operating voltage range [V]
                maxVo = Voc_V(iPlot);
                TotErrorAve(iPlot) = ...
                    sum(abs(CurrentMod(Voltage < maxVo & Voltage > minVo) - ...
                    CurrentExp(Voltage < maxVo & Voltage > minVo)));
                MaxErrorAve(iPlot) = ...
                    max(abs(CurrentMod(Voltage < maxVo & Voltage > minVo) - ...
                    CurrentExp(Voltage < maxVo & Voltage > minVo)));
            end
            
            xlim(app.UIAxes2_2, [minVo 2.5])
            ylim(app.UIAxes2_2, [0 ceil(max(Isc_V)*1000)])
            legend(app.UIAxes2_2, 'Experimental', 'Modeled')
            
            %% Save values in tables
            IterationResults = table();
            IterationMeanResults = table();
            
            IterationResults.Intensity = light_concentration_V;
            IterationResults.Properties.VariableUnits{'Intensity'} = 'sun';
            IterationResults.Isc = Isc_V;
            IterationResults.Properties.VariableUnits{'Isc'} = 'A';
            app.IscnoshadingAEditField.Value = max(Isc_V);
            IterationResults.Voc = Voc_V;
            IterationResults.Properties.VariableUnits{'Voc'} = 'V';
            [~, posJmax] = max(Isc_V);
            app.VocnoshadingVEditField.Value = Voc_V(posJmax);
            IterationResults.Rs = xValues(:,1);
            IterationResults.Properties.VariableUnits{'Rs'} = 'Ohm';
            IterationResults.Rsh = xValues(:,2);
            IterationResults.Properties.VariableUnits{'Rsh'} = 'Ohm';
            IterationResults.n = xValues(:,3);
            IterationResults.Properties.VariableUnits{'n'} = '-';
            IterationResults.MaxError = MaxError;
            IterationResults.Properties.VariableUnits{'MaxError'} = 'A';
            IterationResults.TotError = TotError;
            IterationResults.Properties.VariableUnits{'TotError'} = 'A';
            IterationResults.I0 = I0;
            IterationResults.Properties.VariableUnits{'I0'} = 'A';
            if(~isempty(TempPV))
                IterationResults.TempPV = ones(height(IterationResults),1)*TempPV;
                IterationResults.Properties.VariableUnits{'TempPV'} = 'C';
            end
            
            IterationMeanResults.Rs = mean(IterationResults.Rs);
            IterationMeanResults.Properties.VariableUnits{'Rs'} = 'Ohm';
            app.RsOhmEditField.Value = mean(IterationResults.Rs);
            
            IterationMeanResults.Rsh = ...
                mean(IterationResults.Rsh(IterationResults.Rsh < 10^7)); %Ignore very large values
            IterationMeanResults.Properties.VariableUnits{'Rsh'} = 'Ohm';
            app.RshOhmEditField.Value = mean(IterationResults.Rsh(IterationResults.Rsh < 10^7));
            
            IterationMeanResults.n = mean(IterationResults.n);
            IterationMeanResults.Properties.VariableUnits{'n'} = '-';
            app.nEditField.Value = mean(IterationResults.n);
            
            IterationMeanResults.MaxError = mean(IterationResults.MaxError);
            IterationMeanResults.Properties.VariableUnits{'MaxError'} = 'A';
            
            IterationMeanResults.TotError = mean(IterationResults.TotError);
            IterationMeanResults.Properties.VariableUnits{'TotError'} = 'A';
            
            IterationMeanResults.I0 = mean(IterationResults.I0);
            IterationMeanResults.Properties.VariableUnits{'I0'} = 'A';
            app.I0AEditField.Value = mean(IterationResults.I0);
            
            %% Save to mat file
%             [filename, pathName] = uiputfile({'*.mat', 'All files (*.mat)'},...
%                 'Choose name to save data','');
%             FullFileName = [pathName filename];
%             if(filename ~= 0)
%                 save(FullFileName,'IterationResults','IterationMeanResults')
%             end
            PVparameters = IterationMeanResults;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [FF, Voc, Isc, VmP, ImP, Eff] = calc_FF(~, CV_data, ...
                column_voltage, column_current, flagPrint)
            %Calculates the fill factor for a CV curve of a PV and returns other
            %relevant values, like open circuit voltage [V], short circuit current
            %[A], voltage at maximum power point [V], current at maximum power point
            %[A] and the PV efficiency [%], for the PV efficiency 1 sun illumination is
            %assumed
            Voltage = CV_data(:,column_voltage);
            Current = CV_data(:,column_current);
            PosNeg = Current(2)/abs(Current(2)); %Detect if current positive or neg.
            Voltage = Voltage(PosNeg*Current>0); %Only values above the Voc
            Current = Current(PosNeg*Current>0); %Only values above the Voc
            
            %Power of PV during CV
            mP = abs(Voltage.*Current);
            %Current at maximum power point
            ImP = abs(Current(mP==max(mP)));
            %Voltage at maximum power point
            VmP = Voltage(mP==max(mP));
            %Open Circuit Voltage
            Minimum_Current_point = ...
                abs(Current(:))==min(abs(Current(:)));
            Voc = Voltage(Minimum_Current_point);
            %Short Circuit Current Density
            %Curves usually start at 0V
            Isc = abs(mean(Current(Voltage < 0.05)));
            %Fill factor
            FF = ImP*VmP/Isc/Voc; %decimal value
            %Efficiency at 1 sun illumination (0.1 W/cm2)
            Eff = (VmP*ImP)/(0.1)*100; %[%]
            
            if(flagPrint)
                fprintf('######################\n');
                fprintf('Fill factor: %.2f\n',FF);
                fprintf('OCV: %.2f V\n',Voc);
                fprintf('Isc: %.2f mA cm-2\n',Isc);
                fprintf('Max. Power at %.2f V and %.2f mA cm-2\n',VmP,ImP);
                fprintf('PV efficiency: %.2f%%\n',Eff);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [TimeData, VoltageData, CurrentData] = vectorData(app)
            %% Load data into vectors
            %Number of rows for each file
            vLength = zeros(size(app.RawDataIV,1));
            for iFile = size(app.RawDataIV,1)
                vLength(iFile) = length(app.RawDataIV{iFile,1}.data);
            end
            vLengthTot = sum(vLength); %Total number of rows

            TimeData = zeros(vLengthTot,1);
            TimeColumn = 1;
            for iFile = size(app.RawDataIV,1)
                if(iFile == 1)
                    TimeData(1:vLength(iFile)) = ...
                        app.RawDataIV{iFile,1}.data(:,TimeColumn);
                else
                    TimeData(vLength(iFile-1)+1:vLength(iFile)) = ...
                        app.RawDataIV{iFile,1}.data(:,TimeColumn);
                end
            end
            
            VoltageData = zeros(vLengthTot,1);
            VoltageColumn = 2;
            for iFile = size(app.RawDataIV,1)
                if(iFile == 1)
                    VoltageData(1:vLength(iFile)) = ...
                        app.RawDataIV{iFile,1}.data(:,VoltageColumn);
                else
                    VoltageData(vLength(iFile-1)+1:vLength(iFile)) = ...
                        app.RawDataIV{iFile,1}.data(:,VoltageColumn);
                end
            end

            CurrentData = zeros(vLengthTot,1);
            CurrentColumn = 3;
            for iFile = size(app.RawDataIV,1)
                if(iFile == 1)
                    CurrentData(1:vLength(iFile)) = ...
                        app.RawDataIV{iFile,1}.data(:,CurrentColumn);
                else
                    CurrentData(vLength(iFile-1)+1:vLength(iFile)) = ...
                        app.RawDataIV{iFile,1}.data(:,CurrentColumn);
                end
            end

            %% Remove unnecessary data at beginning and end
            StartTime = app.StartsEditField.Value;
            EndTime = app.EndsEditField.Value;
            VoltageData(TimeData < StartTime | TimeData > EndTime) = [];
            CurrentData(TimeData < StartTime | TimeData > EndTime) = [];
            TimeData(TimeData < StartTime | TimeData > EndTime) = [];

            %% Take only every nth value -> increase calculation speed
            %Reduce the amount of values taken into account
            %Increase speed, use only every nth value of the curve
            nthVal = app.FactortoincreasecalculationspeedEditField.Value;
            if(nthVal < 1), nthVal = 1; end
            TimeData = TimeData( ...
                1:nthVal:floor(length(TimeData)/nthVal)*nthVal);
            VoltageData = VoltageData( ...
                1:nthVal:floor(length(VoltageData)/nthVal)*nthVal);
            CurrentData = CurrentData( ...
                1:nthVal:floor(length(CurrentData)/nthVal)*nthVal);

        end
            
        function DeconvLoss = DeconvoluteConstPVParam(app, TimeData, ...
                VoltageData, CurrentData)
            %This only applies when constant PV values can be assumed
            %Calculates current losses due to PV/EC based on comparison of
            %operating point with IV curves of PV and EC components
            
            if(app.IscnoshadingAEditField.Value == 0)
                warndlg(['Please check if parameters in PV fitting' ...
                    ' tab have been entered'],'Warning');
                return;
            end

            DeconvLoss = table();            

            %% Load parameters such as characteristic PV values
            Rs = app.RsOhmEditField.Value; %Series resistance [Ohm]
            Rsh = app.RshOhmEditField.Value; %Shunt resistance [Ohm]
            n = app.nEditField.Value; %Ideality factor
            %Short circuit current at 1 sun (no shading) [A]
            Isc_initial = app.IscnoshadingAEditField.Value;
            %Open circuit voltage at 1 sun (no shading) [V]
            Voc_initial = app.VocnoshadingVEditField.Value;
            %Standard potential difference for full cell reaction
            E0 = app.E0EditField.Value; %[V]
            TempPV = app.PVtemperatureCEditField.Value + 273.15; %[K]
            
            %% Calculate shading
            %Calculate shading influence for example by vapor condensation
            %at PV front
            ParaCheck = app.ParallelCheck.Value;
            [intensity, shading] = Match_OP_To_PV_Curve(app, ...
                VoltageData, abs(CurrentData), Rs, Rsh, n, Isc_initial, ...
                Voc_initial, TempPV, ParaCheck);
            shading(intensity > 1) = 0;
            intensity(intensity > 1) = 1;
            shading = shading*100; %Decimal to %
            
            %% Calculate losses by shading/EC overpotential
            %Compare operating point current to ideal op. point 
            %(no shading)
            IdealPV_Current = app_calculate_IV_PV(VoltageData, ...
                Isc_initial, Voc_initial, 1, Rs, Rsh, n, TempPV);
            %Compare operating point current to ideal op. point
            %(no overpotential)
            IdealEC_Current = zeros(size(intensity,1),1);
            if app.ParallelCheck.Value
              parforArg = Inf; %Run regular parfor loop
            else
              parforArg = 0; %Runs essentially like a for loop
            end
            parfor (iInt = 1:size(intensity,1), parforArg)
                IdealEC_Current(iInt) = app_calculate_IV_PV(E0, ...
                    Isc_initial, Voc_initial, intensity(iInt), Rs, ...
                    Rsh, n, TempPV);
            end
            
            %Losses defined to be negative currents
            %Current loss by PV shading 
            PVLoss = -(IdealPV_Current - abs(CurrentData));
            %Current loss by EC overpotential
            ECLoss = -(IdealEC_Current - abs(CurrentData));
            %% Calculate average to limit data points in scatter plot to 
            % 500 values
            %Divide total time by 500 to get 500 time steps, then build 
            % averages
            iSteps = 500;
            %TimeInt = TimeData(2) - TimeData(1);
            TimeStep = floor(size(TimeData,1)/iSteps); %in rows
            TimeAvg = zeros(iSteps,1);
            CurrentAvg = zeros(iSteps,1);
            VoltageAvg = zeros(iSteps,1);
            PVLossAvg = zeros(iSteps,1);
            
            %Build average values for Time, Current and Voltage
            for iStep = 1:iSteps
                if(iStep*TimeStep > length(TimeData))
                    CurrentAvg(TimeAvg == 0) = [];
                    VoltageAvg(TimeAvg == 0) = [];
                    TimeAvg(TimeAvg == 0) = [];
                    PVLossAvg(PVLossAvg == 0) = [];
                    break;
                else
                    TimeAvg(iStep) = ...
                        mean(TimeData(((iStep-1)*TimeStep+1): ...
                        iStep*TimeStep,1));
                    CurrentAvg(iStep) = ...
                        mean(CurrentData(((iStep-1)*TimeStep+1): ...
                        iStep*TimeStep,1));
                    VoltageAvg(iStep) = ...
                        mean(VoltageData(((iStep-1)*TimeStep+1): ...
                        iStep*TimeStep,1));
                    PVLossAvg(iStep) = ...
                        mean(PVLoss(((iStep-1)*TimeStep+1): ...
                        iStep*TimeStep,1));
                end
            end
            
            %% Calculate potential current gain (percentage) due to PV
            % cooling
            CurrentGain = PVLoss;
            CurrentGain(CurrentGain < 0) = 0; %[A]
            vFraction = CurrentGain./abs(CurrentData)*100; %[%]
            
            CurrentGainAvg = PVLossAvg;
            CurrentGainAvg(CurrentGainAvg < 0) = 0; %[A]
            vFractionAvg = CurrentGainAvg./abs(CurrentAvg); %[%]
            
            %% Save important values to table
            DeconvLoss.Time = TimeData;
            DeconvLoss.Properties.VariableUnits{'Time'} = 's';
            DeconvLoss.Voltage = VoltageData;
            DeconvLoss.Properties.VariableUnits{'Voltage'} = 'V';
            DeconvLoss.Current = abs(CurrentData);
            DeconvLoss.Properties.VariableUnits{'Current'} = 'A';
            DeconvLoss.PVLoss = PVLoss;
            DeconvLoss.Properties.VariableUnits{'PVLoss'} = 'A';
            DeconvLoss.ECLoss = ECLoss;
            DeconvLoss.Properties.VariableUnits{'ECLoss'} = 'A';
            DeconvLoss.CurrentGain = CurrentGain; %[A]
            DeconvLoss.Properties.VariableUnits{'CurrentGain'}='A';
            DeconvLoss.CurrentFraction = vFraction; %[%]
            DeconvLoss.Properties.VariableUnits{'CurrentFraction'}='%';

            DeconvLoss.Intensity = intensity;
            DeconvLoss.Properties.VariableUnits{'Intensity'} = 'sun';
            DeconvLoss.Shading = shading;
            DeconvLoss.Properties.VariableUnits{'Shading'} = '%';
            
            %For the average values, there will be less rows -> new table
            DeconvLossAvg = table();
            DeconvLossAvg.TimeAvg = TimeAvg;
            DeconvLossAvg.Properties.VariableUnits{'TimeAvg'}='s';
            DeconvLossAvg.CurrentAvg = CurrentAvg;
            DeconvLossAvg.Properties.VariableUnits{'CurrentAvg'}='A';
            DeconvLossAvg.VoltageAvg = VoltageAvg;
            DeconvLossAvg.Properties.VariableUnits{'VoltageAvg'}='V';
            DeconvLossAvg.PVLossAvg = PVLossAvg;
            DeconvLossAvg.Properties.VariableUnits{'PVLossAvg'}='A';
            DeconvLossAvg.CurrentGainAvg = CurrentGainAvg;
            DeconvLossAvg.Properties.VariableUnits{'CurrentGainAvg'}='A';
            DeconvLossAvg.CurrentFractionAvg = vFractionAvg; %[%]
            DeconvLossAvg.Properties.VariableUnits{['Current' ...
                'FractionAvg']}='%';

            SaveDeconv(app); %Save as text file
        end

        function SaveDeconv(app)
            %% Save as *.txt file
            Delimit = app.DelimiterDropDown.Value;
            switch Delimit
                case 'Comma'
                    strDelimit = ',';
                case 'Tab'
                    strDelimit = '\t';
                case 'Space'
                    strDelimit = ' ';
            end

            [filenameDeconv, pathNameDeconv] = ...
                uiputfile({'*.txt', 'All files (*.txt)'},...
                'Choose name to save deconvoluted data');
            FullFileNameDeconv = [pathNameDeconv filenameDeconv];
            if(filenameDeconv ~= 0)
                writetable(app.DeconvLoss, FullFileNameDeconv, ...
                    'Delimiter',strDelimit);
            end
        end
        
        function [intensity, shading] = Match_OP_To_PV_Curve(~, ...
                OpVolt, OpCurrent, Rs, Rsh, n, Isc_initial, ...
                Voc_initial, TempPV, ParaCheck)
            %Calculates the light intensity (shading of PV) by water vapor
            % as it matches to operating point to a correspondent PV I-V
            % curve
            
            %% Find curve for operating point -> determine effective light
            % intensity
            intensity = zeros(size(OpVolt,1),1);
                
            options = optimset('TolFun', 1e-15,'Display','off');
            D = parallel.pool.DataQueue;
            h = waitbar(0, 'Please wait ...'); %Start waitbar
            afterEach(D, @nUpdateWaitbar); %When data sent, update waitbar
            p = 1; %Progress value in waitbar
            
            if ParaCheck
              parforArg = Inf; %Run regular parfor loop
            else
              parforArg = 0; %Runs essentially like a for loop
            end
            parfor (iPoint = 1:size(OpVolt,1), parforArg)
                %Operating Point [V],[A]
                Operating_point = [OpVolt(iPoint) OpCurrent(iPoint)];
                %Start the search at x0 sun
                x0 = Operating_point(2)/Isc_initial;
                %Calc. the diff. of a curve at x light int. and the point
                intensity(iPoint) = fsolve(@(x) ...
                    abs(app_calculate_IV_PV(Operating_point(1), ...
                    Isc_initial, Voc_initial, x, Rs, Rsh, n, TempPV) -...
                    Operating_point(2)), x0, options);
                if ~mod(iPoint, 100), send(D, iPoint); end %Update waitbar
            end
            close(h) %Close waitbar
            
            shading = 1 - intensity;
        
            function nUpdateWaitbar(~)
                waitbar(p/size(OpVolt,1), h); %Update progress
                p = p + 100; %Increase by 100 [unit?]
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function DeconvLoss = ...
                DeconvoluteECFit(app, TimeData, VoltageData, CurrentData)
            %Calculates losses based on fitting an EC curve for every
            %operating point. This fitting process assumes that the EC
            %curve can be estimated by adjusting the available area

            DeconvLoss = table();
            %% Import IV curves
            CurrentUnits = app.CurrentDropDown.Value;
            VoltageUnits = app.VoltageDropDown.Value;
            iHeaderLines = app.HeaderrowsEditField.Value;
            Delimit = app.DelimiterDropDown.Value;
            switch Delimit
                case 'Comma'
                    strDelimit = ',';
                case 'Tab'
                    strDelimit = '\t';
                case 'Space'
                    strDelimit = ' ';
            end

            %Select PV at room temperature file
            [filenamePVBef,pathNamePVBef] = uigetfile({['*.txt; *.csv;' ...
                ' *.mpt'],'All files (*.txt; *.csv; *.mpt)'}, ...
                'Choose PV IV at room temperature','', ...
                'MultiSelect','off');
            if(size(filenamePVBef,2) == 1)
                return; %If no file selected, stop program
            end
            FullfileNamePVBef = [pathNamePVBef filenamePVBef];
            PVBefore = importdata(FullfileNamePVBef, strDelimit, ...
                iHeaderLines);
            PVBefore = PVBefore.data;

            %Select PV at equilibrium temperature file
            [filenamePVAft,pathNamePVAft] = uigetfile({['*.txt; *.csv;' ...
                ' *.mpt'],'All files (*.txt; *.csv; *.mpt)'}, ...
                'Choose PV IV at equilibrium temperature','', ...
                'MultiSelect','off');
            if(size(filenamePVAft,2) == 1)
                return; %If no file selected, stop program
            end
            FullfileNamePVAft = [pathNamePVAft filenamePVAft];
            PVAfter = importdata(FullfileNamePVAft, strDelimit, ...
                iHeaderLines);
            PVAfter = PVAfter.data;

            %Select EC file
            [filenameEC,pathNameEC] = uigetfile({['*.txt; *.csv;' ...
                ' *.mpt'],'All files (*.txt; *.csv; *.mpt)'}, ...
                'Choose EC IV','', 'MultiSelect','off');
            if(size(filenameEC,2) == 1)
                return; %If no file selected, stop program
            end
            FullfileNameEC = [pathNameEC filenameEC];
            EC = importdata(FullfileNameEC, strDelimit, iHeaderLines);
            EC = EC.data;

            %% Convert to SI units
            if(strncmp(VoltageUnits,'mV',2)) %Convert mV to V
                PVBefore(:,1) = PVBefore(:,1)/1000; 
                PVAfter(:,1) = PVAfter(:,1)/1000;
                EC(:,1) = EC(:,1)/1000;
            end
            if(strncmp(CurrentUnits,'mA',2)) %Convert mA to A
                PVBefore(:,2) = PVBefore(:,2)/1000;
                PVAfter(:,2) = PVAfter(:,2)/1000;
                EC(:,2) = EC(:,2)/1000;
            end
            
            %% Load into vectors
            %Measured operating point
            OpPointV = VoltageData;
            OpPointC = abs(CurrentData);
            
            %PV curve measured before operation
            Before_V = PVBefore(:,1); %Voltage is first column
            Before_C = abs(PVBefore(:,2)); %Current is second column
            % Remove values beyond the open circuit voltage
            [~, Voc_Bef, ~, ~, ~, ~] = calc_FF(app, PVBefore, 1, ...
                2, false);
            Before_C(Before_V > Voc_Bef) = [];
            Before_V(Before_V > Voc_Bef) = [];
            %Find current at E0
            %Standard potential difference for full cell reaction
            E0 = app.E0EditField.Value; %[V]
            [~, posBE0] = min(abs(Before_V-E0));
            CurrBeforeE0 = abs(Before_C(posBE0));
            
            %PV curve measured after operation
            After_V = PVAfter(:,1); %Voltage is first column
            After_C = abs(PVAfter(:,2)); %Current is second column
            % Remove values beyond the open circuit voltage
            [~, Voc_Aft, ~, ~, ~, ~] = calc_FF(app, PVAfter, 1, ...
                2, false);
            After_C(After_V > Voc_Aft) = [];
            After_V(After_V > Voc_Aft) = [];
            %Find current at E0
            [~, posAE0] = min(abs(After_V-E0));
            CurrAfterE0 = abs(After_C(posAE0));
            
            %EC curve
            EC_V = EC(:,1); %Voltage is first column
            EC_C = EC(:,2); %Current is second column
            
            %% Calculate active area fraction 
            % (how much is blocked by bubbles)
            
            vAreaFrac = zeros(length(OpPointV),1);
            for iOP = 1:length(OpPointV)
                %Find closest match
                [~, posV] = min(abs(OpPointV(iOP) - EC_V));
                vAreaFrac(iOP) = OpPointC(iOP)/EC_C(posV);
            end
            
            %% Find intersection between EC (with active area fraction)
            % and PV curves
            
            %Int_Bef_V = zeros(length(OpPointV),1);
            Int_Bef_C = zeros(length(OpPointV),1);
            %Int_Aft_V = zeros(length(OpPointV),1);
            Int_Aft_C = zeros(length(OpPointV),1);
            
            D = parallel.pool.DataQueue;
            h = waitbar(0, 'Please wait ...');
            num_files = length(OpPointV);
            % Dummy call to nUpdateWaitbar to initialise
            nUpdateWaitbar(app, num_files, h);
            % Go back to simply calling nUpdateWaitbar with the data
            afterEach(D, @nUpdateWaitbar);
            
            if app.ParallelCheck.Value
              parforArg = Inf; %Run regular parfor loop
            else
              parforArg = 0; %Runs essentially like a for loop
            end
            parfor (iOP2 = 1:length(vAreaFrac), parforArg)
                [~, interBef] = ...
                    app_intersections( ...
                    EC_V,EC_C*vAreaFrac(iOP2),Before_V,Before_C,1);
                if(~isempty(interBef)), Int_Bef_C(iOP2) = interBef; end
                [~, interAft] = ...
                    app_intersections( ...
                    EC_V,EC_C*vAreaFrac(iOP2),After_V,After_C,1);
                if(~isempty(interAft)), Int_Aft_C(iOP2) = interAft; end
                send(D, 1);
            end
            close(h)
            
            %Total distance between intersections [A]
            TotDist = Int_Bef_C - Int_Aft_C;
            Dist_to_Bef = Int_Bef_C - OpPointC; %[A]
            PosP = Dist_to_Bef./TotDist;
            %Current gained due to cooling [A]
            CurrentGain = OpPointC - Int_Aft_C;
            CurrentGain(CurrentGain < 0) = 0;
            
            %% Calculate achieved percentage due to lower PV temperatures
            vFraction = CurrentGain./OpPointC*100; %[%]
            
            %% Calculate losses (diagonal EC curve used to calculate
            % intersection with PVs
            %Set boundaries from 0 to 1
            PosP(PosP > 1) = 1;
            PosP(PosP < 0) = 0;
            
            PVLossPos = Int_Bef_C - OpPointC; %[mA]
            %Current at E0 as a function of the position p (T affects Isc)
            CurrE0 = (1-PosP)*CurrBeforeE0 + PosP*CurrAfterE0;
            ECLossPos = CurrE0 - Int_Bef_C;
            
            %% Write to table
            DeconvLoss.Time = TimeData;
            DeconvLoss.Properties.VariableUnits{'Time'}='s';
            DeconvLoss.Voltage = OpPointV;
            DeconvLoss.Properties.VariableUnits{'Voltage'}='V';
            DeconvLoss.Current = OpPointC;
            DeconvLoss.Properties.VariableUnits{'Current'}='A';
            DeconvLoss.PVLoss = -PVLossPos; %[A]
            DeconvLoss.Properties.VariableUnits{'PVLoss'}='A';
            DeconvLoss.ECLoss = -ECLossPos; %[A]
            DeconvLoss.Properties.VariableUnits{'ECLoss'}='A';
            DeconvLoss.CurrentGain = CurrentGain; %[A]
            DeconvLoss.Properties.VariableUnits{'CurrentGain'}='A';
            DeconvLoss.CurrentFraction = vFraction; %[%]
            DeconvLoss.Properties.VariableUnits{'CurrentFraction'}='%';
            
            DeconvLoss.AreaFraction = vAreaFrac; %[%]
            DeconvLoss.Properties.VariableUnits{'AreaFraction'}='%';

            %% Save as *.txt file
            [filenameDeconv, pathNameDeconv] = ...
                uiputfile({'*.txt', 'All files (*.txt)'},...
                'Choose name to save deconvoluted data', pathNamePVBef);
            FullFileNameDeconv = [pathNameDeconv filenameDeconv];
            if(filenameDeconv ~= 0)
                writetable(DeconvLoss, FullFileNameDeconv, ...
                    'Delimiter',strDelimit);
            end

            function p = nUpdateWaitbar(~, data, h)
                %function to update waitbar during parfor loops
                persistent TOTAL COUNT H
                if nargin == 3
                    % initialisation mode
                    H = h;
                    TOTAL = data;
                    COUNT = 0;
                else
                    % afterEach call, increment COUNT
                    COUNT = 1 + COUNT;
                    p = COUNT / TOTAL;
                    waitbar(p, H);
                end
            end
        end

        function plotDeconv(app, DeconvLoss, PlotUnits)
            %Plots the deconvoluted current losses (and gains)

            %% Adjust units to user input
            TimeData = DeconvLoss.Time;
            if(strncmp(PlotUnits{1,1},'h',1))
                TimeData = TimeData/3600; %s to h
                xAxisLab = 'Time [h]';
            else
                xAxisLab = 'Time [s]';
            end

            CurrentData = DeconvLoss.Current;
            PVLossData = DeconvLoss.PVLoss;
            ECLossData = DeconvLoss.ECLoss;
            CurrentGainData = DeconvLoss.CurrentGain;
            if(strncmp(PlotUnits{3,1},'mA',2))
                CurrentData = CurrentData*1000; %A to mA
                PVLossData = PVLossData*1000; %A to mA
                ECLossData = ECLossData*1000; %A to mA
                CurrentGainData = CurrentGainData*1000; %A to mA
                yAxisLab = 'Current [mA]';
            else
                yAxisLab ='Current [A]';
            end

            %% Plot deconvoluted losses            
            scatter(app.UIAxes_2, TimeData, PVLossData)
            hold(app.UIAxes_2,'on')
            scatter(app.UIAxes_2, TimeData, ECLossData)
            legend(app.UIAxes_2, 'PV losses', 'EC losses', 'Location', ...
                'best')
            app.UIAxes_2.XLabel.String = xAxisLab;
            app.UIAxes_2.YLabel.String = yAxisLab;
            app.UIAxes_2.YLim = [-inf 0];
            hold(app.UIAxes_2,'off')
            
            %% Plot fraction of gained current due to lower PV temperatures
            yyaxis(app.UIAxes_3, 'left')
            cla(app.UIAxes_3)
            aP = area(app.UIAxes_3, TimeData, ...
                [CurrentGainData CurrentData-CurrentGainData],...
                'LineStyle','none');
            aP(1).FaceColor = [0 .45 .74];
            aP(2).FaceColor = [.85 .33 .1];
            set(app.UIAxes_3,{'ycolor'},{[0 0 0]});
            app.UIAxes_3.YLabel.String = yAxisLab;
            yyaxis(app.UIAxes_3, 'right')
            cla(app.UIAxes_3)
            hold(app.UIAxes_3,'on')
            scatter(app.UIAxes_3, TimeData, DeconvLoss.CurrentFraction,...
                'x', 'LineWidth', 1, 'SizeData', 48, 'MarkerEdgeColor', ...
                [.93 .69 .13])
            set(app.UIAxes_3,{'ycolor'},{[.93 .69 .13]});
            app.UIAxes_3.YLabel.String = 'Enabled by cooler T [%]';
            app.UIAxes_3.XLabel.String = xAxisLab;
            legend(app.UIAxes_3, 'Cooler T', 'Total', 'Location', ...
                'best')
            hold(app.UIAxes_3,'off')
        end

    end
    
    methods (Access = private, Static)

    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: 
        % ImporttimedependentcurrentvoltagedataButton
        function ImporttimedependentcurrentvoltagedataButtonPushed(app, event)
            app.RawDataIV = importIV(app);
            app.PlotUnits = {'s';'V';'A'};
            if ~isempty(app.RawDataIV)
                plotIVt(app, app.RawDataIV, app.PlotUnits);
            end
        end

        % Button pushed function: ChangetomAButton
        function ChangetomAButtonPushed(app, event)
            if ~isempty(app.RawDataIV)
                if(contains(app.ChangetomAButton.Text,'mA'))
                    app.PlotUnits{3,1} = 'mA';
                    app.ChangetomAButton.Text = 'Change to A';
                else
                    app.PlotUnits{3,1} = 'A';
                    app.ChangetomAButton.Text = 'Change to mA';
                end
                plotIVt(app, app.RawDataIV, app.PlotUnits);
            end
        end

        % Button pushed function: ChangetohButton
        function ChangetohButtonPushed(app, event)
            if ~isempty(app.RawDataIV)
                if(contains(app.ChangetohButton.Text,'Change to h'))
                    app.PlotUnits{1,1} = 'h';
                    app.ChangetohButton.Text = 'Change to s';
                else
                    app.PlotUnits{1,1} = 's';
                    app.ChangetohButton.Text = 'Change to h';
                end
                plotIVt(app, app.RawDataIV, app.PlotUnits);
            end
        end

        % Button pushed function: StartcalculationButton
        function StartcalculationButtonPushed(app, event)
            if(isempty(app.RawDataIV))
                warndlg(['Please import time-dependent current-voltage' ...
                    ' data first'], 'Warning')
                return;
            end
            [TimeData, VoltageData, CurrentData] = vectorData(app);
            %Are time-dependent EC curves available?
            if(app.TimedependentECIVcurvesavailableCheckBox.Value == 1)%Yes
                %Currently not included in app
                warndlg(['The selected conditions are not ' ...
                            'supported at the moment. Please choose' ...
                            ' different conditions.'], 'Warning')
            else%No
                %Can time-dependent EC curves be estimated?
                if(app.CanECIVcurvesbeestimatedCheckBox.Value == 1)%Yes
                    %Estimated via simple area estimation?
                    if(app.ViasimpleareaestimationCheckBox.Value == 1)%Yes
                        %DOI: 10.1039/D1EE03957A
                        app.DeconvLoss = DeconvoluteECFit(app, TimeData, ...
                            VoltageData, CurrentData);
                    else%No
                        %Currently not included in app
                        warndlg(['The selected conditions are not ' ...
                            'supported at the moment. Please choose' ...
                            ' different conditions.'], 'Warning')
                    end
                else%No
                    %Can PV parameters assumed to be stable?
                    if(app.StablePVparametersCheckBox.Value == 1)%Yes
                        %DOI: 10.1063/1.5142561
                        app.DeconvLoss = DeconvoluteConstPVParam(app, ...
                            TimeData, VoltageData, CurrentData);
                    else%No
                        %Currently not included in app
                        warndlg(['The selected conditions are not ' ...
                            'supported at the moment. Please choose' ...
                            ' different conditions.'], 'Warning')
                    end
                end
            end
            app.PlotUnits_D = {'s';'V';'A'};
            if(~isempty(app.DeconvLoss))
                plotDeconv(app, app.DeconvLoss, app.PlotUnits_D);
            end
        end

        % Button pushed function: EstimatePVparametersButton
        function EstimatePVparametersButtonPushed(app, event)
            app.PVparameters = estimate_PV_Param(app);
%             app.IscnoshadingAEditField.Value = app.PVparameters.Isc1sun;
%             app.VocnoshadingVEditField = app.PVparameters.Voc1sun;
%             app.RsOhmEditField = app.PVparameters.Rs;
%             app.RshOhmEditField = app.PVparameters.Rsh;
%             app.nEditField = app.PVparameters.n;
%             app.I0AEditField = app.PVparameters.I0;
        end

        % Button pushed function: SynccurrentvoltagedataButton
        function SynccurrentvoltagedataButtonPushed(app, event)
            %Import separate current and voltage data and export a synced
            %version with both variables
            %function to import the data
            DelimitV = app.VoltageDropDown_2.Value;
            switch DelimitV
                case 'Comma'
                    strDelimitV = ',';
                case 'Tab'
                    strDelimitV = '\t';
                case 'Space'
                    strDelimitV = ' ';
            end

            DelimitI = app.CurrentDropDown_2.Value;
            switch DelimitI
                case 'Comma'
                    strDelimitI = ',';
                case 'Tab'
                    strDelimitI = '\t';
                case 'Space'
                    strDelimitI = ' ';
            end
            TimeUnits = app.TimeDropDown.Value;
            iHeaderLines = app.HeaderrowsEditField.Value;

            %Select voltage file(s)
            [filenameV,pathNameV] = uigetfile({'*.txt; *.csv; *.mpt',...
                'All files (*.txt; *.csv; *.mpt)'},'Choose voltage data','',...
                'MultiSelect','off');
            if(size(filenameV,2) == 1) %If no file selected, stop function
                return;
            end

            filenameV = cellstr(filenameV);
            FullfileNameTxtV = cellstr(fullfile(pathNameV,filenameV));
            VData = cell(size(filenameV,2),2);

            %Load voltage file(s)
            for iSizeV = 1:size(filenameV,2)
                    Time = app.TimeEditField.Value;
                    Voltage = app.VoltageEditField_2.Value;
                    Columns = [Time Voltage];
                    VData{iSizeV,1} = ...
                        readtable(FullfileNameTxtV{1,iSizeV},...
                        'Delimiter', strDelimitV, 'NumHeaderLines',...
                        iHeaderLines);
                    VData{iSizeV,1} = VData{iSizeV,1}(:,Columns);
                    allVars = 1:width(VData{iSizeV,1});
                    VData{iSizeV,1} = renamevars(...
                        VData{iSizeV,1}, allVars, ["Time","Voltage"]);
                    VData{iSizeV,2} = filenameV{1,iSizeV};
            end
            
            for iSizeV = 1:size(filenameV,2)
                if(strncmp(TimeUnits,'h',1))
                    VData{iSizeV,1}.Time = ...
                        VData{iSizeV,1}.Time*3600; %Convert h to s
                elseif(strncmp(TimeUnits,'min',3))
                    VData{iSizeV,1}.Time = ...
                        VData{iSizeV,1}.Time*60; %Convert min to s
                end
            end

            %Select current file(s)
            [filenameI,pathNameI] = uigetfile({'*.txt; *.csv; *.mpt',...
                'All files (*.txt; *.csv; *.mpt)'},'Choose current data','',...
                'MultiSelect','off');
            if(size(filenameI,2) == 1) %If no file selected, stop function
                return;
            end

            filenameI = cellstr(filenameI);
            FullfileNameTxtI = cellstr(fullfile(pathNameI,filenameI));
            IData = cell(size(filenameI,2),2);

            %Load current file(s)
            for iSizeI = 1:size(filenameI,2)
                    Time = app.TimeEditField.Value;
                    Current = app.CurrentEditField_2.Value;
                    Columns = [Time Current];
                    IData{iSizeI,1} = ...
                        readtable(FullfileNameTxtI{1,iSizeI},...
                        'Delimiter', strDelimitI, 'NumHeaderLines', ...
                        iHeaderLines);
                    IData{iSizeI,1} = IData{iSizeI,1}(:,Columns);
                    allVars = 1:width(IData{iSizeV,1});
                    IData{iSizeV,1} = renamevars(...
                        IData{iSizeV,1}, allVars, ["Time","Current"]);
                    IData{iSizeI,2} = filenameI{1,iSizeI};
            end
            
            for iSizeI = 1:size(filenameI,2)
                if(strncmp(TimeUnits,'h',1))
                    IData{iSizeI,1}.Time = ...
                        IData{iSizeI,1}.Time*3600; %Convert h to s
                elseif(strncmp(TimeUnits,'min',3))
                    IData{iSizeI,1}.Time = ...
                        IData{iSizeI,1}.Time*60; %Convert min to s
                end
            end

            %% Adjust for different log start times
            zerotimeV = ...
                    datenum(app.VoltageEditField.Value,...
                    'yyyy-mm-dd_HH-MM-SS') * 86400; %[s]
            zerotimeI = ...
                    datenum(app.CurrentEditField.Value,...
                    'yyyy-mm-dd_HH-MM-SS') * 86400; %[s]
            timeShift = zerotimeV - zerotimeI; %[s]

            VData{1,1}.Time = VData{1,1}.Time - VData{1,1}.Time(1);
            IData{1,1}.Time = IData{1,1}.Time - IData{1,1}.Time(1) - ...
                timeShift;

            %% Match current with voltage data and create common time axis
            IVData = table();
            TimeIntervV = VData{1,1}.Time(101) - VData{1,1}.Time(100);
            TimeIntervI = IData{1,1}.Time(101) - IData{1,1}.Time(100);
            if(TimeIntervI < TimeIntervV)
                %Current is logged more frequently
                IVData.Time = IData{1,1}.Time;
                %Voltage assignment here
                vVoltage = zeros(length(IData{1,1}.Time),1);
                CurrentTime = IData{1,1}.Time;
                VoltageTime = VData{1,1}.Time;
                Voltage = VData{1,1}.Voltage;
                
                MaxMemory = memory; %See how much memory available
                %Use 90% of available memory to calculate matches
                ChunkSize = floor(MaxMemory.MaxPossibleArrayBytes...
                    *0.9/8/length(CurrentTime));
                
                %Initialize waitbar
                fig1 = uifigure;
                fWait = uiprogressdlg(fig1,'Title','Please wait',...
                    'Message', ...
                    'Estimated remaining time: calculating...', ...
                    'ShowPercentage', 'on', 'Cancelable','on');
                drawnow
                sClock = clock;

                for iTimeI = 1:floor(length(CurrentTime)/ChunkSize)
                    if fWait.CancelRequested %If user cancels
                        close(fig1)
                        return
                    end
                    
                    %Calculate
                    [~, pos] = min(abs(CurrentTime(1+(iTimeI-1)* ...
                        ChunkSize:ChunkSize*iTimeI)'-VoltageTime),[],1);
                    vVoltage(1+(iTimeI-1)*ChunkSize:ChunkSize*iTimeI) = ...
                        Voltage(pos);  

                    %Estimate completion time based on first run of loop
                    if iTimeI == 1
                         is = etime(clock,sClock);
                         esttime = is * ...
                             (floor(length(CurrentTime)/ChunkSize) + 1);
                    end

                    %Update progress of waitbar
                    iProgress = iTimeI/ ...
                        (floor(length(CurrentTime)/ChunkSize)+1);
                    if(iProgress > 1), iProgress = 1; end
                    fWait.Value = iProgress;
                    MinLeft = (esttime-etime(clock,sClock))/60;
                    if(MinLeft < 0), MinLeft = 0; end
                    fWait.Message = ['Estimated remaining time: ' ...
                        num2str(floor(MinLeft), '%02.f') ':'...
                        num2str(mod(MinLeft,1)*60, '%02.f')];
                end

                if(iTimeI > 1) %Only do this if not done in one chunk
                    %Last chunk is smaller
                    [~, pos] = min(abs(CurrentTime(1+(iTimeI)* ...
                        ChunkSize:length(CurrentTime))'-VoltageTime),[],1);
                    vVoltage(1+(iTimeI)*ChunkSize:length(CurrentTime)) =...
                        Voltage(pos);

                    % Update progress of waitbar
                    iProgress = iTimeI/ ...
                        (floor(length(CurrentTime)/ChunkSize)+1);
                    if(iProgress > 1), iProgress = 1; end
                    fWait.Value = iProgress;
                    MinLeft = (esttime-etime(clock,sClock))/60;
                    if(MinLeft < 0), MinLeft = 0; end
                    fWait.Message = ['Estimated remaining time: ' ...
                        num2str(floor(MinLeft), '%02.f') ':'...
                        num2str(mod(MinLeft,1)*60, '%02.f')];
                end
                close(fig1) %Close waitbar

                %Matrix operation is typically faster than parfor
%                 D = parallel.pool.DataQueue;
%                 h = waitbar(0, 'Please wait ...');
%                 num_files = length(CurrentTime);
%                 % Dummy call to nUpdateWaitbar to initialise
%                 nUpdateWaitbar(app, num_files, h);
%                 % Go back to simply calling nUpdateWaitbar with the data
%                 afterEach(D, @nUpdateWaitbar);
% 
%                 parfor iTimeI = 1:length(CurrentTime)
%                     %Find closest time match
%                     [~, pos] = min(abs(CurrentTime(iTimeI)-VoltageTime));
%                     vVoltage(iTimeI) = Voltage(pos); %#ok<PFBNS> 
%                     send(D, 1);
%                 end
%                 close(h)

                IVData.Voltage = vVoltage;
                IVData.Current = IData{1,1}.Current;
            else
                %Equal interval or voltage is logged more frequently
                IVData.Time = VData{1,1}.Time;
                IVData.Voltage = VData{1,1}.Voltage;
                %Current assignment here
                vCurrent = zeros(length(VData{1,1}.Time),1);
                CurrentTime = IData{1,1}.Time;
                VoltageTime = VData{1,1}.Time;
                Current = IData{1,1}.Current;

                MaxMemory = memory; %See how much memory available
                %Use 90% of available memory to calculate matches
                ChunkSize = floor(MaxMemory.MaxPossibleArrayBytes...
                    *0.9/8/length(VoltageTime));
                
                %Initialize waitbar
                fig1 = uifigure;
                fWait = uiprogressdlg(fig1,'Title','Please wait',...
                    'Message', ...
                    'Estimated remaining time: calculating...', ...
                    'ShowPercentage', 'on', 'Cancelable','on');
                drawnow
                sClock = clock;

                for iTimeV = 1:floor(length(VoltageTime)/ChunkSize)
                    if fWait.CancelRequested %If user cancels
                        close(fig1)
                        return
                    end
                    
                    %Calculate
                    [~, pos] = min(abs(VoltageTime(1+(iTimeV-1)* ...
                        ChunkSize:ChunkSize*iTimeV)'-CurrentTime),[],1);
                    vCurrent(1+(iTimeV-1)*ChunkSize:ChunkSize*iTimeV) = ...
                        Current(pos);

                    %Estimate completion time based on first run of loop
                    if iTimeV == 1
                         is = etime(clock,sClock);
                         esttime = is * ...
                             (floor(length(CurrentTime)/ChunkSize) + 1);
                    end

                    %Update progress of waitbar
                    iProgress = iTimeV/ ...
                        (floor(length(CurrentTime)/ChunkSize)+1);
                    if(iProgress > 1), iProgress = 1; end
                    fWait.Value = iProgress;
                    MinLeft = (esttime-etime(clock,sClock))/60;
                    if(MinLeft < 0), MinLeft = 0; end
                    fWait.Message = ['Estimated remaining time: ' ...
                        num2str(floor(MinLeft), '%02.f') ':'...
                        num2str(mod(MinLeft,1)*60, '%02.f')];
                end

                if(iTimeV > 1) %Only do this if not done in one chunk
                    %Last chunk is smaller
                    [~, pos] = min(abs(VoltageTime(1+(iTimeV)* ...
                        ChunkSize:length(VoltageTime))'-CurrentTime),[],1);
                    vCurrent(1+(iTimeV)*ChunkSize:length(VoltageTime)) =...
                        Current(pos);
                    
                    % Update progress of waitbar
                    iProgress = (iTimeV+1)/ ...
                        (floor(length(CurrentTime)/ChunkSize)+1);
                    if(iProgress > 1), iProgress = 1; end
                    fWait.Value = iProgress;
                    MinLeft = (esttime-etime(clock,sClock))/60;
                    if(MinLeft < 0), MinLeft = 0; end
                    fWait.Message = ['Estimated remaining time : ' ...
                        num2str(floor(MinLeft), '%02.f') ':'...
                        num2str(mod(MinLeft,1)*60, '%02.f')];
                end
                close(fig1)
                
                %Matrix operation is typically faster than parfor
%                 D = parallel.pool.DataQueue;
%                 h = waitbar(0, 'Please wait ...');
%                 num_files = length(VoltageTime);
%                 % Dummy call to nUpdateWaitbar to initialise
%                 nUpdateWaitbar(app, num_files, h);
%                 % Go back to simply calling nUpdateWaitbar with the data
%                 afterEach(D, @nUpdateWaitbar);
% 
%                 parfor iTimeV = 1:length(VoltageTime)
%                     %Find closest time match
%                     [~, pos] = min(abs(VoltageTime(iTimeV)-CurrentTime));
%                     vCurrent(iTimeV) = Current(pos); %#ok<PFBNS> 
%                     send(D, 1);
%                 end
%                 close(h)

                IVData.Current = vCurrent;
            end
            
            %% Save as *.txt file
            [filenameSync, pathNameSync] = ...
                uiputfile({'*.txt', 'All files (*.txt)'},...
                'Choose name to save synced data', pathNameV);
            FullFileNameSync = [pathNameSync filenameSync];
            DelimitSync = app.DelimiterDropDown.Value;
            switch DelimitSync
                case 'Comma'
                    strDelimitSy = ',';
                case 'Tab'
                    strDelimitSy = '\t';
                case 'Space'
                    strDelimitSy = ' ';
            end
            if(filenameSync ~= 0)
                writetable(IVData, FullFileNameSync,'Delimiter', ...
                    strDelimitSy);
            end

%             function p = nUpdateWaitbar(~, data, h)
%                 %function to update waitbar during parfor loops
%                 persistent TOTAL COUNT H
%                 if nargin == 3
%                     % initialisation mode
%                     H = h;
%                     TOTAL = data;
%                     COUNT = 0;
%                 else
%                     % afterEach call, increment COUNT
%                     COUNT = 1 + COUNT;
%                     p = COUNT / TOTAL;
%                     waitbar(p, H);
%                 end
%             end
        end

        % Button pushed function: ChangetomAButton_2
        function ChangetomAButton_2Pushed(app, event)
            if ~isempty(app.DeconvLoss)
                if(contains(app.ChangetomAButton_2.Text,'mA'))
                    app.PlotUnits_D{3,1} = 'mA';
                    app.ChangetomAButton_2.Text = 'Change to A';
                else
                    app.PlotUnits_D{3,1} = 'A';
                    app.ChangetomAButton_2.Text = 'Change to mA';
                end
                plotDeconv(app, app.DeconvLoss, app.PlotUnits_D);
            end
        end

        % Button pushed function: ChangetohButton_2
        function ChangetohButton_2Pushed(app, event)
            if ~isempty(app.DeconvLoss)
                if(contains(app.ChangetohButton_2.Text,'Change to h'))
                    app.PlotUnits_D{1,1} = 'h';
                    app.ChangetohButton_2.Text = 'Change to s';
                else
                    app.PlotUnits_D{1,1} = 's';
                    app.ChangetohButton_2.Text = 'Change to h';
                end
                plotDeconv(app, app.DeconvLoss, app.PlotUnits_D);
            end
        end

        % Button pushed function: SavecalculatedresultsButton
        function SavecalculatedresultsButtonPushed(app, event)
            if(istable(app.DeconvLoss))
                SaveDeconv(app); %Save as text file
            else
                warndlg('Please start calculation first.')
            end
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.CATSUIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 3x1 grid
                app.GridLayout.RowHeight = {662, 662, 662};
                app.GridLayout.ColumnWidth = {'1x'};
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = 1;
                app.LeftPanel.Layout.Row = 2;
                app.LeftPanel.Layout.Column = 1;
                app.RightPanel.Layout.Row = 3;
                app.RightPanel.Layout.Column = 1;
            elseif (currentFigureWidth > app.onePanelWidth && currentFigureWidth <= app.twoPanelWidth)
                % Change to a 2x2 grid
                app.GridLayout.RowHeight = {662, 662};
                app.GridLayout.ColumnWidth = {'1x', '1x'};
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = [1,2];
                app.LeftPanel.Layout.Row = 2;
                app.LeftPanel.Layout.Column = 1;
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 2;
            else
                % Change to a 1x3 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {343, '1x', 213};
                app.LeftPanel.Layout.Row = 1;
                app.LeftPanel.Layout.Column = 1;
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = 2;
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 3;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create CATSUIFigure and hide until all components are created
            app.CATSUIFigure = uifigure('Visible', 'off');
            app.CATSUIFigure.AutoResizeChildren = 'off';
            app.CATSUIFigure.Position = [100 100 1020 662];
            app.CATSUIFigure.Name = 'CATS';
            app.CATSUIFigure.Icon = fullfile(pathToMLAPP, 'CATS_resources', 'icon_48.png');
            app.CATSUIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.CATSUIFigure);
            app.GridLayout.ColumnWidth = {343, '1x', 213};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create TabGroup2
            app.TabGroup2 = uitabgroup(app.LeftPanel);
            app.TabGroup2.Position = [9 46 328 610];

            % Create GeneralTab
            app.GeneralTab = uitab(app.TabGroup2);
            app.GeneralTab.Title = 'General';

            % Create ImportsettingsPanel
            app.ImportsettingsPanel = uipanel(app.GeneralTab);
            app.ImportsettingsPanel.TitlePosition = 'centertop';
            app.ImportsettingsPanel.Title = 'Import settings';
            app.ImportsettingsPanel.FontWeight = 'bold';
            app.ImportsettingsPanel.Position = [11 454 304 126];

            % Create DelimiterDropDownLabel
            app.DelimiterDropDownLabel = uilabel(app.ImportsettingsPanel);
            app.DelimiterDropDownLabel.HorizontalAlignment = 'right';
            app.DelimiterDropDownLabel.Position = [10 66 53 22];
            app.DelimiterDropDownLabel.Text = 'Delimiter';

            % Create DelimiterDropDown
            app.DelimiterDropDown = uidropdown(app.ImportsettingsPanel);
            app.DelimiterDropDown.Items = {'Comma', 'Space', 'Tab'};
            app.DelimiterDropDown.Position = [96 66 75 22];
            app.DelimiterDropDown.Value = 'Tab';

            % Create CurrentDropDownLabel
            app.CurrentDropDownLabel = uilabel(app.ImportsettingsPanel);
            app.CurrentDropDownLabel.HorizontalAlignment = 'right';
            app.CurrentDropDownLabel.Position = [2 12 45 22];
            app.CurrentDropDownLabel.Text = 'Current';

            % Create CurrentDropDown
            app.CurrentDropDown = uidropdown(app.ImportsettingsPanel);
            app.CurrentDropDown.Items = {'A', 'mA'};
            app.CurrentDropDown.Position = [50 12 53 22];
            app.CurrentDropDown.Value = 'mA';

            % Create VoltageDropDownLabel
            app.VoltageDropDownLabel = uilabel(app.ImportsettingsPanel);
            app.VoltageDropDownLabel.HorizontalAlignment = 'right';
            app.VoltageDropDownLabel.Position = [107 12 45 22];
            app.VoltageDropDownLabel.Text = 'Voltage';

            % Create VoltageDropDown
            app.VoltageDropDown = uidropdown(app.ImportsettingsPanel);
            app.VoltageDropDown.Items = {'V', 'mV'};
            app.VoltageDropDown.Position = [155 12 50 22];
            app.VoltageDropDown.Value = 'V';

            % Create TimeDropDownLabel
            app.TimeDropDownLabel = uilabel(app.ImportsettingsPanel);
            app.TimeDropDownLabel.HorizontalAlignment = 'right';
            app.TimeDropDownLabel.Position = [211 12 31 22];
            app.TimeDropDownLabel.Text = 'Time';

            % Create TimeDropDown
            app.TimeDropDown = uidropdown(app.ImportsettingsPanel);
            app.TimeDropDown.Items = {'s', 'min', 'h'};
            app.TimeDropDown.Position = [246 12 53 22];
            app.TimeDropDown.Value = 's';

            % Create HeaderrowsEditFieldLabel
            app.HeaderrowsEditFieldLabel = uilabel(app.ImportsettingsPanel);
            app.HeaderrowsEditFieldLabel.HorizontalAlignment = 'right';
            app.HeaderrowsEditFieldLabel.Position = [183 66 74 22];
            app.HeaderrowsEditFieldLabel.Text = 'Header rows';

            % Create HeaderrowsEditField
            app.HeaderrowsEditField = uieditfield(app.ImportsettingsPanel, 'numeric');
            app.HeaderrowsEditField.Position = [264 66 34 22];
            app.HeaderrowsEditField.Value = 1;

            % Create UnitsLabel
            app.UnitsLabel = uilabel(app.ImportsettingsPanel);
            app.UnitsLabel.HorizontalAlignment = 'center';
            app.UnitsLabel.FontWeight = 'bold';
            app.UnitsLabel.Position = [137 36 35 22];
            app.UnitsLabel.Text = 'Units';

            % Create MethodoptionsPanel
            app.MethodoptionsPanel = uipanel(app.GeneralTab);
            app.MethodoptionsPanel.TitlePosition = 'centertop';
            app.MethodoptionsPanel.Title = 'Method options';
            app.MethodoptionsPanel.FontWeight = 'bold';
            app.MethodoptionsPanel.Position = [11 319 304 124];

            % Create TimedependentECIVcurvesavailableCheckBox
            app.TimedependentECIVcurvesavailableCheckBox = uicheckbox(app.MethodoptionsPanel);
            app.TimedependentECIVcurvesavailableCheckBox.Enable = 'off';
            app.TimedependentECIVcurvesavailableCheckBox.Text = {'Time-dependent EC'; 'IV curves available'};
            app.TimedependentECIVcurvesavailableCheckBox.Position = [15 38 145 58];

            % Create CanECIVcurvesbeestimatedCheckBox
            app.CanECIVcurvesbeestimatedCheckBox = uicheckbox(app.MethodoptionsPanel);
            app.CanECIVcurvesbeestimatedCheckBox.Text = {'Can EC IV curves'; 'be estimated?'};
            app.CanECIVcurvesbeestimatedCheckBox.Position = [168 53 134 30];
            app.CanECIVcurvesbeestimatedCheckBox.Value = true;

            % Create ViasimpleareaestimationCheckBox
            app.ViasimpleareaestimationCheckBox = uicheckbox(app.MethodoptionsPanel);
            app.ViasimpleareaestimationCheckBox.Enable = 'off';
            app.ViasimpleareaestimationCheckBox.Text = {'Via simple area '; 'estimation?'};
            app.ViasimpleareaestimationCheckBox.Position = [167 19 132 28];
            app.ViasimpleareaestimationCheckBox.Value = true;

            % Create StablePVparametersCheckBox
            app.StablePVparametersCheckBox = uicheckbox(app.MethodoptionsPanel);
            app.StablePVparametersCheckBox.Enable = 'off';
            app.StablePVparametersCheckBox.Text = 'Stable PV parameters?';
            app.StablePVparametersCheckBox.Position = [16 22 146 22];
            app.StablePVparametersCheckBox.Value = true;

            % Create CalculationoptionsPanel
            app.CalculationoptionsPanel = uipanel(app.GeneralTab);
            app.CalculationoptionsPanel.TitlePosition = 'centertop';
            app.CalculationoptionsPanel.Title = 'Calculation options';
            app.CalculationoptionsPanel.FontWeight = 'bold';
            app.CalculationoptionsPanel.Position = [12 126 303 181];

            % Create StartandendtimeforcurrentvoltageanalysisLabel
            app.StartandendtimeforcurrentvoltageanalysisLabel = uilabel(app.CalculationoptionsPanel);
            app.StartandendtimeforcurrentvoltageanalysisLabel.Position = [30 133 253 22];
            app.StartandendtimeforcurrentvoltageanalysisLabel.Text = 'Start and end time for current-voltage analysis';

            % Create StartsEditFieldLabel
            app.StartsEditFieldLabel = uilabel(app.CalculationoptionsPanel);
            app.StartsEditFieldLabel.HorizontalAlignment = 'right';
            app.StartsEditFieldLabel.Position = [6 100 47 22];
            app.StartsEditFieldLabel.Text = 'Start [s]';

            % Create StartsEditField
            app.StartsEditField = uieditfield(app.CalculationoptionsPanel, 'numeric');
            app.StartsEditField.ValueDisplayFormat = '%11.1f';
            app.StartsEditField.Position = [55 100 91 22];

            % Create EndsEditFieldLabel
            app.EndsEditFieldLabel = uilabel(app.CalculationoptionsPanel);
            app.EndsEditFieldLabel.HorizontalAlignment = 'right';
            app.EndsEditFieldLabel.Position = [154 99 43 22];
            app.EndsEditFieldLabel.Text = 'End [s]';

            % Create EndsEditField
            app.EndsEditField = uieditfield(app.CalculationoptionsPanel, 'numeric');
            app.EndsEditField.ValueDisplayFormat = '%11.1f';
            app.EndsEditField.Position = [199 99 91 22];
            app.EndsEditField.Value = Inf;

            % Create MinimumpotentialdifferenceforfullcellreactionofinterestLabel
            app.MinimumpotentialdifferenceforfullcellreactionofinterestLabel = uilabel(app.CalculationoptionsPanel);
            app.MinimumpotentialdifferenceforfullcellreactionofinterestLabel.HorizontalAlignment = 'right';
            app.MinimumpotentialdifferenceforfullcellreactionofinterestLabel.Position = [16 11 176 28];
            app.MinimumpotentialdifferenceforfullcellreactionofinterestLabel.Text = {'Standard potential difference for'; 'full cell reaction of interest [V]'};

            % Create E0EditField
            app.E0EditField = uieditfield(app.CalculationoptionsPanel, 'numeric');
            app.E0EditField.Position = [201 14 37 22];
            app.E0EditField.Value = 1.33;

            % Create FactortoincreasecalculationspeedEditFieldLabel
            app.FactortoincreasecalculationspeedEditFieldLabel = uilabel(app.CalculationoptionsPanel);
            app.FactortoincreasecalculationspeedEditFieldLabel.HorizontalAlignment = 'right';
            app.FactortoincreasecalculationspeedEditFieldLabel.Position = [91 51 102 28];
            app.FactortoincreasecalculationspeedEditFieldLabel.Text = {'Factor to increase'; 'calculation speed'};

            % Create FactortoincreasecalculationspeedEditField
            app.FactortoincreasecalculationspeedEditField = uieditfield(app.CalculationoptionsPanel, 'numeric');
            app.FactortoincreasecalculationspeedEditField.Limits = [1 10000];
            app.FactortoincreasecalculationspeedEditField.RoundFractionalValues = 'on';
            app.FactortoincreasecalculationspeedEditField.ValueDisplayFormat = '%.0f';
            app.FactortoincreasecalculationspeedEditField.Position = [201 54 37 22];
            app.FactortoincreasecalculationspeedEditField.Value = 10;

            % Create ParallelCheck
            app.ParallelCheck = uicheckbox(app.GeneralTab);
            app.ParallelCheck.Text = 'Enable parallel computing';
            app.ParallelCheck.Position = [13 95 160 22];
            app.ParallelCheck.Value = true;

            % Create PVfittingTab_2
            app.PVfittingTab_2 = uitab(app.TabGroup2);
            app.PVfittingTab_2.Title = 'PV fitting';

            % Create CalculatedparametersPanel
            app.CalculatedparametersPanel = uipanel(app.PVfittingTab_2);
            app.CalculatedparametersPanel.TitlePosition = 'centertop';
            app.CalculatedparametersPanel.Title = 'Calculated parameters';
            app.CalculatedparametersPanel.FontWeight = 'bold';
            app.CalculatedparametersPanel.Position = [3 209 321 221];

            % Create IscnoshadingAEditFieldLabel
            app.IscnoshadingAEditFieldLabel = uilabel(app.CalculatedparametersPanel);
            app.IscnoshadingAEditFieldLabel.HorizontalAlignment = 'right';
            app.IscnoshadingAEditFieldLabel.Position = [26 165 109 22];
            app.IscnoshadingAEditFieldLabel.Text = 'Isc (no shading) [A]';

            % Create IscnoshadingAEditField
            app.IscnoshadingAEditField = uieditfield(app.CalculatedparametersPanel, 'numeric');
            app.IscnoshadingAEditField.Position = [143 165 91 22];

            % Create VocnoshadingVEditFieldLabel
            app.VocnoshadingVEditFieldLabel = uilabel(app.CalculatedparametersPanel);
            app.VocnoshadingVEditFieldLabel.HorizontalAlignment = 'right';
            app.VocnoshadingVEditFieldLabel.Position = [21 134 114 22];
            app.VocnoshadingVEditFieldLabel.Text = 'Voc (no shading) [V]';

            % Create VocnoshadingVEditField
            app.VocnoshadingVEditField = uieditfield(app.CalculatedparametersPanel, 'numeric');
            app.VocnoshadingVEditField.Position = [143 134 91 22];

            % Create RsOhmEditFieldLabel
            app.RsOhmEditFieldLabel = uilabel(app.CalculatedparametersPanel);
            app.RsOhmEditFieldLabel.HorizontalAlignment = 'right';
            app.RsOhmEditFieldLabel.Position = [80 104 56 22];
            app.RsOhmEditFieldLabel.Text = 'Rs [Ohm]';

            % Create RsOhmEditField
            app.RsOhmEditField = uieditfield(app.CalculatedparametersPanel, 'numeric');
            app.RsOhmEditField.Position = [143 104 91 22];

            % Create RshOhmLabel
            app.RshOhmLabel = uilabel(app.CalculatedparametersPanel);
            app.RshOhmLabel.HorizontalAlignment = 'right';
            app.RshOhmLabel.Position = [66 73 72 22];
            app.RshOhmLabel.Text = 'Rsh [Ohm]';

            % Create RshOhmEditField
            app.RshOhmEditField = uieditfield(app.CalculatedparametersPanel, 'numeric');
            app.RshOhmEditField.Position = [143 73 91 22];

            % Create nEditFieldLabel
            app.nEditFieldLabel = uilabel(app.CalculatedparametersPanel);
            app.nEditFieldLabel.HorizontalAlignment = 'right';
            app.nEditFieldLabel.Position = [110 43 25 22];
            app.nEditFieldLabel.Text = 'n';

            % Create nEditField
            app.nEditField = uieditfield(app.CalculatedparametersPanel, 'numeric');
            app.nEditField.Position = [143 43 91 22];

            % Create I0AEditFieldLabel
            app.I0AEditFieldLabel = uilabel(app.CalculatedparametersPanel);
            app.I0AEditFieldLabel.HorizontalAlignment = 'right';
            app.I0AEditFieldLabel.Position = [103 13 34 22];
            app.I0AEditFieldLabel.Text = 'I0 [A]';

            % Create I0AEditField
            app.I0AEditField = uieditfield(app.CalculatedparametersPanel, 'numeric');
            app.I0AEditField.Position = [143 13 91 22];

            % Create ValuesforcalculationPanel
            app.ValuesforcalculationPanel = uipanel(app.PVfittingTab_2);
            app.ValuesforcalculationPanel.TitlePosition = 'centertop';
            app.ValuesforcalculationPanel.Title = 'Values for calculation';
            app.ValuesforcalculationPanel.FontWeight = 'bold';
            app.ValuesforcalculationPanel.Position = [3 442 321 132];

            % Create PVsizecm2EditFieldLabel
            app.PVsizecm2EditFieldLabel = uilabel(app.ValuesforcalculationPanel);
            app.PVsizecm2EditFieldLabel.HorizontalAlignment = 'right';
            app.PVsizecm2EditFieldLabel.Enable = 'off';
            app.PVsizecm2EditFieldLabel.Position = [22 15 79 22];
            app.PVsizecm2EditFieldLabel.Text = 'PV size [cm2]';

            % Create PVsizecm2EditField
            app.PVsizecm2EditField = uieditfield(app.ValuesforcalculationPanel, 'numeric');
            app.PVsizecm2EditField.Enable = 'off';
            app.PVsizecm2EditField.Position = [108 15 28 22];
            app.PVsizecm2EditField.Value = 1;

            % Create PVtemperatureCLabel
            app.PVtemperatureCLabel = uilabel(app.ValuesforcalculationPanel);
            app.PVtemperatureCLabel.HorizontalAlignment = 'right';
            app.PVtemperatureCLabel.Position = [145 15 113 22];
            app.PVtemperatureCLabel.Text = 'PV temperature [°C]';

            % Create PVtemperatureCEditField
            app.PVtemperatureCEditField = uieditfield(app.ValuesforcalculationPanel, 'numeric');
            app.PVtemperatureCEditField.Position = [264 15 45 22];
            app.PVtemperatureCEditField.Value = 25;

            % Create FactortoincreasespeedforPVparameterestimationEditFieldLabel
            app.FactortoincreasespeedforPVparameterestimationEditFieldLabel = uilabel(app.ValuesforcalculationPanel);
            app.FactortoincreasespeedforPVparameterestimationEditFieldLabel.HorizontalAlignment = 'right';
            app.FactortoincreasespeedforPVparameterestimationEditFieldLabel.Position = [7 46 97 56];
            app.FactortoincreasespeedforPVparameterestimationEditFieldLabel.Text = {'Factor to'; 'increase speed'; 'for PV parameter'; 'estimation'};

            % Create SpeedFactor
            app.SpeedFactor = uieditfield(app.ValuesforcalculationPanel, 'numeric');
            app.SpeedFactor.Limits = [1 Inf];
            app.SpeedFactor.RoundFractionalValues = 'on';
            app.SpeedFactor.ValueDisplayFormat = '%.0f';
            app.SpeedFactor.Position = [108 43 28 59];
            app.SpeedFactor.Value = 1;

            % Create MinimumvoltageforPVparameterestimationVLabel
            app.MinimumvoltageforPVparameterestimationVLabel = uilabel(app.ValuesforcalculationPanel);
            app.MinimumvoltageforPVparameterestimationVLabel.HorizontalAlignment = 'right';
            app.MinimumvoltageforPVparameterestimationVLabel.Position = [157 53 97 42];
            app.MinimumvoltageforPVparameterestimationVLabel.Text = {'Minimum voltage'; 'for PV parameter'; 'estimation [V]'};

            % Create MinimumvoltagePVParameter
            app.MinimumvoltagePVParameter = uieditfield(app.ValuesforcalculationPanel, 'numeric');
            app.MinimumvoltagePVParameter.Position = [264 53 45 42];
            app.MinimumvoltagePVParameter.Value = 1.5;

            % Create IVsyncTab
            app.IVsyncTab = uitab(app.TabGroup2);
            app.IVsyncTab.Title = 'IV sync';

            % Create VoltageEditFieldLabel
            app.VoltageEditFieldLabel = uilabel(app.IVsyncTab);
            app.VoltageEditFieldLabel.HorizontalAlignment = 'right';
            app.VoltageEditFieldLabel.Position = [41 519 50 22];
            app.VoltageEditFieldLabel.Text = 'Voltage';

            % Create VoltageEditField
            app.VoltageEditField = uieditfield(app.IVsyncTab, 'text');
            app.VoltageEditField.Position = [104 519 156 22];
            app.VoltageEditField.Value = '2021-06-10_17-22-25';

            % Create CurrentEditFieldLabel
            app.CurrentEditFieldLabel = uilabel(app.IVsyncTab);
            app.CurrentEditFieldLabel.HorizontalAlignment = 'right';
            app.CurrentEditFieldLabel.Position = [42 492 48 22];
            app.CurrentEditFieldLabel.Text = 'Current';

            % Create CurrentEditField
            app.CurrentEditField = uieditfield(app.IVsyncTab, 'text');
            app.CurrentEditField.Position = [104 492 156 22];
            app.CurrentEditField.Value = '2021-06-10_17-22-24.456';

            % Create Starttimesformatyyyymmdd_HHMMSSLabel
            app.Starttimesformatyyyymmdd_HHMMSSLabel = uilabel(app.IVsyncTab);
            app.Starttimesformatyyyymmdd_HHMMSSLabel.HorizontalAlignment = 'center';
            app.Starttimesformatyyyymmdd_HHMMSSLabel.FontWeight = 'bold';
            app.Starttimesformatyyyymmdd_HHMMSSLabel.Position = [22 547 265 22];
            app.Starttimesformatyyyymmdd_HHMMSSLabel.Text = 'Start times (format: yyyy-mm-dd_HH-MM-SS)';

            % Create VoltageEditField_2Label
            app.VoltageEditField_2Label = uilabel(app.IVsyncTab);
            app.VoltageEditField_2Label.HorizontalAlignment = 'right';
            app.VoltageEditField_2Label.Position = [112 365 45 22];
            app.VoltageEditField_2Label.Text = 'Voltage';

            % Create VoltageEditField_2
            app.VoltageEditField_2 = uieditfield(app.IVsyncTab, 'numeric');
            app.VoltageEditField_2.Position = [165 365 41 22];
            app.VoltageEditField_2.Value = 4;

            % Create CurrentEditField_2Label
            app.CurrentEditField_2Label = uilabel(app.IVsyncTab);
            app.CurrentEditField_2Label.HorizontalAlignment = 'right';
            app.CurrentEditField_2Label.Position = [214 365 50 22];
            app.CurrentEditField_2Label.Text = 'Current';

            % Create CurrentEditField_2
            app.CurrentEditField_2 = uieditfield(app.IVsyncTab, 'numeric');
            app.CurrentEditField_2.Position = [271 365 40 22];
            app.CurrentEditField_2.Value = 3;

            % Create TimeEditFieldLabel
            app.TimeEditFieldLabel = uilabel(app.IVsyncTab);
            app.TimeEditFieldLabel.HorizontalAlignment = 'right';
            app.TimeEditFieldLabel.Position = [10 365 31 22];
            app.TimeEditFieldLabel.Text = 'Time';

            % Create TimeEditField
            app.TimeEditField = uieditfield(app.IVsyncTab, 'numeric');
            app.TimeEditField.Position = [47 365 41 22];
            app.TimeEditField.Value = 1;

            % Create CurrentDropDown_2Label
            app.CurrentDropDown_2Label = uilabel(app.IVsyncTab);
            app.CurrentDropDown_2Label.HorizontalAlignment = 'right';
            app.CurrentDropDown_2Label.Position = [169 434 41 22];
            app.CurrentDropDown_2Label.Text = 'Current';

            % Create CurrentDropDown_2
            app.CurrentDropDown_2 = uidropdown(app.IVsyncTab);
            app.CurrentDropDown_2.Items = {'Comma', 'Space', 'Tab'};
            app.CurrentDropDown_2.Position = [214 434 100 22];
            app.CurrentDropDown_2.Value = 'Tab';

            % Create VoltageDropDown_2Label
            app.VoltageDropDown_2Label = uilabel(app.IVsyncTab);
            app.VoltageDropDown_2Label.HorizontalAlignment = 'right';
            app.VoltageDropDown_2Label.Position = [1 434 55 22];
            app.VoltageDropDown_2Label.Text = 'Voltage';

            % Create VoltageDropDown_2
            app.VoltageDropDown_2 = uidropdown(app.IVsyncTab);
            app.VoltageDropDown_2.Items = {'Comma', 'Space', 'Tab'};
            app.VoltageDropDown_2.Position = [60 434 100 22];
            app.VoltageDropDown_2.Value = 'Comma';

            % Create DelimiterLabel
            app.DelimiterLabel = uilabel(app.IVsyncTab);
            app.DelimiterLabel.HorizontalAlignment = 'center';
            app.DelimiterLabel.FontWeight = 'bold';
            app.DelimiterLabel.Position = [135 461 56 22];
            app.DelimiterLabel.Text = 'Delimiter';

            % Create ColumnnumberLabel
            app.ColumnnumberLabel = uilabel(app.IVsyncTab);
            app.ColumnnumberLabel.HorizontalAlignment = 'center';
            app.ColumnnumberLabel.FontWeight = 'bold';
            app.ColumnnumberLabel.Position = [115 395 97 22];
            app.ColumnnumberLabel.Text = 'Column number';

            % Create v01Label
            app.v01Label = uilabel(app.LeftPanel);
            app.v01Label.Position = [11 6 28 22];
            app.v01Label.Text = 'v0.1';

            % Create CenterPanel
            app.CenterPanel = uipanel(app.GridLayout);
            app.CenterPanel.Layout.Row = 1;
            app.CenterPanel.Layout.Column = 2;

            % Create TabGroup
            app.TabGroup = uitabgroup(app.CenterPanel);
            app.TabGroup.Position = [12 6 442 650];

            % Create ImporteddataTab
            app.ImporteddataTab = uitab(app.TabGroup);
            app.ImporteddataTab.Title = 'Imported data';

            % Create UIAxes
            app.UIAxes = uiaxes(app.ImporteddataTab);
            title(app.UIAxes, 'Imported data')
            xlabel(app.UIAxes, 'Time [s]')
            ylabel(app.UIAxes, 'Voltage [V]')
            app.UIAxes.PlotBoxAspectRatio = [2.17910447761194 1 1];
            app.UIAxes.XTickLabelRotation = 0;
            app.UIAxes.YTickLabelRotation = 0;
            app.UIAxes.ZTickLabelRotation = 0;
            app.UIAxes.Box = 'on';
            app.UIAxes.Position = [5 303 428 291];

            % Create ChangetomAButton
            app.ChangetomAButton = uibutton(app.ImporteddataTab, 'push');
            app.ChangetomAButton.ButtonPushedFcn = createCallbackFcn(app, @ChangetomAButtonPushed, true);
            app.ChangetomAButton.FontWeight = 'bold';
            app.ChangetomAButton.Position = [98 265 100 23];
            app.ChangetomAButton.Text = 'Change to mA';

            % Create ChangetohButton
            app.ChangetohButton = uibutton(app.ImporteddataTab, 'push');
            app.ChangetohButton.ButtonPushedFcn = createCallbackFcn(app, @ChangetohButtonPushed, true);
            app.ChangetohButton.FontWeight = 'bold';
            app.ChangetohButton.Position = [247 264 100 23];
            app.ChangetohButton.Text = 'Change to h';

            % Create PVfittingTab
            app.PVfittingTab = uitab(app.TabGroup);
            app.PVfittingTab.Title = 'PV fitting';

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.PVfittingTab);
            title(app.UIAxes2, 'PV fitting individual')
            xlabel(app.UIAxes2, 'Voltage [V]')
            ylabel(app.UIAxes2, 'Current [mA]')
            app.UIAxes2.XTickLabelRotation = 0;
            app.UIAxes2.YTickLabelRotation = 0;
            app.UIAxes2.ZTickLabelRotation = 0;
            app.UIAxes2.Box = 'on';
            app.UIAxes2.Position = [1 332 432 282];

            % Create UIAxes2_2
            app.UIAxes2_2 = uiaxes(app.PVfittingTab);
            title(app.UIAxes2_2, 'PV fitting with averaged parameters')
            xlabel(app.UIAxes2_2, 'Voltage [V]')
            ylabel(app.UIAxes2_2, 'Current [mA]')
            app.UIAxes2_2.XTickLabelRotation = 0;
            app.UIAxes2_2.YTickLabelRotation = 0;
            app.UIAxes2_2.ZTickLabelRotation = 0;
            app.UIAxes2_2.Box = 'on';
            app.UIAxes2_2.Position = [5 43 428 281];

            % Create CalculatedresultsTab
            app.CalculatedresultsTab = uitab(app.TabGroup);
            app.CalculatedresultsTab.Title = 'Calculated results';

            % Create UIAxes_2
            app.UIAxes_2 = uiaxes(app.CalculatedresultsTab);
            title(app.UIAxes_2, 'Calculated losses')
            xlabel(app.UIAxes_2, 'Time [s]')
            ylabel(app.UIAxes_2, 'Current [A]')
            app.UIAxes_2.PlotBoxAspectRatio = [2.17910447761194 1 1];
            app.UIAxes_2.XTickLabelRotation = 0;
            app.UIAxes_2.YTickLabelRotation = 0;
            app.UIAxes_2.ZTickLabelRotation = 0;
            app.UIAxes_2.Box = 'on';
            app.UIAxes_2.Position = [5 335 428 291];

            % Create UIAxes_3
            app.UIAxes_3 = uiaxes(app.CalculatedresultsTab);
            title(app.UIAxes_3, 'Calculated gains')
            xlabel(app.UIAxes_3, 'Time [s]')
            ylabel(app.UIAxes_3, 'Current [A]')
            app.UIAxes_3.PlotBoxAspectRatio = [2.17910447761194 1 1];
            app.UIAxes_3.XTickLabelRotation = 0;
            app.UIAxes_3.YTickLabelRotation = 0;
            app.UIAxes_3.ZTickLabelRotation = 0;
            app.UIAxes_3.Box = 'on';
            app.UIAxes_3.Position = [5 40 428 290];

            % Create ChangetomAButton_2
            app.ChangetomAButton_2 = uibutton(app.CalculatedresultsTab, 'push');
            app.ChangetomAButton_2.ButtonPushedFcn = createCallbackFcn(app, @ChangetomAButton_2Pushed, true);
            app.ChangetomAButton_2.FontWeight = 'bold';
            app.ChangetomAButton_2.Position = [119 10 100 23];
            app.ChangetomAButton_2.Text = 'Change to mA';

            % Create ChangetohButton_2
            app.ChangetohButton_2 = uibutton(app.CalculatedresultsTab, 'push');
            app.ChangetohButton_2.ButtonPushedFcn = createCallbackFcn(app, @ChangetohButton_2Pushed, true);
            app.ChangetohButton_2.FontWeight = 'bold';
            app.ChangetohButton_2.Position = [246 10 100 23];
            app.ChangetohButton_2.Text = 'Change to h';

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 3;

            % Create EstimatePVparametersButton
            app.EstimatePVparametersButton = uibutton(app.RightPanel, 'push');
            app.EstimatePVparametersButton.ButtonPushedFcn = createCallbackFcn(app, @EstimatePVparametersButtonPushed, true);
            app.EstimatePVparametersButton.FontWeight = 'bold';
            app.EstimatePVparametersButton.Position = [18 461 178 50];
            app.EstimatePVparametersButton.Text = 'Estimate PV parameters';

            % Create StartcalculationButton
            app.StartcalculationButton = uibutton(app.RightPanel, 'push');
            app.StartcalculationButton.ButtonPushedFcn = createCallbackFcn(app, @StartcalculationButtonPushed, true);
            app.StartcalculationButton.FontWeight = 'bold';
            app.StartcalculationButton.Position = [18 403 178 50];
            app.StartcalculationButton.Text = 'Start calculation';

            % Create ImporttimedependentcurrentvoltagedataButton
            app.ImporttimedependentcurrentvoltagedataButton = uibutton(app.RightPanel, 'push');
            app.ImporttimedependentcurrentvoltagedataButton.ButtonPushedFcn = createCallbackFcn(app, @ImporttimedependentcurrentvoltagedataButtonPushed, true);
            app.ImporttimedependentcurrentvoltagedataButton.FontWeight = 'bold';
            app.ImporttimedependentcurrentvoltagedataButton.Position = [18 518 178 60];
            app.ImporttimedependentcurrentvoltagedataButton.Text = {'Import time-dependent'; 'current-voltage data'};

            % Create SynccurrentvoltagedataButton
            app.SynccurrentvoltagedataButton = uibutton(app.RightPanel, 'push');
            app.SynccurrentvoltagedataButton.ButtonPushedFcn = createCallbackFcn(app, @SynccurrentvoltagedataButtonPushed, true);
            app.SynccurrentvoltagedataButton.FontWeight = 'bold';
            app.SynccurrentvoltagedataButton.Position = [18 586 178 55];
            app.SynccurrentvoltagedataButton.Text = 'Sync current-voltage data';

            % Create SavecalculatedresultsButton
            app.SavecalculatedresultsButton = uibutton(app.RightPanel, 'push');
            app.SavecalculatedresultsButton.ButtonPushedFcn = createCallbackFcn(app, @SavecalculatedresultsButtonPushed, true);
            app.SavecalculatedresultsButton.FontWeight = 'bold';
            app.SavecalculatedresultsButton.Position = [18 344 178 50];
            app.SavecalculatedresultsButton.Text = 'Save calculated results';

            % Show the figure after all components are created
            app.CATSUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = CATS_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.CATSUIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.CATSUIFigure)
        end
    end
end