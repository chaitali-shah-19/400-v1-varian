classdef ExperimentFunctions_SINGLE_FID_20240909 < handle
    
    properties
        
        interfaceSigGen
        interfaceSigGen2 %for the second SG386 source
        interfaceBNC
        RawData
        AvgData
        Std_dev
        delay_init_and_end_seq
        statuswbexp =0; %will be used as a handle to close the waitbar in case of stop scan
        wf
        com
        mw2
        mw
        arb1
        arb2
        arb3
        arb4
        arb5
        arb6
        arb7
        arb8
        arb9
        arb10
        arb11
        arb12
        arb_mag
        pw1
        pw2
        pw3
        pw4
        pw5
        pw6
        pw7
        pw8
        laser
        gpib_LockInAmp
        att
        ramp
        psu
        psu2
        tabor
        srs
        laserdome
        swabian
        zaber
    end
    
    methods
        
        %%%%% Run experiment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function RunExperiment(obj,handles)
            import PulseStreamer.*;
            warning('off','all');
            handles.Imaginghandles.ImagingFunctions.interfaceDataAcq.hasAborted = 0;    
            err = 0;

            %% Test Sage write
            Sage_write('1');
            obj.com.StopBuffer(1);
            
            helper_scan=1;
            for aux=1:1:handles.ExperimentalScan.vary_points(1)  %loop over first scan dimension
                %                         pause(3);
                disp(['Running experiment number:' num2str(aux)])
                if ~handles.Imaginghandles.ImagingFunctions.interfaceDataAcq.hasAborted                    
                    %% SETUP AWG + SRS
                    if length(handles.Array_PSeq{1}.Channels)>=39
                        if handles.Array_PSeq{1}.Channels(39).Enable && ~handles.Imaginghandles.ImagingFunctions.interfaceDataAcq.hasAborted
                            
                            %AWG parameters
                            awg_center_freq=handles.Array_PSeq{helper_scan,aux}.Channels(39).Frequency;
                            awg_bw_freq=handles.Array_PSeq{helper_scan,aux}.Channels(39).Phase;
                            awg_amp=handles.Array_PSeq{helper_scan,aux}.Channels(39).Amplitude;
                            sweep_freq=handles.Array_PSeq{helper_scan,aux}.Channels(39).PhaseQ;
                            sweep_sigma=handles.Array_PSeq{helper_scan,aux}.Channels(39).FreqIQ;
                            symm=handles.Array_PSeq{helper_scan,aux}.Channels(39).Phasemod;
                            %SRS parameters
                            srs_freq=handles.Array_PSeq{helper_scan,aux}.Channels(39).FreqmodI;
                            srs_amp=handles.Array_PSeq{helper_scan,aux}.Channels(39).FreqmodQ;
                            
                            Sage_write(['2,',num2str(awg_center_freq),',',num2str(awg_bw_freq),',',num2str(awg_amp),',',...
                                num2str(sweep_freq),',',num2str(sweep_sigma),',',num2str(symm),',',num2str(srs_freq),',',num2str(srs_amp)]);
                            disp('Wrote Sage write 2: Make MW waveform');
                            pause(1.5);
                        end
                    end
                    %% SETUP NMR pulse parameters
                    
                    if length(handles.Array_PSeq{1}.Channels)>=38
                        if handles.Array_PSeq{1}.Channels(38).Enable && ~handles.Imaginghandles.ImagingFunctions.interfaceDataAcq.hasAborted
                            if aux<handles.ExperimentalScan.vary_points(1)
                                pw=handles.Array_PSeq{helper_scan,aux+1}.Channels(38).Frequency;
                                PW2=handles.Array_PSeq{helper_scan,aux+1}.Channels(38).FreqmodQ;
                                tacq=handles.Array_PSeq{helper_scan,aux+1}.Channels(38).Phase;
                                sequence_type=handles.Array_PSeq{helper_scan,aux+1}.Channels(38).PhaseQ;
                                gain=handles.Array_PSeq{helper_scan,aux+1}.Channels(38).FreqIQ;
                                tof=handles.Array_PSeq{helper_scan,aux+1}.Channels(38).AmpIQ;
                            else
                                pw=42;
                                tacq=24;
                                gain=26;
                                PW2= 42;
                                tof=5.1; %3.7
                                sequence_type=2;
                            end
                            
                            if sequence_type == 2  %spin lock sequence type, can uncomment others later
                                pw_current=handles.Array_PSeq{helper_scan,aux}.Channels(38).Frequency;
                                PW2_current=handles.Array_PSeq{helper_scan,aux}.Channels(38).FreqmodQ;
                                tacq_current=handles.Array_PSeq{helper_scan,aux}.Channels(38).Phase;
                                gain_current=handles.Array_PSeq{helper_scan,aux}.Channels(38).FreqIQ;
                                tof_current=handles.Array_PSeq{helper_scan,aux}.Channels(38).AmpIQ;
                                [nc,Tmax]= make_pw_nmr_8_SINGLEFID20240909(pw_current,tacq_current,gain_current,PW2_current,tof_current);
                                pause(0.1);
                                %for next cycle
                                make_pw_nmr_8_SINGLEFID20240909(pw,tacq,gain,PW2,tof);
                                Tmax=Tmax*14*1e-3; %Tmax * number of pulses in sequence pulsed_spinlock
                                nc=nc*14; %nc * number of aquire windows in pulsed_spinlock = total number of pulses
                                
                                %Insted of sage_write_2, we do following
                                pw = pw_current*1e-6;
                                PW2 = PW2_current*1e-6;
                                numberOfPulses_total = nc;
                                tacq = tacq_current;
                                numberOfPulses= floor(numberOfPulses_total/Tmax); %in 1 second %will be 1 for FID
                                loops=Tmax;
                                
                                Sage_write(['3,',num2str(pw),',',num2str(PW2),',',num2str(numberOfPulses_total),',',num2str(tacq),',',...
                                    num2str(numberOfPulses),',',num2str(loops),',',num2str(Tmax)]);
                                disp('Wrote Sage write 3: Setting up receiver');
                            else
                                disp('Sequence type not accepted\n')
                            end
                            pause(0.2);
                        end
                    end
                    
                    %%EXTRA PAUSE
                    pause(0.5);
                    
                    if handles.Imaginghandles.ImagingFunctions.interfaceDataAcq.hasAborted==1
                        disp('You have aborted!');
                        obj.com.StopBuffer(1); obj.com.StopBuffer(3); obj.com.StopBuffer(7); obj.com.StopBuffer(6); obj.com.StopBuffer(9);
                        obj.com.KillAll();
                        handles.ExperimentFunctions.swabian.reset(); 
                        break;
                    else
                        
                        %% Initialize actuator position for very first motion
                        if handles.Array_PSeq{1}.Channels(12).Enable && aux==1 && ~isempty(handles.Array_PSeq{helper_scan,aux}.Channels(12).Frequency)
                            start_pos=obj.com.GetFPosition(obj.com.ACSC_AXIS_0);
                            Position=handles.Array_PSeq{helper_scan,aux}.Channels(12).Frequency;
                            obj.com.SetVelocity(obj.com.ACSC_AXIS_0,5); %If the shuttler isn't in the right position at the start of the experiment, it will pull it back to the right position with this velocity. Note, previously set at 50 *upside down face emoji*
                            obj.com.SetAcceleration(obj.com.ACSC_AXIS_0,50);
                            obj.com.SetDeceleration(obj.com.ACSC_AXIS_0,50);
                            obj.com.SetJerk(obj.com.ACSC_AXIS_0,500);
                            
                            obj.com.ToPoint(obj.com.ACSC_AMF_WAIT, obj.com.ACSC_AXIS_0, -Position); %MOVE
                            obj.com.Go(obj.com.ACSC_AXIS_0);
                            displacement=abs(abs(start_pos)-abs(Position));
                            disp('Setting initial motion parameters');
                            if displacement>2
                                disp('Shuttler not in place. MOVING!');
                                obj.com.WaitMotionEnd (obj.com.ACSC_AXIS_0, 1000*10*displacement/50);
                            end
                        end
                        
                        if handles.Array_PSeq{1}.Channels(12).Enable
                            if ~isempty(handles.Array_PSeq{helper_scan,aux}.Channels(12).Frequency)
                                disp('Getting motion parameters from NVautomizer');
                                Position=handles.Array_PSeq{helper_scan,aux}.Channels(12).Frequency;
                                Velocity=handles.Array_PSeq{helper_scan,aux}.Channels(12).Amplitude;
                                Accn=handles.Array_PSeq{helper_scan,aux}.Channels(12).Phase;
                                Jerk=handles.Array_PSeq{helper_scan,aux}.Channels(12).AmpIQ;
                                c_position=handles.Array_PSeq{helper_scan,aux}.Channels(12).FreqIQ;
                                Wait_time=handles.Array_PSeq{helper_scan,aux}.Channels(12).FreqmodQ;
                                T1delay_time=handles.Array_PSeq{helper_scan,aux}.Channels(24).Frequency; %for wait time experiments

                                %% Used for acquistion delay
                                disp('Writing motion parameters from NVautomizer');
                                obj.com.SetVelocity(obj.com.ACSC_AXIS_0,Velocity);
                                obj.com.SetAcceleration(obj.com.ACSC_AXIS_0,Accn);
                                obj.com.SetDeceleration(obj.com.ACSC_AXIS_0,Accn);
                                obj.com.SetJerk(obj.com.ACSC_AXIS_0,Jerk);
                                
                                obj.com.WriteVariable(T1delay_time*1e3, 'V2', obj.com.ACSC_NONE); %write delay time to Motion Manager
                                obj.com.WriteVariable(Wait_time*1e3, 'V1', obj.com.ACSC_NONE); %need to do 1e3 because ACS times are in ms
                                obj.com.WriteVariable(-Position, 'V0', obj.com.ACSC_NONE);
                                
                                
                                if handles.Array_PSeq{1}.Channels(12).Enable && ~handles.Imaginghandles.ImagingFunctions.interfaceDataAcq.hasAborted
                                    obj.com.ToPoint(obj.com.ACSC_AMF_WAIT, obj.com.ACSC_AXIS_0,-c_position);  %go up to coil position
                                    disp('Run Buffer 1'); %changed from buffer 3 to 1 to run delay exp
                                    obj.com.RunBuffer(1);
                                end
                            end
                        end
                        
                        
                        %% turn laser on during exp. %%
                        if handles.Array_PSeq{1}.Channels(10).Enable
                            if ~isempty(handles.Array_PSeq{helper_scan,aux}.Channels(10).Frequency)
                                disp('Getting laser times from NVautomizer');
                                laser_on_time=handles.Array_PSeq{helper_scan,aux}.Channels(10).Frequency;
                                laser_delay_time=handles.Array_PSeq{helper_scan,aux}.Channels(10).Phase;
                            end
                        end


                        %% Main Sequence
                        if ~err && ~handles.Imaginghandles.ImagingFunctions.interfaceDataAcq.hasAborted
                            disp('Load sequence to Swabian');
                            [stPB,durPB,err] = obj.ListOfStatesAndDurationsPB(handles.Array_PSeq{helper_scan,aux},handles);
                            [pattern0, pattern1, pattern2, pattern3, pattern4,...
                                pattern5, pattern6, pattern7] = make_swabian_pattern (stPB, durPB);

                            sequence = Sequence();
                            sequence.setDigital(0, pattern0);
                            sequence.setDigital(1, pattern1);
                            sequence.setDigital(2, pattern2);
                            sequence.setDigital(3, pattern3);
                            sequence.setDigital(4, pattern4);
                            sequence.setDigital(5, pattern5);
                            sequence.setDigital(6, pattern6);
                            sequence.setDigital(7, pattern7);
                            sequence = sequence.repeat(1);
                            n_runs =1;
                            finalState = OutputState([],0,0);

                        end
                        
                        if err
                            disp('You have aborted!');
                            handles.Imaginghandles.ImagingFunctions.interfaceDataAcq.hasAborted = 1;
                            break;
                        end
                        
                        % END Load Sequence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        if handles.Imaginghandles.ImagingFunctions.interfaceDataAcq.hasAborted==1
                            disp('You have aborted!');
                            obj.com.StopBuffer(1); obj.com.StopBuffer(3); obj.com.StopBuffer(7);obj.com.StopBuffer(6);
                            obj.com.KillAll();
                            handles.ExperimentFunctions.swabian.reset(); 
                        else
                            %%%% MAIN FIRE %%%%%%%%%%%%
                            disp('Firing Swabian');

                            handles.ExperimentFunctions.swabian.stream(sequence, n_runs, finalState); %will fire Swabian
                            
                            if handles.Array_PSeq{1}.Channels(12).Enable % For motion
                                if ~handles.Imaginghandles.ImagingFunctions.interfaceDataAcq.hasAborted
                                    pause(Wait_time+3); %WAIT TIME OF 90 SECONDS + 3 SECONDS DELAY TO ENSURE WE DONT CUT MW OFF EARLY
                                    Sage_write(['4']);
                                    disp('Wrote Sage write 4: Move down and read the data in');
%                                     fmirror = handles.ExperimentFunctions.zaber.detectDevices();
%                                     fmirror_LRpos = fmirror(1).moveAbsolute(-22803); %Absolute Coil Pos [data]
%                                     fmirror_UDpos = fmirror(2).moveAbsolute(-5910); %Absolute Coil Pos [data]
%                                     fprintf('Laser mirror LR moved to magnet pos:\n %.2f.\n', fmirror(1).getPosition());
%                                     fprintf('Laser mirror UD moved to magnet pos:\n %.2f.\n', fmirror(2).getPosition());
                                    pause(1.03*sum(durPB)-Wait_time-3);
                                end
                                
                                handles.ExperimentFunctions.swabian.reset();
                                if handles.Imaginghandles.ImagingFunctions.interfaceDataAcq.hasAborted==1
                                    disp('You have aborted!');
                                    obj.com.StopBuffer(1); obj.com.StopBuffer(3); obj.com.StopBuffer(7);obj.com.StopBuffer(6);
                                    obj.com.KillAll();
                                    handles.ExperimentFunctions.swabian.reset();
                                end

                                %Sage_write(['4']);
%                                 disp('Wrote Sage write 4: Move down and read the data in');
                                
                                t = timer;
                                t.StartFcn = @(~,~) disp('Laser sequence complete, waiting at bottom of magnet');
                                t.TimerFcn = @(~,~) disp('Full NMR pulse sequence complete!');
                                t.StopFcn =  @(~,~) stop(t);
                                t.StartDelay =10+30-laser_on_time-laser_delay_time; %180s exp+30+T1delay_time %100=70s exp +30s wait
                                %t.StartDelay =10; %100
                                t.ExecutionMode = 'singleshot';
                                start(t);
                                wait(t);
% 
%                                 fmirror_LRpos = fmirror(1).moveAbsolute(-23530); %Absolute Hyperpolarize Pos [data]
%                                 fmirror_UDpos = fmirror(2).moveAbsolute(-5810); %Absolute Hyperpolarize Pos [data]
%                                 fprintf('Laser mirror LR moved to hyperpolarize pos:\n %.2f.\n', fmirror(1).getPosition());
%                                 fprintf('Laser mirror UD moved to hyperpolarize pos:\n %.2f.\n', fmirror(2).getPosition());
                            end
                        end
                    end
                    
                    if handles.Imaginghandles.ImagingFunctions.interfaceDataAcq.hasAborted==1
                        obj.com.StopBuffer(1); obj.com.StopBuffer(3); obj.com.StopBuffer(7);obj.com.StopBuffer(6);
                        obj.com.KillAll();
                        handles.ExperimentFunctions.swabian.reset(); 
                        break;
                    end
                end
                
                %handles.ExperimentFunctions.swabian.reset(); 
                %% Wait for shuttler to return up
                pause(0.5);
                t2 = timer;
                t2.StartFcn = @(~,~) disp('Wait for the shuttler to move back up');
                t2.TimerFcn = @(~,~) disp('Shuttler should be up');
                t2.StopFcn =  @(~,~) stop(t2);
                t2.StartDelay = 5; % (put 200s here to match when shuttler reaches top - 20240427) shuttle time is 130s for velocity 5 plus additional 600s wait time to accomodate temperature stabilization of 5 min experiments
                t2.ExecutionMode = 'singleshot';
                start(t2);
                wait(t2);
                
%                 pause(150);
%                 pause(30);
                obj.com.StopBuffer(1);
              
                disp(['Experiment ' num2str(aux) ' complete']);
                disp('--------------------------------------');
                
            end %END loop over first scan dimension end aux
        end
        
        
        function AbortRun(obj,handles)
            
            obj.com.KillAll();
            obj.com.StopBuffer(1); obj.com.StopBuffer(3); obj.com.StopBuffer(7);obj.com.StopBuffer(6);
            handles.ExperimentFunctions.swabian.reset(); 
            
            handles.Imaginghandles.ImagingFunctions.interfaceDataAcq.hasAborted = 1;
            handles.Tracker.hasAborted=1;
            
        end
        
        %%%%% Sequencing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function [err] = LoadSequenceToExperimentPB(obj,PS,rep,handles)
            
            [stPB,durPB,err] = obj.ListOfStatesAndDurationsPB(PS,handles);
            if ~err
                handles.Imaginghandles.ImagingFunctions.interfacePulseGen.sendSequence(stPB,durPB,rep);
            else
                return;
            end
            
        end
        
        function UpdateAvgDataPlot(obj,handles,is_first)
            
            cla(handles.axes_AvgData);
            colors = ['k' 'r' 'b' 'g' 'm' 'y']; % enough for 6 rises of SPD Acq
            if isfield(handles,'scan_nonlinear') && handles.scan_nonlinear==1
                x=handles.nonlinear_data.x;
            else
                x = linspace(handles.ExperimentalScan.vary_begin(1),handles.ExperimentalScan.vary_end(1),handles.ExperimentalScan.vary_points(1));
            end
            
            for j=1:1:length(obj.AvgData)
                y = obj.AvgData{j};
                hold(handles.axes_AvgData,'on');
                line(j)=plot(x,y,['.-' colors(j)],'Parent',handles.axes_AvgData,'MarkerSize',15,'LineWidth',0.5);
                hold(handles.axes_AvgData,'on');
                plot(x,smooth(y),['-' colors(j)],'Parent',handles.axes_AvgData,'LineWidth',1.5);
                
                title(handles.axes_AvgData,datestr(now,'yyyy-mm-dd-HHMMSS'));
                if get(handles.checkbox_errorbar,'Value')
                    errorbar(x, y, obj.Std_dev{j}, '.','Parent',handles.axes_AvgData)
                end
            end
            grid(handles.axes_AvgData, 'on');
            hold(handles.axes_AvgData,'off');
            
            if is_first == 1
                
                a = {'Rise 1', 'Rise 2', 'Rise 3', 'Rise 4', 'Rise 5', 'Rise 6'};
                legend(handles.axes_AvgData,line,a(1:1:length(obj.AvgData)));
                xlabel(handles.axes_AvgData,handles.ExperimentalScan.vary_prop{1});
                ylabel(handles.axes_AvgData,'kcps');
                title(handles.axes_AvgData,datestr(now,'yyyy-mm-dd-HHMMSS'));
            end
            grid(handles.axes_AvgData, 'on');
            drawnow();
            
        end
        
        
        function PinesDataPlot(obj,handles,is_first)
            
            cla(handles.axes_AvgData); %clear axis
            colors = ['k' 'r' 'b' 'g' 'm' 'y']; % enough for 6 rises of SPD Acq
            x = linspace(handles.ExperimentalScan.vary_begin(1),handles.ExperimentalScan.vary_end(1),handles.ExperimentalScan.vary_points(1));
            set(handles.axes_AvgData,'xlim',[min(x) max(x)]);
            diff_x=diff(x);
            delta_x=diff_x(1);
            for j=1:1:length(obj.AvgData)
                y = obj.AvgData{j};
                hold(handles.axes_AvgData,'on');
                line(1)=plot(x(1:2:end)+delta_x,y(1:2:end),['.-' colors(1)],'Parent',handles.axes_AvgData,'MarkerSize',15,'LineWidth',0.5);
                line(2)=plot(x(2:2:end),y(2:2:end),['.-' colors(2)],'Parent',handles.axes_AvgData,'MarkerSize',15,'LineWidth',0.5);
            end
            grid(handles.axes_AvgData, 'on');
            hold(handles.axes_AvgData,'off');
            
            if is_first == 1
                
                a = {'Rise 1', 'Rise 2', 'Rise 3', 'Rise 4', 'Rise 5', 'Rise 6'};
                legend(handles.axes_AvgData,line,a(1:1:length(obj.AvgData)));
                xlabel(handles.axes_AvgData,handles.ExperimentalScan.vary_prop{1});
                ylabel(handles.axes_AvgData,'kcps');
                title(handles.axes_AvgData,datestr(now,'yyyy-mm-dd-HHMMSS'));
            end
            grid(handles.axes_AvgData, 'on');
            drawnow();
            
        end
        
        
        function UpdateSingleDataPlot(obj,handles,is_first)
            
            cla(handles.axes_RawData);
            colors = ['k' 'r' 'b' 'g' 'm' 'y']; % enough for 6 rises of SPD Acq
            if isfield(handles,'scan_nonlinear') && handles.scan_nonlinear==1
                x=handles.nonlinear_data.x;
            else
                x = linspace(handles.ExperimentalScan.vary_begin(1),handles.ExperimentalScan.vary_end(1),handles.ExperimentalScan.vary_points(1));
            end
            for j=1:1:length(obj.RawData)
                y = obj.RawData{j};
                hold(handles.axes_RawData,'on');
                plot(x,y,['.-' colors(j)],'Parent',handles.axes_RawData);
                title(handles.axes_AvgData,datestr(now,'yyyy-mm-dd-HHMMSS'));
            end
            hold(handles.axes_RawData,'off');
            
            if is_first == 1
                a = {'Rise 1', 'Rise 2', 'Rise 3', 'Rise 4', 'Rise 5', 'Rise 6'};
                legend(handles.axes_RawData,a(1:1:length(obj.RawData)));
                xlabel(handles.axes_AvgData,{handles.ExperimentalScan.vary_prop{1}});
                ylabel(handles.axes_AvgData,'kcps');
                title(handles.axes_AvgData,datestr(now,'yyyy-mm-dd-HHMMSS'));
            end
            grid(handles.axes_AvgData,'on');
            drawnow();
            
        end
        
        function LoadExpFromScan(obj,handles)
            
            
            
            set(handles.sliderExpScan,'Visible','off');
            set(handles.text_vary2_param,'Visible','off');
            set(handles.slider_each_avg,'Visible','off');
            set(handles.avg_nb,'Visible','off');
            exps = get(handles.popup_experiments,'String');
            selectedExp= exps{get(handles.popup_experiments,'Value')};
            
            if length(exps)>1
                
                fp = getpref('nv','SavedExpDirectory');
                if strcmp(selectedExp,'Current Experiment')
                    selectedExp= exps{get(handles.popup_experiments,'Value')+1}; %loads the first exp in the list, that corresponds to the displayed 'Current Experiment'
                end
                
                SavedExp = load(fullfile(fp,selectedExp));
                
                
                handles.ExperimentalScanDisplayed = SavedExp.Scan;
                set(handles.edit_NoteExp, 'String', handles.ExperimentalScanDisplayed.Notes)
                set(handles.edit_NoteExp,'Enable','on');
                set(handles.button_NoteExp,'Enable','on');
                
                if length(handles.ExperimentalScanDisplayed.vary_points) == 2 %2d scan
                    set(handles.sliderExpScan,'Visible','on');
                    set(handles.text_vary2_param,'Visible','on');
                    set(handles.sliderExpScan,'Max', handles.ExperimentalScanDisplayed.vary_points(2));
                    set(handles.sliderExpScan,'Min', 1);
                    set(handles.sliderExpScan,'Value', 1);
                    set(handles.text_vary2_param,'String',sprintf('%s = %s',handles.ExperimentalScanDisplayed.vary_prop{2},num2str(handles.ExperimentalScanDisplayed.vary_begin(2))));
                end
                
                if ~isempty(handles.ExperimentalScanDisplayed.scan_nonlinear)
                    if handles.ExperimentalScanDisplayed.scan_nonlinear~=0
                        handles.scan_nonlinear=handles.ExperimentalScanDisplayed.scan_nonlinear;
                        handles.nonlinear_data=handles.ExperimentalScanDisplayed.nonlinear_data;
                        set(handles.checkbox18,'Value',1);
                    else
                        set(handles.checkbox18,'Value',0);
                        
                    end
                    
                    
                end
                
                
                set(handles.uipanel90,'Visible','on');
                
                set(handles.text_scan_avg, 'String', handles.ExperimentalScanDisplayed.Averages);
                set(handles.text_scan_reps, 'String', handles.ExperimentalScanDisplayed.Repetitions);
                set(handles.text_sequence,'String',handles.ExperimentalScanDisplayed.Sequence);
                set(handles.displayed_power,'String',sprintf('Power = %0.3f +- %0.3f mW',handles.ExperimentalScanDisplayed.Laserpower(1),handles.ExperimentalScanDisplayed.Laserpower(2)));
                set(handles.text_name_displayed_seq,'String',handles.ExperimentalScanDisplayed.SequenceName);
                
                
                obj.ExpPlot(handles,handles.ExperimentalScanDisplayed,1);
                obj.EachAvgPlot(handles,handles.ExperimentalScanDisplayed,1,1);
                
                
                %load variable data to table_show_float
                a = size(handles.ExperimentalScanDisplayed.Variable_values);
                tabledata = cell(a(2), 6);
                
                for p=1:1:a(2)  %number of lines in matrix
                    tabledata{p,1} = handles.ExperimentalScanDisplayed.Variable_values{p}.name;
                    
                    %needs something here to display correctly
                    if length(handles.ExperimentalScanDisplayed.vary_prop) == 2
                        j_end = 2;
                    else
                        j_end = 1;
                    end
                    
                    
                    for j=1:1:j_end
                        if strcmp(handles.ExperimentalScanDisplayed.Variable_values{p}.name, handles.ExperimentalScanDisplayed.vary_prop{j})
                            tabledata{p,3} = true;
                            tabledata{p,4} = handles.ExperimentalScanDisplayed.vary_begin(j);
                            tabledata{p,5} = handles.ExperimentalScanDisplayed.vary_step_size(j);
                            tabledata{p,6} = handles.ExperimentalScanDisplayed.vary_end(j);
                        else
                            tabledata{p,2} = handles.ExperimentalScanDisplayed.Variable_values{p}.value;
                        end
                    end
                end
                
                set(handles.table_show_float,'data',tabledata);
                
                %load variable data to table_show_float
                a = size(handles.ExperimentalScanDisplayed.Bool_values);
                
                if a(2) ~= 0
                    tabledata = cell(a(2), 2);
                    for p=1:1:a(2)  %number of lines in matrix
                        tabledata{p,2} = handles.ExperimentalScanDisplayed.Bool_values{p}.name;
                        if handles.ExperimentalScanDisplayed.Bool_values{p}.value
                            tabledata{p,1} = true;
                        else
                            tabledata{p,1} = false;
                        end
                    end
                    
                    set(handles.table_show_bool,'data',tabledata);
                    
                end
                
                set(handles.button_NoteExp, 'Enable','on');
                set(handles.edit_NoteExp, 'Enable','on');
                
                %plot single averages
                if handles.ExperimentalScanDisplayed.Averages > 1
                    set(handles.slider_each_avg,'Max', handles.ExperimentalScanDisplayed.Averages);
                    set(handles.slider_each_avg,'Min', 1);
                    set(handles.slider_each_avg,'Value', 1);
                    set(handles.avg_nb,'String','1');
                    set(handles.slider_each_avg,'Visible','on');
                    set(handles.avg_nb,'Visible','on');
                end
                
                
                % Update handles structure
                gobj = findall(0, 'Name', 'Experiment');
                guidata(gobj, handles);
                
            end
            
        end
        
    end
    
    methods (Static)
        
        function UpdateExpScanList(handles)
            
            fp = 'C:\QEG2\21-C\SavedExperiments';%getpref('nv','SavedExpDirectory');
            files = dir(fp);
            
            datenums = [];
            found_filenames = {};
            % sort files by modifying date desc
            for k=1:length(files)
                % takes the date from the filename
                r = regexp(files(k).name, '(\d{4}-\d{2}-\d{2}-\d{6})','tokens');
                if ~isempty(r) %remove dirs '.' and '..', and other files
                    datenums(end+1) = datenum(r{1}, 'yyyy-mm-dd-HHMMSS');
                    found_filenames{end+1} = files(k).name;
                end
                
            end
            if ~isempty(datenums)
                [~,dind] = sort(datenums,'descend');
            end;
            
            fn{1} = 'Current Experiment';
            
            for k=1:length(datenums),
                
                
                fn{end+1} = found_filenames{dind(k)};
                
                
            end
            set(handles.popup_experiments,'String',fn);
            
            gobj = findall(0,'Name','Experiment');
            guidata(gobj,handles);
            
        end
        
        function SaveExp(handles)
            
            fp = getpref('nv','SavedExpDirectory');
            
            if length(handles.ExperimentalScan.vary_prop) == 2 %2d scan
                handles.num = 2;
                text = ['-vary-',handles.ExperimentalScan.vary_prop{1}, '-',handles.ExperimentalScan.vary_prop{2},'-'];
                
            else %1d scan
                handles.num = 1;
                text = ['-vary-' handles.ExperimentalScan.vary_prop{1} '-'];
            end
            
            a = datestr(now,'yyyy-mm-dd-HHMMSS');
            fn = [num2str(handles.num),'DExp-seq-' handles.ExperimentalScan.SequenceName(1:1:end-4),text,a];
            
            fullFN = fullfile(fp,fn);
            
            Scan = handles.ExperimentalScan;
            
            Scan.DateTime = a;
            
            if ~isempty(Scan.ExperimentData),
                save(fullFN,'Scan');
            else
                uiwait(warndlg({'No data found for current experiment. Saving aborted.'}));
                return;
            end
            
            %send_email(fn,'',[fullFN '.mat']);
            
            gobj = findall(0,'Name','Experiment');
            guidata(gobj,handles);
            handles = guidata(gobj);
            
        end
        
        function ExpPlot(handles,scan,sliderValue)
            
            cla(handles.axes_AvgData,'reset')
            colors = ['k' 'r' 'b' 'g' 'm' 'y']; % enough for 6 rises of SPD Acq
            if ~isempty(scan.scan_nonlinear) && scan.scan_nonlinear==1
                x=scan.nonlinear_data.x;
            else
                x = linspace(scan.vary_begin(1),scan.vary_end(1),scan.vary_points(1));
            end
            %x = linspace(scan.vary_begin(1),scan.vary_end(1),scan.vary_points(1));
            for j=1:1:length(scan.ExperimentData{sliderValue})
                y = scan.ExperimentData{sliderValue}{j};
                %plot(x,y,['.-' colors(j)],'Parent',handles.axes_AvgData,'MarkerSize',15);
                line(j)=plot(x,y,['.-' colors(j)],'Parent',handles.axes_AvgData,'MarkerSize',15,'LineWidth',0.5);
                hold(handles.axes_AvgData,'on');
                plot(x,smooth(y),['-' colors(j)],'Parent',handles.axes_AvgData,'LineWidth',1.5);
                if get(handles.checkbox_errorbar, 'Value')
                    errorbar(x, y, scan.ExperimentDataError{sliderValue}{j}, '.','Parent',handles.axes_AvgData)
                end
            end
            hold(handles.axes_AvgData,'off');
            a = {'Rise 1', 'Rise 2', 'Rise 3', 'Rise 4', 'Rise 5', 'Rise 6'};
            legend(handles.axes_AvgData,line,a(1:1:length(scan.ExperimentData{sliderValue})));
            xlabel(handles.axes_AvgData,scan.vary_prop{1});
            ylabel(handles.axes_AvgData,'kcps');
            title(handles.axes_AvgData,scan.DateTime);
            grid(handles.axes_AvgData, 'on');
            drawnow();
        end
        
        %         function ExpPlot_spectrum(handles,fidname,FTtype)
        %             fname=fidname;
        %             [y, SizeTD2, SizeTD1] = Qfidread(fname, 17000 ,1, 5, 1);
        %             MatrixOut = matNMRFT1D(y, FTtype);
        %             cla(handles.axes_AvgData,'reset')
        %             %colors = ['k' 'r' 'b' 'g' 'm' 'y']; % enough for 6 rises of SPD Acq
        %                 plot(abs(MatrixOut),['.-' 'b'],'Parent',handles.axes_AvgData,'MarkerSize',15,'LineWidth',0.5);
        %                 hold(handles.axes_AvgData,'on');
        %             hold(handles.axes_AvgData,'off');
        %             xlabel(handles.axes_AvgData,'frequency');
        %             ylabel(handles.axes_AvgData,'A');
        %             title(handles.axes_AvgData,scan.DateTime);
        %             grid(handles.axes_AvgData, 'on');
        %             drawnow();
        %              %plot(abs(MatrixOut));
        %         end
        
        function EachAvgPlot(handles,scan,sliderValue,twodsliderValue)
            
            colors = ['k' 'r' 'b' 'g' 'm' 'y']; % enough for 6 rises of SPD Acq
            if ~isempty(scan.scan_nonlinear) && scan.scan_nonlinear==1
                x=scan.nonlinear_data.x;
            else
                x = linspace(scan.vary_begin(1),scan.vary_end(1),scan.vary_points(1));
            end
            %x = linspace(scan.vary_begin(1),scan.vary_end(1),scan.vary_points(1));
            for j=1:1:length(scan.ExperimentDataEachAvg{twodsliderValue}{sliderValue})
                y = scan.ExperimentDataEachAvg{twodsliderValue}{sliderValue}{j};
                plot(x,y,['.-' colors(j)],'Parent',handles.axes_RawData);
                hold(handles.axes_RawData,'on');
                if get(handles.checkbox_errorbar_each_avg, 'Value')
                    errorbar(x, y, scan.ExperimentDataErrorEachAvg{twodsliderValue}{sliderValue}{j}, '.','Parent',handles.axes_RawData)
                end
            end
            hold(handles.axes_RawData,'off');
            a = {'Rise 1', 'Rise 2', 'Rise 3', 'Rise 4', 'Rise 5', 'Rise 6'};
            legend(handles.axes_RawData,a(1:1:length(scan.ExperimentDataEachAvg{twodsliderValue}{sliderValue})));
            xlabel(handles.axes_RawData,scan.vary_prop{1});
            ylabel(handles.axes_RawData,'kcps');
            title(handles.axes_RawData,datestr(now,'yyyy-mm-dd-HHMMSS'));
            grid(handles.axes_RawData, 'on');
            drawnow();
        end
        
        
        
        function [IShape,AOMShape,MWShape,QShape,SPDShape,CShape,NShape,NullShape,err] = ListOfStatesAndDurationsAWG(PulseSequence,handles)
            
            err = 0;
            % This function takes the PulseSequence and transforms it into a format
            % that is readable by the AWG (exclusively)
            
            % get max and min time for this sequence
            mintime = PulseSequence.GetMinRiseTime();
            maxtime = PulseSequence.time_pointer;
            
            num_points = int64((maxtime-mintime)*handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate);
            
            NullShape = zeros(1,num_points);
            IShape = zeros(1,num_points);
            AOMShape = zeros(1,num_points);
            MWShape = zeros(1,num_points);
            QShape = zeros(1,num_points);
            SPDShape = zeros(1,num_points);
            CShape = zeros(1,num_points);
            NShape = zeros(1,num_points);
            
            %AWG states must be a multiple of the clock frequency
            for hlp=1:1:length(PulseSequence.Channels)
                if PulseSequence.Channels(hlp).Enable
                    
                    if ~isequal(mod(round(PulseSequence.Channels(hlp).RiseDurations*1e11),round(1e11/handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate)),zeros(1,length(PulseSequence.Channels(hlp).RiseDurations)))
                        % without round won't work; multiplying by 1e11
                        err = 1;
                        uiwait(warndlg('Sequence durations not a multiple of AWG sample rate. Abort. (Or: you enabled a channel that you are not using)'));
                        return;
                    end
                end
            end
            
            if PulseSequence.Channels(1).Enable  %pulsed mode for laser
                
                for p=1:1:length(PulseSequence.Channels(1).RiseTimes)
                    AOMShape(int32(1+(PulseSequence.Channels(1).RiseTimes(p)+mintime)*handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate):1:int32((PulseSequence.Channels(1).RiseTimes(p)+mintime+PulseSequence.Channels(1).RiseDurations(p))*handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate)) = 1;
                end
            end
            
            if PulseSequence.Channels(2).Enable  %pulsed mode for mw
                
                for p=1:1:length(PulseSequence.Channels(2).RiseTimes)
                    MWShape(int32(1+(PulseSequence.Channels(2).RiseTimes(p)+mintime)*handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate):1:int32((PulseSequence.Channels(2).RiseTimes(p)+mintime+PulseSequence.Channels(2).RiseDurations(p))*handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate)) = 1;
                end
                
                n = round(100e9*(1/handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate/2))/100e9:round(100e9*(1/handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate))/100e9:round(100e9*((maxtime-mintime)-(1/handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate/2)))/100e9;
                IShape =  sin(2*pi*PulseSequence.Channels(2).FreqIQ*n + PulseSequence.Channels(2).Phasemod).*PulseSequence.Channels(2).Ampmod + PulseSequence.Channels(2).FreqmodI;
                QShape =  cos(2*pi*PulseSequence.Channels(2).FreqIQ*n + PulseSequence.Channels(2).Phasemod).*PulseSequence.Channels(2).Ampmod + PulseSequence.Channels(2).FreqmodQ; %or sin(x + pi/2)
                
            end
            
            if PulseSequence.Channels(3).Enable  %pulsed mode for SPD
                
                for p=1:1:length(PulseSequence.Channels(3).RiseTimes)
                    SPDShape(int32(1+(PulseSequence.Channels(3).RiseTimes(p)+mintime)*handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate):1:int32((PulseSequence.Channels(3).RiseTimes(p)+mintime+PulseSequence.Channels(3).RiseDurations(p))*handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate)) = 1;
                end
                
            end
            
            if PulseSequence.Channels(4).Enable  %pulsed mode for C
                
                n = (1/handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate/2):(1/handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate):((maxtime-mintime)-(1/handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate/2));
                CShape = sin(2*pi*PulseSequence.Channels(4).Frequency*n + PulseSequence.Channels(4).Phasemod).*PulseSequence.Channels(4).Ampmod + PulseSequence.Channels(4).FreqmodI;
            end
            
            if PulseSequence.Channels(5).Enable %pulsed mode for N
                %Masashi commented 20130409;
                %n = (1/handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate/2):(1/handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate):((maxtime-mintime)-(1/handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate/2));
                n = round(100e9*(1/handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate/2))/100e9:round(100e9*(1/handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate))/100e9:round(100e9*((maxtime-mintime)-(1/handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate/2)))/100e9;
                NShape = sin(2*pi*PulseSequence.Channels(5).Frequency*n + PulseSequence.Channels(5).Phasemod).*PulseSequence.Channels(5).Ampmod + PulseSequence.Channels(5).FreqmodI;
            end
            
            %Make sure first and last states of the array are 0, otherwise markers will be "on"
            %at undesirable places
            
            %put a zero-shape in the beginning and end
            additional_zeros = handles.ExperimentFunctions.delay_init_and_end_seq*handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate; %ALWAYS integer bc min sampl rate 1e7
            
            NullShape = [zeros(1,additional_zeros) NullShape zeros(1,additional_zeros)];
            IShape = [zeros(1,additional_zeros) IShape zeros(1,additional_zeros)];
            AOMShape =  [zeros(1,additional_zeros) AOMShape zeros(1,additional_zeros)];
            MWShape =  [zeros(1,additional_zeros) MWShape zeros(1,additional_zeros)];
            QShape =  [zeros(1,additional_zeros) QShape zeros(1,additional_zeros)];
            SPDShape =  [zeros(1,additional_zeros) SPDShape zeros(1,additional_zeros)];
            CShape =  [zeros(1,additional_zeros) CShape zeros(1,additional_zeros)];
            NShape = [zeros(1,additional_zeros) NShape zeros(1,additional_zeros)];
            
%             figure();
%             subplot(2,2,1);
%             plot(MWShape);
%             subplot(2,2,2);
%             plot(QShape);
%             subplot(2,2,3);
%             plot(CShape);
%             subplot(2,2,4);
%             plot(NShape);
        end
        
        function [statesPB,durationsPB,err] = ListOfStatesAndDurationsPB(PulseSequence,handles)
            
            err = 0;
            
            % This function takes the PulseSequence and transforms it into a format
            % that is readable by the Pulse Blaster. For the Pulse Blaster,
            % this format consists of the binary value representing the states of the four lines, plus
            % the duration of each state in s.
            
            %Start all channels at off
            %Everytime there is an event related to a channel, bitflip operation
            
            %initialize to sth
            Big_matrix_times = [];
            Big_matrix_indices = [];
            
            current_state = 0;
            
            for k=1:1:length(PulseSequence.Channels)
                
                if PulseSequence.Channels(k).Enable  %pulsed mode
                    
                    Final_times = PulseSequence.Channels(k).RiseTimes + PulseSequence.Channels(k).RiseDurations;
                    
                    Unsorted_times = [PulseSequence.Channels(k).RiseTimes Final_times];
                    %Times = [to^1, tf^1, to^2, tf^2, ... to^last, tf^last]
                    Times = sort(Unsorted_times);
                    Channel_ind = k*ones(1,length(Times));
                    
                    Big_matrix_times = [Big_matrix_times Times];
                    Big_matrix_indices = [Big_matrix_indices Channel_ind];
                    
                    
                    
                    %else, or off, does not need to do anything
                    
                end
                
            end
            
            % Treat PB data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            [Sorted_big_matrix_times, Ind] = sort(Big_matrix_times);
            Sorted_big_matrix_indices = Big_matrix_indices(Ind);
            
            Sorted_list_states = [];
            
            for k=1:1:(length(Sorted_big_matrix_times)-1)
                
                new_state = SingleBitFlip(current_state, Sorted_big_matrix_indices(k));
                
                Sorted_list_states = [Sorted_list_states new_state];
                
                current_state = new_state;
                
            end
            
            Final_sorted_big_matrix_times = diff(Sorted_big_matrix_times);
            Final_sorted_list_states = Sorted_list_states;
            
            A = find(Final_sorted_big_matrix_times > 1e-12);    %~= 0 gives a numerical bug
            Final_sorted_list_states = Final_sorted_list_states(A);
            Final_sorted_big_matrix_times = Final_sorted_big_matrix_times(A);
            
            % first and last states are 0
            Final_sorted_big_matrix_times = [handles.ExperimentFunctions.delay_init_and_end_seq Final_sorted_big_matrix_times handles.ExperimentFunctions.delay_init_and_end_seq];
            Final_sorted_list_states = [0 Final_sorted_list_states 0];
            
            durationsPB = Final_sorted_big_matrix_times;
            statesPB = Final_sorted_list_states;
            
            %PB states must be a multiple of the clock frequency
            if ~isequal(mod(round(durationsPB*1e11),round(1e11/handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate)),zeros(1,length(durationsPB)))
                % without round won't work; multiplying by 1e11 as in
                % send_sequence in pulse blaster class
                err = 1;
                uiwait(warndlg('Sequence durations not a multiple of Pulse Blaster sample rate. Abort.'));
                return;
            end
            
            %PB has a latency of 5 clock cycles;
            latency = 5*(1/handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SampleRate);
            if ~isempty(find(durationsPB<latency-10*eps, 1)) %10eps for numerical value
                %err = 1;
                err_latency=1
                % uiwait(warndlg('Pulse Blaster latency critera not met. Abort.'));
                return;
            end
            
            
            % END Treat PB data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
        end
        
        function TrackCenterPB(handles)
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.StartProgramming();
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.PBInstruction(1,handles.Imaginghandles.ImagingFunctions.interfacePulseGen.INST_CONTINUE, 0, 750,handles.Imaginghandles.ImagingFunctions.interfacePulseGen.ON);
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.PBInstruction(1,handles.Imaginghandles.ImagingFunctions.interfacePulseGen.INST_BRANCH, 0, 150,handles.Imaginghandles.ImagingFunctions.interfacePulseGen.ON);
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.StopProgramming();
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.PBStart();
            handles.Imaginghandles.ImagingFunctions.TrackCenter(handles.Imaginghandles);
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.PBStop();
        end
        
        function TrackCenterAWG(handles)
            
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.writeToSocket('AWGCONTROL:RMODE CONTINUOUS');
            
            NullShape = zeros(6*1e4,1); %length here is arbitrary
            OnShape = [zeros(1,1)' ones(6*1e4-2,1)' zeros(1,1)']';
            
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.clear_waveforms();
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.create_waveform('POWMEAS',NullShape,OnShape,NullShape);
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.setSourceWaveForm(1,'POWMEAS');
            
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SetChannelOn(1);
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.AWGStart();
            go_on = 0;
            while ~go_on
                go_on = handles.Imaginghandles.ImagingFunctions.interfacePulseGen.writeReadToSocket('AWGControl:RSTate?');
            end
            
            handles.Imaginghandles.ImagingFunctions.TrackCenter(handles.Imaginghandles);
            
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SetChannelOff(1);
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.AWGStop();
        end
        
        function [power_array_in_V] = PowerMeasurementPB(handles)
            
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.StartProgramming();
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.PBInstruction(1,handles.Imaginghandles.ImagingFunctions.interfacePulseGen.INST_CONTINUE, 0, 750,handles.Imaginghandles.ImagingFunctions.interfacePulseGen.ON);
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.PBInstruction(1,handles.Imaginghandles.ImagingFunctions.interfacePulseGen.INST_BRANCH, 0, 150,handles.Imaginghandles.ImagingFunctions.interfacePulseGen.ON);
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.StopProgramming();
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.PBStart();
            [laser_power_in_V,std_laser_power_in_V,power_array_in_V] = handles.Imaginghandles.ImagingFunctions.MonitorPower();
            pause(1);
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.PBStop();
            
            [laser_power,std_laser_power] = PhotodiodeConversionVtomW(laser_power_in_V,std_laser_power_in_V);
            set(handles.text_power,'String',sprintf('Power = %0.3f +- %0.3f mW',laser_power,std_laser_power));
            
        end
        
        function [power_array_in_V] = PowerMeasurementAWG(handles)
            
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.writeToSocket('AWGCONTROL:RMODE CONTINUOUS');
            
            NullShape = zeros(6*1e4,1); %length here is arbitrary
            OnShape = [zeros(1,1)' ones(6*1e4-2,1)' zeros(1,1)']';
            
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.clear_waveforms();
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.create_waveform('POWMEAS',NullShape,OnShape,NullShape);
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.setSourceWaveForm(1,'POWMEAS');
            
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SetChannelOn(1);
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.AWGStart();
            while ~go_on
                go_on = handles.Imaginghandles.ImagingFunctions.interfacePulseGen.writeReadToSocket('AWGControl:RSTate?');
            end
            
            [laser_power_in_V,std_laser_power_in_V,power_array_in_V] = handles.Imaginghandles.ImagingFunctions.MonitorPower();
            handles.Imaginghandles.ImagingFunctions.interfacePulseGen.SetChannelOff(1);
            
            [laser_power,std_laser_power] = PhotodiodeConversionVtomW(laser_power_in_V,std_laser_power_in_V);
            set(handles.text_power,'String',sprintf('Power = %0.3f +- %0.3f mW',laser_power,std_laser_power));
        end
        
        
        
        
    end
    
end