%% Clear everything
clear;
close;

%% Set defaults Vars
savealldata=false;
savesinglechunk=false;
savemultwind=false;
sampleRate = 1000e6; %sampling rate for digitizing signal
sampleRateDAC = 9e9; %sampling rate for producing chirp
adcDualChanMode = 1;
% fullScaleMilliVolts =1000;
trigSource = 1; % 1 = external-trigger
dacChanInd = 1; % 1 = chan 1 and 2 = chan 2
adcChanInd = 0; % ADC Channel 1
measurementTimeSeconds = 7; %Integer
delay = 0.00000173; % dead time


% remoteAddr = '192.168.1.2'; % old computer
% remoteAddr = '192.168.10.5'; % new computer
remoteAddr = '192.168.0.122'; % new computer
remotePort = 2020;
localPort = 9090;

off = 0;
on = 1;
pfunc = ProteusFunctions;

dll_path = 'C:\\Windows\\System32\\TEPAdmin.dll';

cType = "DLL";  %"LAN" or "DLL"

paranoia_level = 2;

if cType == "LAN"
    try
        connStr = strcat('TCPIP::',connStr,'::5025::SOCKET');
        inst = TEProteusInst(connStr, paranoia_level);
        
        res = inst.Connect();
        assert (res == true);
    catch ME
        rethrow(ME)
    end   
else
    asm = NET.addAssembly(dll_path);

    import TaborElec.Proteus.CLI.*
    import TaborElec.Proteus.CLI.Admin.*
    import System.*
    
    admin = CProteusAdmin(@OnLoggerEvent);
    rc = admin.Open();
    assert(rc == 0);   
    
    try
        slotIds = admin.GetSlotIds();
        numSlots = length(size(slotIds));
        assert(numSlots > 0);
        
        % If there are multiple slots, let the user select one ..
        sId = slotIds(1);
        if numSlots > 1
            fprintf('\n%d slots were found\n', numSlots);
            for n = 1:numSlots
                sId = slotIds(n);
                slotInfo = admin.GetSlotInfo(sId);
                if ~slotInfo.IsSlotInUse
                    modelName = slotInfo.ModelName;
                    if slotInfo.IsDummySlot
                        fprintf(' * Slot Number:%d Model %s [Dummy Slot].\n', sId, modelName);
                    else
                        fprintf(' * Slot Number:%d Model %s.\n', sId, modelName);
                    end
                end
            end
            pause(0.1);
            choice = 8 %input('Enter SlotId ');
            fprintf('\n');
            sId = uint32(choice);
        end
        
        % Connect to the selected instrument ..
        should_reset = true;
        inst = admin.OpenInstrument(sId, should_reset);
        instId = inst.InstrId;
        
    catch ME
        admin.Close();
        rethrow(ME) 
    end    
end
    
    % ---------------------------------------------------------------------
    % Setup instrument
    % ---------------------------------------------------------------------
    
    res = inst.SendScpi('*IDN?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nConnected to ''%s''\n', netStrToStr(res.RespStr));
    
    pause(0.01);
    
    res = inst.SendScpi('*CLS'); % clear
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi('*RST'); % reset
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi(':FREQ:RAST 9E9'); % set sample clock
    assert(res.ErrCode == 0);
    fprintf('Sample clock set\n')
    
    fprintf('Reset complete\n');
    
    % ---------------------------------------------------------------------
    % ADC Config
    % ---------------------------------------------------------------------

res = inst.SendScpi('INST:CHAN 1');% num2str(dacChanInd)]); %
assert(res.ErrCode == 0);
inst.SendScpi(':DIG:DDC:MODE COMP');
inst.SendScpi(':DIG:DDC:CFR1 0.0');
inst.SendScpi(':DIG:DDC:PHAS1 90.0');
fprintf('ADC Configured\n');
    
    %% Measurement Loop
    u2 = udp(remoteAddr, 'RemotePort', remotePort, 'LocalPort', localPort);
    try
        fopen(u2);
        connect=on;
        fprintf('Server connected and started\n');
    catch
        fclose(u2);
        connect=off;
        fprintf('Server not connected\n');
    end
    
    % Loop to repeatedly wait for messages and send replies
    % Break or Ctrl+C to get out of loop
    while ( connect==on )
        % Wait for message
        while(u2.BytesAvailable == 0)
            % If no bytes in u2 buffer, wait 10ms then check again
            pause(0.01);
        end
        %         cmdBytes = fread(u2);
        readBytes = fscanf(u2)
        dataBytes=1;counter=1;
        while ~isempty(readBytes)
            [tempBytes,readBytes]=strtok(readBytes,',');
            cmdBytes(counter)=str2num(tempBytes);
            counter=counter+1;
        end
        %dataBytes=strtok(dataBytes,',');
        cmdByte=cmdBytes(1);
        
        % 1 - Init
        % 2 - Make MW and play
        % 3 - Setup Receive
        % 4 - MAIN READOUT CODE
        % 5 - Disconnect Intrument
        % 6 - Cleanup
        % 7 - Play MW chirp waveform
        
        %measure_time=str2double(dataBytes); %in sec
        disp(['Command byte is:' cmdByte]);
        switch cmdByte
            case 1 % Init

                disp('test worked');
                
            case 2 % Make MW and play

                
                awg_center_freq = cmdBytes(2);
                awg_bw_freq = cmdBytes(3);
                awg_amp = cmdBytes(4);
                sweep_freq = cmdBytes(5);
                sweep_sigma = cmdBytes(6);
                symm = cmdBytes(7);
                srs_freq = cmdBytes(8);
                srs_amp = cmdBytes(9);
                


                sweepers_en=0;
                bits=16;
                segment=1;
                rampTime = 1/sweep_freq;
                fCenter = awg_center_freq - srs_freq;
                fStart = fCenter - 0.5*awg_bw_freq;
                fStop = fCenter + 0.5*awg_bw_freq;
                dt = 1/sampleRateDAC;
              
                t =  0:1/sampleRateDAC:rampTime;
                dacWave = chirp(t,fStart,rampTime,fStop);
                seglenTrunk = (floor(length(dacWave)/64))*64;
                dacWave = dacWave(1:seglenTrunk);

                maxSig = max(dacWave);
                verticalScale = ((2^bits/2)-1);
                vertScaled = (dacWave/ maxSig) * verticalScale;
                dacSignal = uint16(vertScaled + verticalScale);
                dacSignal = typecast(dacSignal, 'uint8');
                fprintf('waveform length - ');
                fprintf(num2str(length(dacSignal)));
                fprintf('\n') ;
               % plot(dacSignal(1:1000))



                % Define segment 1
                res = inst.SendScpi(strcat(':TRAC:DEF 1,',num2str(length(dacSignal)))); % define memory location 1 of length dacSignal
                assert(res.ErrCode == 0);
                

                % select segmen 1 as the the programmable segment
                res = inst.SendScpi(':TRAC:SEL 1');
                assert(res.ErrCode == 0);

                res = inst.SendScpi(':TRAC:FORM U16');
                assert(res.ErrCode == 0);

                % Download the binary data to segment 1
                %res = inst.WriteBinaryData(':TRAC:DATA 0,#',dacSignal);% fliplr(dacSignal) to reverse sweep; dacSignal to set sweep back
                res = inst.WriteBinaryData(strcat(':TRAC:DATA 0,#4',num2str(length(dacSignal)),','),dacSignal);% fliplr(dacSignal) to reverse sweep; dacSignal to set sweep back
                
                assert(res.ErrCode == 0);

                srs_freq_str = [':SOUR:NCO:CFR1  ' sprintf('%0.2e', srs_freq)];
                res = inst.SendScpi(srs_freq_str);
                assert(res.ErrCode == 0);
                srs_freq_str = [':SOUR:NCO:CFR2 ' sprintf('%0.2e', srs_freq)];
                res = inst.SendScpi(srs_freq_str);
                assert(res.ErrCode == 0);

                res=inst.SendScpi(':NCO:SIXD1 ON');
                assert(res.ErrCode == 0);

                res = inst.SendScpi(':SOUR:MODE DUC'); % IQ MODULATOR --
                assert(res.ErrCode == 0);
                %THINK OF A MIXER, BUT DIGITAL

                sampleRateDAC_str = [':FREQ:RAST ' sprintf('%0.2e', sampleRateDAC)];
                res = inst.SendScpi(sampleRateDAC_str); % set sample clock
                assert(res.ErrCode == 0);


               
                % Play MW chirp waveform
                % ---------------------------------------------------------------------
                % Play segment 1 in channel 1
                % ---------------------------------------------------------------------
                res = inst.SendScpi('INST:CHAN 1'); % plays MW out of DAC channel 1
                assert(res.ErrCode == 0);

                res = inst.SendScpi(':SOUR:FUNC:MODE:SEG 1');
                assert(res.ErrCode == 0);

                %amp_str = [':SOUR:VOLT ' sprintf('%0.2f', awg_amp)];
                res = inst.SendScpi(':SOUR:VOLT MAX');
                %res = inst.SendScpi(amp_str);
                
                assert(res.ErrCode == 0);

                res = inst.SendScpi(':OUTP ON');
                assert(res.ErrCode == 0);

                fprintf('Waveform generated and playing\n');
              
            case 3 % Setup Receive
                sampleRate = 1000e6;
                adcDualChanMode = 1;
                fullScaleMilliVolts =1000;
                adcChanInd = 0; % ADC Channel 1 %%2
                trigSource = 1; % 1 = external-trigger

                res = inst.SendScpi(':DIG:MODE DUAL');
                assert(res.ErrCode == 0);

                fprintf(strcat([':DIG:FREQ ',num2str(sampleRate), '\n']));
                cmd = [':DIG:FREQ ' num2str(sampleRate)];
                res = inst.SendScpi(cmd);
                assert(res.ErrCode == 0);

                % Enable capturing data from channel 1
                res = inst.SendScpi(':DIG:CHAN:SEL 1'); 
                assert(res.ErrCode == 0);

                res = inst.SendScpi(':DIG:CHAN:STATE ENAB');
                % Select internal or external as start-capturing trigger:
                %res = inst.SendScpi(':DIG:TRIG:SOURCE CPU');

                res = inst.SendScpi(':DIG:TRIG:SOURCE EXT');
                assert(res.ErrCode == 0);

                %more external capture commands?
                %Set type to edge trigger
                res = inst.SendScpi(':DIG:TRIG:TYPE EDGE');
                assert(res.ErrCode == 0);

                %Set trigger to falling edge
                res = inst.SendScpi(':DIG:TRIG:SLOP NEG');
                assert(res.ErrCode == 0);

                %Set trigger threshold to 0.3 V
                %Check that LEV1 vs LEV2 is correct and
                %that 0.3 is pos/neg
                res = inst.SendScpi('DIG:TRIG:LEV1 0.3');  % AWG will record anything above 0.3V
                assert(res.ErrCode == 0);

                %Set external trigger delay
                %300 value is 4.1e-6, test for 400 system and re-enter correct val
                res = inst.SendScpi('DIG:TRIG:DEL:EXT 10.5e-06');%IG:TRIG:DEL:EXT 5.5e-06'
                assert(res.ErrCode == 0);
                fprintf('Instrument setup complete and ready to aquire\n');

                pw = cmdBytes(2);
                PW2 = cmdBytes(3);
                numberOfPulses_total = cmdBytes(4);
                tacq = cmdBytes(5);
                numberOfPulses = cmdBytes(6);
                loops = cmdBytes(7);
                Tmax = cmdBytes(8);

                %readLen = (tacq+4)*2^(-20)/(2^(-30))-96; % copied from digTest3 as comparison only!!!!
                readLen = (32+4)*2^(-20)/4/16/(2^(-30))-96; % must be divisible by 96, HARDCODED TACQ TO BE 32!!!! DONT FORGET TO SWITCH BACK!!!!!
                %readLen=readLen*2; %to acquire 16 us of Varian acq window
                %If above line is commented, we are only acquiring 8us window
                readLen=readLen/2; %this will acquire only 4us windows

                netArray = NET.createArray('System.UInt16', 4*readLen*numberOfPulses); %total array -- all memory needed
                % sets the memory size of Proteus
                cmd = [':DIG:ACQuire:FRAM:DEF ' num2str(numberOfPulses*loops) ',' num2str(2*readLen)]; % if we set this to just readLen here, would it only capture 16 us?
                res = inst.SendScpi(cmd);
                assert(res.ErrCode == 0);
                fprintf('Waiting for Shuttle\n');

                pause(0.2);

            case 4 % MAIN READOUT CODE
                % Set this to capture all data
                res = inst.SendScpi('DIG:ACQ:CAPT:ALL');
                assert(res.ErrCode == 0);

                % Other code has it set the number
                % of frames to capture
                cmd = [':DIG:ACQuire:FRAM:CAPT 1,' num2str(numberOfPulses*loops)]; %hard-coded in firstIndex=1
                %setting full data as total # frames
                res = inst.SendScpi(cmd);
                assert(res.ErrCode == 0);

                % Close instrument in case improperly closed from prior run
                res = inst.SendScpi(':DIG:INIT OFF');
                assert(res.ErrCode == 0);
                % Starts the digitizer to wait for
                % capture
                res = inst.SendScpi(':DIG:INIT ON');
                assert(res.ErrCode == 0);
    
                fprintf('Waiting for trigger \n');
                rTrg=0;
                fprintf(['pulses expected =', num2str([numberOfPulses*loops]), '\n']);
                %pause(192); 
                while rTrg < numberOfPulses*loops %set to total # triggers
                    %res = inst.SendScpi(':DIG:TRIG:IMM');
                    res = inst.SendScpi(":DIG:ACQuire:FRAM:STATus?");
                    assert(res.ErrCode == 0);
                    respFields= split(convertCharsToStrings(char(res.RespStr)),',');
                    %fprintf(1, '\nADC Status  ''%s''\n', char(res.RespStr));
                    paramlen = size(respFields);
                    if paramlen(1) >=4
                        rTrg = str2num(respFields(4)); %the 4th element of respFields is numberOfPulses
                    end
                    pause(2);
                end
                fprintf(1, '\nADC Status  ''%s''\n', char(res.RespStr));

                % Turn proteus acquisition off before sending data to computer
                res = inst.SendScpi(':DIG:INIT OFF');
                assert(res.ErrCode == 0);

                % Choose what to read (only the frame-data without the header in this example)
                res = inst.SendScpi(':DIG:DATA:TYPE FRAM');
                assert(res.ErrCode == 0);

                fprintf('Transfering acquired data to computer....\n')

                pulseAmp = [];
                relPhase = [];
                %TESTpulses = [];

                for n = 1:loops
                    fprintf('Start Read %d .... ', n);
                    firstIndex = ((n-1)*numberOfPulses)+1;
                    tic
                    % Select frames starting at next index
                    cmd = [':DIG:DATA:SEL ' num2str(firstIndex) ',' num2str(numberOfPulses)];
                    res = inst.SendScpi(cmd);
                    assert(res.ErrCode == 0);

                    rc = inst.ReadMultipleAdcFrames(adcChanInd, firstIndex, numberOfPulses, netArray); %this is where the device reads
                    assert(rc == 0);
                    samples = double(netArray); %get the data
                    samples = samples(1:2:length(samples)); 
                    %
                    samples = samples(1:2:length(samples));
                    fprintf('Read %d Complete\n', n);
                    toc

                    tic
                    %delete mem
                    fprintf('Clear mem %d .... ', n);
                    rc =inst.WipeMultipleAdcFrames(adcChanInd, firstIndex, numberOfPulses, 0);
                    assert(rc == 0);
                    %                                                     cmd = [':DIG:ACQ:ZERO ' num2str(firstIndex) ',' num2str(numberOfPulses) ',' num2str(0)];
                    %                                                     res = inst.SendScpi(cmd);
                    assert(res.ErrCode == 0);
                    fprintf('Clear mem %d Complete\n', n);
                    toc

                    tic
                    fprintf('Starting iteration %d data processing....', n);
                    pulses = reshape(samples, [], numberOfPulses); % reshape samples into a more useful array (2 dimensions)

                   
                    clear samples;

                    % Do we want to include this zero padding before above data processing?
                    power_of_2 = floor(log2(readLen)); % this is normally 15 for 32us captures
                    padded_len= 2^(power_of_2) ;%2^15;
                    % We already essentially do the next few lines
                    dF = sampleRate/16/padded_len; %set the discretization of freq in terms of sampleRate
                    f = -sampleRate/16/2:dF:sampleRate/16/2-dF; %sampleRate sets the 'bandwidth'
                    %
                    [v,b1]=min(abs(f-20.00706e6)); %picks out the 20MHz component
                    [v,b2]=min(abs(f+20.00706e6));

                    % This truncates data but we don't use?
                    eff_read=100:readLen-100;
                    cyclesPoints = 50;


                    for i = 1:numberOfPulses
                        pulse = pulses(:, i);
                        if n == 1
                            if i == numberOfPulses%before i == 2
                                figure(4);clf;
                                plot(pulse);
                                xlabel('Number of Points [a.u.]');
                                ylabel('Signal[a.u.]');
                            end
                        end
                        pulseMean = mean(pulse);
                        pulseDC = pulse - pulseMean; % remove DC
                        X = fftshift(fft(pulseDC,padded_len)); % perform FFT and zero-pad to the padded_len
                        linFFT = (abs(X)/readLen);
                        amp=(linFFT(b1) + linFFT(b2));
                        phase1=angle(X(b1));
                        phase2=angle(X(b2));
                        %phase=phase1-phase2; %why does this exist?
                        pulseAmp(i+(numberOfPulses*(n-1))) = amp;
                        relPhase(i+(numberOfPulses*(n-1))) = phase1;
                        %pulseAmp = [pulseAmp, amp];
                        %relPhase = [relPhase, phase1];
                    end
                    %save(['X:\d'] ,'pulses');

                    clear pulses;
                    fprintf('Data processing iteration %d complete!\n', n);
                    toc
                end

                ivec=1:length(pulseAmp);
                delay2 = 0.000003; % dead time the unknown one, this is actually rof3 -Ozgur aka this is probably something we have to change?
                time_cycle = ((PW2)+delay2+(tacq+2+4+10)*1e-6)*1; %this may be set differently as we plan to trigger on falling edge, time to do a pulse-acquire
                % using pulsed_spinlock_altskip, we take two time cycles to collect a data point
                time_axis=time_cycle.*ivec;
                %drop first point
                time_axis(1)=[];pulseAmp(1)=[];relPhase(1)=[]; %why only drop the first point? why drop it at all?
                a = datestr(now,'yyyy-mm-dd-HHMMSS');
                fn = sprintf([a,'_Proteus400']);
                % Save data
                try
                    figure(1);clf;
                    scatter(time_axis,pulseAmp,'r','.');
                    xlabel('Time [s]');
                    ylabel('Total Signal [a.u.]');
                    set(gca,'ylim',[0 max(pulseAmp)*1.1]);
                    figure(12);clf;
                    plot(time_axis,relPhase);
                    xlabel('Time [s]');
                    ylabel('Relative Phase');
                catch
                    disp('Plot error occured');
                end

                fprintf('Writing data to Quantum Stream:.....\n');
                save(['X:\' fn],'pulseAmp','time_axis','relPhase'); %maybe it works
                fprintf('Save complete\n');

            
                fprintf('Complete ready for next aquisition!\n');

                %for relay expts/T1(B)
                %pause(40); %40 sec laser pumping

           
                % Turn OFF AWG and SRS
              %  if length(handles.Array_PSeq{1}.Channels)>=39
                %    if handles.Array_PSeq{1}.Channels(39).Enable
                        res = inst.SendScpi(':OUTP OFF');
                        assert(res.ErrCode == 0);
                       % rc = admin.CloseInstrument(instId);
                       % admin.Close();
                        fprintf(1, 'END\n');
                   % end
               % end
                pause(2);


            case 5 % Shutdown disconnect instr

                connect = 0;

                
        end % end switch
        
    end % end while
    
    res = inst.SendScpi(':SYST:ERR?');
    fprintf(1, '\nEnd - server stopped!! \nInstrument Error Status = %s\n', netStrToStr(res.RespStr));
    
    if cType == "LAN"
    inst.Disconnect();
    else
        admin.CloseInstrument(instId);    
        admin.Close();
    end     
clear inst;
clear;

%fclose();

% Function netStrToStr
function str = netStrToStr(netStr)
    try
        str = convertCharsToStrings(char(netStr));
    catch
        str = '';
    end
end

%% Function - On Logger-Event
function OnLoggerEvent(~, e)
    try
        % print only:
        % TaborElec.Proteus.CLI.LogLevel.Critical (1)
        % TaborElec.Proteus.CLI.LogLevel.Error    (2)
        % TaborElec.Proteus.CLI.LogLevel.Warning  (3)
        if int32(e.Level) >= 1 &&  int32(e.Level) <= 3
            msg = netStrToStr(e.Message);
            if strlength(msg) > 0
                fprintf(2, '\n ~ %s\n', msg(1));
            end
        end
        System.Diagnostics.Trace.WriteLine(e.Message);
    catch ME
        rethrow(ME)
    end
end

% case 6 % Program MW chirp waveform
%                 
%                 % ---------------------------------------------------------------------
%                 % Download chirp waveform of 1024 points to segment 1 of channel 1
%                 % ---------------------------------------------------------------------                
%                 global sampleRateDAC
%                 sampleRateDAC = 9e9
%                 awg_center_freq = cmdBytes(2);
%                 awg_bw_freq = cmdBytes(3);
%                 awg_amp = cmdBytes(4);
%                 sweep_freq = cmdBytes(5);
%                 sweep_sigma = cmdBytes(6);
%                 symm = cmdBytes(7);
% %                 srs_freq = cmdBytes(8); % should be 0.3625e9
%                 srs_freq = 0.3625e9; % new value for good chirp
%                 
%                 srs_amp = cmdBytes(9);
%                 pol_times = 60
%                 starting_pol = '-';
% 
% %                 pol_times = [cmdBytes(10) cmdBytes(11) cmdBytes(12) cmdBytes(13) cmdBytes(14) cmdBytes(15)];
% %                 pol_times = nonzeros(pol_times);
% %                 starting_pol_sign = cmdBytes(16);
% %                 
% %                 if starting_pol_sign == 1
% %                     starting_pol = '+';
% %                 elseif starting_pol_sign ==-1
% %                     starting_pol = '-';
% %                 else
% %                     disp('starting polarization value is wrong (only +1 or -1)')
% %                     break
% %                 end
%                 
%        
% %                 awg_center_freq = 3.775e9;
% %                 awg_bw_freq = 24e6;
% %                 awg_amp = 1.2;
% %                 sweep_freq = 750;
% %                 sweep_sigma = 0.1;
% %                 symm = 0;
% %                 srs_freq = 0.3625e9; % new value for good chirp
% %                 srs_amp = 3;
% %                 pol_times = [t1 t2 t3];
% %                 pol_times = nonzeros(pol_times);
% %                 starting_pol_sign = 1;
% %                 
%                 inst.SendScpi("*CLS")
%                 inst.SendScpi("*RST")
%                 res = inst.SendScpi(['INST:CHAN ' num2str(dacChanInd)]); % select channel 2
%                 assert(res.ErrCode == 0);
%                 
%                 % check if its the good variable bc im not sure
%                 inst.SendScpi([':FREQ:RAST ' num2str(sampleRateDAC)]);
%                 assert(res.ErrCode == 0);
%                 
%                 inst.SendScpi(":INIT:CONT OFF");
%                 assert(res.ErrCode == 0);
%                 
%                 inst.SendScpi(":INIT:CONT ON");
%                 assert(res.ErrCode == 0);
%                                 
%                 rampTime = 1/sweep_freq;
%                 fCenter = awg_center_freq - srs_freq;
%                 fStart = fCenter - 0.5*awg_bw_freq;
%                 disp(['fstart = ' num2str(fStart)]);
%                 fStop = fCenter + 0.5*awg_bw_freq;
%                 dt = 1/sampleRateDAC;
%                 
%                 chirps{1}.dacSignal = makeChirp(sampleRateDAC, rampTime, dt, fStart, fStop, bits);   
%                 chirps{2}.dacSignal = fliplr(chirps{1}.dacSignal);
%                 chirps{1}.segLen = length(chirps{1}.dacSignal);
%                 chirps{2}.segLen = length(chirps{2}.dacSignal);
%                 chirps{1}.segm = 1;
%                 chirps{2}.segm = 2;
%                 chirps{1}.srs_freq = 0.3625e9;
%                 chirps{2}.srs_freq = 0.3625e9;
%                 fprintf('waveform length - ');
%                 fprintf(num2str(length(chirps{1}.dacSignal)));
%                 fprintf('\n') ;
% 
% 
% task_list = build_tasktable(inst,pol_times,chirps,sampleRateDAC,'first sign',starting_pol);
% 
% % % Play seg 1
% % res = inst.SendScpi(':SOUR:FUNC:MODE:SEGM 1');
% % assert(res.ErrCode == 0);
%  
% res = inst.SendScpi(':SOUR:VOLT MAX');
% assert(res.ErrCode == 0);
% 
% res = inst.SendScpi(':OUTP ON');
% assert(res.ErrCode == 0);
% 
% inst.SendScpi(':TRIG:ACTIVE:SEL TRG2');
% inst.SendScpi(':TRIG:LEV 1.0');
% inst.SendScpi(':TRIG:ACTIVE:STAT ON');
% 
% % create the task table
% create_task_table(inst,task_list, sampleRateDAC);
% 
% % write task table
% inst.SendScpi(':TASK:COMP:WRITE');
% fprintf(1, 'SEQUENCE CREATED!\n');
% 
% fprintf(1, 'SETTING AWG OUTPUT\n');
% 
% inst.SendScpi(':SOUR:FUNC:MODE TASK');
% 
% res = inst.SendScpi(':OUTP ON');
% assert(res.ErrCode == 0);
%                 
%                 
% % ---------------------------------------------------------------------
% % Play segment 1 in channel 1
% % ---------------------------------------------------------------------
% 
% res = inst.SendScpi(['INST:CHAN ' num2str(dacChanInd)]); % select channel 2
% assert(res.ErrCode == 0);
% 
% res = inst.SendScpi(':SOUR:FUNC:MODE TASK');
% assert(res.ErrCode == 0);
% 
% 
% amp_str = ':SOUR:VOLT MAX';
% res = inst.SendScpi(amp_str);
% assert(res.ErrCode == 0);
% 
% res = inst.SendScpi(':OUTP ON');
% assert(res.ErrCode == 0);
% 
% fprintf('Waveform generated and playing\n');
% 
% function task_list = build_tasktable(inst,pol_times,chirps,sampleRateDAC,varargin)
% % function that create the task liste before creating the task table
% % input: - pol_times : list that contains the different DNP time for each
% % steps. 
% %          - chirps: 2 dimension structure that contains the information
% %          for each chirp, with positive and negative DNP respectively
% %         
% %          - sampleRateDAC: sample rate for the creation of the mw
% %
% % options : -  'first sign' could be '+' or '-' to stqrt with positive or
% % negative DNP. positive DNP is the default value
% %
% %               - 'test' could be 'on' or 'off' depending on if the chirp
% %               structure contains real chirps (on) or arbitrary wave form
% %               for the test (off)
%     
%     test = 'off';
%     chirp_num = 1;
%     if isempty(varargin) == 0
%         for i = length(varargin)/2
%             if strcmp(varargin(i),'first sign') == 1
%                 if strcmp(varargin(i+1),'+') ==1
%                     chirp_num = 1;
%                 elseif strcmp(varargin(i+1),'-')
%                     chirp_num = 2;
%                 end
%             elseif strcmp(varargin(i),'test') == 1
%                 test = varargin(i+1);
%             end
%         end
%     end
%     
%     for i = 1:2
%              % Define segment i 
%         res = inst.SendScpi([':TRAC:DEF ',num2str(chirps{i}.segm),',',num2str(chirps{i}.segLen)]); % define memory location 1 of length dacSignal
%         assert(res.ErrCode == 0);
% 
%         % select segment as the the programmable segment
%         res = inst.SendScpi(sprintf(':TRAC:SEL %d',chirps{i}.segm));
%         assert(res.ErrCode == 0); 
% 
% %         % Download the binary data to segment
% %         res = inst.WriteBinaryData(':TRAC:DATA 0,', chirps{i}.dacSignal);
% %         assert(res.ErrCode == 0);
%         
%         % Download the binary data to segment
%         prefix = ':TRAC:DATA 0,';
%         
%         global bits
%         if (bits==16)
%             myWfm = uint16(chirps{i}.dacSignal);
%             myWfm = typecast(myWfm, 'uint8');
%         else
%             myWfm = uint8(chirps{i}.dacSignal);
%         end
%         
%         res = inst.WriteBinaryData(prefix, myWfm);
%         
%         if strcmp(test,'off') == 1
%             srs_freq_str = [':SOUR:NCO:CFR1 ' sprintf('%0.2e',  sampleRateDAC - chirps{i}.srs_freq)]; %srs_freq
%             res = inst.SendScpi(srs_freq_str);
%             assert(res.ErrCode == 0);
% 
%            inst.SendScpi(':NCO:SIXD1 ON');
% 
%             rc = inst.SendScpi(':SOUR:MODE DUC');
%             assert(rc.ErrCode == 0);
%             
%             try
%                 sampleRateDAC_str = [':FREQ:RAST ' sprintf('%0.2e', sampleRateDAC)];
%                 res = inst.SendScpi(sampleRateDAC_str); % set sample clock
%                 assert(res.ErrCode == 0);
%                 catch
%                 sampleRateDAC_str = [':FREQ:RAST ' sprintf('%0.2e', sampleRateDAC)];
%                 res = inst.SendScpi(sampleRateDAC_str); % set sample clock
%                 assert(res.ErrCode == 0);
%             end
%     end
%     
%     for i = 1:length(pol_times)
%         task_list{i}.len = chirps{chirp_num}.segLen;
%         task_list{i}.segm = chirps{chirp_num}.segm;
%         task_list{i}.dacSignal = chirps{chirp_num}.dacSignal;
%         task_list{i}.time = pol_times(i);
%         if chirp_num == 1 
%             chirp_num = 2;
%         else
%             chirp_num = 1;
%         end
%     end
%     end
% end
% 
% function dacWav = makeChirp(sampleRateDAC, rampTime, dt, fStart, fStop, bits)            
% 
%     t = 0:1/sampleRateDAC:rampTime;
%     dacWave = chirp(t,fStart,rampTime,fStop);
%     seglenTrunk = (floor(length(dacWave)/ 64))*64;
%     dacWave = dacWave(1:seglenTrunk);
%     dacWav = ampScale(bits, dacWave);
% 
% end
% 
% % Scale to FSD
% function dacSignal = ampScale(bits, rawSignal)
%  
%   maxSig = max(rawSignal);
%   verticalScale = ((2^bits)/2)-1;
% 
%   vertScaled = (rawSignal / maxSig) * verticalScale;
%   dacSignal = (vertScaled + verticalScale);
%   %plot(dacSignal);
% 
% %   if bits > 8
% %       dacSignal16 = [];
% %       sigLen = length(dacSignal);
% %       k=1;
% %       for j = 1:2:sigLen*2;
% %         dacSignal16(j) = bitand(dacSignal(k), 255);
% %         dacSignal16(j+1) = bitshift(dacSignal(k),-8);
% %         k = k + 1;
% %       end
% %       dacSignal = dacSignal16;
% %   end
% end
