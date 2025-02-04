function [nc,Tmax]=make_pw_nmr_10(pw,tacq,gain,PW2,tof,nc1T)

    %[fid, msg] = fopen(['Z:/vnmrsys/maclib/ga_save'], 'w');
     [fid, msg] = fopen(['\\192.168.0.219\vnmrsys\maclib\ga_save'], 'w');

    
    if fid == -1
        error(msg);
    end
    fprintf(fid, '%s \n');
    str = verbatim;
    %{
"macro for acquisition and subsequent saving of the FID"
"rewriting of command name (ga --> ga2, modified command as new ga)?"
    
seqfil = 'pi_half_interrupt_psl'

$counter = counter
$counter_ct = counter_ct
$pslabel = pslabel
$date = '_'
shell('date +%Y-%m-%d_%H.%M.%OS; cat'):$date
$1 = $date+'_'+$pslabel

//svf('/home/vnmr1/vnmrsys/data/' + $1)    
//svfj('/home/vnmr1/vnmrsys/data/' + $1)
    
if (end_flag = 'y') then
	counter_ct = 1	
	end_flag = 'n'
endif

if ($counter_ct = $counter) then
	counter_ct = 1
	end_flag = 'y'
endif
    
%}
 sw=125e4;
%  sw=700000;
% neco=80;
dtacq=(1/sw)*1e6;
% tacq=32;
% tacq=32;
neco=tacq/dtacq;
acq_time=neco*dtacq;
np=2*262144;
% rof2_rem=mod(pw+acq_time+ddrtc,dtacq);
% rof2=dtacq*2 - rof2_rem;
rof1=10;
rof2=2; %changed from 2 
rof3=3e-6;%17e-6; %34e-6; 
tdead=rof1+rof2+rof3+4;
%nc=9999;%



% Tmax=350;
% Tmax=200;
% Tmax=1000;
%Tmax=3000/4;
%Tmax=5000/4;
%Tmax=1000/4;
%Tmax=50;
Tmax=70000/14; %changed to 16ms, denominator here should be the number of pulses within the brackets in pulse sequence 
%pw2 = 60;%32.5; %pw;%63; %30*9/2; %40*9/2; % 198; % 45; % 477; % pi is estimated to occur at pw = 106us
%tof = 4;
%tof=3.70; %17.46
tof=-tof*1000;
% pi/2 is estimated to occur at pw = 42us
% pw2 = 18; %first pi/2 pulse 
 %Tmax=500;
%nc_max=floor(Tmax*1e3/(neco*dtacq + PW2 + tdead));
nc_max=floor(Tmax*1e3/(neco*dtacq + PW2 + tdead)); %1e3 for denominator us to ms unit conversion

nc1T = nc1T*1000/14;
nc1 = floor(nc1T*1e3/(neco*dtacq + PW2 + tdead));

tpwr = 42;
% pw=30;
%nc=500;
nc=nc_max; 
%nc1=floor(228200/14);
nc2=nc-nc1;
%nc=0; %added
%nc=4;

% nc=1694;



fprintf(fid,'%s \n',str);
fprintf(fid,'%s \n',['tpwr =' num2str(tpwr)]);
fprintf(fid,'%s \n',['pw2 =' num2str(PW2)]);
fprintf(fid,'%s \n',['pw =' num2str(pw)]);
fprintf(fid,'%s \n',['nc =' num2str(nc)]);
fprintf(fid,'%s \n',['nc1 =' num2str(nc1)]);
fprintf(fid,'%s \n',['nc2 =' num2str(nc2)]);
fprintf(fid,'%s \n',['sw =' num2str(sw)]);
% fprintf(fid,'%s \n',['at =' num2str(at)]);
fprintf(fid,'%s \n',['neco =' num2str(neco)]);
fprintf(fid,'%s \n',['gain =' num2str(gain)]);
fprintf(fid,'%s \n',['rof1 =' num2str(rof1)]);
fprintf(fid,'%s \n',['rof2 =' num2str(rof2)]);
fprintf(fid,'%s \n',['rof3 =' num2str(rof3)]);
fprintf(fid,'%s \n',['tof =' num2str(tof)]);
         
   str = verbatim;
    %{
au
wbs('save_current_wft')
wexp('ga_loop')
%}
   fprintf(fid,'%s \n',str);
    
   fclose(fid);
   
end