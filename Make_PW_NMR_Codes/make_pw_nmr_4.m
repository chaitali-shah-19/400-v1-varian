function make_pw_nmr_4(pw)

    [fid, msg] = fopen(['Z:/vnmrsys/maclib/ga_save'], 'w');
    if fid == 1
        error(msg);
    end
    fprintf(fid, '%s \n');
    str = verbatim;
    %{
"macro for acquisition and subsequent saving of the FID"
"rewriting of command name (ga --> ga2, modified command as new ga)?"
    
seqfil = 'cpmgexp'

$counter = counter
$counter_ct = counter_ct
$pslabel = pslabel
$date = '_'
shell('date +%Y-%m-%d_%H.%M.%OS; cat'):$date
$1 = $date+'_'+$pslabel

folder_name = '/home/ashok/FID_spectra/'+ $1	
shell('mkdir '+ folder_name)


if (end_flag = 'y') then
	counter_ct = 1	
	end_flag = 'n'
endif

if ($counter_ct = $counter) then
	counter_ct = 1
	end_flag = 'y'
endif
    
%}
%sw=125e4;
%sw=625000;
%sw=3e7;
% neco=80;
np=2*262144;
Tmax=500;
at = Tmax/1000;
sw=np/(2*at);


% valid sw values are given by: 
% sw = 10^7*(F*2^(E+1)+2^(2+5*(E-1)))^-1
% for E = [1 2 3]
% and for F = [0 1 2 .. 25 27 29] when E = 1
% and for F = [0 1 2 .. 255] when E = 2
% and for F = [0 1 2 .. 2869] when E = 3
% and additionally sw = 5e6
valid_sw = [5e6];
f_array_e_1 = [1:25 27 29]; 
f_array_e_2 = [1:255];
f_array_e_3 = [1:2869];

e_value = 1;
for f = f_array_e_1
    sw_calc = 10^7*(f*2^(e_value+1)+2^(2+5*(e_value-1)))^-1;
    valid_sw = [valid_sw sw_calc];
end

e_value = 2;
for f = f_array_e_2
    sw_calc = 10^7*(f*2^(e_value+1)+2^(2+5*(e_value-1)))^-1;
    valid_sw = [valid_sw sw_calc];
end

e_value = 3;
for f = f_array_e_3
    sw_calc = 10^7*(f*2^(e_value+1)+2^(2+5*(e_value-1)))^-1;
    valid_sw = [valid_sw sw_calc];
end



disp(['requested sw:' num2str(sw)]);
disp(['requested Tmax:' num2str(Tmax)]);

[val,idx]=min(abs(valid_sw-sw));
sw=valid_sw(idx);


at=np/(2*sw);
Tmax = at *1000;
disp(['adjusted sw:',num2str(sw)]);
disp(['adjusted Tmax:',num2str(Tmax)]);
disp(['number of points:',num2str(np)]);

dtacq=(1/sw)*1e6
tacq=32;
neco=tacq/dtacq
acq_time=neco*dtacq;

% rof2_rem=mod(pw+acq_time+ddrtc,dtacq);
% rof2=dtacq*2 - rof2_rem;
rof1=2;
rof2=2;
tdead=rof1+rof2+2;
nc=934;



nc_max=floor(Tmax*1e3/(neco*dtacq + pw+tdead));
nc=nc_max;

% nc=1694;

fprintf(fid,'%s \n',str);

fprintf(fid,'%s \n',['pw =' num2str(pw)]);
fprintf(fid,'%s \n',['nc =' num2str(nc)]);
fprintf(fid,'%s \n',['pw2 =' num2str(9*21.28)]);
%fprintf(fid,'%s \n',['np =' num2str(np)]);
fprintf(fid,'%s \n',['sw =' num2str(sw)]);
fprintf(fid,'%s \n',['at =' num2str(Tmax)]);
fprintf(fid,'%s \n',['neco =' num2str(neco)]);
fprintf(fid,'%s \n',['rof1 =' num2str(rof1)]);
fprintf(fid,'%s \n',['rof2 =' num2str(rof2)]);
% fprintf(fid,'%s \n',['rof3 =' num2str(rof3)]);
         
   str = verbatim;
    %{
au
wbs('save_current_wft')
wexp('ga_loop')
%}
   fprintf(fid,'%s \n',str);
    
   fclose(fid);
   
end