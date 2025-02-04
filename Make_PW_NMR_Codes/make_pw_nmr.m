function make_pw_nmr(pw,d3,loop_num)

    [fid, msg] = fopen(['Z:/vnmrsys/maclib/ga_save'], 'w');
    if fid == 1
        error(msg);
    end
    fprintf(fid, '%s \n');
    str = verbatim;
    %{
"macro for acquisition and subsequent saving of the FID"
"rewriting of command name (ga --> ga2, modified command as new ga)?"

seqfil = 'test_sequence'
    
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
    
sw 0000

at = 300e-6
np = 150
nf = loop_number
proc = 'ft'   
    
%}
  fprintf(fid,'%s \n',str);
  

   fprintf(fid,'%s \n',['pw =' num2str(pw)]);
   fprintf(fid,'%s \n',['d3 =' num2str(d3)]);
   fprintf(fid,'%s \n',['loop_number =' num2str(loop_num)]);
   fprintf(fid,'%s \n',['nf =' num2str(loop_num)]); 
   fprintf(fid,'%s \n',['cf = 1' ]); 
   
   str = verbatim;
    %{
au
werr('when_error')
wbs('save_current_wft')
wexp('ga_loop')
%}
   fprintf(fid,'%s \n',str);

   fclose(fid);
   
end