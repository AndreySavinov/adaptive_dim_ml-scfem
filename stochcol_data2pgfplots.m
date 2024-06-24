%STOCHCOL_DATA2PGFPLOTS creates .dat files to be read by pgfplots latex package
% 
% Provide the name of the file.
% 
% This file will be populated with the data produced at the end of
% running the adaptive code.
% The data have to be loaded in the workspace.
%
% The default data file created contains the following columns:
% - iterations                          (iter)
% - total number of dofs                (dofs)
% - total (direct) error estimate       (error_d)
% - spatial (direct) error estimate     (error_s_d)
% - parametric (direct) error estimate  (error_p_d)
% - effectivity indices                 (eff_ind_d)
%
%   TIFISS scriptfile: AB; 25 Novemeber 2022
% Copyright (c) 2018 A. Bespalov, L. Rocchi

% File name
  filename = 'sl_sc_test.dat';

% Create data structure
  data = [(1:iter)', ...
          (dof(1:iter))',                 ... 
          (error_d_iter(1:iter))',                  ...
          (err_s_d_iter(1:iter))',                  ...
          (err_p_d_iter(1:iter))',                  ...
          (eff_ind_d_iter(1:iter))'];
      
% Open file
  fid = fopen(filename,'w');
  if fid == -1, error('Error opening file!\n'); end        
% Save data
  fprintf(fid,'iter\t dofs\t error_d\t error_s_d\t error_p_d\t eff_ind_d\n');
  fprintf(fid,'%d\t %d\t %.9f\t %.9f\t %.9f\t %.9f\n',data');
% Close file
  fclose(fid);

% end scriptfile