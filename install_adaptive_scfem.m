%%% Run it to fit functionality to your system
if ispc
    system('copy .\stochcol_testproblems_pc.m  .\stochcol_testproblems.m');
    system('copy .\goafem\goafem_stochcol_adaptive_global_settings_pc.m  .\goafem\goafem_stochcol_adaptive_global_settings.m');
elseif ismac || isunix
    !/bin/cp   ./stochcol_testproblems_unix.m  ./stochcol_testproblems.m
    !/bin/cp   ./goafem/goafem_stochcol_adaptive_global_settings_unix.m  ./goafem/goafem_stochcol_adaptive_global_settings.m
else 
    error('System is not supported');
end