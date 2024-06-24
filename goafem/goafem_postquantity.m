% GOAFEM_POSTQUANTITY resizes quantities to remove superfluous iterations
%
% TR; 14 July 2022
err_p_iter_u(iter+1:end) = [];
err_s_iter_u(iter+1:end) = []; 
error_iter_u(iter+1:end) = [];
err_p_iter_z(iter+1:end) = [];
err_s_iter_z(iter+1:end) = [];
error_iter_z(iter+1:end) = [];
error_iter(iter+1:end) = [];

B_iter(iter+1:end) = [];
Goal_unc_iter(iter+1:end) = [];
Goal_iter(iter+1:end) = [];

err_p_dir_iter_u(iter+1:end) = [];
err_s_dir_iter_u(iter+1:end) = [];
error_dir_iter_u(iter+1:end) = [];
err_p_dir_iter_z(iter+1:end) = [];
err_s_dir_iter_z(iter+1:end) = [];
error_dir_iter_z(iter+1:end) = [];
error_dir_iter(iter+1:end) = [];

sols_u_iter(iter+1:end) = [];
sols_z_iter(iter+1:end) = [];
paras_sg_iter(iter+1:end) = [];
paras_fem_iter(iter+1:end) = [];
dof(iter+1:end) = [];