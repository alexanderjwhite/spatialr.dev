for lam_param in `echo 0 0.01 0.03 0.08 0.22 0.60 1.67 4.64 12.92 35.94 100`; do
	for comp_param in `echo none u v`; do 
		sed "s/comp_param=comp_param/comp_param=\"$comp_param\"/g" sim_spatial_alex.R|sed "s/lam_param=lam_param/lam_param=$lam_param/g" >sim_spatial_alex.$comp_param.$lam_param.R; 
		sed "s/sub.R/sim_spatial_alex.$comp_param.$lam_param.R/g" sub.sh > sim_spatial_alex.$comp_param.$lam_param.R.sh; 
		sbatch sim_spatial_alex.$comp_param.$lam_param.R.sh;
	done
done

