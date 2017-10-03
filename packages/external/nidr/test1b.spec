KEYWORD12 variables {N_vam3(start,0,stop)} TAG "Variables:top"
	[ beta_uncertain INTEGER > 0 {N_vam(intz,numBetaUncVars)}
	  alphas ALIAS buv_alphas REALLIST {N_vam(RealLb,betaUncAlphas)}
	  betas ALIAS buv_betas REALLIST {N_vam(RealLb,betaUncBetas)}
	  lower_bounds ALIAS buv_lower_bounds REALLIST {N_vam(RealLd,betaUncLowerBnds)}
	  upper_bounds ALIAS buv_upper_bounds REALLIST {N_vam(RealLd,betaUncUpperBnds)}
	  [ descriptors ALIAS buv_descriptors STRINGLIST {N_vae(ulbl,UncVar_beta)} ]
	  ]
